function generate_data(n, seed, density)
  function random_conflicts()
    E = Vector{Set}()
    D = 0
    while D < density
      D = 2*length(E)/(n*(n-1))
      new_edge = Set(rand(1:n,2))
      if !∈(new_edge,E) && length(new_edge) == 2
        push!(E,new_edge)
      end
    end
    return E
  end

  Random.seed!(seed)
  to_print = []
  push!(to_print,"InstanceID")
  push!(to_print,"$(seed)_$(n)")
  push!(to_print,"NumObjects")
  push!(to_print,n)
  push!(to_print,"Probabilities")
  p = rand(0.5:0.01:1.0,n)
  for i in 1:n
    push!(to_print,p[i])
  end
  push!(to_print,"Rewards")
  c = rand(1:1:1000,n)
  for i in 1:n
    push!(to_print,-c[i])
  end
  push!(to_print,"Conflicts")
  conflicts = random_conflicts()
  for (i,j) in conflicts
    push!(to_print,"$i $j")
  end
  writedlm(string(dirname(@__FILE__),"/../data/Smuggler/$(seed)_$(n)_$(density).txt"),to_print)
end


#-----------------------------------------------------------------------------------

function read_data(instance::AbstractString)
  file = readdlm(string(dirname(@__FILE__),"/../data/Smuggler/",instance))
  @info "read instance $(file[2,1])"
  data = Data()
  data.constraint_uncertainty = false
  data.n = file[4,1]
  line = findfirst(file.=="Probabilities")[1]
  data.p = file[line+1:line+data.n,1]
  line = findfirst(file.=="Rewards")[1]
  data.hc = file[line+1:line+data.n,1]
  line = findfirst(file.=="Conflicts")[1]
  data.conflicts = Set.([file[i,1:2] for i in (line+1):size(file)[1]])
  data.compatibility = trues(data.n,data.n)
  for (i,j) in data.conflicts
    data.compatibility[i,j] = false
    data.compatibility[j,i] = false
  end
  data.E = Vector{Vector{Int64}}(undef, data.n)
  for i in 1:data.n
    data.E[i] = findall(x->x!=1,data.compatibility[i,:])
  end
  return data
end

#-----------------------------------------------------------------------------------

mutable struct Data
  n::Int
  M::Int
  E::Vector{Vector{Int64}}
  hc::Vector{Float64}           # nominal rewards
  p::Vector{Float64}            # probabilities
  conflicts::Vector{Set}
  compatibility::Matrix{Bool}
  constraint_uncertainty::Bool
  nX_before_dom::Int
  Data() = new()
end

#-----------------------------------------------------------------------------------

function h(x)
  return prod(data.p[l] for l in x)
end

#-----------------------------------------------------------------------------------

function generate_solutions(ub)
  function gen_sol(x::BitArray,i::Int64,ℓ::Int64,admissible::BitArray)
    time_elapsed = time() - starting_time
    if time_elapsed < TIME_LIM
      selected = findall(x)
      c = sum((1+ρ)*data.hc[selected])
      p = prod(data.p[selected])
      if p*c ≤ ub # weak and quick optimality check
        worstcase = create_model(TIME_LIM)
        @timeit to "build model" begin
          @variable(worstcase, (1+ρ)*data.hc[i] ≤ ξ[i in selected] ≤ (1-ρ)*data.hc[i] )
          @constraint(worstcase, sum( ((ξ[i] - data.hc[i])/data.hc[i])^2 for i in selected) ≤ Ω)
          @objective(worstcase, Min, sum(ξ[i] for i in selected)) # we are being optimistic here
        end
        @timeit to "optimize model" optimize!(worstcase)
        if p*objective_value(worstcase) ≤ ub # filter by optimality
          push!(X,copy(x))
        end
      end
      if i < data.n
        i = i+1
        x[i] = false
        gen_sol(x,i,ℓ,admissible)
        if ℓ < B
          x[i] = true
          if admissible[i] # filter by conflicts constraints
            ℓ += 1
            p *= data.p[i]
            c += (1+ρ)*data.hc[i]
            next_admissible = copy(admissible)
            next_admissible[i+1:end] = admissible[i+1:end].* data.compatibility[i,i+1:end]
            gen_sol(x,i,ℓ,next_admissible)
          end
        end
      end
    end
  end
  starting_time = time()
  X = Vector{BitArray{1}}() # set of all feasible vectors for the knapsack
  x = falses(data.n)
  i = 0
  ℓ = 0
  admissible = trues(data.n)
  println("Computing X")
  gen_sol(x,i,ℓ,admissible)
  println("|X|=$(length(X))")
  data.nX_before_dom = length(X)

  # Next we perform dominance over a full ellipsoid set because the comparison can be made much faster
  non_dom = trues(length(X))
  @timeit to "dominance" begin
    p = zeros(length(X))
    μ = Vector{Dict{Int64,Float64}}(undef,length(X))
    μ_square = Vector{Dict{Int64,Float64}}(undef,length(X))
    sum_μ = Vector{Float64}(undef,length(X))
    sqrtΩ = sqrt(Ω)
    for i in 1:length(X)
      μ[i] = Dict{Int64,Float64}()
      μ_square[i] = Dict{Int64,Float64}()
      x_i = findall(X[i])
      p[i] = prod(data.p[x_i])
      for l in x_i
        μ[i][l] = p[i]*data.hc[l]
        μ_square[i][l] = (μ[i][l])^2
      end
      sum_μ[i] = sum(μ[i][l] for l in x_i)
    end
    for i in 1:length(X)
      x_i = findall(X[i])
      for j in 1:length(X)
        sumsquares = 0
        @timeit to "get xj" x_j = findall(X[j])
        @timeit to "define sumsquare" begin
          for l in x_j
            if X[i] == 1
              sumsquares += (μ[i][l] - μ[j][l])^2
            else
              sumsquares += μ_square[j][l]
            end
          end
          for l in x_i
            if X[j] == 0
              sumsquares += μ_square[i][l]
            end
          end
        end
        # now we maximize μ^T ξ over the ellipsoid, performaning the change of variable η_i = (ξ_i - ĉ_i)/ĉ_i
        constant_term = sum_μ[j] - sum_μ[i]
        norm_term = sqrtΩ*sqrt(sumsquares)
        # We remove solution i if there exists a solution j that is cheaper than solution i over the full elliposid.
        if constant_term + norm_term < - ϵ
          non_dom[i] = false
          break
        end
      end
    end
  end
  X = X[findall(non_dom)]
  println("After dominance check: |X|=$(length(X))")
  return X
end

#-----------------------------------------------------------------------------------

function compute_basic_ub(x)
  worstcase = create_model(TIME_LIM)
  selected = findall(x .> 0.99)
  @variable(worstcase, (1+ρ)*data.hc[i] ≤ ξ[i in selected] ≤ (1-ρ)*data.hc[i] )
  @constraint(worstcase, sum( ((ξ[i] - data.hc[i])/data.hc[i])^2 for i in selected) ≤ Ω )
  @objective(worstcase, Max, sum(ξ[i] for i in selected))
  optimize!(worstcase)
  cost = objective_value(worstcase)
  for i in selected
    cost *= data.p[i]
  end
  return cost
end


#-----------------------------------------------------------------------------------

function generate_neighbors(sol,aa)
  X = Vector{BitArray{1}}()
  for k in 1:size(sol,1)
    x = falses(data.n)
    for j in N
      if sol[k,j]>0 x[j]=true end
    end
    push!(X,x)
    sol_items=findall(x)
    for i in sol_items
      admissible=trues(data.n)
      others = setdiff(sol_items,i)
      for o in others
        admissible[findall(x->x==false,data.compatibility[o,:])].=false
      end
      for iprime in setdiff(findall(admissible),sol_items)
        y = copy(x)
        y[i]=false
        y[iprime]=true
        push!(X,y)
      end
    end
  end
  return X
end

#-----------------------------------------------------------------------------------

"Add the constraints forming X"
function add_constraints_X(model::Model, x::Array{VariableRef,1})
  @constraint(model,[i in 1:data.n],length(data.E[i])*x[i]+sum(x[j] for j in data.E[i])≤ length(data.E[i]))
  @constraint(model,sum(x[i] for i in 1:data.n) ≤ B)
end

#-----------------------------------------------------------------------------------

"Create aa and cc with k rows. Notice aa is empty for this application."
function populate_scenarios()
  cc = Array{Float64,2}(undef,1,data.n)
  aa = Vector{Array{Float64,2}}(undef,1)
  aa[1] = Array{Float64,2}(undef,1,data.n)
  for l in N cc[1,l] = data.hc[l] end
  # for i = 1:k-1
  #   new_cost = [data.hc[l] for l in N] #This might be improved by adding other scenarios
  #   cc = [cc; new_cost']
  # end
  # Notice the first argument returned is aa, which is empty for this application
  return [], cc
end

#-----------------------------------------------------------------------------------

"This function is non-trivial when the uncertainty set is defined through an extended formulation."
function add_scenario(ξ)
  return ξ, []
end

#-----------------------------------------------------------------------------------

function build_and_solve_separation(x,v)
  time_remaining = TIME_LIM - (time()-STARTTIME)
  # Compute the worst-case cost:
  # max{η : ξ \in Ξ, η ≤ f(ξ,x^k), ∀ k in [k]}
  worstcase = create_model(TIME_LIM)
  @variable(worstcase, η)
  @variable(worstcase, (1+ρ)*data.hc[i] ≤ ξ[i in N] ≤ (1-ρ)*data.hc[i])
  @constraint(worstcase, sum(((ξ[i] - data.hc[i])/data.hc[i])^2 for i in N) ≤ Ω )
  @constraint(worstcase, sep_cost[k in K], η ≤ sum(ξ[i]*x[k,i] for i in N)*prod(data.p[i]^x[k,i] for i in N))
  @objective(worstcase, Max, η)
  @timeit to "separation" optimize!(worstcase)
  if termination_status(worstcase) == MOI.OPTIMAL
    return value.(ξ),objective_value(worstcase)
  else
    return [],-Inf
  end
end

#-----------------------------------------------------------------------------------

"The heuristic constraint generation requires a starting upper bound. We provide a simple greedy bound instead
of the true static solution, which would require to solve a MINLP."

function solve_static()
  function callback(cb_data)
      selected = findall([ callback_value(cb_data,x[i]) > 0.9 for i in 1:data.n])
      # compute the worst ξ ∈ Ξ
      worstcase = create_model(TIME_LIM)
      @variable(worstcase, (1+ρ)*data.hc[i] ≤ ξ[i in selected] ≤ (1-ρ)*data.hc[i] )
      @constraint(worstcase, sum( ((ξ[i] - data.hc[i])/data.hc[i])^2 for i in selected) ≤ Ω )
      @objective(worstcase, Max, sum(ξ[i] for i in selected))
      optimize!(worstcase)

      #@warn worstcase

      #@warn "$(objective_value(worstcase)) ? > ? $(callback_value(cb_data,ω) + ϵ)"

      if objective_value(worstcase) > callback_value(cb_data,ω) + ϵ
        # add a cut related to an optimal scenario ξ
        con = @build_constraint( ω ≥ sum(value(ξ[i]) * x[i] for i in selected) + sum(data.hc[i] * x[i] for i in setdiff(1:data.n,selected)) )
        MOI.submit(static, MOI.LazyConstraint(cb_data), con)
        ncuts += 1
      end
   end

  ncuts = 0
  static = create_model(TIME_LIM)
  @variable(static, x[i in N], Bin,container=Array)
  @variable(static, sum(data.hc) ≤ ω ≤ 0 )
  add_constraints_X(static,x)
  @objective(static, Min, ω)
  MOI.set(static, MOI.LazyConstraintCallback(), callback)
  optimize!(static)
  probability = prod(data.p[findall(value.(x) .> 0.9)])
  #@warn probability
  @info "static value is $(probability*objective_value(static)), obtained by generating $ncuts cuts"

  "
  selected = []
  for index in sortperm(data.hc)
    length(selected) == L && break
    # check if index′ compatible with all selected items
    if all([ data.compatibility[index, index′] for index′ in selected ])
      push!(selected,index)
    end
  end

  probability = prod(data.p[selected])

  # compute the worst ξ ∈ Ξ
  worstcase = create_model(TIME_LIM)
  @variable(worstcase, (1-ρ)*data.hc[i] ≤ ξ[i in selected] ≤ (1+ρ)*data.hc[i] )
  @constraint(worstcase, sum( (ξ[i] - data.hc[i])^2 for i in selected) ≤ Ω )
  @objective(worstcase, Max, sum(ξ[i] for i in selected))
  optimize!(worstcase)"

  return probability*objective_value(static)
end
