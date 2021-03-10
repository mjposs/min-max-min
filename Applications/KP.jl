function generate_data(n, seed, M, density, constraint_uncertainty)
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

  function diamond_uniform(n)
      v = rand(n-1)
      v = [0; sort(v); 1]
      w = [v[i]-v[i-1] for i in 2:(n+1)]
      invert = rand(n)
      for i in 1:n if invert[i] ≤ 0.5 w[i] = -w[i] end end
      return w
  end

  Random.seed!(seed)
  to_print = []
  push!(to_print,"InstanceID")
  push!(to_print,"$(seed)_$(n)_$M")
  push!(to_print,"NumObjects")
  push!(to_print,n)
  push!(to_print,"NumFactors")
  push!(to_print,M)
  push!(to_print,"Weights")
  c = 100*rand(n)
  for i in 1:n
    push!(to_print,c[i])
  end
  push!(to_print,"Budget")
  push!(to_print,sum(c)/2)
  push!(to_print,"Profits")
  for i in 1:n
    push!(to_print,-c[i]/5)
  end
  push!(to_print,"FactorLoading_Φ")
  Φ = Matrix{Float64}(undef, 0, M)
  for i in 1:n
    Φ = [Φ; diamond_uniform(M)']
  end
  for i in 1:n
    push!(to_print,Φ[i,:])
  end
  if constraint_uncertainty
    push!(to_print,"FactorLoading_Ψ")
    Ψ = Matrix{Float64}(undef, 0, M)
    for i in 1:n
      Ψ = [Ψ; diamond_uniform(M)']
    end
    for i in 1:n
      push!(to_print,Φ[i,:])
    end
  end
  push!(to_print,"Conflicts")
  conflicts = random_conflicts()
  for (i,j) in conflicts
    push!(to_print,"$i $j")
  end
  if constraint_uncertainty writedlm(string(dirname(@__FILE__),"/../data/KPneg_weightuncertainty/$(seed)_$(n)_$(M)_$(density).txt"),to_print)
  else writedlm(string(dirname(@__FILE__),"/../data/KPneg/$(seed)_$(n)_$(M)_$(density).txt"),to_print)
  end
end


#-----------------------------------------------------------------------------------

function read_data(instance::AbstractString, constraint_uncertainty::Bool)
  if constraint_uncertainty file = readdlm(string(dirname(@__FILE__),"/../data/KPneg_weightuncertainty/",instance))
  else file = readdlm(string(dirname(@__FILE__),"/../data/KP/",instance))
  end
  @info "read instance $(file[2,1])"
  data = Data()
  data.n = file[4,1]
  data.M = file[6,1]
  data.w = file[8:8+data.n-1,1]
  data.B = Vector{Float64}(undef,1)
  data.B[1] = file[8+data.n+1,1]
  data.hc = file[11+data.n:11+2*data.n-1,1]
  data.Φ = file[12+2*data.n:11+3*data.n,1:4]
  if constraint_uncertainty
    data.ncons = 1
    data.Ψ = file[13+3*data.n:12+4*data.n,1:4]
    data.constraint_uncertainty = true
  else
    data.ncons = 0
    data.constraint_uncertainty = false
  end
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
  B::Vector{Float64} # vector of uncertain rhs. contains one element in this application.
  E::Vector{Vector{Int64}}
  w::Vector{Float64}
  hc::Vector{Float64}
  Φ::Matrix{Float64}
  Ψ::Matrix{Float64}
  conflicts::Vector{Set}
  compatibility::Matrix{Bool}
  ncons::Int #number of uncertain constraints. equal to 1 in this application.
  constraint_uncertainty::Bool
  Data() = new()
end

#-----------------------------------------------------------------------------------

function generate_solutions(ub)
  function gen_sol(x::BitArray,i::Int64,w::Float64,c::Float64,tΦ::Array{Float64,1},tΨ::Array{Float64,1},admissible::BitArray)
    time_elapsed = time() - starting_time
    if time_elapsed < TIME_LIM
      if i == data.n
        if c - sum(abs.(tΦ)) ≤ ub # filter by optimality
          push!(X,copy(x))
        end
      else
        i = i+1
        x[i] = false
        gen_sol(x,i,w,c,tΦ,tΨ,admissible)
        x[i] = true
        if admissible[i] # filter by conflicts constraints
          if w - sum(abs.(tΨ)) ≤ data.B[1] # filter by capacity constraint
            tΦ += 0.5 * data.Φ[i,:] .* data.hc[i]
            if data.constraint_uncertainty
              tΨ += 0.5 * data.Ψ[i,:] .* data.w[i]
            end
            c += data.hc[i]
            w += data.w[i]
            next_admissible = copy(admissible)
            next_admissible[i+1:end] = admissible[i+1:end].* data.compatibility[i,i+1:end]
            gen_sol(x,i,w,c,tΦ,tΨ,next_admissible)
          end
        end
      end
    end
  end
  starting_time = time()
  X = Vector{BitArray{1}}() # set of all feasible vectors for the knapsack
  x = falses(data.n)
  i = 0
  w = 0.0
  c = 0.0
  tΦ = zeros(data.M)
  tΨ = zeros(data.M)
  admissible = trues(data.n)
  println("Computing X")
  gen_sol(x,i,w,c,tΦ,tΨ,admissible)
  time_elapsed = time() - starting_time
  println("|X|=$(length(X))")
  non_dom = trues(length(X)) # given x ∈ X, if ∃ y ∈ X : y ⊃ x, then x is dominated by y
  for i in 1:length(X)
   x = X[i]
   for y in X
     if ⊋(findall(y),findall(x))
       non_dom[i] = false
       break
     end
   end
  end
  X = X[findall(non_dom)]
  println("After dominance check: |X|=$(length(X))")
  return X
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
        if !data.constraint_uncertainty
          # check the knapsack constraint is satisfied
          if data.B[1] - sum(data.w[s] for s in sol_items)+data.w[i]-data.w[iprime]>0
            y[i]=false
            y[iprime]=true
            push!(X,y)
          end
        else
          # add only solutions that satisfy the knapsack constraint for all scenarios
          # generated so far
          if data.B[1] - maximum(sum(aa[1][:,s] for s in sol_items) + aa[1][:,i] - aa[1][:,iprime])>0
            y[i]=false
            y[iprime]=true
            push!(X,y)
          end
        end
      end
    end
  end
  return X
end

#-----------------------------------------------------------------------------------

"Add the constraints forming X"
function add_constraints_X(model::Model, x::Array{VariableRef,1})
  if !data.constraint_uncertainty  # the knapsack constraint is in X
    @constraint(model, sum(data.w[i]*x[i] for i in N) ≤ data.B[1])
  end
  #@constraint(model,[(i,j) in data.conflicts],x[i]+x[j] ≤ 1)
  @constraint(model,[i in 1:data.n],length(data.E[i])*x[i]+sum(x[j] for j in data.E[i])≤ length(data.E[i]))
  #@constraint(model,[i in 1:length(data.conflicts)],x[collect(data.conflicts[i])[1]]+x[collect(data.conflicts[i])[2]] <= 1)
end

#-----------------------------------------------------------------------------------

"Create aa and cc with k rows"
function populate_scenarios()
  cc = Array{Float64,2}(undef,1,data.n)
  aa = Vector{Array{Float64,2}}(undef,1)
  aa[1] = Array{Float64,2}(undef,1,data.n)
  for l in N cc[1,l] = data.hc[l] end
  for i = 1:k-1
    new_cost = [data.hc[l] for l in N] #This might be improved by adding other scenarios
    cc = [cc; new_cost']
  end
  for l in N aa[1][1,l] = data.w[l] end
  for i = 1:k-1
    new_weight = [data.w[l] for l in N] #This might be improved by adding other scenarios
    aa[1] = [aa[1]; new_weight']
  end
  return aa, cc
end

#-----------------------------------------------------------------------------------

function add_scenario(δ)
  new_cost = [(1+sum(data.Φ[i,j]*δ[j] for j in 1:data.M)/2)data.hc[i] for i in N]
  if data.constraint_uncertainty
    new_weights = [ [(1+sum(data.Ψ[i,j]*δ[j] for j in 1:data.M)/2)data.w[i] for i in N] for cons in 1:1 ] # only one constraint is considered
  else
    new_weights = undef
  end
  return new_cost, new_weights
end

#-----------------------------------------------------------------------------------

function build_and_solve_separation(x)
  feasible = true
  if data.constraint_uncertainty
    # Search for a violated constraint. Here this is an LP but could be generalized
    # to more constraints using the indicator-constraints proposed by Subramaniam et al.
    time_remaining = TIME_LIM - TimerOutputs.tottime(to)/10^9
    separate = create_model(max(0,time_remaining))
    @variable(separate, z)
    @variable(separate, -1 ≤ ξ[1:data.M] ≤ 1)
    @constraint(separate, sep_cost[k in K], z <= sum((1+sum(data.Ψ[i,j]*ξ[j] for j in 1:data.M)/2)data.w[i]*x[k,i] for i in N))
    @objective(separate, Max, z)
    @timeit to "constraint" optimize!(separate)
    if termination_status(separate)==MOI.OPTIMAL
      feasible = objective_value(separate) < data.B[1]
    end
  end
  if !data.constraint_uncertainty || feasible
    # Compute the worst-case only if x is feasible
    time_remaining = TIME_LIM - TimerOutputs.tottime(to)/10^9
    separate = create_model(max(0,time_remaining))
    @variable(separate, z)
    @variable(separate, -1 ≤ ξ[1:data.M] ≤ 1)
    @constraint(separate, sep_cost[k in K], z <= sum((1+sum(data.Φ[i,j]*ξ[j] for j in 1:data.M)/2)data.hc[i]*x[k,i] for i in N))
    @objective(separate, Max, z)
    @timeit to "objective" optimize!(separate)
  end
  if termination_status(separate)==MOI.OPTIMAL
    return value.(ξ),objective_value(separate),feasible
  else
    return [],-Inf, false
  end
end

#-----------------------------------------------------------------------------------

"This is the local search heuristic proposed in André B. Chassein, Marc Goerigk, Jannis Kurtz, Michael Poss:
Faster algorithms for min-max-min robustness for combinatorial problems with budgeted uncertainty. Eur. J. Oper. Res. 279(2): 308-319 (2019)"

function heuristic_dualization()
  N = 1:data.n
  heursol = +Inf
  alphas = zeros(k)
  for j in K
    alphas[j] = (2*j)/((k)*(k+1))
  end
  timexstep = 0
  timealphastep = 0
  while true
    xstep = create_model(max(0,TIME_LIM-timexstep-timealphastep))
    @variable(xstep, 0 ≤ u[j in 1:data.M])
    @variable(xstep, 0 ≤ v[j in 1:data.M])
    @variable(xstep, x[k in K, i in N], Bin, container=Array)
    @objective(xstep,  Min, sum(sum(data.hc[i]*x[k,i]*alphas[k] for i in N) for k in K) + sum(v[j]+u[j] for j in 1:data.M))
    @constraint(xstep, dualconst[j in 1:data.M], u[j] - v[j] == sum(sum(data.Φ[i,j]*data.hc[i]/2*x[k,i] for i in N)*alphas[k] for k in K))
    for k in K add_constraints_X(xstep,x[k,:]) end
    timexstep += @elapsed optimize!(xstep)
    if primal_status(xstep) == MOI.FEASIBLE_POINT
      global xsol = value.(x)
      if heursol <= objective_value(xstep) break end
      heursol = objective_value(xstep)
    end
    αstep = create_model(max(0,TIME_LIM-timexstep-timealphastep))
    @variable(αstep, 0 ≤ u[j in 1:data.M])
    @variable(αstep, 0 ≤ v[j in 1:data.M])
    @variable(αstep, 0 ≤ α[k in K])
    @objective(αstep, Min, sum(sum(data.hc[i]*xsol[k,i]*α[k] for i in N) for k in K) + sum(v[j]+u[j] for j in 1:data.M))
    @constraint(αstep, dualconst[j in 1:data.M], u[j] - v[j] == sum(sum(data.Φ[i,j]*data.hc[i]/2*xsol[k,i] for i in N)*α[k] for k in K))
    @constraint(αstep, sum(α[j] for j in K) == 1)
    timealphastep += @elapsed optimize!(αstep)
    if primal_status(αstep) == MOI.FEASIBLE_POINT
      alphas=value.(α)
    end
    if timexstep+timealphastep > TIME_LIM break end
  end
  return heursol,xsol
end

#-----------------------------------------------------------------------------------

"This is the monolithic reformulation from Grani Adiwena Hanasusanto, Daniel Kuhn, Wolfram Wiesemann:
K-Adaptability in Two-Stage Robust Binary Programming. Oper. Res. 63(4): 877-891 (2015)"

function exact_dualization()
  dualized = create_model(TIME_LIM)
  @variable(dualized, 0 ≤ u[j in 1:data.M])
  @variable(dualized, 0 ≤ v[j in 1:data.M])
  @variable(dualized, 0 ≤ α[k in K] ≤ 1)
  @variable(dualized, 0 ≤ ρ[k in K,i in N] ≤ 1)
  @variable(dualized, x[k in K, i in N], Bin,container=Array)
  for k in K  add_constraints_X(dualized,x[k,:]) end
  @constraint(dualized, sum(α[k] for k in K) == 1)
  @constraint(dualized,[j in 1:data.M], u[j] - v[j] == sum(sum(data.Φ[i,j]*data.hc[i]/2*ρ[k,i] for i in N) for k in K))
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≤ x[k,i])
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≤ α[k])
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≥ x[k,i]+α[k]-1)
  @objective(dualized,Min, sum(sum(data.hc[i]*ρ[k,i] for i in N) for k in K) + sum(u[j]+v[j] for j in 1:data.M))
  optimize!(dualized)
  return objective_value(dualized),MOI.get(dualized, MOI.RelativeGap()), value.(x)
end

#-----------------------------------------------------------------------------------

"The heuristic constraint generation requires a starting upper bound, obtained through a classical
dualized model."

function solve_static()
  dualized = create_model(TIME_LIM)
  @variable(dualized, 0 ≤ u[j in 1:data.M])
  @variable(dualized, 0 ≤ v[j in 1:data.M])
  @variable(dualized, x[i in N], Bin,container=Array)
  add_constraints_X(dualized,x)
  @constraint(dualized,[j in 1:data.M], u[j] - v[j] == sum(data.Φ[i,j]*data.hc[i]/2*x[i] for i in N))
  if data.constraint_uncertainty
    @variable(dualized, 0 ≤ u_weight[j in 1:data.M])
    @variable(dualized, 0 ≤ v_weight[j in 1:data.M])
    @constraint(dualized,[j in 1:data.M], u_weight[j] - v_weight[j] == sum(data.Ψ[i,j]*data.w[i]/2*x[i] for i in N))
    @constraint(dualized, sum(data.w[i]*x[i] for i in N) + sum(u_weight[j]+v_weight[j] for j in 1:data.M) ≤ data.B[1])
  end
  @objective(dualized,Min, sum(data.hc[i]*x[i] for i in N) + sum(u[j]+v[j] for j in 1:data.M))
  optimize!(dualized)
  return objective_value(dualized)
end
