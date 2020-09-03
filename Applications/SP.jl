function read_data(instance::AbstractString)
  data = Data()
  file = readdlm(string(dirname(@__FILE__),"/../data/SP/",instance))
  println("********** Read instance $(file[2,1]) ****************")
  data.nV = file[4,1]
  data.n = file[6,1]
  coord = file[8:(7+data.nV),1:2]
  data.A = file[(9+data.nV):(8+data.nV+data.n),1:data.nV]
  data.s = file[10+data.nV+data.n,1]
  data.t = file[12+data.nV+data.n,1]
  data.b=zeros(Int64,data.nV)
  for i in 1:data.nV
    data.b[i]=0
    if i==data.s
      data.b[i]=1
    elseif i==data.t
      data.b[i]=-1
    end
  end
  data.from = Vector{Int}(undef,data.n)
  data.to = Vector{Int}(undef,data.n)
  for i in 1:data.n
    data.from[i] = findall(data.A[i,1:data.nV] .== 1)[1]
    data.to[i] = findall(data.A[i,1:data.nV] .== -1)[1]
  end
  data.δ⁺ = Array{Set,1}(undef,data.nV)
  data.δ⁻ = Array{Set,1}(undef,data.nV)
  for i in 1:data.nV
    data.δ⁺[i] = Set()
    data.δ⁻[i] = Set()
  end
  for i in 1:data.n
    push!(data.δ⁺[data.from[i]],i)
    push!(data.δ⁻[data.to[i]],i)
  end
  data.d = ones(data.n)
  data.hc = ones(data.n)
  for i in 1:data.n
    data.hc[i] = norm(coord[data.from[i],1:2]-coord[data.to[i],1:2])
    data.d[i] = data.hc[i]/2
  end
  return data
end

#-----------------------------------------------------------------------------------

mutable struct Data
  n::Int
  hc::Vector{Float64} # nominal cost
  d::Vector{Float64} # deviation
  from::Vector{Int}
  to::Vector{Int}
  δ⁺::Array{Set,1} #outgoing arcs
  δ⁻::Array{Set,1} #incoming arcs
  nV::Int
  A::Matrix{Int}
  s::Int
  t::Int
  b::Vector{Int} #RHS
  Data() = new()
end

#-----------------------------------------------------------------------------------

function get_all_SP(sp,ub)
  function getAllPathsUtil(arc::Int, u::Int, t::Int, visited::BitArray{1}, path::BitArray{1}, cost::Float64)
    # Mark the current node as visited and store in path
    if cost + sp[u] < ub
      visited[u] = true
      if arc > 0
        path[arc] = true
      end
      # If current vertex is same as destination, then store current path
      if u == t
        push!(X,copy(path))
      else
        for new_arc in data.δ⁺[u]
          if !visited[data.to[new_arc]]
            cost = cost + data.hc[new_arc]
            getAllPathsUtil(new_arc, data.to[new_arc], t, visited, path, cost)
            cost = cost - data.hc[new_arc]
          end
        end
      end
      # Remove current vertex from path[] and mark it as unvisited
      if arc > 0
        path[arc] = false
      end
      visited[u] = false
      return X
    end
  end
  X = Vector{BitArray{1}}()
  path = falses(data.n)
  cost = 0.0
  # Mark all the vertices as not visited
  visited = falses(data.nV)
  # Call the recursive helper function to print all paths
  getAllPathsUtil(0, data.s, data.t, visited, path, cost)
  return X
end

#-----------------------------------------------------------------------------------

function generate_solutions(ub)
  g = DiGraph(data.nV)
  for i in 1:data.n
    add_edge!(g, data.from[i], data.to[i])
  end
  dist = zeros(data.nV,data.nV)
  for i in 1:data.n
    dist[data.from[i],data.to[i]] = data.hc[i]
  end
  sp = Vector(undef,data.nV)
  for v in 1:data.nV
    all_dist = dijkstra_shortest_paths(g, v, dist).dists
    sp[v] = all_dist[data.t]
  end
  print("Computing X ")
  @time X_bool = get_all_SP(sp,ub) # set of all feasible vectors for theSP
  println("|X|=$(length(X_bool))")
  return X_bool
end

#-----------------------------------------------------------------------------------

function add_constraints_X(model::Model, x::Array{VariableRef,1})
  @constraint(model, [n in 1:data.nV], sum(x[i] for i in data.δ⁺[n])-sum(x[i] for i in data.δ⁻[n]) ≥ data.b[n])
end

#-----------------------------------------------------------------------------------

"Create cc with k rows"
function populate_scenarios()
  cc = Array{Float64,2}(undef,1,data.n)
  for l in N cc[1,l] = data.hc[l] end
  for i = 1:k-1
    delta = zeros(data.n)
    for j=1:Γ
      delta[i+j]=1
    end
    new_cost = [data.hc[l] + data.d[l]*delta[l] for l in N]
    cc = [cc; new_cost']
  end
  return cc
end

#-----------------------------------------------------------------------------------

function add_scenario(δ)
  return [data.hc[l] + data.d[l]*δ[l] for l in N]
end

#-----------------------------------------------------------------------------------

function build_and_solve_separation(x)
  time_remaining = TIME_LIM-(TIME_GENSOL + TIME_HEURISTIC + TIME_BS + TIME_LAZY_MIP)
  separate = create_model(max(0,time_remaining))
  @variable(separate, z)
  @variable(separate, 0 ≤ δ[N] ≤ 1)
  @constraint(separate, sep_cost[i in K], z <= sum((data.hc[l] + data.d[l]*δ[l])*x[i,l] for l in N))
  @constraint(separate, sum(δ[l] for l in N) <= Γ)
  @objective(separate, Max, z)
  optimize!(separate)
  if termination_status(separate)==MOI.OPTIMAL
    return value.(δ),objective_value(separate)
  else
    return [],-Inf
  end
end

#-----------------------------------------------------------------------------------

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
    @variable(xstep, 0 ≤ w)
    @variable(xstep, 0 ≤ v[i in N])
    @variable(xstep, x[k in K, i in N], Bin, container=Array)
    @objective(xstep,  Min, sum(sum(data.hc[i]*x[k,i]*alphas[k] for i in N) for k in K) + Γ*w + sum(v[i] for i in N))
    @constraint(xstep, dualconst[i in N], w + v[i] ≥ sum(data.d[i]*x[k,i]*alphas[k] for k in K))
    for k in K add_constraints_X(xstep,x[k,:]) end
    timexstep += @elapsed optimize!(xstep)
    if primal_status(xstep) == MOI.FEASIBLE_POINT
      global xsol = JuMP.value.(x)
      if heursol <= JuMP.objective_value(xstep) break end
      heursol = JuMP.objective_value(xstep)
    end
    αstep = create_model(max(0,TIME_LIM-timexstep-timealphastep))
    @variable(αstep, 0 ≤ w)
    @variable(αstep, 0 ≤ α[k in K])
    @variable(αstep, 0 ≤ v[i in N])
    @objective(αstep, Min, sum(sum(data.hc[i]*xsol[k,i]*α[k] for i in N) for k in K) + Γ*w + sum(v[i] for i in N))
    @constraint(αstep, dualconst[i in N], w + v[i] ≥ sum(data.d[i]*xsol[k,i]*α[k] for k in K))
    @constraint(αstep, sum(α[j] for j in K) == 1)
    optimize!(αstep)
    timealphastep += @elapsed optimize!(αstep)
    if primal_status(αstep) == MOI.FEASIBLE_POINT
      alphas=JuMP.value.(α)
    end
    if timexstep+timealphastep > TIME_LIM break end
  end
  return heursol,xsol
end

#-----------------------------------------------------------------------------------

function exact_dualization()
  dualized = create_model(TIME_LIM)
  @variable(dualized, 0 ≤ w)
  @variable(dualized, 0 ≤ v[i in N])
  @variable(dualized, 0 ≤ α[k in K] ≤ 1)
  @variable(dualized, 0 ≤ ρ[k in K,i in N] ≤ 1)
  @variable(dualized, x[k in K, i in N], Bin,container=Array)
  for k in K  add_constraints_X(dualized,x[k,:]) end
  @constraint(dualized, sum(α[k] for k in K) == 1)
  @constraint(dualized,[i in N], w + v[i] ≥ sum(data.d[i]ρ[k,i] for k in K))
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≤ x[k,i])
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≤ α[k])
  @constraint(dualized,[i in N,k in K], ρ[k,i] ≥ x[k,i]+α[k]-1)
  @objective(dualized,Min, sum(sum(data.hc[i]*ρ[k,i] for i in N) for k in K) + Γ*w + sum(v[i] for i in N))
  optimize!(dualized)
  return JuMP.objective_value(dualized), MOI.get(dualized, MOI.RelativeGap()), value.(x)
end
