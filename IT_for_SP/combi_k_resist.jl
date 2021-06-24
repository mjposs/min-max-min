# Compatible with JuMP v0.18

mutable struct Stat
  count1::Int64
  count2::Int64
  count3::Int64
  maxmintime::Float64
  inter_computation_time::Float64
end

function compute_robust_sol_KP(data, Γ, heur, problem)::Float64
  opt = sum(data.hc) + sum(data.d)
  print("Preprocessing: solving robust models ... ")
  if problem == "KP"
    d = copy(data.d)
    push!(d,0.0)
    for dev in d
      model = create_model()
      @variable(model, x[1:data.n], Bin)
      @constraint(model, knp, sum(data.w[i] * x[i] for i in 1:data.n) ≥ data.C)
      @objective(model, Min, sum( ( data.hc[i] + max(0, data.d[i]-dev) )*x[i] for i in 1:data.n ) )
      print("$dev ")
      optimize!(model)
      profit = Γ * dev + objective_value(model)
      if profit < opt
        opt = profit
      end
    end
  else
    d = copy(data.d)
    push!(d,0.0)
    for dev in d
      g = DiGraph(data.nV)
      for i in 1:data.n
        add_edge!(g, data.from[i], data.to[i])
      end
      dist = zeros(data.nV,data.nV)
      for i in 1:data.n
        dist[data.from[i],data.to[i]] = data.hc[i] + max(0, data.d[i] - dev)
      end
      all_dist = dijkstra_shortest_paths(g, data.s, dist).dists

      profit = Γ * dev + all_dist[data.t]
      if profit < opt
        opt = profit
      end
    end
  end

  println("")
  println("Robust solution is $opt")
  println("Alternating heuristic is $heur")

  if heur < opt
    return convert(Float64, heur)
  else
    return convert(Float64, opt)
  end
end

function generate_solutions_KP(data, Γ, solution, heur)
  # compute the optimal solution of the min-max problem
  function gen_sol(x::BitArray,i::Int64,w::Int64,c::Float64,opt::Float64)
    if time() - STARTTIME > TIME_LIMIT
      return 0
    end
    if (i == data.n)
      push!(X,copy(x))
    else
      i = i+1

      x[i] = true
      if c + data.hc[i] < heur # filter by optimality
        c = c + data.hc[i]
        gen_sol(x,i,w,c,heur)
        c = c - data.hc[i]
      end

      x[i] = false
      if w - data.w[i] >= data.C  # filter by feasibility
        w = w - data.w[i]
        gen_sol(x,i,w,c,opt)
      end
    end
  end

  X = Vector{BitArray{1}}() # set of all feasible vectors for the knapsack
  x = falses(data.n)
  i = 0
  w = sum(data.w)
  c = 0.0

  gen_sol(x,i,w,c,heur)

  println("|X|=$(size(X)[1])")
  solution.n_sol = size(X)[1]

  #Y = zeros(0,data.n) # set of maximal feasible vectors for the knapsack
  #for i in 1:size(X)[1]
  #   y = X[i,1:data.n]
  #  dominated = ( (y'data.hc)[1] < opt) # remove immediately solutions whith deterministic cost below the robust optimal solution
  #  for j in 1:size(X)[1]
  #    if dominated
  #      break
  #    end
  #    x = X[j,1:data.n]
  #    dominated = dominate(y,x)
  #  end
  #  if !dominated
  #    Y = [Y; y']
  #  end
  #end
  #println("|Y|=$(size(Y)[1])")
  #solution.n_sol_red = size(Y)[1]

return X
end


function generate_solutions_SP_BB(data, Γ, solution, heur)
  m = length(data.hc)
  n = data.nV
  model = create_model()
  @variable(model, 0 <= f[1:m] <= 1)
  @constraint(model, equilibrium[v in setdiff(1:n,[data.s,data.t])], sum(f[i] for i in data.δ⁺[v]) - sum(f[i] for i in data.δ⁻[v]) == 0 )
  @constraint(model, sum(f[i] for i in data.δ⁺[data.s]) == 1)
  @objective(model, Min, sum(data.hc[i]*f[i] for i in 1:m))

  #print(model)

  Tree = [ Vector{BitArray}() for i in 1:(m+1) ]
  push!(Tree[1],BitArray(0))
  nnodes = 1
  X_bool = Vector{BitArray{1}}()
  starttime = time()

  while true && (time() - starttime < TIME_LIMIT)
    availablenodes = [ i for i in 1:(m+1) if !isempty(Tree[i]) ]
    deepest = 0
    if isempty(availablenodes)
      break
    else
      deepest = maximum(availablenodes)
    end
    currentnode = pop!(Tree[deepest])

    for i in 1:m
      set_lower_bound(f[i], 0.0)
      set_upper_bound(f[i], 1.0)
    end
    for i in 1:length(currentnode)
      if currentnode[i]
        set_lower_bound(f[i], 1.0)
      else
        set_upper_bound(f[i], 0.0)
      end
    end
    optimize!(model)

    relax = objective_value(model)
    if relax < heur
      if deepest == m + 1
        push!(X_bool,copy(currentnode))
      else
        leftbranch = copy(currentnode)
        push!(leftbranch, false)
        push!(Tree[deepest+1],leftbranch)

        rightbranch = copy(currentnode)
        push!(rightbranch, true)
        push!(Tree[deepest+1],rightbranch)

        nnodes = nnodes + 2
      end
    end
  end

  println("** Enumerating paths require $nnodes nodes")
  return X_bool
end

function generate_solutions_SP(data, Γ, solution, heur)
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
  @time X_bool = get_all_paths(data, heur, sp) # set of all feasible vectors for theSP
  println("|X|=$(length(X_bool))")

  #X = BitArray(length(X_bool),data.n)
  #print("creating matrix X ")
  #@time X = [ X_bool[s][i] for s in 1:length(X_bool), i in 1:data.n ]

return X_bool
end


function algo_combi_k_resist(heur,data)

    function populate(types_of_deviations,k,typeod)
        if k == K
            if length(typeod) > 0
                push!(types_of_deviations,typeod)
            end
        else
            k = k+1
            populate(types_of_deviations,k,typeod)
            typeod2 = copy(typeod)
            push!(typeod2,k)
            populate(types_of_deviations,k,typeod2)
        end
    end

    function remove(to_remove)
        all_but_s = setdiff(1:length(costs),to_remove)
        X = X[all_but_s]
        X_set = X_set[all_but_s]
        costs = costs[all_but_s]
        costs_dev = costs_dev[all_but_s]
        costs_Γ = costs_Γ[all_but_s,:]
        max_dev_index = max_dev_index[all_but_s]
    end

    function next(k::Int,next_dev::Int,available_dev::Vector{Float64})::Int
        next_index = next_dev + 1
        if next_index > length(max_dev_index[X_ω[ω[k]][σ[k]]])
            return 0
        elseif available_dev[max_dev_index[X_ω[ω[k]][σ[k]]][next_index]] < 1 - 1/granularity + 0.00001
            return next_index
        else
            next(k,next_index,available_dev)
        end
    end

    function maxmin_discrete(deviations::Array{Int,1},γ::Int)
        ############## body of the function where the actual maxmin is computed
        i = length(deviations)
        cost = Vector(undef,K)
        if i == TOD
            for kk in 1:K
                cost[kk] = costs[X_ω[ω[kk]][σ[kk]]]
            end
            for tod in 1:TOD
                for dev in 1:deviations[tod]
                    for kk in types_of_deviations[tod]
                        #println(types_of_deviations[tod])
                        #println(shared_dev)
                        cost[kk] = cost[kk] + shared_dev[tod,dev]
                    end
                end
            end
            mincost = minimum(cost)
            if mincost > maxmincost
                maxmincost = mincost
            end
            ################ the recursive part to build the possible vectors of deviations ##########
        elseif length(deviations) == TOD - 1
            push!(deviations,min(γ,bounds[i + 1]))
            maxmin_discrete(deviations,0)
        else
            for gamma in 0:min(γ,bounds[i + 1])
                deviations2 = copy(deviations)
                push!(deviations2,gamma)
                maxmin_discrete(deviations2,γ-gamma)
            end
        end
    end

    function maxmin_convex(ω::Vector{Int},σ::Vector{Int})
        dev_indexes = Array{Int,1}()
        for i in N
            for k in 1:K
                if X[X_ω[ω[k]][σ[k]]][i]
                    push!(dev_indexes,i)
                    break
                end
            end
        end

        I = 1:length(dev_indexes)
        Ik = Array{Set{Int},1}(K)
        for k in 1:K
            Ik[k] = Set{Int}()
            for i in 1:length(dev_indexes)
                if X[X_ω[ω[k]][σ[k]]][dev_indexes[i]]
                    push!(Ik[k],i)
                end
            end
        end

        slave = create_model()
        #slave = Model()

        @variable(slave, Ω ≥ 0)
        @variable(slave, c[I] ≥ 0)
        @variable(slave, 0 ≤ δ[I] ≤ 1)

        @constraint(slave, sum(δ[i] for i in I) ≤ Γ)
        @constraint(slave, [k in 1:K], Ω ≤ costs[X_ω[ω[k]][σ[k]]] + sum(data.d[dev_indexes[i]]*δ[i] for i in Ik[k]))
        @objective(slave, Max, Ω)

        optimize!(slave)

        return objective_value(slave)
    end

    function next_continuous(k::Int,current_dev::Int,dev_left::Array{Float64,1})
        next_index = current_dev + 1
        if dev_left[max_dev_index[X_ω[ω[k]][σ[k]]][next_index]] > 0.0001
            return next_index
        else
            next_continuous(k,next_index,dev_left)
        end
    end

    function enumerate_resist(k::Int,statistics::Stat) # enumerate all k-solutions between first and last
        ω_first[k] = Γ*granularity
        if k > 1
            ω_first[k] = ω[k-1]
        end
        ω_last[k] = 1
        if k == K
            ω_last[k] = max( 1,Γ*granularity + 1 - sum(ω[1:K-1]) )
        end

        #println("Before the ω loop for k=$k: ω_first[k]=$(ω_first[k]) and ω_last[k]=$(ω_last[k])")
        for ω[k] in ω_first[k]:-1:ω_last[k]
            #println("ω[$k] = $(ω[k])")
            σ_first[k] = 1
            if k > 1
                if ω[k] == ω[k-1]
                    σ_first[k] = σ[k-1]+1
                end
            end
            σ_last[k] = length(X_ω[ω[k]])

            #println("Before the σ loop for k=$k: σ_first[k]=$(σ_first[k]) and σ_last[k]=$(σ_last[k])")
            for σ[k] in σ_first[k]:σ_last[k]
                #println("σ       $k    $(σ_first[k])     $(σ[k])   $(σ_last[k])    $(length(X_ω[ω[k]]))")
                #println("σ[$k]=$(σ[k])   $BREAK")

                k==1 && BREAK && break
                time() - STARTTIME > TIME_LIMIT && break

                if k == 1 && discrete # allocate the last deviations for the heuristic
                    for γ in 1:rem(Γ,K)
                        costs_dev[X_ω[ω[1]][σ[1]]] = costs_dev[X_ω[ω[1]][σ[1]]] + data.d[max_dev_index[X_ω[ω[1]][σ[1]]][γ+div(Γ,K)]]
                    end
                end

                if k < K
                    enumerate_resist(k+1,statistics)

                    ################ STOP CALLING ITSELF, BODY OF THE FUNCTION
                else
                    start = time()
                    statistics.count1 = statistics.count1 + 1
                    if minimum( [ costs_dev[X_ω[ω[kk]][σ[kk]]] for kk in 1:K ] ) < minmaxmincost
                        # cheap lower bound on the max-min
                        statistics.count2 = statistics.count2 + 1

                        # computing the true max-min is expensive because it requires the interesection between the k solutions
                        # Thus, we compute a greedy lower bound on the max-min
                        c = [costs[X_ω[ω[k]][σ[k]]] for k in 1:K]
                        #if discrete
                        next_dev = fill(1,K)
                        available_dev = zeros(data.n)
                        granularity = 1 # remove this to experiment stronger (but slower) LB in the convex case
                        #number_deviations = [ 0 for kk in 1:K ]
                        for γ in 1:Γ*granularity
                          if sum(next_dev) > 0
                            k_min = 0
                            c_min = 999999999999
                            for kk in 1:K
                              if next_dev[kk] > 0 && c[kk] < c_min
                                k_min = kk
                                c_min = c[kk]
                              end
                            end
                            index = max_dev_index[X_ω[ω[k_min]][σ[k_min]]][next_dev[k_min]]
                            available_dev[index] = available_dev[index] + 1/granularity
                            for kk in 1:K
                              if X[X_ω[ω[kk]][σ[kk]]][index]
                                c[kk] = c[kk] + 1/granularity * data.d[index]
                                if available_dev[index] + 1/granularity > 1 - 1/granularity + 0.00001
                                  #number_deviations[kk] = number_deviations[kk] + 1
                                  next_dev[kk] = next(kk,next_dev[kk],available_dev)
                                end
                              end
                            end
                          else
                            break
                          end
                        end

                        if minimum(c) <= minmaxmincost # compute the true max-min only if the LB is low enough
                            if discrete
                                statistics.count3 = statistics.count3 + 1
                                start_inter_computation_time = time()
                                for tod in 1:TOD
                                    shared_set = []
                                    for k in types_of_deviations[tod]
                                        if length(shared_set) == 0
                                            shared_set = X_set[X_ω[ω[k]][σ[k]]]
                                        else
                                            shared_set = intersect(shared_set,X_set[X_ω[ω[k]][σ[k]]])
                                        end
                                    end
                                    for k in setdiff(1:K,types_of_deviations[tod])
                                        shared_set = setdiff(shared_set,X_set[X_ω[ω[k]][σ[k]]])
                                    end
                                    shared_dev[tod,1:length(shared_set)] = sort([ data.d[i] for i in shared_set ], rev=true)
                                    if (length(types_of_deviations[tod]) == 1)
                                        bounds[tod] = min( length(shared_set) , ω[types_of_deviations[tod][1]] )
                                    else
                                        bounds[tod] = min( length(shared_set) , maximum( [ω[types_of_deviations[tod][kk]] for kk in 1:length(types_of_deviations[tod])] ) )
                                    end
                                end
                                statistics.inter_computation_time = statistics.inter_computation_time + time() - start_inter_computation_time
                                maxmincost = 0
                                startmaxmintime = time()
                                maxmin_discrete(Vector{Int}(),Γ)
                                statistics.maxmintime = statistics.maxmintime + time() - startmaxmintime
                            elseif convex
                                statistics.count3 = statistics.count3 + 1
                                startmaxmintime = time()
                                maxmincost = maxmin_convex(ω,σ)
                                statistics.maxmintime = statistics.maxmintime + time() - startmaxmintime
                            end
                            if maxmincost < minmaxmincost
                                minmaxmincost = maxmincost
                                BREAK = true # BREAK ALL FOR-LOOPS
                            end
                        end
                    end
                end
                if k == 1 && discrete
                    for γ in 1:rem(Γ,K) # desallocate the last deviations for the heuristic
                        costs_dev[X_ω[ω[1]][σ[1]]] = costs_dev[X_ω[ω[1]][σ[1]]] - data.d[max_dev_index[X_ω[ω[1]][σ[1]]][γ+div(Γ,K)]]
                    end
                end
                if k == 1 && !BREAK
                    push!(to_remove,X_ω[ω[1]][σ[1]])
                end
            end
        end
    end

    ############################### some initialization ############################
    maxmincost = 0
    minmaxmincost = sum(data.hc) + sum(data.d)
    heur = compute_robust_sol_KP(data, Γ, heur, problem)
    if minmaxmincost > heur
        minmaxmincost = heur
    end
    sol = solution()
    LP = false
    N = 1:data.n
    if problem == "KP"
        X = generate_solutions_KP(data, Γ, sol, heur)
    elseif problem == "SP"
        X = generate_solutions_SP(data, Γ, sol, heur)
    end
    START = time()
    r = size(X)[1]
    number_initial_sol = r
    S = 1:r
    X_time = time() - STARTTIME
    cheap_order = sortperm(data.d,rev=true)
    original_order = invperm(cheap_order)
    data.d = data.d[cheap_order]
    data.hc = data.hc[cheap_order]

    BREAK = false
    X_set, X_ω, costs, max_dev_index, costs_dev, costs_Γ = 0,0,0,0,0,0
    types_of_deviations, TOD, shared_dev, bounds, pre_computed_bounds, to_remove, prep_time, = 0,0,0,0,0,0,0
    statistics = Stat(0,0,0,0.0,0.0)

    try
      print("Changing order of items in X ")
      @time X = [ X[s][cheap_order] for s in S ]
      print("Define X_set ")
      @time X_set = [ [i for i in N if X[s][i] ] for s in S ]
      print("Define costs ")
      @time costs = [ sum(data.hc[X_set[s]]) for s in S ]
      cheap_order = sortperm(costs)
      print("Reorganize X ")
      @time X = X[cheap_order]
      print("Reorganize X_set ")
      @time X_set = X_set[cheap_order]
      print("Compute costs ")
      @time costs = [ costs[s] for s in cheap_order ]
      print("Compute max_dev_index ")
      @time max_dev_index = [ X_set[s][sortperm(data.d[X_set[s]],rev=true)] for s in S ]
      print("Compute costs_dev ")
      @time costs_dev = [ costs[s] + sum(data.d[max_dev_index[s][1:min(div(Γ,K),length(max_dev_index[s]))]]) for s in S ]
      @time if convex
          for s in S
              if length(max_dev_index[s]) >= div(Γ,K)+1
                  costs_dev[s] = costs_dev[s] + rem(Γ,K)/K * data.d[max_dev_index[s][div(Γ,K)+1]]
              end
          end
      end
      costs_Γ = Matrix{Float64}(undef,r,granularity*Γ)
      print("Compute costs_Γ ")
      @time for s in S, g in 1:granularity*Γ
          γ = g/granularity
          floor_γ = convert(Int,floor(γ))
          frac_γ = γ - floor_γ
          costs_Γ[s,g] = costs[s] + sum(data.d[max_dev_index[s][1:min(floor_γ,length(max_dev_index[s]))]])
          if length(max_dev_index[s]) >= 1+floor_γ
            costs_Γ[s,g] = costs_Γ[s,g] + frac_γ * data.d[max_dev_index[s][1+floor_γ]]
          end
      end

      ############################### start the enumeration algorithm #####################################

      types_of_deviations = Vector{OrderedSet{Int}}()
      populate(types_of_deviations,0,OrderedSet{Int}())
      TOD = length(types_of_deviations)
      shared_dev = Matrix(undef,TOD,data.n)
      bounds = Vector{Int}(undef,TOD)
      pre_computed_bounds = [ length(X_set[s]) for s in 1:r ]
      to_remove = []

      prep_time = time() - START
      ############### Repeat until loop ################
      println("Starting enumeration part ...")
      X_ω = Vector{Vector{Int}}(undef,granularity*Γ)
      while true
          r = size(X)[1]
          S = 1:r
          for g in 1:granularity*Γ
              X_ω[g] = Vector{Int}()
          end
          γ = Vector{Int}(undef,r)
          for s in S
              if costs[s] > minmaxmincost
                  γ[s] = 0
              else
                  #println("$minmaxmincost      $(costs_Γ[s,:])")
                  γ[s] = minimum( [ g for g in 1:granularity*Γ if costs_Γ[s,g] >= minmaxmincost ] ) # Can be sped-up if need be?
                  push!(X_ω[γ[s]],s)
              end
          end
          println("Current best cost is $minmaxmincost")
          #for g in 1:granularity*Γ
          #    println("$g  $(length(X_ω[g]))")
          #end

          BREAK = false
          @time enumerate_resist(1,statistics)

          for s in S
              if γ[s] == 0
                  push!(to_remove,s)
              end
          end
          remove(to_remove)
          to_remove = []

          !BREAK && break
          println(time() - STARTTIME)
      end
    catch error
      println(error)
       if isa(error, OutOfMemoryError)
        statistics.maxmintime = -1
        statistics.inter_computation_time = -1
      end
    end
    ############# End of the repeat loop ###################

    println("optimal solution = $minmaxmincost")

    return (time() - STARTTIME, minmaxmincost, number_initial_sol, X_time, prep_time, statistics.count1, statistics.count2, statistics.count3, statistics.maxmintime, statistics.inter_computation_time)
end
