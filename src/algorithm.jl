" The following initializes the master problem and runs a loop iterating between adding scenarios
and computing a solution for the current subset of scenarios. It is described in Algorithm 1 of the paper."

function scenario_generation()
  if heuristicmode == "DUAL"
    global TIME_HEURISTIC = @elapsed heur,incumbent = heuristic_dualization()
  elseif heuristicmode == "HEURCG"
    global TIME_HEURISTIC = @elapsed lb, heur, gap, lm, incumbent = heuristic_scenario_generation()
  end
  ub = heur
  println("Heuristic value is $ub")
  global TIME_GENSOL = @elapsed solutions = generate_solutions(ub)
  global NUMBER_GENSOL = size(solutions,1)
  for kk in (size(solutions,1)+1):k push!(solutions,solutions[1]) end # add useless solutions in case s < k. Could be avoided by adapting the code but these instances are easy anyway
  S = 1:size(solutions,1)
  cc = populate_scenarios()
  new_scenario_added = true
  x_sol_array = Vector{Vector{Float64}}()
  sol_items = Vector{Vector{Int}}()
  fill_solutions_bit(solutions, x_sol_array, sol_items)
  lb = minimum([sum([x_sol_array[s][l]*cc[1,l] for l in N]) for s in S]) # starting lb: min_x∈X u^Tx for some u ∈ U
  M = 1:size(cc,1)
  d = Array{Float64,2}(undef,size(solutions,1),size(cc,1))
  global NUMITER = 0 # count the number of iterations of the scenarios generation loop
  gap = 100
  if TIME_HEURISTIC + TIME_GENSOL < TIME_LIM
    println("IT LB  UB  gap time")
    while new_scenario_added
      δ = []
      global NUMITER += 1
      d = update_d(x_sol_array,cc,d,S,M)
      global TIME_BS += @elapsed x, lb = binary_search(x_sol_array,sol_items, lb, ub, min(gap,1),d,cc,δ)
      if length(x) > 0
        global TIME_LAZY_MIP += @elapsed δ, sep_value = build_and_solve_separation(x)
        if sep_value < ub  && sep_value > -Inf
          ub = round(sep_value, digits = 4)
          incumbent = x
        end
      end
      gap = 100*round((ub - lb)/abs(ub), digits = 4)
      if gap > ϵ && length(δ)>0
          new_cost = add_scenario(δ)
          M = 1:(length(M)+1)
          cc = [cc; new_cost']
      else
        new_scenario_added = false
      end
      TT = TIME_GENSOL + TIME_HEURISTIC + TIME_BS + TIME_LAZY_MIP
      println("$(NUMITER) $(round(lb, digits = 4)) $ub $gap $TT") #print output
      if TT ≥ TIME_LIM
        break
      end
    end
  end
  return lb, ub, gap, length(M), incumbent, heur
end

#-----------------------------------------------------------------------------------
"Heuristic variant of the previous function, desribed in Algorithm 3 of the paper."

function heuristic_scenario_generation()
    time_ns = 0
    cc = populate_scenarios()
    δ = []
    time_ns += @elapsed detsol = find_new_solution(cc[1,:],Inf,[],[],false,false,TIME_LIM)
    incumbent = Array{Int64}(undef, k, n)
    for k in K incumbent[k,:] = round_solution(detsol) end
    δ, ub = build_and_solve_separation(incumbent)
    println("Heuristic value is $ub")
    global TIME_GENSOL = @elapsed solutions = generate_neighbors(round_solution(incumbent),ub)
    global NUMBER_GENSOL = size(solutions,1)
    for kk in (size(solutions,1)+1):k push!(solutions,solutions[1]) end # add useless solutions in case s < k. Could be avoided by adapting the code but these instances are easy anyway
    S = 1:size(solutions,1)
    new_scenario_added = true
    x_sol_array = Vector{Vector{Float64}}()
    sol_items = Vector{Vector{Int}}()
    fill_solutions_bit(solutions, x_sol_array, sol_items)
    lb = sum([incumbent[1,i]*cc[1,i] for i in N])
    M = 1:size(cc,1)
    global NUMITER = 0 # count the number of iterations of the scenarios generation loop
    gap = 100*round((ub - lb)/abs(ub), digits = 4)
    if TIME_HEURISTIC + TIME_GENSOL < TIME_LIM
        println("IT LB  UB  gap time")
        while new_scenario_added
            global NUMITER += 1
            dtime = @elapsed d = calculate_d(x_sol_array,cc,S,M)
            println("recalculate d time: $dtime")
            global TIME_BS += @elapsed x, lb = binary_search(x_sol_array,sol_items, lb, ub, min(gap,1),d,cc,δ)
            if length(x) > 0
              global TIME_LAZY_MIP += @elapsed δ, sep_value = build_and_solve_separation(x)
              if sep_value < ub  && sep_value > -Inf
                ub = round(sep_value, digits = 4)
                incumbent = x
              end
            end
            gap = 100*round((ub - lb)/abs(ub), digits = 4)
            if gap > ϵ && length(δ)>0
                new_cost = add_scenario(δ)
                cc = [cc; new_cost']
                M = 1:size(cc,1)
            else
              new_scenario_added = false
            end
            new_costs = [sum([x_sol_array[s][l]*cc[length(M),l] for l in N]) for s in S]
            d = [d new_costs]
            UB = minimum(d[:,end])
            if new_scenario_added
              time_ns += @elapsed new_solution = relax_and_fix(cc[end,:],UB)
              println(findall(x->x!=0,new_solution))
              println("find new sol time: $time_ns")
              timesoladd = @elapsed solutions_to_add = generate_neighbors(reshape(round_solution(new_solution),length(new_solution),1)',ub)
              println("add solution time: $timesoladd")
              for sol in solutions_to_add
                if !(sol in solutions)
                  NUMBER_GENSOL += 1
                  push!(solutions,sol)
                end
              end
              x_sol_array = Vector{Vector{Float64}}()
              sol_items = Vector{Vector{Int}}()
              fill_solutions_bit(solutions, x_sol_array, sol_items)
              S = 1:size(solutions,1)
            end
            TT = TIME_GENSOL + TIME_HEURISTIC + TIME_BS + TIME_LAZY_MIP + time_ns
            println("$(NUMITER) $(round(lb, digits = 4)) $ub $gap $TT") #print output
            if TT ≥ TIME_LIM
              break
            end
        end
    end
    return lb, ub, gap, length(M), incumbent
end

#-----------------------------------------------------------------------------------

"The following computes an optimal solution for the current set of scenarios using
a binary search algorithm based on the optimal solution cost, in line with the radius
formulation used for the p-center problem. It is described in Algorithm 2 of the paper."

function binary_search(x_sol_array,sol_items,lb,ub,ϵ_BS,d,cc,δ)
  S = 1:length(x_sol_array)
  M = 1:size(cc,1)
  last = false # equal to true if this is the last iteration: the threshold has been reached but a feasible solution is still needed
  ITER_BS = 0
  global TIME_SCP = 0
  x = []
  while (ub - lb)/abs(ub) > ϵ_BS/100 || last || ITER_BS == 0 # the absolute value is needed to handle negative objectives
    ITER_BS += 1
    if TIME_SCP ≥ TIME_LIM-(TIME_GENSOL + TIME_HEURISTIC + TIME_BS + TIME_LAZY_MIP)
      break
    end
    if last
      r = ub # set r large enough to find a feasble solution
    else
      r = (lb+ub)/2 # otherwise binary search classical update of r
    end
    TIME_SCP += @elapsed cov = [ 2*(d[s,m]-r)/abs(d[s,m]+r) ≤ ϵ_p for s in S, m in M]
    infeasible = false
    if minimum([length(findall(cov[:,m])) for m in M]) == 0
      infeasible = true
    else
      TIME_SCP += @elapsed SS = calculate_nondom_solutions(cov,last,M,S)
      TIME_SCP += @elapsed MM = calculate_nondom_scenarios(cov,M,S)
      TIME_SCP += @elapsed objval,solutions_selected = build_and_solve_covering_MIP(SS,MM,cov,TIME_SCP)
      infeasible = objval > 0.9
    end
    if !infeasible
      ub = r
      minarray=[]
      for kk in 1:k-1
        if kk == 1
          minarray = min.(d[solutions_selected[kk],:],d[solutions_selected[kk+1],:])
        elseif kk ≥ 2
          minarray = min.(minarray,d[solutions_selected[kk+1],:])
        end
      end
      value = maximum(minarray)
      if value ≤ ub
        ub = value
      end
      if last
        x = compute_x_values(S, sol_items, solutions_selected)
        break
      end
    else
      if last break end
      lb = r
    end
    if (ub - lb)/abs(ub) ≤ ϵ_BS/100
      last = true
    end
  end
  return x, lb
end

#-----------------------------------------------------------------------------------

function calculate_nondom_solutions(cov,last,M,S)
  sizes = [length(findall(cov[s,:])) for s in S]
  map = [BitSet([j for j in M if cov[s, j]]) for s in S]
  discarded = falses(length(S))
  for s in 1:(length(S)+1-k)
    if length(S)-length(discarded)==k break end
    if discarded[s] continue end
    for ss in 1:(length(S)+1-k)
      if s == ss continue
      elseif discarded[ss] continue
      elseif issubset(map[ss], map[s])
          discarded[ss] = true
          break
      elseif issubset(map[s], map[ss])
          discarded[s] = true
          break
      end
    end
  end
  SS = [s for s in S if !discarded[s]]
  return SS
end

#-----------------------------------------------------------------------------------

function calculate_nondom_scenarios(cov,M,S)
  sizes = [length(findall(cov[:,m])) for m in M]
  map = [BitSet([s for s in S if cov[s, j]]) for j in M]
  discarded = falses(length(M))
  for m in 1:(length(M)+1-k)
      if discarded[m] continue end
      for mm in 1:(length(M)+1-k)
        if m == mm continue
        elseif sizes[m] < sizes[mm] continue
        elseif discarded[mm] continue
        elseif ⊆(map[mm],map[m])
          discarded[m] = true
          break
        elseif ⊆(map[m],map[mm])
          discarded[mm] = true
          break
        end
      end
  end
  MM = [m for m in M if !discarded[m]]
  return MM
end

#-----------------------------------------------------------------------------------

" For a given radius, test whether the associated covering problem is feasible using
a MILP formulation with slack variables."

function build_and_solve_covering_MIP(SS,MM,cov,TIME_SCP)
  time_remaining = TIME_LIM-(TIME_GENSOL + TIME_HEURISTIC + TIME_BS + TIME_LAZY_MIP + TIME_SCP)
  cover = create_model(time_remaining)
  @variable(cover, y[SS], Bin)
  @variable(cover, slack[MM] ≥ 0)
  @constraint(cover, covering[m in MM], sum(y[s] for s in SS if cov[s,m]) >= 1-slack[m]) # covering constraints
  @constraint(cover, lessthank, sum(y[s] for s in SS) == k) # cardinality constraint
  @objective(cover, Min, sum(slack[m] for m in MM))
  optimize!(cover)
  sol_used = []
  if termination_status(cover) == MOI.OPTIMAL
    for s in SS
      if value.(y)[s] > 0.9
        push!(sol_used,s)
      end
    end
    return objective_value(cover),sol_used
  else
    return Inf,sol_used
  end
end

#-----------------------------------------------------------------------------------

function update_d(x_sol_array,cc,d,S,M)
    if NUMITER == 1
      for s in S, m in M
        d[s,m] = sum([x_sol_array[s][l]*cc[m,l] for l in N])
      end
    else
      new_costs = [sum([x_sol_array[s][l]*cc[length(M),l] for l in N]) for s in S]
      d = [d new_costs]
    end
  return d
end

#-----------------------------------------------------------------------------------

function calculate_d(x_sol_array,cc,S,M)
  d = Array{Float64,2}(undef,size(x_sol_array,1),size(cc,1))
  for s in S, m in M
    d[s,m] = sum([x_sol_array[s][l]*cc[m,l] for l in N])
  end
  return d
end

#-----------------------------------------------------------------------------------

function find_new_solution(c_vec,UB,fixedzero,fixedone,relax,bound,TIME)
  sol_model=create_model(TIME)
  if relax
    @variable(sol_model, 0<=x[i in N]<=1, container=Array)
  else
    @variable(sol_model, x[i in N], Bin, container=Array)
  end
  add_constraints_X(sol_model,x)
  @constraint(sol_model,[i in fixedzero],x[i]==0)
  @constraint(sol_model,[i in fixedone],x[i]==1)
  if bound @constraint(sol_model,sum([x[i]*c_vec[i] for i in N])<=UB) end
  @objective(sol_model,Min,sum([x[i]*c_vec[i] for i in N]))
  optimize!(sol_model)
  return value.(x)
end

#-----------------------------------------------------------------------------------

function relax_and_fix(cc,UB)
  fixedzero = []
  fixedone = []
  time_remaining = TIME_LIM-(TIME_HEURISTIC+TIME_GENSOL+TIME_BS+TIME_LAZY_MIP)
  time_relax = @elapsed solution = find_new_solution(cc,UB,fixedzero,fixedone,true,true,time_remaining)
  for i in N
    if solution[i]>=1-ϵ_p
      push!(fixedone,i)
    elseif solution[i]<=ϵ_p
      push!(fixedzero,i)
    end
  end
  solution = find_new_solution(cc,UB,fixedzero,fixedone,false,false,(time_remaining-time_relax))
  return solution
end


#-----------------------------------------------------------------------------------

function compute_x_values(S, sol_items, solutions_used)::Array{Int64,2}
    x = zeros(k, data.n)
    k_ = 0
    for s in solutions_used
        k_ += 1
        for l in sol_items[S[s]]
            x[k_ ,l] = 1
        end
    end
    return x
end

#-----------------------------------------------------------------------------------

function fill_solutions_bit(solutions, x_sol_array, sol_items)
    S = 1:size(solutions,1)
    for s in S
        sol = zeros(data.n)
        for i in N
            if solutions[s][i]
              sol[i] = 1
            end
        end
        push!(x_sol_array, sol)
        items = Vector{Int}()
        for l in N
            if sol[l] == 1
                push!(items, l)
            end
        end
        push!(sol_items, items)
    end
end

#-----------------------------------------------------------------------------------

function round_solution(sol)
  for i=1:length(sol)
    if sol[i]<=0.00001
      global sol[i] = 0
    elseif sol[i] >=0.99999
      global sol[i]= 1
    end
  end
  return sol
end
