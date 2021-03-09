" The following initializes the master problem and runs a loop iterating between adding scenarios
and computing a solution for the current subset of scenarios. It is described in Algorithms 1 and 3 of the paper.
  It is assumed throughout that the number of uncertain constraints is equal to data.ncons."

function scenario_generation(heuristic)
  aa, cc = populate_scenarios()
  if !heuristic
    if heuristicmode == "DUAL"
      @timeit to "heuristic" heur,incumbent = heuristic_dualization()
    elseif heuristicmode == "HEURCG"
      @timeit to "heuristic" begin
        lb, heur, gap, lm, incumbent = scenario_generation(true)
      end
    end
    ub = heur
    @info "exact mode"
    println("Heuristic value is $ub")
    @timeit to "gen_sol" solutions = generate_solutions(ub)
  else
    # HEURISTIC MODE: start with one solution and its neighbours
    @info "heuristic mode"
    detsol = find_new_solution(cc[1,:],aa[:][1,:],Inf,[],[],false,false,TIME_LIM)
    incumbent = Array{Int64}(undef, k, n)
    for k in K incumbent[k,:] = round_solution(detsol) end
    @timeit to "gen_sol" solutions = generate_neighbors(round_solution(incumbent),aa)
    @timeit to "static ub" ub = solve_static()
  end
  for kk in (size(solutions,1)+1):k push!(solutions,solutions[1]) end # add useless solutions in case s < k. Could be avoided by adapting the code but these instances are easy anyway
  S = 1:size(solutions,1)
  new_scenario_added = true
  x_sol_array = Vector{Vector{Float64}}()
  sol_items = Vector{Vector{Int}}()
  fill_solutions_bit(solutions, x_sol_array, sol_items)
  lb = minimum(sum(x_sol_array[s][l]*cc[1,l] for l in N) for s in S) # starting lb: min_x∈X u^Tx for some u ∈ U
  M = 1:size(cc,1)
  d = Matrix{Float64}(undef,size(solutions,1),size(cc,1)) # cost of each solution in each scenario
  a_cov = Vector{Matrix{Bool}}() # constraint satisfaction of each solution in each scenario
  for cons in 1:data.ncons
    push!(a_cov,Matrix{Bool}(undef,size(solutions,1),size(cc,1)))
  end
  NUMITER = 0 # count the number of iterations of the scenarios generation loop
  gap = 100
  if TimerOutputs.tottime(to)/10^9 < TIME_LIM
    println("IT LB  UB  gap time")
    while new_scenario_added
      δ = []
      NUMITER += 1
      if heuristic
        # The following might be optimized by updating instead of creating from scratch
        # but I don't think this is bottleneck
        @timeit to "update d" d = calculate_d(x_sol_array,cc,S,M)
        if data.ncons > 0
          @timeit to "update a_cov" a_cov = calculate_a(x_sol_array,aa,S,M)
        end
      else
        @timeit to "update d" d = update_d(x_sol_array,cc,d,S,M,NUMITER)
        if data.ncons > 0
          @timeit to "update a_cov" a_cov = update_a(x_sol_array,aa,a_cov,S,M,NUMITER)
        end
      end
      @timeit to "binary search" x, lb = binary_search(x_sol_array,sol_items, lb, ub, min(gap,1),d,a_cov,cc,δ)
      if length(x) > 0
        @timeit to "separation" δ, sep_value, feasible = build_and_solve_separation(x)
        if feasible && sep_value < ub  && sep_value > -Inf
          ub = round(sep_value, digits = 4)
          incumbent = x
        end
      end
      gap = 100*round((ub - lb)/abs(ub), digits = 4)
      if gap > ϵ && length(δ)>0
          new_cost, new_weights = add_scenario(δ)
          M = 1:(length(M)+1)
          cc = [cc; new_cost']
          for cons in 1:data.ncons
            aa[cons] = [aa[cons]; new_weights[cons]']
          end
      else
        new_scenario_added = false
      end
      if heuristic && new_scenario_added
        # HEURISTIC MODE: add a good solution for the last scenario and its neighbors along with the scenario
        @timeit to "relax and fix" new_solution = relax_and_fix(cc[end,:],aa[:][end,:],ub)
        @debug findall(x->x!=0,new_solution)
        @timeit to "generate neighbors" solutions_to_add = generate_neighbors(reshape(round_solution(new_solution),length(new_solution),1)',aa)
        for sol in solutions_to_add
          if !(sol in solutions)
            push!(solutions,sol)
          end
        end
        x_sol_array = Vector{Vector{Float64}}()
        sol_items = Vector{Vector{Int}}()
        fill_solutions_bit(solutions, x_sol_array, sol_items)
        S = 1:size(solutions,1)
      end
      TotalTime = round(TimerOutputs.tottime(to)/10^9, digits = 2)
      LowerBound = round(lb, digits = 4)
      println("$NUMITER $LowerBound $ub $gap $TotalTime") #print output
      TotalTime > TIME_LIM && break
    end
  end
  return lb, ub, gap, length(M), incumbent
end

#-----------------------------------------------------------------------------------

"The following computes an optimal solution for the current set of scenarios using
a binary search algorithm based on the optimal solution cost, in line with the radius
formulation used for the p-center problem. It is described in Algorithm 2 of the paper."

function binary_search(x_sol_array,sol_items,lb,ub,ϵ_BS,d,a_cov,cc,δ)
  S = 1:length(x_sol_array)
  M = 1:size(cc,1)
  last = false # equal to true if this is the last iteration: the threshold has been reached but a feasible solution is still needed
  ITER_BS = 0
  x = []
  while (ub - lb)/abs(ub) > ϵ_BS/100 || last || ITER_BS == 0 # the absolute value is needed to handle negative objectives
    ITER_BS += 1
    TimerOutputs.tottime(to)/10^9 > TIME_LIM && break
    if last
      r = ub # set r large enough to find a feasble solution
    else
      r = (lb+ub)/2 # otherwise binary search classical update of r
    end
    @timeit to "build cov" cov = [ 2*(d[s,m]-r)/abs(d[s,m]+r) ≤ ϵ_p for s in S, m in M ]
    for cons in 1:data.ncons
      @timeit to "update cov" cov = cov .& a_cov[cons] # this might be sped-up by avoiding repeating 0 elements and using sparse representations
    end
    infeasible = false
    if minimum(length(findall(cov[:,m])) for m in M) == 0
      infeasible = true
    else
      @timeit to "dominance solutions" SS = calculate_nondom_solutions(cov,last,M,S)
      @timeit to "dominance scenarios" MM = calculate_nondom_scenarios(cov,M,S)
      @timeit to "covcering" objval,solutions_selected = build_and_solve_covering_MIP(SS,MM,cov)
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
      last && break
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

function build_and_solve_covering_MIP(SS,MM,cov)
  time_remaining = TIME_LIM - TimerOutputs.tottime(to)/10^9
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

function update_d(x_sol_array,cc,d,S,M,NUMITER)
    if NUMITER == 1
      for s in S, m in M
        d[s,m] = sum(x_sol_array[s][l]*cc[m,l] for l in N)
      end
    else
      new_costs = [sum(x_sol_array[s][l]*cc[length(M),l] for l in N) for s in S]
      d = [d new_costs]
    end
  return d
end

#-----------------------------------------------------------------------------------

function update_a(x_sol_array,aa,a_cov,S,M,NUMITER)
    if NUMITER == 1
      for cons in data.ncons, s in S, m in M
        a_cov[cons][s,m] = sum(x_sol_array[s][l]*aa[cons][m,l] for l in N) ≤ data.B[cons]
      end
    else
      for cons in data.ncons, s in S, m in M
        new_cov = [sum(x_sol_array[s][l]*aa[cons][length(M),l] for l in N) ≤ data.B[cons] for s in S]
        a_cov[cons] = [a_cov[cons] new_cov]
      end
    end
  return a_cov
end

#-----------------------------------------------------------------------------------

function calculate_d(x_sol_array,cc,S,M)
  d = Matrix{Float64}(undef,length(S),length(M))
  for s in S, m in M
    d[s,m] = sum(x_sol_array[s][l]*cc[m,l] for l in N)
  end
  return d
end

#-----------------------------------------------------------------------------------

function calculate_a(x_sol_array,aa,S,M)
  a_cov = Vector{Matrix{Bool}}()
  for cons in 1:data.ncons
    push!(a_cov,Matrix{Bool}(undef,length(S),length(M)))
    for s in S, m in M
      a_cov[cons][s,m] = sum(x_sol_array[s][l]*aa[cons][m,l] for l in N) ≤ data.B[cons]
    end
  end
  return a_cov
end

#-----------------------------------------------------------------------------------

function find_new_solution(c_vec,a_vec,UB,fixedzero,fixedone,relax,bound,TIME)
  sol_model = create_model(TIME)
  if relax
    @variable(sol_model, 0<=x[i in N]<=1, container=Array)
  else
    @variable(sol_model, x[i in N], Bin, container=Array)
  end
  add_constraints_X(sol_model,x)
  if data.ncons > 0
    @constraint(sol_model, [cons in 1:data.ncons], sum(a_vec[cons][i]*x[i] for i in N) ≤ data.B[cons] )
  end
  @constraint(sol_model,[i in fixedzero],x[i]==0)
  @constraint(sol_model,[i in fixedone],x[i]==1)
  if bound @constraint(sol_model,sum(x[i]*c_vec[i] for i in N) ≤ UB) end
  @objective(sol_model,Min,sum(x[i]*c_vec[i] for i in N))
  optimize!(sol_model)
  return value.(x)
end

#-----------------------------------------------------------------------------------

function relax_and_fix(cc,aa,UB)
  fixedzero = []
  fixedone = []
  time_remaining = TIME_LIM-TimerOutputs.tottime(to)/10^9
  time_relax = @elapsed solution = find_new_solution(cc,aa,UB,fixedzero,fixedone,true,true,time_remaining)
  for i in N
    if solution[i]>=1-ϵ_p
      push!(fixedone,i)
    elseif solution[i]<=ϵ_p
      push!(fixedzero,i)
    end
  end
  solution = find_new_solution(cc,aa,UB,fixedzero,fixedone,false,false,(time_remaining-time_relax))
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
