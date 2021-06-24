function rcg()
  LB = -Inf
  UB = Inf

  N = 1:data.n
  master = create_model(TIME_LIM)
  @variable(master, x[k in K, i in N], Bin, container=Array)
  @variable(master, η ≥ LB)
  y = Vector{Vector{VariableRef}}()
  w = Vector{Matrix{VariableRef}}()
  for k in K add_constraints_X(master,x[k,:]) end
  #symmetry breaking constraints
  @constraint(master, [k in 1:(length(K)-1)], sum(i*x[k,i] for i in N) + 1 ≤ sum(i*x[k+1,i] for i in N))
  @objective(master, Min, η)

  incumbent = []
  it = 0
  new_cost = populate_scenarios()[2]
  push!(y, @variable(master, [k in K], base_name="y[$it]", Bin))
  push!(w, @variable(master, [k in K, i in N], base_name="w[$it]", Bin))
  @constraint(master, sum(y[it+1][k] for k in K) == 1)
  @constraint(master, [k in K, i in N], w[it+1][k,i] ≥ x[k,i] + y[it+1][k] -1)
  @constraint(master, [k in K, i in N], w[it+1][k,i] ≤ x[k,i])
  @constraint(master, [k in K, i in N], w[it+1][k,i] ≤ y[it+1][k])
  @constraint(master, η ≥ sum(new_cost[i]*w[it+1][k,i] for k in K, i in N))
  gap = 100
  new_scenario_added = true
  while new_scenario_added && time()-STARTTIME < TIME_LIM
    it += 1

    optimize!(master)
    # the second argument is the cost of the solution for the current set of scenarios. It is only
    # needed when constraint uncertainty arise
    if termination_status(master) == MOI.OPTIMAL
      LB = objective_value(master)
      ξ, objective_slave, feasible = build_and_solve_separation(value.(x),0)
      if length(ξ)>0 && objective_slave < UB
        UB = objective_slave
        incumbent = value.(x)
      end
    else
      break
    end
    gap = 100*round((UB - LB)/abs(UB), digits = 4)
    if gap > ϵ && length(ξ)>0
      (new_cost, new_weight) = add_scenario(ξ)
      push!(y, @variable(master, [k in K], base_name="y[$it]", Bin))
      push!(w, @variable(master, [k in K, i in N], base_name="w[$it]", Bin))
      @constraint(master, sum(y[it+1][k] for k in K) == 1)
      @constraint(master, [k in K, i in N], w[it+1][k,i] ≥ x[k,i] + y[it+1][k] -1)
      @constraint(master, [k in K, i in N], w[it+1][k,i] ≤ x[k,i])
      @constraint(master, [k in K, i in N], w[it+1][k,i] ≤ y[it+1][k])
      @constraint(master, η ≥ sum(new_cost[i]*w[it+1][k,i] for k in K, i in N))
    else new_scenario_added = false
    end
    println("$LB  $UB $gap $it")
  end
  return LB, UB, incumbent, gap, it
end
