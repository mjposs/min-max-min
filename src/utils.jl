using JuMP, CPLEX, LightGraphs, DelimitedFiles, LinearAlgebra, Random, TimerOutputs

include("../src/algorithm.jl")

const ϵ_p = 10^-5          #precision epsilon
const TIME_LIM = 7200
const Optimizer = "CPLEX"
#const GUROBI_ENV = Gurobi.Env()
const to = TimerOutput()

#-----------------------------------------------------------------------------------

function create_model(TIME_LIM)
  if Optimizer=="CPLEX"
    M = Model(CPLEX.Optimizer)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_THREADS"), 1)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_SCRIND"), 0)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_MIPDISPLAY"), 1)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_TILIM"), max(0,TIME_LIM))
  elseif Optimizer=="Gurobi"
    M = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    set_optimizer_attribute(M, "OutputFlag", 0)
    set_optimizer_attribute(M, "TimeLimit", TIME_LIM)
    set_optimizer_attribute(M, "Threads", 1)
    #M = Model(with_optimizer(Gurobi.Optimizer,GUROBI_ENV,OutputFlag= 0,TimeLimit=max(0,TIME_LIM), Threads=1))
  end
  return M
end
