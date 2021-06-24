using JuMP
using CPLEX
#using Gurobi
using LightGraphs
using DelimitedFiles
using LinearAlgebra
using Random
using TimerOutputs

include("../src/algorithm.jl")
include("../src/RCG.jl")

const ϵ = 10^-5
const TIME_LIM = 7200
#const Optimizer = "CPLEX"
const Optimizer = "CPLEX"
#const GUROBI_ENV = Gurobi.Env()
const to = TimerOutput()
"For non-linear objective function, different choices arise for the objective of the local search step.
Set the following to LINEAR in case of KP or SP"
#const local_search_mode = "PRODUCT"
const local_search_mode = "LINEAR"

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
    set_optimizer_attribute(M, "TimeLimit", max(0,TIME_LIM))
    set_optimizer_attribute(M, "Threads", 1)
    #M = Model(with_optimizer(Gurobi.Optimizer,GUROBI_ENV,OutputFlag= 0,TimeLimit=max(0,TIME_LIM), Threads=1))
  end
  return M
end
