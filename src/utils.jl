using JuMP, Gurobi, LightGraphs, DelimitedFiles, LinearAlgebra, Random

include("../src/algorithm.jl")

TIME_BS = 0 #time spent in the solution of master
TIME_LAZY_MIP = 0  #time spent in the separation routine
TIME_HEURISTIC = 0 #time spent in the heuristic
TIME_GENSOL = 0
NUMITER = 0        #number of iterations of the extended static algorithm
const ϵ_p = 10^-5          #precision epsilon
const TIME_LIM = 7200
const Optimizer = "Gurobi"
const GUROBI_ENV = Gurobi.Env()

#-----------------------------------------------------------------------------------

function create_model(TIME_LIM)
  if Optimizer=="CPLEX"
    M = Model(with_optimizer(CPLEX.Optimizer))
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
