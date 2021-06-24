using JuMP,LightGraphs,DataStructures,DelimitedFiles,LinearAlgebra
using Gurobi

include("data_MMM.jl")
include("graph.jl")
include("combi_k_resist.jl")

#-----------------------------------------------------------------------------------

function create_model()
  if Optimizer=="CPLEX"
    M = Model(with_optimizer(CPLEX.Optimizer))
    MOI.set(M, MOI.RawParameter("CPX_PARAM_THREADS"), 1)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_SCRIND"), 0)
    MOI.set(M, MOI.RawParameter("CPX_PARAM_MIPDISPLAY"), 1)
    #MOI.set(M, MOI.RawParameter("CPX_PARAM_TILIM"), max(0,TIME_LIM))
  elseif Optimizer=="Gurobi"
    M = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    set_optimizer_attribute(M, "OutputFlag", 0)
    #set_optimizer_attribute(M, "TimeLimit", TIME_LIM)
    set_optimizer_attribute(M, "Threads", 1)
    #M = Model(with_optimizer(Gurobi.Optimizer,GUROBI_ENV,OutputFlag= 0,TimeLimit=max(0,TIME_LIM), Threads=1))
  end
  return M
end

#-----------------------------------------------------------------------------------

const Optimizer = "Gurobi"
const GUROBI_ENV = Gurobi.Env()

TIME_LIMIT = 7200
problem = ("SP")
separation = "binary"
count1, count2, count3 = 0, 0, 0 # used to assess the quality of the upper bounds on the objective function
minmaxmincost = 0
LB = 0
maxmintime, inter_computation_time, test_time = 0, 0, 0

appfolder = dirname(@__FILE__)

heuristic = Dict()
for K in [2,3,4], Γ in [3,6], n in 20:5:30
  file = readdlm("$(dirname(@__FILE__))/heuristic_results/mip_vs_heuristic_n_$(n)_B_$(Γ)_K_$(K).txt")
  for i in 1:100
      heuristic[K,n,Γ,i] = file[i,2]
  end
end

K = parse(Int, ARGS[1])
Γ = parse(Int, ARGS[2])
discrete = ARGS[3] == "discrete"
convex = ARGS[3] == "convex"
#convex = true
#outputfilename2 = "res/SP_resistLP_$(ARGS[3])_B_$(Γ)_K_$(K).csv"
outputfilename2 = "res/SP_resistLP_convex_B_$(Γ)_K_$(K).csv"
#outputfilename = "res/SP_$(ARGS[3])_B_$(Γ)_K_$(K).csv"
n1,n2 = 20, 30
i1, i2 = 1, 100
if length(ARGS) > 3
    n1 = parse(Int, ARGS[4])
    n2 = parse(Int, ARGS[5])
    i1 = parse(Int, ARGS[6])
    i2 = parse(Int, ARGS[7])
end
granularity = 1 # If greater than 1, means fractional deviations are considered in X_ω, e.g. = 2 => deviations of 0.5 are considered. Keep this equal to 1 for the discrete model, otherwise not valid
if convex
    granularity = 10
else
    granularity = 1
end
#@info "Solving $(ARGS[3]) problems with K=$K and Γ=$Γ, granularity=$granularity"
@info "Solving convex problems with K=$K and Γ=$Γ, granularity=$granularity"
σ = Vector{Int}(undef,K)
ω = Vector{Int}(undef,K)
σ_first = Vector{Int}(undef,K)
ω_first = Vector{Int}(undef,K)
σ_last = Vector{Int}(undef,K)
ω_last = Vector{Int}(undef,K)
@warn "warming up"
STARTTIME = time()
data = read_dataMMM_SP("/../data/SP/1_20.txt")
algo_combi_k_resist(heuristic[2,20,3,1],data)
@warn "main loop"
for n in n1:5:n2, i in i1:i2
    GC.gc()
    instance = "/../data/SP/$(i)_$(n).txt"
    global data = read_dataMMM_SP(instance)
    global starttime = time()
    global STARTTIME = time()


    (runtime, opt, r, X_time, prep_time, count1, count2, count3, maxmintime, inter_computation_time) = algo_combi_k_resist(heuristic[K,n,Γ,i],data)
    solved = runtime < TIME_LIMIT
    if !solved
        runtime = TIME_LIMIT
    end
    output = "$solved $n $i $opt $r $count1 $count2 $count3 $runtime $(round(100*maxmintime/runtime)) $(round(100*inter_computation_time/runtime)) $(round(100*X_time/runtime)) $(round(100*prep_time/runtime))"
    @info output
    outputfile = open(outputfilename2, "a")
    println(outputfile, output)
    close(outputfile)
end
