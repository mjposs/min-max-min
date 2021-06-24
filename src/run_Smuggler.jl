include("utils.jl")
include("../Applications/Smuggler.jl")

#ENV["JULIA_DEBUG"] = Main
# for n in 20:20:100
#     for seed in 1:10
#         generate_data(n, seed, 0.5)
#     end
# end

n = 20
k = 2
inst_id = 1
density = 0.5
heuristicmode = "HEURCG"
const constraint_uncertainty = false

instance = "$(inst_id)_$(n)_$(density).txt"
data = read_data(instance)
K = 1:k
N = 1:data.n
B = 5
Ω = 1
ρ = 0.5


# The following excutions avoid warming-up issues
@info "Executing functions to warm up ..."
oldstd = stdout
redirect_stdout(open("/dev/null", "w"))
reset_timer!(to::TimerOutput)
STARTTIME = time()
scenario_generation(false,true)
redirect_stdout(oldstd) # recover original stdout
@info "warm-up finished!"

#Now comes the instance that is tested, passed by the arguments
n = parse(Int, ARGS[1])
k = parse(Int, ARGS[2])
inst_id = parse(Int, ARGS[3])

instance = "$(inst_id)_$(n)_$(density).txt"
outputfilename = "/../results/Smuggler_scenario_generation_$(n)_$(k)_$(B)_$(density).csv"
instance_being_solved = "Solving instance $inst_id with k=$k and density=$density and ϵ=$ϵ ..."
@info instance_being_solved
data = read_data(instance)
K = 1:k
N = 1:data.n

reset_timer!(to::TimerOutput)
STARTTIME = time()
runtime = @elapsed lb, ub, solution, gap, NUMBER_GENSC, NUMITER, NUMBER_GENSOL = scenario_generation(false,false)
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $NUMBER_GENSC $(data.nX_before_dom) $NUMBER_GENSOL $ϵ $runtime $NUMITER $lb $ub $gap $(TimerOutputs.time(to["heuristic"])/10^9) $(TimerOutputs.time(to["gen_sol"])/10^9)"
println(output)
println()
println(outputfile, output)
close(outputfile)
show(to)
