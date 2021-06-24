include("utils.jl")
include("../Applications/SP.jl")

n = 20
k = 2
Γ = 3
inst_id = 1
heuristicmode = "DUAL"

instance = "$(inst_id)_$n.txt"
data = read_data(instance)
K = 1:k
N = 1:data.n

# The following excutions avoid warming-up issues
const STARTTIME = time()
@info "Executing functions to warm up ..."
oldstd = stdout
redirect_stdout(open("/dev/null", "w"))
scenario_generation(false,false)
redirect_stdout(oldstd) # recover original stdout
@info "warm-up finished!"

# Now comes the instance that is tested, passed by the arguments
n = parse(Int, ARGS[1])
k = parse(Int, ARGS[2])
inst_id = parse(Int, ARGS[3])
Γ = parse(Int, ARGS[4])

instance = "$(inst_id)_$n.txt"
@info "Solving instance $inst_id with k=$k and Γ=$Γ and ϵ=$ϵ ..."
data = read_data(instance)
K = 1:k
N = 1:data.n

outputfilename = "/../results/SP_scenario_generation_$(n)_$(k)_$(Γ).csv"
reset_timer!(to::TimerOutput)
STARTTIME = time()
runtime = @elapsed lb, ub, solution, gap, NUMBER_GENSC, NUMITER, NUMBER_GENSOL = scenario_generation(false,true)
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $Γ $NUMBER_GENSC $NUMBER_GENSOL $ϵ $runtime $NUMITER $lb $ub $gap $(TimerOutputs.time(to["heuristic"])/10^9) $(TimerOutputs.time(to["gen_sol"])/10^9)"
println(output)
println()
println(outputfile, output)
close(outputfile)
show(to)

# outputfilename = "/../results/SP_rc_generation_$(n)_$(k)_$(Γ).csv"
# reset_timer!(to::TimerOutput)
# STARTTIME = time()
# runtime = @elapsed lb, ub, incumbent, gap, it = rcg()
# outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
# output = "$n $inst_id $k $Γ $runtime $it $lb $ub $gap"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)
# show(to)
