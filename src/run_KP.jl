include("utils.jl")
include("../Applications/KP.jl")

n = 20
k = 2
inst_id = 1
density = 0.5
M = 4

const heuristicmode = "HEURCG"
#const heuristicmode = "DUAL"
const formulation = "AGG"
#const formulation = "CONF"
const constraint_uncertainty = false

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
data = read_data(instance,constraint_uncertainty)
K = 1:k
N = 1:data.n

# The following executions avoid warming-up issues
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
constraint_uncertainty = parse(Bool, ARGS[4])
formulation = ARGS[5]
M = 4

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
if constraint_uncertainty
    suffix = "_WU_"
else
    suffix = "_"
end
instance_being_solved = "Solving instance $inst_id with k=$k and M=$M and density=$density and ϵ=$ϵ "
if constraint_uncertainty
    instance_being_solved *="and constraint uncertainty ..."
else
    instance_being_solved *="and deterministic constraint ..."
end
@info instance_being_solved
data = read_data(instance,constraint_uncertainty)
K = 1:k
N = 1:data.n

outputfilename = "/../results/KP$(suffix)scenario_generation_$(formulation)_$(n)_$(k).csv"
reset_timer!(to::TimerOutput)
STARTTIME = time()
runtime = @elapsed lb, ub, solution, gap, NUMBER_GENSC, NUMITER, NUMBER_GENSOL = scenario_generation(false,true)
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $NUMBER_GENSC $NUMBER_GENSOL $ϵ $runtime $NUMITER $lb $ub $gap $(TimerOutputs.time(to["heuristic"])/10^9) $(TimerOutputs.time(to["gen_sol"])/10^9)"
println(output)
println()
println(outputfile, output)
close(outputfile)
show(to)

# outputfilename = "/../results/KP$(suffix)rc_generation_$(n)_$(k)_$(density).csv"
# reset_timer!(to::TimerOutput)
# STARTTIME = time()
# runtime = @elapsed lb, ub, incumbent, gap, it = rcg()
# outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
# output = "$n $inst_id $k $M $runtime $it $lb $ub $gap"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)
# show(to)
