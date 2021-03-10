include("utils.jl")
include("../Applications/KP.jl")

n = 100
k = 2
inst_id = 1
density = 0.5
ϵ = 0.005
M = 4
heuristicmode = "HEURCG"
constraint_uncertainty = false

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
data = read_data(instance,constraint_uncertainty)
K = 1:k
N = 1:data.n

# The following excutions avoid warming-up issues
@info "Executing functions to warm up ..."
oldstd = stdout
redirect_stdout(open("/dev/null", "w"))
scenario_generation(false)
#heuristic_dualization()
#exact_dualization()
redirect_stdout(oldstd) # recover original stdout
@info "warm-up finished!"

# Now comes the instance that is tested, passed by the arguments
#n = parse(Int, ARGS[1])
#k = parse(Int, ARGS[2])
#inst_id = parse(Int, ARGS[3])
#density = parse(Float64,ARGS[4])
#ϵ = parse(Float64, ARGS[5])
M = 4

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
if constraint_uncertainty
    suffix = "_WU_"
else
    suffix = ""
end
outputfilename = "/../results/KP$(suffix)scenario_generation_$(n)_$(k)_$(density).csv"
outputfilenameheur = "/../results/KP$(suffix)scenario_generation_heur_$(n)_$(k)_$(density).csv"
outputfilenamedual = "/../results/KP$(suffix)scenario_generation_dual_$(n)_$(k)_$(density).csv"
outputfilenamedualheur = "/../results/KP$(suffix)scenario_generation_dual_heur_$(n)_$(k)_$(density).csv"
instance_being_solved = "Solving instance $inst_id with k=$k and M=$M and density=$density and ϵ=$ϵ "
if constraint_uncertainty
    instance_being_solved +="and constraint uncertainty ..."
else
    instance_being_solved +="and deterministic constraint ..."
end
@info instance_being_solved
data = read_data(instance,constraint_uncertainty)
K = 1:k
N = 1:data.n

# runtime = @elapsed ub, gap, incumbent = exact_dualization() #
# outputfile = open(string(dirname(@__FILE__),outputfilenamedual), "a")
# output = "$n $inst_id $k $M $ub $runtime $gap"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)
#
reset_timer!(to::TimerOutput)
runtime = @elapsed lb, ub, gap, M, NUMITER, NUMBER_GENSOL = scenario_generation(false)
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $M $NUMBER_GENSOL $ϵ $runtime $NUMITER $lb $ub $gap"
println(output)
println()
println(outputfile, output)
close(outputfile)
show(to)
#
# runtime = @elapsed ub,incumbent = heuristic_dualization()
# outputfile = open(string(dirname(@__FILE__),outputfilenamedualheur), "a")
# output = "$n $inst_id $k $M $runtime $ub"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)

# runtime = @elapsed lb, ub, gap, nscenarios, incumbent = heuristic_scenario_generation()   #
# outputfile = open(string(dirname(@__FILE__),outputfilenameheur), "a")
# output = "$n $inst_id $k $M $NUMBER_GENSOL $nscenarios $ϵ $runtime $NUMITER $lb $ub $gap"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)
