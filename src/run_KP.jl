include("utils.jl")
include("../Applications/KP.jl")

n = 20
k = 2
inst_id = 1
density = 0.5
ϵ = 0.005
M = 4
heuristicmode = "HEURCG"

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
outputfilename = "/../result/KPscenario_generation_$(n)_$(k)_$(density).csv"
outputfilenameheur = "/../result/KPscenario_generation_heur_$(n)_$(k)_$(density).csv"
outputfilenamedual = "/../result/KPscenario_generation_dual_$(n)_$(k)_$(density).csv"
data = read_data(instance)
K = 1:k
N = 1:data.n

# The following excutions avoid warming-up issues
@elapsed ub, incumbent = heuristic_dualization()
@elapsed lb, gap, incumbent = exact_dualization()
@elapsed lb, ub, gap, nscenarios, incumbent, heur = scenario_generation()
@elapsed lb, ub, gap, nscenarios, incumbent = heuristic_scenario_generation()

# Now comes the instance that is tested, passed by the arguments
n = parse(Int, ARGS[1])
k = parse(Int, ARGS[2])
inst_id = parse(Int, ARGS[3])
density = parse(Float64,ARGS[4])
ϵ = parse(Float64, ARGS[5])
M = 4

instance = "$(inst_id)_$(n)_$(M)_$(density).txt"
outputfilename = "/../results/KPscenario_generation_$(n)_$(k)_$(density).csv"
outputfilenameheur = "/../results/KPscenario_generation_heur_$(n)_$(k)_$(density).csv"
outputfilenamedual = "/../results/KPscenario_generation_dual_$(n)_$(k)_$(density).csv"
outputfilenamedualheur = "/../results/KPscenario_generation_dual_heur_$(n)_$(k)_$(density).csv"
println("Solving instance $inst_id with k=$k and M=$M and density=$density and ϵ=$ϵ ...")
data = read_data(instance)
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
runtime = @elapsed lb, ub, gap, nscenarios, incumbent, heur = scenario_generation()   #
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $M $NUMBER_GENSOL $nscenarios $ϵ $runtime $NUMITER $lb $ub $gap $(TIME_HEURISTIC) $(TIME_GENSOL)"
println(output)
println()
println(outputfile, output)
close(outputfile)
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
