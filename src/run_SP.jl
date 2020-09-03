include("utils.jl")
include("../Applications/SP.jl")

n = 20
k = 2
Γ = 3
inst_id = 1
ϵ = 0.005
heuristicmode = "DUAL"

instance = "$(inst_id)_$n.txt"
outputfilename = "/../result/SP_$(n)_$(k)_$(Γ).csv"
outputfilenamedual = "/../result/SP_dual_$(n)_$(k)_$(Γ).csv"
#println("Solving instance $inst_id with k=$k and Γ=$Γ and ϵ=$ϵ ...")
data = read_data(instance)
K = 1:k
N = 1:data.n

# The following excutions avoid warming-up issues
time = @elapsed lb, gap, incumbent = exact_dualization()
time = @elapsed lb, ub, gap, nscenarios, incumbent, heur = scenario_generation()

# Now comes the instance that is tested, passed by the arguments
n = parse(Int, ARGS[1])
k = parse(Int, ARGS[2])
Γ = parse(Int, ARGS[3])
inst_id = parse(Int, ARGS[4])
ϵ = parse(Float64, ARGS[5])

instance = "$(inst_id)_$n.txt"
outputfilename = "/../results/SP_$(n)_$(k)_$(Γ).csv"
outputfilenamedual = "/../results/SP_dual_$(n)_$(k)_$(Γ).csv"
println("Solving instance $inst_id with k=$k and Γ=$Γ and ϵ=$ϵ ...")
data = read_data(instance)
K = 1:k
N = 1:data.n

# time = @elapsed lb, gap, incumbent = exact_dualization() #
# outputfile = open(string(dirname(@__FILE__),outputfilenamedual), "a")
# output = "$n $inst_id $k $Γ $lb $time $gap"
# println(output)
# println()
# println(outputfile, output)
# close(outputfile)

time = @elapsed lb, ub, gap, nscenarios, incumbent, heur = scenario_generation()   #
outputfile = open(string(dirname(@__FILE__),outputfilename), "a")
output = "$n $inst_id $k $Γ $NUMBER_GENSOL $nscenarios $ϵ $time $NUMITER $lb $ub $gap $(TIME_HEURISTIC) $(TIME_GENSOL)"
println(output)
println()
println(outputfile, output)
close(outputfile)
