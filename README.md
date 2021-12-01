# min-max-min

This repository contains the algorithms used to solve min-sup-min robust optimization problems studied in the paper

*Ayşe Nur Arslan, Michael Poss, and Marco Silva: Min-sup-min robust combinatorial optimization with few recourse solutions. Accepted for publication in INFORMS Journal on Computing. Available at https://hal.archives-ouvertes.fr/hal-02939356/document*

The following algorithms are available:
* The monholithic reformulation from [HKW15](https://doi.org/10.1287/opre.2015.1392 "K-Adaptability in Two-Stage Robust Binary Programming"), see function `exact_dualization()` located in each application file
* The local search heuristic from [CGKP19](https://doi.org/10.1016/j.ejor.2019.05.045 "Faster algorithms for min-max-min robustness for combinatorial problems with budgeted uncertainty"), see function `heuristic_dualization()` located in `Applications/SP.jl`
* The scenario generation algorithm from [GKP20](https://doi.org/10.1016/j.dam.2020.07.011 "Min-max-min robustness for combinatorial problems with discrete budgeted uncertainty."), see function `rcg()` located in `src/RCG.jl`
* The scenario generation algorithm (exact and heuristic) described in Algorithm 1 of the paper, see function `scenario_generation()` located in `src/algorithm.jl`
* The iterative algorithm (IT) from [CGKP19](https://doi.org/10.1016/j.ejor.2019.05.045 "Faster algorithms for min-max-min robustness for combinatorial problems with budgeted uncertainty."), see function `algo_combi_k_resist()` located in `IT_for_SP/combi_k_resist.jl`

## Guide

The code currently contains three applications: 
* The shortest path problem (SP) first introduced by [HKW15](https://doi.org/10.1287/opre.2015.1392 "K-Adaptability in Two-Stage Robust Binary Programming")
* The knapsack problem with conflicts (KP), with and without constraint uncertainty, inspired by the capital budgeting from [HKW15](https://doi.org/10.1287/opre.2015.1392 "K-Adaptability in Two-Stage Robust Binary Programming")
* The non-linear smuggler problem (Smuggler), inspired by [GR20](https://doi.org/10.1016/j.dam.2019.08.012 "On the complexity and approximation of the maximum expected value all-or-nothing subset.") 

Additional applications can be added by creating the corresponding files. To test one of the applications, unzip the corresponding data files, and execute the corresponding run file with julia, passing the instance details as parameters.
