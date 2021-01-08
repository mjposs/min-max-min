# min-max-min

This repository contains the algorithms used to solve min-max-min robust optimization problems studied in the paper

*Ayşe Nur Arslan, Michael Poss, and Marco Silva: Min-max-min robust combinatorial optimization with few recourse solutions. Available at https://hal.archives-ouvertes.fr/hal-02939356/document*

Four algorithms are available:
* The monholithic reformulation from [HKW15](https://doi.org/10.1287/opre.2015.1392 "K-Adaptability in Two-Stage Robust Binary Programming"), see function `exact_dualization()`
* The local search heuristic from [CGKP19](https://doi.org/10.1016/j.ejor.2019.05.045 "Faster algorithms for min-max-min robustness for combinatorial problems with budgeted uncertainty"), see function `heuristic_dualization()`
* The scenario generation algorithm described in Algorithm 1 of the paper, see function `scenario_generation()`
* The heuristic variant described in Algorithm 3 of the paper, see function `heuristic_scenario_generation()`

## Guide

The code currently contains two applications: the shortest path problem (SP) and the knapsack problem with conflicts (KP). Additional applications can be added by creating the corresponding files. To test one of the two applications, unzip the corresponding data files, and execute the corresponding run file with julia, passing the instance details as parameters. Notice that the current version of the code does not contain the functions in SP.jl to apply `heuristic_scenario_generation()` to SP.
