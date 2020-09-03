# min-max-min

This repository contains the algorithms used to solve min-max-min robust optimization problems studied in the paper

*Ayşe Nur Arslan, Michael Poss, and Marco Silva: Min-max-min robust combinatorial optimization with few recourse solutions. Available at XXXX hal reference missing XXXX*

Four algorithms are available:
* The monholithic reformulation from [HKW15](https://doi.org/10.1287/opre.2015.1392 "K-Adaptability in Two-Stage Robust Binary Programming")
* The local search heuristic from [CGKP19](https://doi.org/10.1016/j.ejor.2019.05.045 "Faster algorithms for min-max-min robustness for combinatorial problems with budgeted uncertainty")
* The scenario generation algorithm described in Algorithm 1 of the paper.
* The heuristic variant described in Algorithm 3 of the paper.

## Guide

The code currently contains two applications: the shortest path problem (SP) and the knapsack problem with conflicts (KP). Additional applications can be added by creating the appropriate file. To test one of the two applications, unzip the corresponding data files, and execute the corresponding run file with julia, passing the instance details as parameters.
