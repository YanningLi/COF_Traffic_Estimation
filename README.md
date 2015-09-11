Yanning Li, Sep 01, 2015
emlynlyn@gmail.com

## Traffic estimation based on convex optimizatoin framework.

This code can be used for traffic estimation on a network (connection, merge, diverge), given the boudnary conditions. 

# Features:
1. This code supports uneven discretization of both time and space. 

2. The matlab scripts will be too complicated in a large problem. Hence, those code are completely rewritten in classes. 

3. This code changed all indexing to dictionaries using MATLAB struct feature.

# Remarks:
1. This version only supports initial and boundary conditins.

2. This version only support three types of uncontrolled junctions: connection, merge, diverge.

3. A more complte version is under development which supports controlled ramps, model predictive control, internal and density conditions, and various modifications to improve the readability and reliability of the code.


# Structure of the code:
script: This script defines the network and application. Run as main function.

---- initNetwork: the class that contains all necessary topology information and data for a network.

---- optProgram: the class that builds an optimization program and output solutions.

---- postSolution: this class visualize the solution. It also iteratively update the discretization grid to resolve the discretization error.

# Installation:
1. Install CPLEX.

2. Add the following path to MATLAB:
---- (CPLEX installation folder)/IBM/ILOG/CPLEX_Stuidio1262/cplex/matlab
---- (this toolbox folder)
---- (this toolbox folder)/COF_toolbox/matlab
---- (this toolbox folder)/UCBerkeley_LWR_solver/matlab_package/source

Run the examples:
  Simply run the example scripts in MATLAB.
Three exampels are provided, respectively
- fwd_sim_one_link.m: the estimation of traffic density on one link.
- fwd_sim_two_links.m: the estimation of traffic density on two links (a connection).
- fwd_sim_merge.m: the estimation of traffic density on a merge network (two links in and one link out).
