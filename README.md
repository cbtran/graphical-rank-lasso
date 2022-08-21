# Graphical Rank Lasso - code
This repo contains code to reproduce all numerical results from our paper [A Completely Tuning-Free and Robust Approach to Sparse Precision Matrix Estimation](https://proceedings.mlr.press/v162/tran22b.html).

## R files

1. `rank-lasso-source-cpp.R`: Source file to fit linear regression model using Rank Lasso. This file requires C++ source files from folder `src/`

2. `rank-lasso-graph-source.R`: Source file to generate graphs and fit model. 

3. `rank-lasso-graph-reproduce.R`: Code to run simulation studies.

4. `rank-lasso-graph-realdata.R`: Code to run real data example.
