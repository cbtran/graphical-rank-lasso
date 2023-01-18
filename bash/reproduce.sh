#!/bin/bash

num_sim=100
n=100
graph_type_list=("random" "band" "cluster")
d_list=(25 50 100 200 400)
df_list=(3 5 10 "Inf")

for graph in ${graph_type_list[@]}
do
  for d in ${d_list[@]}
  do
    for graph in ${graph_type_list[@]}
    do
        for df in ${df_list[@]}
        do
            Rscript ./rank-lasso-graph-reproduce.R --num_sim=$num_sim --n=$n --d=$d --graph_type=$graph --df=$df
        done
    done
  done
done