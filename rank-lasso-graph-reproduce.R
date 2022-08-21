library(knitr)
library(glmnet)
library(robsel)
library(igraph)
library(mltools)
library(huge)
library(flare)
library(mvtnorm)
library(caret)
library(parallel)
library(gridExtra)
numCores <- detectCores()
setwd("~/Desktop/rank-lasso/code")
source("./rank-lasso-source-cpp.R")
source("./rank-lasso-graph-source.R")
fileFolder= "./random-df10/"
num_sim <- 50
n <- 100
################### random Graph ###################################
#################### d = 25 ###################
d = 25
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta
d25.random.data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  d25.random.data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
}
set.seed(num_sim + 10)
d25.random.validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
#Rank Lasso NS
d25.random.graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  d25.random.graphical.rank.lasso.all[[i]] <- rank_lasso_graph(d25.random.data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
d25.random.graphical.rank.lasso.icov <- lapply(d25.random.graphical.rank.lasso.all, get, x="omega")
d25.random.graphical.rank.mcp.icov <- lapply(d25.random.graphical.rank.lasso.all, get, x="omega_scad")
# 
# #TIGER
d25.random.tiger.lambda <- sqrt(log(25) / n)
d25.random.tiger.fit <- lapply(1:num_sim, function(i)
  huge(d25.random.data.list[[i]], lambda=d25.random.tiger.lambda, method = "tiger", sym="and"))
d25.random.tiger.icov <- sapply(d25.random.tiger.fit, get, x="icov")
# 
# #CLIME
d25.random.clime.icov <- lapply(1:num_sim, function(i)
  graph_validation(d25.random.validation.set, sugm(d25.random.data.list[[i]], method = "clime", sym="and",
                                                   nlambda=10, lambda.min.ratio=0.35)))

# #glasso
d25.random.glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(d25.random.validation.set, huge(d25.random.data.list[[i]], method = "glasso",
                                                   nlambda=10, lambda.min.ratio=0.35)))

d25.random.table <- cbind(model="random", d=d, 
                          rbind(
                            graph_est_and_selection(d25.random.clime.icov, true_graph, true_omega, num_sim,
                                                    "CLIME"),
                            graph_est_and_selection(d25.random.glasso.icov, true_graph, true_omega, num_sim,
                                                    "GLasso"),
                            graph_est_and_selection(d25.random.tiger.icov, true_graph, true_omega, num_sim, 
                                                    "TIGER"),
                            graph_est_and_selection(d25.random.graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                    "gRankLasso"),
                            graph_est_and_selection(d25.random.graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                    "gRankMCP")
                          ))


d25.random.table
write.csv(d25.random.table, paste0(fileFolder, "d25-random-df10.csv"))
save(n, d, num_sim, g, d25.random.data.list, d25.random.validation.set,
     d25.random.graphical.rank.lasso.all,
     d25.random.tiger.fit, 
     d25.random.clime.icov, d25.random.glasso.icov,
     file=paste0(fileFolder,"d25-random-df10.RData"))

#################### d = 50 ###################
d = 50
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta
d50.random.data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  d50.random.data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
}
set.seed(num_sim + 10)
d50.random.validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)

# 
# #TIGER
d50.random.tiger.lambda <- sqrt(log(50) / n)
d50.random.tiger.fit <- lapply(1:num_sim, function(i)
  huge(d50.random.data.list[[i]], lambda=d50.random.tiger.lambda, method = "tiger", sym="and"))
d50.random.tiger.icov <- sapply(d50.random.tiger.fit, get, x="icov")
#
# #CLIME
d50.random.clime.icov <- lapply(1:num_sim, function(i)
  graph_validation(d50.random.validation.set, sugm(d50.random.data.list[[i]], method = "clime", sym="and",
                                                   nlambda=10, lambda.min.ratio=0.35)))

# #glasso
d50.random.glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(d50.random.validation.set, huge(d50.random.data.list[[i]], method = "glasso",
                                                   nlambda=10, lambda.min.ratio=0.35)))

#Rank Lasso NS
d50.random.graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  d50.random.graphical.rank.lasso.all[[i]] <- rank_lasso_graph(d50.random.data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
d50.random.graphical.rank.lasso.icov <- lapply(d50.random.graphical.rank.lasso.all, get, x="omega")
d50.random.graphical.rank.mcp.icov <- lapply(d50.random.graphical.rank.lasso.all, get, x="omega_scad")
##############


d50.random.table <- cbind(model="random", d=d, 
                          rbind(
                            graph_est_and_selection(d50.random.clime.icov, true_graph, true_omega, num_sim,
                                                    "CLIME"),
                            graph_est_and_selection(d50.random.glasso.icov, true_graph, true_omega, num_sim,
                                                    "GLasso"),
                            graph_est_and_selection(d50.random.tiger.icov, true_graph, true_omega, num_sim, 
                                                    "TIGER"),
                            graph_est_and_selection(d50.random.graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                    "gRankLasso"),
                            graph_est_and_selection(d50.random.graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                    "gRankMCP")
                          ))


d50.random.table
write.csv(d50.random.table, paste0(fileFolder, "d50-random-df10.csv"))
save(n, d, num_sim, g, d50.random.data.list, d50.random.validation.set,
     d50.random.graphical.rank.lasso.all,
     d50.random.tiger.fit,
     d50.random.clime.icov, d50.random.glasso.icov,
     file=paste0(fileFolder,"d50-random-df10.RData"))


#################### d = 100 ###################
d = 100
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta
d100.random.data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  d100.random.data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
}

set.seed(num_sim + 10)
d100.random.validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
# 
# #TIGER
d100.random.tiger.lambda <- sqrt(log(100) / n)
d100.random.tiger.fit <- lapply(1:num_sim, function(i)
  huge(d100.random.data.list[[i]], lambda=d100.random.tiger.lambda, method = "tiger", sym="and"))
d100.random.tiger.icov <- sapply(d100.random.tiger.fit, get, x="icov")
#
# #CLIME
d100.random.clime.icov <- lapply(1:num_sim, function(i)
  graph_validation(d100.random.validation.set, sugm(d100.random.data.list[[i]], method = "clime", sym="and",
                                                   nlambda=10, lambda.min.ratio=0.35)))

# #glasso
d100.random.glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(d100.random.validation.set, huge(d100.random.data.list[[i]], method = "glasso",
                                                   nlambda=10, lambda.min.ratio=0.35)))

#Rank Lasso NS
d100.random.graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  d100.random.graphical.rank.lasso.all[[i]] <- rank_lasso_graph(d100.random.data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
d100.random.graphical.rank.lasso.icov <- lapply(d100.random.graphical.rank.lasso.all, get, x="omega")
d100.random.graphical.rank.mcp.icov <- lapply(d100.random.graphical.rank.lasso.all, get, x="omega_scad")

###
d100.random.table <- cbind(model="random", d=d, 
                          rbind(
                            graph_est_and_selection(d100.random.clime.icov, true_graph, true_omega, num_sim,
                                                    "CLIME"),
                            graph_est_and_selection(d100.random.glasso.icov, true_graph, true_omega, num_sim,
                                                    "GLasso"),
                            graph_est_and_selection(d100.random.tiger.icov, true_graph, true_omega, num_sim, 
                                                    "TIGER"),
                            graph_est_and_selection(d100.random.graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                    "gRankLasso"),
                            graph_est_and_selection(d100.random.graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                    "gRankMCP")
                          ))


d100.random.table
write.csv(d100.random.table, paste0(fileFolder, "d100-random-df10.csv"))
save(n, d, num_sim, g, d100.random.data.list, d100.random.validation.set,
     d100.random.graphical.rank.lasso.all,
     d100.random.tiger.fit,
     d100.random.clime.icov, d100.random.glasso.icov,
     file=paste0(fileFolder,"d100-random-df10.RData"))

#################### d = 200 ###################
d = 200
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta
d200.random.data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  d200.random.data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
}

set.seed(num_sim + 10)
d200.random.validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)

# 
# #TIGER
d200.random.tiger.lambda <- sqrt(log(200) / n)
d200.random.tiger.fit <- lapply(1:num_sim, function(i)
  huge(d200.random.data.list[[i]], lambda=d200.random.tiger.lambda, method = "tiger", sym="and"))
d200.random.tiger.icov <- sapply(d200.random.tiger.fit, get, x="icov")
#
# #CLIME
d200.random.clime.icov <- lapply(1:num_sim, function(i)
  graph_validation(d200.random.validation.set, sugm(d200.random.data.list[[i]], method = "clime", sym="and",
                                                    nlambda=10, lambda.min.ratio=0.35)))

# #glasso
d200.random.glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(d200.random.validation.set, huge(d200.random.data.list[[i]], method = "glasso",
                                                    nlambda=10, lambda.min.ratio=0.35)))
#Rank Lasso NS
d200.random.graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  d200.random.graphical.rank.lasso.all[[i]] <- rank_lasso_graph(d200.random.data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
d200.random.graphical.rank.lasso.icov <- lapply(d200.random.graphical.rank.lasso.all, get, x="omega")
d200.random.graphical.rank.mcp.icov <- lapply(d200.random.graphical.rank.lasso.all, get, x="omega_scad")

###


d200.random.table <- cbind(model="random", d=d, 
                           rbind(
                             graph_est_and_selection(d200.random.clime.icov, true_graph, true_omega, num_sim,
                                                     "CLIME"),
                             graph_est_and_selection(d200.random.glasso.icov, true_graph, true_omega, num_sim,
                                                     "GLasso"),
                             graph_est_and_selection(d200.random.tiger.icov, true_graph, true_omega, num_sim, 
                                                     "TIGER"),
                             graph_est_and_selection(d200.random.graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                     "gRankLasso"),
                             graph_est_and_selection(d200.random.graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                     "gRankMCP")
                           ))


d200.random.table
write.csv(d200.random.table, paste0(fileFolder, "d200-random-df10.csv"))
save(n, d, num_sim, g, d200.random.data.list, d200.random.validation.set,
     d200.random.graphical.rank.lasso.all,
     d200.random.tiger.fit,
     d200.random.clime.icov, d200.random.glasso.icov,
     file=paste0(fileFolder,"d200-random-df10.RData"))



#################### d = 400 ###################
d = 400
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta
d400.random.data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  d400.random.data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)
}

set.seed(num_sim + 10)
d400.random.validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=10)

# #TIGER
d400.random.tiger.lambda <- sqrt(log(400) / n)
d400.random.tiger.fit <- lapply(1:num_sim, function(i)
  huge(d400.random.data.list[[i]], lambda=d400.random.tiger.lambda, method = "tiger", sym="and"))
d400.random.tiger.icov <- sapply(d400.random.tiger.fit, get, x="icov")
#
# #CLIME
d400.random.clime.icov <- mclapply(1:num_sim, function(i)
  graph_validation(d400.random.validation.set, sugm(d400.random.data.list[[i]], method = "clime", sym="and",
                                                    nlambda=10, lambda.min.ratio=0.35)), mc.cores = 16)

# #glasso
d400.random.glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(d400.random.validation.set, huge(d400.random.data.list[[i]], method = "glasso",
                                                    nlambda=10, lambda.min.ratio=0.35)))

#Rank Lasso NS
d400.random.graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  d400.random.graphical.rank.lasso.all[[i]] <- rank_lasso_graph(d400.random.data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
d400.random.graphical.rank.lasso.icov <- lapply(d400.random.graphical.rank.lasso.all, get, x="omega")
d400.random.graphical.rank.mcp.icov <- lapply(d400.random.graphical.rank.lasso.all, get, x="omega_scad")
# 

d400.random.table <- cbind(model="random", d=d, 
                           rbind(
                             graph_est_and_selection(d400.random.clime.icov, true_graph, true_omega, num_sim,
                                                     "CLIME"),
                             graph_est_and_selection(d400.random.glasso.icov, true_graph, true_omega, num_sim,
                                                     "GLasso"),
                             graph_est_and_selection(d400.random.tiger.icov, true_graph, true_omega, num_sim, 
                                                     "TIGER"),
                             graph_est_and_selection(d400.random.graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                     "gRankLasso"),
                             graph_est_and_selection(d400.random.graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                     "gRankMCP")
                           ))


d400.random.table
write.csv(d400.random.table, paste0(fileFolder, "d400-random-df10.csv"))
save(n, d, num_sim, g, d400.random.data.list, d400.random.validation.set,
     d400.random.graphical.rank.lasso.all,
     d400.random.tiger.fit,
     d400.random.clime.icov, d400.random.glasso.icov,
     file=paste0(fileFolder,"d400-random-df10.RData"))














