suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(robsel))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(qgraph))
suppressPackageStartupMessages(library(huge))
suppressPackageStartupMessages(library(flare))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(mltools))
numCores <- detectCores()
source("./rank-lasso-source-cpp.R")
source("./rank-lasso-graph-source.R")

args <- R.utils::commandArgs(asValues=TRUE)
print(args)
num_sim <- as.integer(args$num_sim)
n <- as.integer(args$n)
graph_type <- args$graph_type
d <- as.integer(args$d)
if (args$df == "Inf") {
    df = Inf
} else {
    df <- as.integer(args$df)
}
fileFolder <-  paste0("./results/d", d, "-", graph_type, "-df", df)
dir.create(fileFolder, recursive = TRUE)

# Generate graph
set.seed(1)
g <- generator(d=d, graph="random", prob=0.05, tau=1.5, v=0.3)
true_omega = g$omega
true_sigma = g$sigma
true_graph = g$theta

# Generate data
data.list <- list()
for (i in 1:num_sim) {
  set.seed(i)
  data.list[[i]] <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=df)
}
set.seed(num_sim + 10)
validation.set <- mvtnorm::rmvt(n=n, sigma=true_sigma, df=df)

# Graphical Rank Lasso
graphical.rank.lasso.all <- list()
for (i in 1:num_sim) {
  graphical.rank.lasso.all[[i]] <- rank_lasso_graph(data.list[[i]])
  cat(paste0("done sim:", i, "\n"))
}
graphical.rank.lasso.icov <- lapply(graphical.rank.lasso.all, get, x="omega")
graphical.rank.mcp.icov <- lapply(graphical.rank.lasso.all, get, x="omega_scad")


# TIGER

tiger.lambda <- sqrt(log(d) / n)
tiger.fit <- lapply(1:num_sim, function(i)
  huge(data.list[[i]], lambda=tiger.lambda, method = "tiger", sym="and"))
tiger.icov <- sapply(tiger.fit, get, x="icov")


# CLIME
clime.icov <- lapply(1:num_sim, function(i)
  graph_validation(validation.set, sugm(data.list[[i]], method = "clime", sym="and",
                                                   nlambda=10, lambda.min.ratio=0.35)))

# Glasso
glasso.icov <- lapply(1:num_sim, function(i)
  graph_validation(validation.set, huge(data.list[[i]], method = "glasso",
                                                   nlambda=10, lambda.min.ratio=0.35)))

table <- cbind(model="random", d=d, 
                          rbind(
                            graph_est_and_selection(clime.icov, true_graph, true_omega, num_sim,
                                                    "CLIME"),
                            graph_est_and_selection(glasso.icov, true_graph, true_omega, num_sim,
                                                    "GLasso"),
                            graph_est_and_selection(tiger.icov, true_graph, true_omega, num_sim, 
                                                    "TIGER"),
                            graph_est_and_selection(graphical.rank.lasso.icov, true_graph, true_omega, num_sim, 
                                                    "gRankLasso"),
                            graph_est_and_selection(graphical.rank.mcp.icov, true_graph, true_omega, num_sim,
                                                    "gRankMCP")
                          ))


write.csv(table, paste0(fileFolder, "/d", d, "-", graph_type, "-df", df, ".csv"))
save(n, d, num_sim, g, data.list, validation.set,
     graphical.rank.lasso.all,
     tiger.fit, 
     clime.icov, 
     glasso.icov,
     file=paste0(fileFolder, "/d", d, "-", graph_type, "-df", df, ".RData"))