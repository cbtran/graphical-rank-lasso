library(knitr)
library(glmnet)
library(robsel)
library(igraph)
library(mltools)
library(BDgraph)
library(CVglasso)
library(flare)
library(mvtnorm)
library(caret)
library(parallel)
library(gridExtra)
numCores <- detectCores()
source("./rank-lasso-source-cpp.R")
source("./rank-lasso-graph-source.R")
fileFolder= "./real-data/"
data( geneExpression ) 
df <- scale(geneExpression, center = T, scale = T)
d <- dim(geneExpression)[2]
n <- dim(geneExpression)[1]
sample.bdmcmc <- bdgraph(data = geneExpression, method = "ggm", algorithm = "bdmcmc", iter = 60000, 
                         g.prior = 0.1, burnin = 30000, save = T)
sample.bdmcmc.summary <- summary(sample.bdmcmc)
sample_selected_graph <- sample.bdmcmc.summary$selected_g
sample_selected_plinks <-  sample.bdmcmc.summary$p_links
sample_select_significant <- BDgraph::select(sample.bdmcmc, cut = 0.6)
sample_select_significant <- sample_select_significant+t(sample_select_significant)
diag(sample_select_significant) = 1
graph.index <- upper.tri(sample_selected_graph)

#Glasso
glasso.fit <- huge(df, lambda.min.ratio = 0.35, nlambda = 10, method="glasso")
glasso.cv <- glasso.chau.cv(glasso.fit)
glasso.cv
glasso.graph <- (glasso.fit$icov[[glasso.cv$opt.idx]] !=0)
sum(glasso.graph[graph.index] != 0)
sum((glasso.graph[graph.index] != 0 ) & sample_select_significant[graph.index] != 0)


#TIGER
tiger <- huge(df, lambda = sqrt(log(d) / n), method = "tiger", sym = "and")
tiger.graph <- (tiger$icov[[1]]!=0)
image(Matrix(tiger.graph))
sum(tiger.icov[graph.index] != 0)
sum((tiger.icov[graph.index] != 0 ) & sample_select_significant[graph.index]!=0)

#Rank Lasso
real.data.rank.fit <- rank_lasso_graph(df)
real.data.rank.omega <- real.data.rank.fit$omega
real.data.rank.mcp <- real.data.rank.fit$omega_scad
real.data.rank.omega.graph <- (real.data.rank.omega!=0)
diag(real.data.rank.omega.graph) = 1
real.data.rank.mcp.graph <- (real.data.rank.mcp!=0)
diag(real.data.rank.mcp.graph) = 1
#Rank
sum(real.data.rank.omega.graph[graph.index] != 0 ) 
sum((real.data.rank.omega.graph[graph.index] != 0 ) & sample_select_significant[graph.index]!=0)

#MCP
sum(real.data.rank.mcp.graph[graph.index] != 0 ) 
sum((real.data.rank.mcp.graph[graph.index] != 0 ) & sample_select_significant[graph.index]!=0)