suppressPackageStartupMessages(library(qgraph))
rank_lasso_graph <- function(data, solver="gurobi", m=FALSE, verbose=TRUE, parallel=TRUE) {
  n <- dim(data)[1]
  p <- dim(data)[2]

  omega <- matrix(0, p, p)
  omega_scad <- matrix(0, p, p)
  if (verbose) {
    cat("Solving Rank Lasso Graph\n")
  }

  if (parallel) {
    numCores <- detectCores()
    graph.fits <- mclapply(1:p, function(i) rank_lasso_solver(x=data[,-i], y=data[,i], solver=solver, verbose=FALSE), mc.cores = numCores)
  } else {
    graph.fits <- lapply(1:p, function(i) rank_lasso_solver(x=data[,-i], y=data[,i], solver=solver, verbose=FALSE))                                                              
  }
  for (i in 1:p) {
    res <- graph.fits[[i]]
    omega[i,i] <- 1/res$mse
    omega[-i,i] <- -1/res$mse * res$beta
    
    omega_scad[i,i] <- 1/res$mse_scad
    omega_scad[-i,i] <- -1/res$mse_scad * res$beta_scad

    
  }

  omega <- pmin(omega, t(omega))
  omega_scad <- pmin(omega_scad, t(omega_scad))
  returns <- list()
  returns$omega <- omega
  returns$omega_scad <- omega_scad

  return(returns)
}


graph_rev_perf <- function(graph_list, true_graph) {
  num_sims <- length(graph_list)
  grp_ind <- upper.tri(graph_list[[1]])
  graph_pred <- lapply(graph_list, function(x) factor((x!=0)[grp_ind], levels=c("FALSE", "TRUE")))
  graph_ref <- factor((true_graph!=0)[grp_ind], levels=c("FALSE", "TRUE"))
  confMat_list <- lapply(1:num_sims, function(x) caret::confusionMatrix(graph_pred[[x]], graph_ref, positive = "TRUE"))
  TP <- c()
  FP <- c()
  True_P <- c()
  True_N <- c()
  mcc <- c()
  for (i in 1:num_sims) {
    TP <- c(TP, confMat_list[[i]]$table[2,2])
    FP <- c(FP, confMat_list[[i]]$table[2,1])
    True_P <- c(True_P, sum(confMat_list[[i]]$table[,2]))
    True_N <- c(True_N, sum(confMat_list[[i]]$table[,1]))
    mcc <- c(mcc, mcc(graph_pred[[i]], graph_ref))
  }

  FDR <- FP/pmax((FP + TP), 1)
  TPR <- TP/(True_P)
  FPR <- FP/(True_N)
  return(list(TP = TP, FP = FP, FDR = FDR, MCC=mcc, TPR=TPR, FPR=FPR))
}

graph_est_and_selection <- function(graph_list, true_graph, true_omega, num_sim, method) {
  graph.L1norm <- sapply(1:num_sim, function(i) norm(graph_list[[i]] - true_omega, type = "1"))
  graph.L2norm <- sapply(1:num_sim, function(i) norm(graph_list[[i]] - true_omega, type = "2"))
  graph.Fnorm <- sapply(1:num_sim, function(i) norm(graph_list[[i]] - true_omega, type = "F"))
  graph.perf <- graph_rev_perf(graph_list, true_graph)
  return.dataframe <- data.frame(
    method=method,
    "L1.norm.Mean" = mean(graph.L1norm),
    "L1.norm.SD" =  sd(graph.L1norm),
    "L2.norm.Mean" = mean(graph.L2norm),
    "L2.norm.SD" = sd(graph.L2norm),
    "F.norm.Mean" = mean(graph.Fnorm),
    "F.norm.SD" = sd(graph.Fnorm),
    "TPR" = mean(graph.perf$TPR),
    "FPR" = mean(graph.perf$FPR),
    "FDR" = mean(graph.perf$FDR),
    "MCC.Mean" = mean(graph.perf$MCC),
    "MCC.SD" = sd(graph.perf$MCC)
  )
}

# Graph generator. Adapted from huge package on CRAN (https://cran.r-project.org/web/packages/huge/index.html). 
# Modified according to the graph described in our paper.
generator <- function (d, graph, v = 0.3, tau=1.5, u = 0.1, 
                       g = NULL, prob = NULL, verbose=FALSE)
{
  if(is.null(g)){
    g = 1
  }
  gcinfo(FALSE)
  if (verbose) 
    cat(graph, "graph structure....")
  if (graph == "cluster") {
      g = ceiling(d/20)
  }
  if (graph == "band") {
      g = 3
  }
  if (graph == "random") {
    if (is.null(prob)) 
      prob = min(1, 3/d)
    prob = sqrt(prob/2) * (prob < 0.5) + (1 - sqrt(0.5 - 0.5 * prob)) * (prob >= 0.5)
  }

  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small + 1
  g.list = c(rep(n.small, g.small), rep(n.large, g.large))
  g.ind = rep(c(1:g), g.list)
  rm(g.large, g.small, n.small, n.large, g.list)
  gc()
  theta = matrix(0, d, d)
  
  if (graph == "band") {
    for (i in 1:g) {
      diag(theta[1:(d - i), (1 + i):d]) = 1
      diag(theta[(1 + i):d, 1:(d - 1)]) = 1
    }
  }
  
  if (graph == "cluster") {
    for (i in 1:g) {
      tmp = which(g.ind == i)
      tmp2 = matrix(runif(length(tmp)^2, 0, 0.5), length(tmp), length(tmp))
      tmp2 = tmp2 + t(tmp2)
      theta[tmp, tmp][tmp2 < prob] = 1
      rm(tmp, tmp2)
      gc()
    }
  }

  if (graph == "random") {
    tmp = matrix(runif(d^2, 0, 0.5), d, d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    rm(tmp)
    gc()
  }
  
  diag(theta) = 0
  omega = theta * v
  d_mat = diag(d)
  diag(d_mat)[(d/2+1):d] = tau
  diag(omega) = abs(min(eigen(omega)$values)) + u + 0.1
  omega = d_mat %*% omega %*% d_mat
  sigma = solve(omega)

  sim = list(sigma = sigma, omega = omega, theta = theta, graph.type = graph)
  
  return(sim)
}


graph_validation <- function(validation.set, model.fit) {
  loss.vec <- rep(NA, length(model.fit$lambda))
  
  #Fit model on validation set
  if (model.fit$method=="clime") {
    validation.fit <- sugm(validation.set, method = "clime", sym="and", lambda = model.fit$lambda)
  } else if (model.fit$method == "glasso") {
    validation.fit <- huge(validation.set, method = "glasso", lambda = model.fit$lambda)
  }
  
  #Compute loss on validation set
  for (i in 1:length(model.fit$lambda)){
    loss.vec[i] <- graph_likelihood_loss(cov(validation.fit$data), validation.fit$icov[[i]])
  }
  
  #Select best
  best.lam.ind <- which.min(loss.vec)
  return(model.fit$icov[[best.lam.ind]])
}

# Likelihood loss. Adapted from flare package on CRAN (https://cran.r-project.org/web/packages/flare/index.html). 
graph_likelihood_loss <- function(sigma, omega) {
  ot <- as.numeric(unlist(determinant(omega)))
  if (ot[2]<=0) {
    warning("Precision matrix estimate is not positive definite!")
  }
  tmp <- (sum(diag(sigma%*%omega))  - ot[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}
