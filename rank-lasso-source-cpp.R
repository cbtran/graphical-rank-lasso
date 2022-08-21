library(Rcpp)
#library(Rglpk)
library(gurobi)
#library(glasso)
#library(lpSolve)
library(rqPen)
library(pracma)
library(huge)
sourceCpp("src/pairwise_diff_vec.cpp")
sourceCpp("src/constraints_mat.cpp")

rank_lasso_solver <- function(x, y, alpha0=0.1, c=1.01, B=500, solver=c("gurobi", "ecos", "glpk", "lpsolve"), sym=FALSE, m=FALSE, verbose=FALSE) {
  n <- dim(x)[1]
  p <- dim(x)[2]  
    if (sym & !m) {
        sym_mess <- "symmetric kernel: i =/= j"
    } else if (!sym & !m) {
        sym_mess <- "asymmetric kernel: i < j"
    } else if (!sym & m>0) {
        sym_mess <- "asymmetric kernel: i < j with subsampling"
    }
    mess <- paste0("Solving Rank Lasso with ", sym_mess, ", and solver: ", solver, "\n")
    if (verbose) {
        cat(mess)
    }

    num_terms = n*(n-1)
    sub_sample_indx <- 1:2
    
    #Compute lambda
    lambdlist <- c()
    for (i in 1:B) {
        set.seed(i)
        tau <- sample(1:n, size=n, replace = F)
        e <- 2*tau - (n+1)
        S_n <- -2/(n*(n-1)) * (t(x) %*% e)
        S_n <- c * norm(S_n, type = "I")
        lambdlist <- c(lambdlist, S_n)
    }
    lambdlist <- sort(lambdlist)
    lambd <- lambdlist[as.integer(B*(1-alpha0))]

    #Rank lasso inequality signs
    if ((sym == FALSE) & !m) {
      num_terms = n*(n-1)/2
    } else if (m > 0) {
      num_terms = n*m
    }
    
    if (solver=="gurobi") {
      f.dir <- c(rep("=", num_terms), rep(">=", 2*p))
    } else {
      f.dir <- c(rep("==", num_terms), rep(">=", 2*p))
    }
    
    result <- "NON-OPTIMAL"
    model <- list()
    while(result != "OPTIMAL") {
      
      if (m > 0) {
        sub_sample_indx <- sort(sample(1:(n*(n-1)/2), size=num_terms, replace=FALSE))
      }
      
    #Rank Lasso
    #Rank lasso obj function
    f.obj <- c(rep(1/num_terms, 2*num_terms), rep(lambd,p), rep(0,p))
    
    #Rank lasso constraints
    if (solver=="lpsolve") {
        constr_mat_list <- create_constr_mat_dense(x, num_terms, sym, sub_sample_indx)
        f.con <- constr_mat_list$constr_mat
        
    } else {
        constr_mat_list <- create_constr_mat(x, num_terms, sym, sub_sample_indx)
        f.con <- constr_mat_list$constr_mat
    }

    # Set right hand side coefficients
    f.rhs.lwr <- c(pairwise_vec_cpp(y, num_terms, sym, sub_sample_indx), rep(0, 2*p))
    
    
    # Set bound constraints
    lb <- c(rep(0, 2*num_terms), rep(0, p), rep(-Inf, p))
    ub <- c(rep(Inf, 2*num_terms), rep(Inf, p), rep(Inf, p))
    
    #Solve
      if (solver=="glpk") {
        solutions <- Rglpk_solve_LP(obj=f.obj, mat=f.con, dir=f.dir, rhs=f.rhs.lwr, 
                                    bounds=list(lower=list(ind=(2*num_terms+p+1):(2*num_terms+2*p), val=rep(-Inf, p)),
                                                upper=list(ind=(2*num_terms+p+1):(2*num_terms+2*p), val=rep(Inf, p))), 
                                    max=F)$solution
      } else if (solver=="lpsolve") {
        solutions <- lp("min", f.obj, f.con, f.dir, f.rhs.lwr)$solution
      } else if (solver=="gurobi") {
        
        model$A          <- f.con
        model$obj        <- f.obj
        model$modelsense <- 'min'
        model$rhs        <- f.rhs.lwr
        model$sense      <- f.dir
        model$lb         <- lb
        
        fit <- gurobi(model, params = list(OutputFlag=0))
        result <- fit$status
        
        if (result != "OPTIMAL") {
          next
        }
        solutions <- fit$x
        
      }
    }
    

    
    beta_rank <- tail(solutions, p)
    returns <- list()
    returns$mse <- mean((x%*%beta_rank - y)^2)
    returns$beta <- beta_rank
    
    #Second-stage
    #Compute lambda
    n.eta <- 5
    eta.list <- logseq(lambd, lambd*1.02, n=n.eta)
    lambda.list <- list()
    for (i in 1:n.eta){
      lambda.list[[i]] <- mcp_deriv(abs(beta_rank), eta.list[i], 3)
    }

    # #Solve
    solutions.list <- list()
    for (i in 1:n.eta){
      f.obj <- c(rep(1/num_terms, 2*num_terms), lambda.list[[i]], rep(0,p))
      model$obj        <- f.obj
      solutions_scad <- gurobi(model, params = list(OutputFlag=0))$x
      solutions.list[[i]] <- tail(solutions_scad, p)
    }

    X_tilde <- pairwise_mat_cpp(x, n*(n-1), sym=T, 1:2)
    y_tilde <- pairwise_vec_cpp(y, n*(n-1), sym=T, 1:2)
    hbic_vec <- rep(NA, n.eta)
    for (i in 1:n.eta) {
      hbic_vec[i] <- HBIC(x_tilde = X_tilde, y_tilde = y_tilde,
                          eta = eta.list[i], beta_eta = solutions.list[[i]], n, p)
    }
    beta <- solutions.list[[which.min(hbic_vec)]]
    returns$beta_scad <- beta
    returns$mse_scad <- mean((x%*%beta - y)^2)
    returns$eta <- eta.list
    returns$hbic_vec <- hbic_vec
    returns$eta.min <- which.min(hbic_vec)
    rm(model)
    return(returns)
}

HBIC <- function(x_tilde, y_tilde, eta, beta_eta, n, p) {
  num_coef <- sum(beta_eta != 0)
  log(sum(abs(x_tilde%*%beta_eta - y_tilde))) + num_coef * log(log(n)) / n * log(p)
}
