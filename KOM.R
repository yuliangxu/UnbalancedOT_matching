




# rm(list = ls())

# set.seed(17)
# install.packages("quadprog")
# install.packages("ipw")
# install.packages("sandwich")
library(gurobi)
#library(ipw)
library(sandwich)
library(sbw)
library(quadprog)
library(CBPS)


################################################################################
#FUNCTIONS SBW

sbw_result <- function(dta,degree,BAL_TOLS){
  
  if(degree==1){
    data_frame = as.data.frame(cbind(dta$Y, dta$Tr, dta$Z1, dta$Z2))
    names(data_frame) = c("y", "t", "Z1", "Z2")
    
    # Treatment indicator
    t_ind = "t"
    
    # Moment covariates
    bal_covs =c("Z1", "Z2")
    
    # Moment tolerances
    bal_tols = BAL_TOLS
    
    # Whether the moment tolerances are expressed in standard deviations
    bal_tols_sd = TRUE
    
    # Here, the 0's in t_ind are weighted to "represent" the 1's and estimate the average treatment on the treated
    target = "all"
    
    # Here, the "ell-2" norm is used to minimize the variance of the weights
    l_norm = "l_2"
    
    # Minimum value of the weights
    w_min = 0
    
    # Here, the weights are constrained to add up to one
    normalize = 1
    
    # Solver
    solver = "gurobi"
    
    # Output display
    display = 0
    
    # Find optimal weights
    #recovering the covariate structure of both the treated and control units
    #in the weighted CONTROL units
    out0 = sbw(data_frame, t_ind, bal_covs, bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display)
    
    # data_frame_te = as.data.frame(cbind(dta$Y, dta$Tr, dta$Z1, dta$Z1))
    # names(data_frame_te) = c("y", "t", "Z1", "Z2")
    # sbw(data_frame_te, t_ind, c("Z1","Z2"), bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display)$data_frame_weights$weights
    
    # data_sbw_final <- data.frame(dta$Tr,out0$data_frame_weights$weights)
    # colnames(data_sbw_final) <- c("Tr","sbw")
    data_sbw_final <- data.frame(rep(NA,1))
    colnames(data_sbw_final) <- c("sbw")
    if (out0$status=="OPTIMAL"){
      data_sbw_final <- data.frame(dta$Tr,out0$data_frame_weights$weights)
      colnames(data_sbw_final) <- c("Tr","sbw")
    }
    
    #recovering the covariate structure of both the treated and control units
    #in the weighted TREATED units
    t1 <- 1-dta$Tr
    
    data_frame2 = as.data.frame(cbind(dta$Y, t1, dta$Z1, dta$Z2))
    names(data_frame2) = c("y", "t1", "Z1", "Z2")
    data_frame2 <- data_frame2[order(data_frame2$t1),]
    
    out1 = sbw(data_frame2, "t1", bal_covs, bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display)
    #out1$data_frame_weights$weights
    
    #out1$data_frame_weights$weights[which(data_sbw_final$sbw==0)]
    #head(data_sbw_final)
    
    # data_temp <- as.data.frame(cbind(data_frame2$t1,out1$data_frame_weights$weights))
    # data_temp <- data_temp[order(data_temp$V1,decreasing = TRUE),]
    #data_sbw_final$sbw[which(data_sbw_final$Tr==1)] <- out1$data_frame_weights$weights[which(out1$data_frame_weights$weights!=0)]
    
    if (out1$status=="OPTIMAL" ){
      data_sbw_final$sbw[which(data_sbw_final$Tr==1)] <- out1$data_frame_weights$weights[which(out1$data_frame_weights$weights!=0)]
    }
    
    sbw_counter <- 0
    sbw_counter2 <- 0
    if (out0$status!="OPTIMAL" | out1$status!="OPTIMAL"){ sbw_counter <- 1}
    
  }#end if degree 1
  
  if(degree==2){
    data_frame = as.data.frame(cbind(dta$Y, dta$Tr, dta$Z1, dta$Z2,
                                     dta$Z1*dta$Z2,
                                     dta$Z1^2,dta$Z2^2,
                                     dta$Z1^3,dta$Z2^3,
                                     dta$Z1*dta$Z2^2,dta$Z2*dta$Z1^2))
    names(data_frame) = c("y", "t", "Z1", "Z2",
                          "Z1Z2",
                          "Z1_2","Z2_2",
                          "Z1_3","Z2_3",
                          "Z1Z2_2","Z2Z1_2")
    
    
    # Treatment indicator
    t_ind = "t"
    
    # Moment covariates
    bal_covs =c("Z1", "Z2",
                "Z1Z2",
                "Z1_2","Z2_2",
                "Z1_3","Z2_3",
                "Z1Z2_2","Z2Z1_2")
    
    # Moment tolerances
    bal_tols = BAL_TOLS
    
    # Whether the moment tolerances are expressed in standard deviations
    bal_tols_sd = TRUE
    
    # Here, the 0's in t_ind are weighted to "represent" the 1's and estimate the average treatment on the treated
    target = "all"
    
    # Here, the "ell-2" norm is used to minimize the variance of the weights
    l_norm = "l_2"
    
    # Minimum value of the weights
    w_min = 0
    
    # Here, the weights are constrained to add up to one
    normalize = 1
    
    # Solver
    solver = "gurobi"
    
    # Output display
    display = 0
    
    # Find optimal weights
    #recovering the covariate structure of both the treated and control units
    #in the weighted CONTROL units
    out0 = try(sbw(data_frame, t_ind, bal_covs, bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display))
    print(out0$status)
    
    data_sbw_final <- data.frame(rep(NA,1))
    colnames(data_sbw_final) <- c("sbw")
    
    if (out0$status=="OPTIMAL"){
      data_sbw_final <- data.frame(dta$Tr,out0$data_frame_weights$weights)
      colnames(data_sbw_final) <- c("Tr","sbw")
    }
    
    # data_frame_te = as.data.frame(cbind(dta$Y, dta$Tr, dta$Z1, dta$Z1))
    # names(data_frame_te) = c("y", "t", "Z1", "Z2")
    # sbw(data_frame_te, t_ind, c("Z1","Z2"), bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display)$data_frame_weights$weights
    
    #recovering the covariate structure of both the treated and control units
    #in the weighted TREATED units
    t1 <- 1-dta$Tr
    
    data_frame2 = as.data.frame(cbind(dta$Y, t1, dta$Z1, dta$Z2,
                                      dta$Z1*dta$Z2,
                                      dta$Z1^2,dta$Z2^2,
                                      dta$Z1^3,dta$Z2^3,
                                      dta$Z1*dta$Z2^2,dta$Z2*dta$Z1^2))
    names(data_frame2) = c("y", "t1", "Z1", "Z2",
                           "Z1Z2",
                           "Z1_2","Z2_2",
                           "Z1_3","Z2_3",
                           "Z1Z2_2","Z2Z1_2")
    
    data_frame2 <- data_frame2[order(data_frame2$t1),]
    
    out1 = try(sbw(data_frame2, "t1", bal_covs, bal_tols, bal_tols_sd, target, l_norm, w_min, normalize, solver, display))
    #out1$data_frame_weights$weights
    print(out1$status)
    
    #out1$data_frame_weights$weights[which(data_sbw_final$sbw==0)]
    #head(data_sbw_final)
    
    
    
    # data_temp <- as.data.frame(cbind(data_frame2$t1,out1$data_frame_weights$weights))
    # data_temp <- data_temp[order(data_temp$V1,decreasing = TRUE),]
    if (out1$status=="OPTIMAL" ){
      data_sbw_final$sbw[which(data_sbw_final$Tr==1)] <- out1$data_frame_weights$weights[which(out1$data_frame_weights$weights!=0)]
    }
    
    sbw_counter <- 0
    sbw_counter2 <- 0
    if (out0$status!="OPTIMAL" | out1$status!="OPTIMAL"){ sbw_counter2 <- 1}
    
    
    
  }#end if degree 2
  
  
  
  sbw <- data_sbw_final$sbw
  
  return(list(sbw=sbw,sbw_counter=sbw_counter,sbw_counter2=sbw_counter2))
}




################################################################################
#FUNCTIONS GPML


#Gaussian kernel for minimization
Kgauss <- function(vect,x){
  #print("Gaussian kernel")
  x <- scale(x)
  Scale <- vect[2]^2
  B <-  exp( -(Scale) * as.matrix(dist(x, upper=T, diag=T))^2 )
  return(B)
}


Klinear <-  function(vect,x){
  xs <- scale(x)
  offset <- vect[2]
  B <- (xs%*%t(xs) + offset)
  return(B)
}

Kpoly <-  function(vect,x){
  xs <- scale(x)
  offset <- vect[2]
  B <- (xs%*%t(xs) + offset)^3
  return(B)
}


#Gaussian kernel for analysis
KGram <- function(vect,X,type_kernel){
  
  if(type_kernel==0){
    K <- Kmatern(vect,X)
  }
  if(type_kernel==1){
    K <- Kgauss(vect,X)
  }
  if(type_kernel==2){
    K <- Klinear(vect,X)
  }
  if(type_kernel==3){
    K <- Kpoly(vect,X)
  }
  return(K)
  
}

#minus log likelihood
minus_log_likelihood <- function(vect,X,Y,nte,type_kernel){
  sn2 <- vect[3]^2
  gamma2 <- vect[1]^2
  #scale <- 1
  if(type_kernel==0){ K_plus <- Kmatern(vect,X) + sn2*diag(nte) }
  if(type_kernel==1){ K_plus <- gamma2*Kgauss(vect,X) + sn2*diag(nte) }
  if(type_kernel==2){ K_plus <- gamma2*Klinear(vect,X) + sn2*diag(nte) }
  if(type_kernel==3){ K_plus <- gamma2*Kpoly(vect,X) + sn2*diag(nte) }
  
  K_plus_inv <- try(solve(K_plus, tol = 1e-21))
  
  # print(paste("yx:","K_plus_inv = ",K_plus_inv))
  # print(paste("yx:","class(K_plus_inv) = ",class(K_plus_inv)))
  
  # if(class(K_plus_inv) != "try-error"){
  #   z <- determinant(K_plus, logarithm=TRUE)
  #   K_plus_log_det <- as.numeric((z$sign*z$modulus)) # log-determinant of K_plus
  #   out <- 0.5 * ( t(Y) %*% K_plus_inv %*% Y ) + 0.5 * K_plus_log_det + (nte/2)*log(2*pi)
  # }
  
  # yx - modified
  if(!("try-error" %in% class(K_plus_inv))){
    z <- determinant(K_plus, logarithm=TRUE)
    K_plus_log_det <- as.numeric((z$sign*z$modulus)) # log-determinant of K_plus
    out <- 0.5 * ( t(Y) %*% K_plus_inv %*% Y ) + 0.5 * K_plus_log_det + (nte/2)*log(2*pi)
  }
  
  return(out)
}



#######################################################
# KOW
#######################################################



KOM <- function(dta,Type_kernel,hpara,GPML){
  
  print("yx --  begin KOM")
  X <- as.matrix(data.frame(dta$Z1,dta$Z2))
  
  y <- dta$Y
  
  TT <- dta$Tr
  
  t1 <-as.integer(TT)
  t0 <-as.integer((1-TT))
  
  if(GPML == 1){
    
    ht <- dta$Tr
    Xtemp <- data.frame(X,y,Tr)
    Xsample <- Xtemp
    #Xsample <- Xtemp[sample(nrow(Xtemp), n), ]
    
    X1 <- Xsample[which(Xsample$Tr==1),]
    X0 <- Xsample[which(Xsample$Tr==0),]
    
    y1 <- y[which(Xsample$Tr==1)]
    y0 <- y[which(Xsample$Tr==0)]
    
    n1 <- length(which(Xsample$Tr==1))
    n0 <- length(which(Xsample$Tr==0))
    
    # GPML optimization --- "L-BFGS-B" algorithm
    print("GPML Optimization")
    
    # GPML optimization --- "L-BFGS-B" algorithm
    
    tol <- 1e-08
    
    system_time_GPML0 <- system.time(res.optim2_0 <- try(optim(par=c(1, 1, 1),
                                                               fn=minus_log_likelihood,
                                                               method=c("L-BFGS-B"),
                                                               lower = rep(tol,3),
                                                               hessian=TRUE,
                                                               control=list(trace=0, maxit=1000),
                                                               X=X0[,1:(dim(X0)[2]-2)],
                                                               Y=y0,
                                                               nte=n0,
                                                               type_kernel=Type_kernel)) )
    
    if(class(res.optim2_0) == "try-error"){
      res.optim2_0 <- list(par=c(1,1,1))
    }
    #res.optim2_0$par
    
    # #
    
    # Xsample <- Xtemp
    #
    # X1 <- Xsample[which(Xsample$Tr==1),]
    #
    # y1 <- y[which(Xsample$Tr==1)]
    #
    # n1 <- length(which(Xsample$Tr==1))
    
    system_time_GPML1 <- system.time( res.optim2_1 <- try(optim(par=c(1, 1, 1),
                                                                fn=minus_log_likelihood,
                                                                method=c("L-BFGS-B"),
                                                                lower = rep(tol,3),
                                                                hessian=TRUE,
                                                                control=list(trace=0, maxit=1000),
                                                                X=X1[,1:(dim(X1)[2]-2)],
                                                                Y=y1,
                                                                nte=n1,
                                                                type_kernel=Type_kernel)) )
    
    #res.optim2_1$par
    if(class(res.optim2_1) == "try-error"){
      res.optim2_1 <- list(par=c(1,1,1))
    }
    
    if ((class(res.optim2_0) != "try-error") & (class(res.optim2_1) != "try-error")){
      res.optim2_0 <- list(par=c(res.optim2_0$par[1],res.optim2_0$par[2],res.optim2_0$par[3]))
      res.optim2_1 <- list(par=c(res.optim2_1$par[1],res.optim2_1$par[2],res.optim2_1$par[3]))
    }
  }#end if GPML
  
  if(GPML == 0){
    #lambda <- 10
    #print(paste("Hyperparameters are fixed: ", hyp1,hyp2,hyp3))
    
    res.optim2_0 <- list(par=c(hyp1[1],hyp2[1],hyp3[1]))
    res.optim2_1 <- list(par=c(hyp1[2],hyp2[2],hyp3[2]))
    #print(res.optim2_1)
  }
  
  #compute K
  GPML_run <- 0
  system_time_matrices <- system.time({
    K1 <- KGram(res.optim2_1$par,X,Type_kernel)
    K0 <- KGram(res.optim2_0$par,X,Type_kernel)
    
    gc()
    
    print("Building matrices")
    onesn <- rep(1/n,n)
    
    #Quadratic part
    # I1 K1 I1 + I0 K0 I0
    I1KI1 <- outer(t1, t1)*K1
    I0KI0 <- outer(t0, t0)*K0
    
    #Linear part
    # en^T K1 I1 + en^T K0 I0
    #K1 I1, K0 I0
    KI1 <- diag(t1)%*%K1
    KI0 <- diag(t0)%*%K0
    
    # en^T K1 I1, en^T K0 I0
    onesnKI1 <- t(onesn)%*%KI1
    onesnKI0 <- t(onesn)%*%KI0
    
    rm(list=c("K1","K0","KI1","KI0"))
    gc()
    
    ##########################################################
    # SATE
    #
    
    tol <- 1e-08
    
    
    lambda1 <- res.optim2_1$par[3]^2
    lambda0 <- res.optim2_0$par[3]^2
    
    # lambda1 <- lamda0 <- 0.001
    
    lambda <- lambda1*diag(t1) + lambda0*diag(t0)
    
    #Update Q
    Kcircl <- 0.5*( I1KI1 + I0KI0 + lambda )
    
    rm(list=c("I1KI1","I0KI0"))
    gc()
    
    #Update c
    #Gurobi:
    en_Kdiam <- -(onesnKI1 + onesnKI0)
    
  })#end system.time
  
  rm(list = c("onesnKI1","onesnKI0"))
  
  print("Solving QP")
  
  model <- list()
  model$A          <- matrix(c(((t1)/n),((t0)/n)), nrow=2, byrow=T)
  model$rhs        <- c(1,1)
  model$modelsense <- "min"
  model$Q          <- Kcircl
  model$obj        <- en_Kdiam
  model$sense      <- c("=")
  model$lb <- rep(tol,n)
  model$vtypes <- "C"
  
  
  params <- list(Presolve=2,OutputFlag=0,QCPDual=0)
  
  system_time_gurobi <- system.time(resg <- try(gurobi(model,params)))
  res <- resg
  
  if (class(res) != "try-error"){
    result <- list(res=res,kow=res$x, res.optim2_1_par = res.optim2_1$par,res.optim2_0_par = res.optim2_0$par, lambda0 = lambda0, lambda1 = lambda1,
                   system_time_GPML0=system_time_GPML0, system_time_GPML1=system_time_GPML1, system_time_gurobi=system_time_gurobi, system_time_matrices=system_time_matrices)
  }
  return(result)
}




