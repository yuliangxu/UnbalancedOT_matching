gc(reset = TRUE)

sim = 10 # 80 100 120
n0 = 100*sim
n1 = 10*sim

# for more treated than control
# n0 = 10*sim
# n1 = 100*sim

d = 2
n = n1 + n0
n_rep = 50

linear_bool = T
covariate_case = 3
library(ATEHonest)

n_method = 6
all_result_CATT = matrix(NA,ncol=n_method,nrow = n_rep)
all_result_ATE = matrix(NA,ncol=n_method,nrow = n_rep)
all_result_time = matrix(NA,ncol=n_method-1,nrow = n_rep)

for(iter in 1:n_rep){
  print(paste("iter = ",iter))
  # data generation ---------------------------------------------------------
  
  idx0 = sort(sample(1:n,n0))
  idx = rep(1,n)
  idx[idx0] = 0
  idx1 = which(idx == 1)
  # # case 1 (less overlap)
  # if(covariate_case == 1){
  #   X = matrix(NA,nrow = n, ncol = d)
  #   X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 2),ncol = d),
  #                    matrix(rnorm(n0/2*d, mean = -1, sd = 1),ncol = d))
  #   X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1, sd = 2),ncol = d),
  #                    matrix(rnorm(n1/2*d, mean = 1, sd = 1),ncol = d))
  #   # par(mfrow = c(1,1))
  #   plot_cov(X,idx0,idx1,main="case 1 (less overlap)")
  # }
  
  # case 3 (much less overlap)
  if(covariate_case == 3){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 2),ncol = d),
                     matrix(rnorm(n0/2*d, mean = -1, sd = 1),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 2, sd = 2),ncol = d),
                     matrix(rnorm(n1/2*d, mean = 1.5, sd = 1),ncol = d))
    # par(mfrow = c(1,1))
    plot_cov(X,idx0,idx1,main="case 3 (less overlap)")
  }
  
  
  # case 2 (more overlap)
  if(covariate_case == 2){
    X = matrix(NA,nrow = n, ncol = d)
    X[idx0,] = rbind(matrix(rnorm(n0/2*d, mean = -1, sd = 0.5),ncol = d),
                     matrix(rnorm(n0/2*d, mean = 1/2, sd = 0.5),ncol = d))
    X[idx1,] = rbind(matrix(rnorm(n1/2*d, mean = 1/2, sd = 0.5),ncol = d),
                     matrix(rnorm(n1/2*d, mean = -1/2, sd = 0.5),ncol = d))
    # par(mfrow = c(1,1))
    # plot_cov(X,idx0,idx1,main="case 2 (more overlap)")
  }
  
  # plot_cov(X,idx0,idx1,main="cov dist")
  
  if(d==8){
    Y0_all = -1 - apply(X[,1:5],1,sum) * apply(X[,6:8],1,sum) + apply(X[,1:5],1,sum)+X[,4] + rnorm(n,sd = 1)
    Y1_all = 2 + 2*apply(X[,1:4],1,sum) + apply(X[,5:8],1,sum)  + rnorm(n,sd = .5)
  }
  if(d == 4){
    Y0_all = -1 - apply(X[,1:2],1,sum) * apply(X[,3:4],1,sum)+ rnorm(n,sd = 1)
    Y1_all = 2 + 2*apply(X[,1:2],1,sum) + apply(X[,3:4],1,sum)  + rnorm(n,sd = .5)
  }
  if(d == 2){
    Y0_all = -1 - X[,1] * X[,2]  + rnorm(n,sd = 1)
    Y1_all = 2 + 2*X[,1] + X[,2] + rnorm(n,sd = .5)
  }
  
  Y0_group = cbind(Y0_all[idx0], Y1_all[idx0])
  Y1_group = cbind(Y0_all[idx1], Y1_all[idx1])
  
  
  X0 = X[idx0,]
  X1 = X[idx1,]
  Y0 = Y0_all[idx0]
  Y1 = Y1_all[idx1]
  
  Y_obs = rep(NA,n)
  Y_obs[idx0] = Y0_all[idx0]
  Y_obs[idx1] = Y1_all[idx1]
  
  D = rep(NA,n)
  D[idx0] = FALSE
  D[idx1] = TRUE
  
  
  true_CATT = mean(Y1_all[idx1]) - mean(Y0_all[idx1])
  
  true_ATE = mean(Y1_all) - mean(Y0_all)
  
  raw_ATE = mean(Y_obs[idx1]) - mean(Y_obs[idx0])
  
  print(paste("true_CATT=",true_CATT))
  
  # dt <- NSWexper[c(1:25, 421:445), ]
  # Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
  # D0 <- distMat(dt[, 2:10], Ahalf, method="manhattan", dt$treated)
  # ## Distance matrix for variance estimation
  # DM <- distMat(dt[, 2:10], Ahalf, method="manhattan")
  # c1 <- ATTOptEstimate(y=dt$re78, d=dt$treated, D0=D0, C=1, DM=DM,
  #                      opt.criterion="RMSE")
  # print(c1)
  
  # ATEHonest ---------------------------------------------------------------
  t0 = Sys.time()
  D0 <- distMat(X,  Ahalf = diag(ncol(X)), method="manhattan", D)

  DM <- distMat(X, chol(solve(cov(X))), method="euclidean")

  c1 <- ATTOptEstimate(y=Y_obs, d=D, D0=D0, C=1, DM=DM, opt.criterion="RMSE")
  print(c1)
  ATEHonest_CATT = c1$e["att"]
  t1 = Sys.time()
  ATEHonest_time = difftime(t1,t0,units = "secs")
  
  
  # amlinear ----------------------------------------------------------------
  
  # library(devtools)
  # tar_file <- "./glmnet_2.0-18.tar.gz"  # Specify the path to your tarball file
  # destination_dir <- "./"  # Specify the directory where you want to extract the tarball
  # untar(tar_file, exdir = destination_dir)
  # detach("package:glmnet", unload = TRUE)
  # remove.packages("glmnet")
  # devtools::install_version("glmnet", "2.0-18")
  # install_github("davidahirshberg/amlinear")
  # library(amlinear)
  # library(CVXR)
  # library(testthat)
  # library(glmnet)
  
  # source("./amlinear/cvtype.R")
  # tau.hat = average_partial_effect(X, Y_obs, 1*D)
  # amlinar = tau.hat[1]
  # 
  # cbind(amlinar, ATEHonest_CATT, true_ATE, true_CATT)
  
  # n = 300; p = 600; nclust = 10
  # beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))
  # clust.alpha = rep(c(0.3, 0.7), nclust/2)
  # cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
  # cluster = sample.int(nclust, n, replace = TRUE)
  # X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
  # W = rbeta(n, clust.alpha[cluster], 1 - clust.alpha[cluster])
  # Y = X %*% beta + rnorm(n, 0, 1) + 2 * W * X[,2]
  # cv.kfold = 4
  # tau.hat = average_partial_effect(X, Y, W)
  # 
  # print( paste("true tau:", round(mean(2 * cluster.center[,2]), 2) ))
  # print(paste("point estimate:", round(tau.hat[1], 2)))
  # print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))
  
  # library(amlinear)
  # library(glmnet)
  # library(CVXR)
  # tau.hat = average_partial_effect(X, Y_obs, 1*D)
  # print(paste("true tau:", round(mean(2 * cluster.center[,2]), 2)))
  # print(paste("point estimate:", round(tau.hat[1], 2)))
  # print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))
  
  

  # KOM ---------------------------------------------------------------------

  t0 = Sys.time()
  source("./KOM/KOM.R")
  dta = data.frame(cbind(Y=Y_obs, Tr = D, Z1 = X[,1],Z2 = X[,2]))
  dta <- dta[order(dta$Tr,decreasing=TRUE),]
  Type_kernel = 2
  hpara <- c(1,1,1)
  GPML <- 1
  resultGPML  <- KOM(dta,Type_kernel,hpara,GPML)
  print(resultGPML$res$statu)
  res_lm_GPML <- lm(Y ~ Tr, data = dta, weights=resultGPML$kow)
  coefGPML    <- res_lm_GPML$coefficients[2]
  KOM_ATE = coefGPML
  t1 = Sys.time()
  KOM_time = difftime(t1,t0,units = "secs")
  
  # UOT ---------------------------------------------------------------------
  
  t0 = Sys.time()
  a1 = rep(1/n1,n1)
  a0 = rep(1/n0,n0)
  center0 = apply(X0, 2, mean)
  center1 = apply(X1, 2, mean)
  
  M = matrix(NA,n0,n1)
  for(i in 1:n0){
    for(j in 1:n1){
      M[i,j] = sqrt(sum((X0[i,] - X1[j,])^2))
    }
  }
  M = M/max(M)
  M = t(M)
  library(reticulate)
  use_python("~/Library/r-miniconda-arm64/bin/python", required = T)
  ot <- import("ot")
  reg_m_kl = reg =  1e-3/2
  maxiter = 2000
  ubG = ot$sinkhorn_unbalanced(a1, a0, M, reg,reg_m_kl, div='kl',numItermax=as.integer(maxiter))
  # summary(apply(ubG,1,sum))
  Y0_ot = rep(NA,n)
  Y0_ot[idx0] = Y0
  Y0_ot[idx1] = (ubG %*% Y0) / apply(ubG,1,sum)
  Y1_ot = rep(NA,n)
  Y1_ot[idx0] = (t(ubG) %*% Y1 ) / apply(ubG,2,sum)
  Y1_ot[idx1] = Y1
  
  UOT_ATE = mean( Y1_ot,na.rm = T) - mean(Y0_ot,na.rm = T)
  UOT_CATT = mean( Y1_ot[idx1],na.rm = T) - mean(Y0_ot[idx1],na.rm = T)
  t1 = Sys.time()
  UOT_time =  difftime(t1,t0,units = "secs")
  
  # KNN ---------------------------------------------------------------------
  
  t0 = Sys.time()
  source("./R/thresh_ot_func.R")
  # KNN (k=3)
  knn3 = get_knn(M,3)
  Y0_knn3 = rep(NA,n); Y1_knn3 = rep(NA,n)
  Y0_knn3[idx0] = Y0; Y1_knn3[idx1] = Y1;
  Y0_knn3[idx1] = apply(knn3$X1_knn,1,function(x){mean(Y0[x])})
  Y1_knn3[idx0] = apply(knn3$X0_knn,1,function(x){mean(Y1[x])})
  
  knn3_CATT = mean( Y1_knn3[idx1],na.rm = T) - mean(Y0_knn3[idx1],na.rm = T)
  knn3_ATE = mean( Y1_knn3,na.rm = T) - mean(Y0_knn3,na.rm = T)
  t1 = Sys.time()
  knn3_time =  difftime(t1,t0,units = "secs")
  
  # KNN (k=1)
  t0 = Sys.time()
  knn1 = get_knn(M,1)
  Y0_knn1 = rep(NA,n); Y1_knn1 = rep(NA,n)
  Y0_knn1[idx0] = Y0; Y1_knn1[idx1] = Y1;
  Y0_knn1[idx1] = apply(knn1$X1_knn,2,function(x){mean(Y0[x])})
  Y1_knn1[idx0] = apply(knn1$X0_knn,2,function(x){mean(Y1[x])})
  knn1_CATT = mean( Y1_knn1[idx1],na.rm = T) - mean(Y0_knn1[idx1],na.rm = T)
  knn1_ATE = mean( Y1_knn1,na.rm = T) - mean(Y0_knn1,na.rm = T)
  t1 = Sys.time()
  knn1_time =  difftime(t1,t0,units = "secs")
  
  # summary -----------------------------------------------------------------
  
  result = cbind(true_CATT, UOT_CATT,ATEHonest_CATT, knn1_CATT, knn3_CATT,raw_ATE)
  n_method = dim(result)[2]
  out_nrep = result[1,2:n_method] - result[1,1]
  
  result_ATE = cbind(true_ATE, UOT_ATE,KOM_ATE, knn1_ATE, knn3_ATE,raw_ATE)
  
  all_result_CATT[iter,] = result
  all_result_ATE[iter,] = result_ATE
  all_result_time[iter,] = cbind(KOM_time, UOT_time, ATEHonest_time,knn1_time, knn3_time)

}

out = list(all_result_CATT = all_result_CATT,
           all_result_ATE = all_result_ATE)
out_time = all_result_time
saveRDS(list(out=out, out_time=out_time),"./result/compare2_case3_d2_nsim10.rds")

sum_stats = function(all_result, n_method){
  list(bias = apply(abs(all_result[,2:n_method]- all_result[,1]),2,mean) ,
       Var = apply(all_result[,2:n_method],2,var) ,
       Mse = apply(all_result[,2:n_method],2,function(x){mean((x - all_result[,1])^2)})
  )
}

all_sum_stats = lapply(out,sum_stats,n_method)
out_CATT = with(all_sum_stats$all_result_CATT,rbind(bias,Var,Mse))
colnames(out_CATT) = colnames(cbind(UOT_CATT,ATEHonest_CATT, knn1_CATT, knn3_CATT,raw_ATE))
View(round(out_CATT, digits = 3))


out_ATE = with(all_sum_stats$all_result_ATE,rbind(bias,Var,Mse))

out_time_mean = apply(all_result_time,2,mean)

colnames(out_CATT) = colnames(cbind(UOT_CATT,ATEHonest_CATT, knn1_CATT, knn3_CATT,raw_ATE))
colnames(out_ATE) = colnames( cbind( UOT_ATE,KOM_ATE, knn1_ATE, knn3_ATE,raw_ATE))
names(out_time_mean) = colnames( cbind(KOM_time, UOT_time, ATEHonest_time,knn1_time, knn3_time))

View(round(out_CATT, digits = 3))
View(round(out_ATE, digits = 3))
View(round(t(as.matrix(out_time_mean)), digits = 3))

# # Load ggplot2
# library(ggplot2)
# 
# # Create data
# data <- data.frame(
#   name=colnames(result) ,
#   value=result[1,]
# )
# ggplot(data, aes(x=name, y=value)) +
#   geom_bar(stat = "identity", width=0.2) +
#   geom_hline(yintercept = true_CATT, color = "red", linetype = "dashed") +
#   geom_hline(yintercept = ATEHonest_CATT, color = "blue", linetype = "dashed")
# print(result)
