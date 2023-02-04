path_ATO = "ATO.R"
get_mse_all = function(case,sim, path=path_ATO){
  e <- new.env()
  e$covariate_case =  case
  e$sim = sim
  source(path, local = e)
  
  cbind(e$bias , e$VAR, e$MSE)
}


ATE_all = NULL
all_sim = c(2,3,5,10,30,50)
n_sim = length(all_sim)
all_result = array(dim = c(11,3,length(all_sim)))
i=0
for(sim in all_sim){
  print(paste("sim = ",sim))
  ATE_mse = get_mse_all(case=4,sim)
  i = i+1
  all_result[,,i] = ATE_mse
}
beepr::beep()
saveRDS(all_result,paste("result/all_result_d2_linear_","case1_",".rds",sep="")) # create your own result folder


for(d in c(2,3,4)){
  for(case in 1:4){
    for(linear in c("linear","nonlinear")){
      print(paste("plot:n_plot_case",case,",d=",d,",",linear))
      
      all_result = readRDS(paste("result/all_result_d",d,"_",linear,"_","case",case,"_",".rds",sep=""))
      stats_name = c("bias","var","mse")
      
      # label = colnames(output_all)
      label = output = c("ATE_TV",
                         "ATE_TVow",
                         "ATE_OT", 
                         "ATE_OT2", 
                         "ATE_OT3", 
                         "ATE_knn3", 
                         "ATE_knn1",
                         "ATE_ow",
                         "ATE_dr", 
                         "ATE_dr_oracle",
                         "ATE_ps")
      
      stats = 1
      dim(all_result)[1]
      which_method = c(1,3,4,5,6,7,8)
      label[which_method]
      
      
      n_methods = length(which_method)
      df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                       method=rep(label[which_method], each=n_sim))
      library(ggplot2)
      p1 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method)) + guides(colour = "none")+
        ggtitle(stats_name[stats]) + xlab("n1/10, n0/100") + ylab(stats_name[stats])
      
      stats = 2
      df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                       method=rep(label[which_method], each=n_sim))
      p2 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method)) + guides(colour = "none")+ 
        ggtitle(stats_name[stats])+ xlab("n1/10, n0/100") + ylab(stats_name[stats])
      
      stats = 3
      df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                       method=rep(label[which_method], each=n_sim) )
      p3 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method))+
        theme(legend.position = 'right') +
        ggtitle(stats_name[stats])+ xlab("n1/10, n0/100") + ylab(stats_name[stats])
      
      
      g_all = egg::ggarrange(p1, p2, p3,nrow = 1, byrow = F,top=paste("case",case,",d=",d,",",linear,sep="") )
      ggsave(paste("./plot/n_plot_case",case,",d=",d,",",linear,".pdf",sep=""),plot = g_all,device="pdf",width = 9, height = 3, unit="in")
      
    }
  }
}

all_result = readRDS(paste("result/all_result_d4_linear_","case1_",".rds",sep=""))
# par(mfrow = c(1,3))
stats_name = c("bias","var","mse")

# label = colnames(output_all)
label = output = c("ATE_TV",
                   "ATE_TVow",
                   "ATE_OT", 
                   "ATE_OT2", 
                   "ATE_OT3", 
                   "ATE_knn3", 
                   "ATE_knn1",
                   "ATE_ow",
                   "ATE_dr", 
                   "ATE_dr_oracle",
                   "ATE_ps")

stats = 1
dim(all_result)[1]
which_method = c(1,3,4,5,6,7,8)
label[which_method]


n_methods = length(which_method)
df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                 method=rep(label[which_method], each=n_sim))
library(ggplot2)
p1 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method)) + guides(colour = "none")+
  ggtitle(stats_name[stats]) + xlab("n1/10, n0/100") + ylab(stats_name[stats])

stats = 2
df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                 method=rep(label[which_method], each=n_sim))
p2 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method)) + guides(colour = "none")+ 
  ggtitle(stats_name[stats])+ xlab("n1/10, n0/100") + ylab(stats_name[stats])

stats = 3
df <- data.frame(x=rep(all_sim, n_methods), val=c(t(all_result[which_method,stats,])), 
                 method=rep(label[which_method], each=n_sim) )
p3 = ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=method))+
  theme(legend.position = 'right') +
  ggtitle(stats_name[stats])+ xlab("n1/10, n0/100") + ylab(stats_name[stats])


egg::ggarrange(p1, p2, p3,nrow = 1, byrow = F,top="case1,d=4,linear")

