library(parallel)
library(foreach)
library(doParallel)

path_simus=""

source(paste0(path_simus,"/PK_PD_param.R"))
##########################################################################

t=seq(0,650,by=0.1)

# n=100000
# Nsim=1000
# no_cores=40


n=20
Nsim=10
no_cores=2


N_pack=n/Nsim
seed=25011995
vect_packages=c("deSolve","MASS")

cl <- makeCluster(no_cores)
registerDoParallel(cl)

Rmax <- foreach(ip=1:N_pack,.combine=rbind,.inorder=FALSE,.packages = vect_packages) %dopar% {  
  res=matrix(NA,nrow=Nsim*length(t_doses),ncol=nrow(DL)+2)
  colnames(res)=c("id","admin",paste0("DL",1:nrow(DL)))
  
  for(i in 1:Nsim){
    
    cat(ip,"/",N_pack,"-",i,"/",Nsim,"\n")
    
    set.seed(seed+(ip-1)*Nsim+i)
    res[(i-1)*length(t_doses)+1:length(t_doses),]=cbind(rep((ip-1)*Nsim+i,length(t_doses)),1:length(t_doses),
                                                        PK_PD_simu_Rmax_local(t=t,
                                                                              doses=DL,
                                                                              t_doses=t_doses,
                                                                              t_inf=t_inf,
                                                                              PK_model=PK_model,
                                                                              PK_param_pop=PK_param_pop,
                                                                              PD_model=PD_model,
                                                                              PD_param_pop=PD_param_pop,
                                                                              PD_initial_values=PD_initial_values))
    
  }
  
  
  res
  
} 


#save(Rmax, file = paste0(path_simus,"/scenario/Rmax/Rmax_local_test.Rdata"))


