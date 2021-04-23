#######################################################################################################
#######################################################################################################
############################################Simulation Run#############################################
#######################################################################################################
#######################################################################################################
library(dfcrm)
##################Scenarios definition
logit=function(x) {
  log(x/(1-x))
}

invlogit=function (x) {
  1/(1 + exp(-x))
}


pTox_calculation_CRS=function(threshold,w_alpha,Rmax_global){
  res=NULL
  
  if(w_alpha==0){
    res=colMeans(Rmax_global>=threshold,na.rm=T)
  } else{
    res=colMeans(1-pnorm((log(threshold)-log(as.matrix(Rmax_global)))/w_alpha),na.rm=T)
  }
  return(as.numeric(res))
}

pTox_conditional=function(doses=c(1,5,10,rep(20,4)),
                          param=list(
                            beta0,
                            beta1,
                            beta2
                          )){
  
  with(param,{
    
    cum_doses_1=c(0,cumsum(doses)[-length(doses)])
    
    cond_prob=invlogit(
      beta0+beta1*doses+beta2*cum_doses_1
    )
    
    return(cond_prob)
    
  })
}

pTox_conditional_log=function(doses=c(1,5,10,rep(20,4)),
                              param=list(
                                beta0,
                                beta1,
                                beta2
                              )){
  
  with(param,{
    
    cum_doses_1=c(0,cumsum(doses)[-length(doses)])
    
    cond_prob=invlogit(
      beta0+beta1*log(doses)+beta2*log(cum_doses_1+1)
    )
    
    return(cond_prob)
    
  })
}

pTox_oDLT_sc=function(DL,pTox_cond_fct,param_cond_fct,gamma){
  pTox_oDLT_CRS_cond=t(apply(DL,1,pTox_cond_fct,
                             param_cond_fct))
  
  pTox_oDLT_noCRS_cond=pTox_oDLT_CRS_cond/gamma
  
  return(list(
    pTox_oDLT_CRS_cond,pTox_oDLT_noCRS_cond
  ))
  
}

pTox_tot=function(threshold=700, 
                  w_alpha=0.1, 
                  Rmax_global, 
                  DL=DL_kg, 
                  pTox_cond_fct,
                  param_cond_fct,
                  gamma){
  
  pTox_CRS=pTox_calculation_CRS(threshold,w_alpha,Rmax_global)
  
  pTox_oDLT_cond=pTox_oDLT_sc(DL,pTox_cond_fct,param_cond_fct,gamma)
  
  pTox_oDLT_CRS_cond=pTox_oDLT_cond[[1]]
  
  pTox_oDLT_noCRS_cond=pTox_oDLT_cond[[2]]
  
  pTox_oDLT_CRS=apply(pTox_oDLT_CRS_cond,1,function(x) 1-prod(1-x))
  
  pTox_oDLT_noCRS=apply(pTox_oDLT_noCRS_cond,1,function(x) 1-prod(1-x))
  
  return(
    list(
      pTox_CRS=pTox_CRS,
      pTox_oDLT=as.numeric(
        pTox_oDLT_CRS*pTox_CRS+pTox_oDLT_noCRS*(1-pTox_CRS)
      ),
      pTox_DLT=as.numeric(
        1-(1-pTox_oDLT_noCRS)*(1-pTox_CRS)
      ),
      pTox_oDLT_CRS=pTox_oDLT_CRS,
      pTox_oDLT_noCRS=pTox_oDLT_noCRS,
      pTox_oDLT_CRS_cond=pTox_oDLT_CRS_cond,
      pTox_oDLT_noCRS_cond=pTox_oDLT_noCRS_cond
      
    )
  )
}

##Function to optimise the value of the threshold
optim_scenario=function(threshold=700,
                        Rmax_global,
                        DL=DL_kg,
                        w_alpha=0.1,
                        pTox_cond_fct,
                        param_cond_fct,
                        gamma=2,
                        target_DLT=0.3,index_DLT=4){
  
  total_pTox=pTox_tot(
    threshold, 
    w_alpha, 
    Rmax_global, 
    DL, 
    pTox_cond_fct,
    param_cond_fct,
    gamma
  )
  
  return((total_pTox$pTox_DLT[index_DLT]-target_DLT)^2)
  
}


######################################################################################################
#########################################Simulation parameters########################################
######################################################################################################
path_simus="C:/Users/Emma/Dropbox/These/Redaction/joint_model_paper/Review1/code"

load(paste0(path_simus,"/scenario/Rmax/Rmax_local.Rdata"))
Rmax_global=aggregate(Rmax[,grep("DL",colnames(Rmax))],by=list(Rmax[,"id"]), max,na.rm=T)[,-1]
rm(Rmax)

DL_kg=matrix(data=c(
  1,5,10,rep(20,4),
  1,5,10,rep(25,4),
  1,5,10,rep(30,4),
  1,5,10,rep(35,4),
  1,5,10,rep(40,4),
  1,5,10,rep(45,4),
  1,5,10,rep(50,4),
  5,10,25,rep(75,4),
  10,25,50,rep(100,4)
),ncol=7,byrow=T)


load(paste0(path_simus,"/scenario/Rmax/Rmax_local_prediction.Rdata"))

Rmax_global_predict=aggregate(Rmax[,grep("DL",colnames(Rmax))],by=list(Rmax[,"id"]), max,na.rm=T)[,-1]
rm(Rmax)
DL_kg_predict=matrix(data=c(
  5,10,rep(30,5),
  1,5,10,30,rep(60,3)
),ncol=7,byrow=T)

#################Création des scénarios
scenario=list()
scenario_predict=list()

###############Sc1 classique from set 1

id_sc=1
w_alpha=0.1
beta0=-8.5
beta1=1.3
beta2=0.25
gamma=2
target_DLT=0.3
index_DLT=4
id_DL=c(1,2,3,6,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(300,1500),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)

pTox_sc_predict=pTox_tot(threshold=scenario[[id_sc]]$threshold, w_alpha=scenario[[id_sc]]$w_alpha,
                         Rmax_global=Rmax_global_predict,
                         DL=DL_kg_predict,
                         pTox_cond_fct=pTox_conditional_log,
                         param_cond_fct=list(beta0=scenario[[id_sc]]$beta0,beta1=scenario[[id_sc]]$beta1,beta2=scenario[[id_sc]]$beta2),
                         gamma=scenario[[id_sc]]$gamma)

scenario_predict[[id_sc]]=list(
  id_DL=1:nrow(DL_kg_predict),
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc_predict$pTox_CRS,
  pTox_oDLT=pTox_sc_predict$pTox_oDLT,
  pTox_DLT=pTox_sc_predict$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc_predict$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc_predict$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc_predict$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc_predict$pTox_oDLT_noCRS_cond
)

###############Set 1 gamma=8

id_sc=2
w_alpha=0.1
beta0=-8.3
beta1=1.6
beta2=0.2
gamma=8
target_DLT=0.3
index_DLT=4
id_DL=c(1,2,3,6,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(300,1500),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)

pTox_sc_predict=pTox_tot(threshold=scenario[[id_sc]]$threshold, w_alpha=scenario[[id_sc]]$w_alpha,
                         Rmax_global=Rmax_global_predict,
                         DL=DL_kg_predict,
                         pTox_cond_fct=pTox_conditional_log,
                         param_cond_fct=list(beta0=scenario[[id_sc]]$beta0,beta1=scenario[[id_sc]]$beta1,beta2=scenario[[id_sc]]$beta2),
                         gamma=scenario[[id_sc]]$gamma)

scenario_predict[[id_sc]]=list(
  id_DL=1:nrow(DL_kg_predict),
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc_predict$pTox_CRS,
  pTox_oDLT=pTox_sc_predict$pTox_oDLT,
  pTox_DLT=pTox_sc_predict$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc_predict$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc_predict$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc_predict$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc_predict$pTox_oDLT_noCRS_cond
)


###############Sc1 classique from set 1, more oDLT

id_sc=3
w_alpha=0.1
beta0=-7.5
beta1=1.4
beta2=0
gamma=2
target_DLT=0.3
index_DLT=4
id_DL=c(1,2,3,6,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(300,2000),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)

pTox_sc_predict=pTox_tot(threshold=scenario[[id_sc]]$threshold, w_alpha=scenario[[id_sc]]$w_alpha,
                         Rmax_global=Rmax_global_predict,
                         DL=DL_kg_predict,
                         pTox_cond_fct=pTox_conditional_log,
                         param_cond_fct=list(beta0=scenario[[id_sc]]$beta0,beta1=scenario[[id_sc]]$beta1,beta2=scenario[[id_sc]]$beta2),
                         gamma=scenario[[id_sc]]$gamma)

scenario_predict[[id_sc]]=list(
  id_DL=1:nrow(DL_kg_predict),
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc_predict$pTox_CRS,
  pTox_oDLT=pTox_sc_predict$pTox_oDLT,
  pTox_DLT=pTox_sc_predict$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc_predict$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc_predict$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc_predict$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc_predict$pTox_oDLT_noCRS_cond
)

###############Sc1 classique from set 1, more CRS

id_sc=4
w_alpha=0.1
beta0=-9.8
beta1=1.4
beta2=0.35
gamma=2
target_DLT=0.3
index_DLT=4
id_DL=c(1,2,3,6,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(300,2000),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)

pTox_sc_predict=pTox_tot(threshold=scenario[[id_sc]]$threshold, w_alpha=scenario[[id_sc]]$w_alpha,
                         Rmax_global=Rmax_global_predict,
                         DL=DL_kg_predict,
                         pTox_cond_fct=pTox_conditional_log,
                         param_cond_fct=list(beta0=scenario[[id_sc]]$beta0,beta1=scenario[[id_sc]]$beta1,beta2=scenario[[id_sc]]$beta2),
                         gamma=scenario[[id_sc]]$gamma)

scenario_predict[[id_sc]]=list(
  id_DL=1:nrow(DL_kg_predict),
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc_predict$pTox_CRS,
  pTox_oDLT=pTox_sc_predict$pTox_oDLT,
  pTox_DLT=pTox_sc_predict$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc_predict$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc_predict$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc_predict$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc_predict$pTox_oDLT_noCRS_cond
)


###############Set 1, MTD=6

id_sc=5
w_alpha=0.1
beta0=-14.8
beta1=1.5
beta2=1
gamma=2
target_DLT=0.3
index_DLT=6
id_DL=c(1,2,3,6,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(300,2000),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)


###############Set 1 new sc MTD s2
id_sc=6
w_alpha=0.1
beta0=-9.9
beta1=1.6
beta2=0.5
gamma=2
target_DLT=0.3
index_DLT=2
id_DL=c(1,3,5,7,8,9)

wm=getprior(0.05,0.3,4,6) 

op1=optimize(optim_scenario,interval = c(200,1500),
             Rmax_global=Rmax_global[,id_DL],
             DL=DL_kg[id_DL,],
             w_alpha=w_alpha,
             pTox_cond_fct=pTox_conditional_log,
             param_cond_fct=list(beta0,beta1,beta2),
             gamma=gamma,
             target_DLT=target_DLT,index_DLT=index_DLT)

pTox_sc=pTox_tot(threshold=op1$minimum, w_alpha=w_alpha,
                 Rmax_global=Rmax_global[,id_DL],
                 DL=DL_kg[id_DL,],
                 pTox_cond_fct=pTox_conditional_log,
                 param_cond_fct=list(beta0,beta1,beta2),
                 gamma=gamma)

pTox_sc

scenario[[id_sc]]=list(
  id_DL=id_DL,
  threshold=op1$minimum,
  w_alpha=w_alpha,
  beta0=beta0,
  beta1=beta1,
  beta2=beta2,
  gamma=gamma,
  target=target_DLT,
  pTox_CRS=pTox_sc$pTox_CRS,
  pTox_oDLT=pTox_sc$pTox_oDLT,
  pTox_DLT=pTox_sc$pTox_DLT,
  pTox_oDLT_CRS=pTox_sc$pTox_oDLT_CRS,
  pTox_oDLT_noCRS=pTox_sc$pTox_oDLT_noCRS,
  pTox_oDLT_CRS_cond=pTox_sc$pTox_oDLT_CRS_cond,
  pTox_oDLT_noCRS_cond=pTox_sc$pTox_oDLT_noCRS_cond,
  wm=wm
)



###########################################

# save(scenario,file=paste0(path_simus,"/scenario/scenario.Rdata"))
# save(scenario_predict,file=paste0(path_simus,"/scenario/scenario_predict.Rdata"))
