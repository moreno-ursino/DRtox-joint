#######################################################################################################
#######################################################################################################
############################################PK/PD functions############################################
#######################################################################################################
#######################################################################################################


############################################################################
##################################Packages##################################
############################################################################

library(deSolve) #To resolve differential equations
library(MASS) 

############################################################################
################################PK functions################################
############################################################################

c_infusion_1cpt_multiple=function(t,doses=c(100,100,100),t_doses=c(0,5,10),
                                  t_inf=c(1,1,1),
                                  PK_param_individual=list(Cl=0.2,V=2)){
  ##1 compartment PK model with multiple infusion administrations
  with(PK_param_individual,{
    k=Cl/V
    
    f=function(t1){
      
      doses_cur=doses[(t1-t_doses)>=0]
      t_doses_cur=t_doses[(t1-t_doses)>=0]
      t_inf_cur=t_inf[(t1-t_doses)>=0]
      
      n_doses=length(doses_cur)
      
      return(
        ifelse(
          (t1-t_doses_cur[n_doses])<=t_inf_cur[n_doses],
          sum(doses_cur[-n_doses]/(t_inf_cur[-n_doses]*k*V)*(1-exp(-k*t_inf_cur[-n_doses]))*
                exp(-k*(t1-t_doses_cur[-n_doses]-t_inf_cur[-n_doses])) ) + 
            doses_cur[n_doses]/(t_inf_cur[n_doses]*k*V)*(1-exp(-k*(t1-t_doses_cur[n_doses]))) ,
          sum(doses_cur/(t_inf_cur*k*V)*(1-exp(-k*t_inf_cur))*exp(-k*(t1-t_doses_cur-t_inf_cur)) )
        )
      )
      
    }
    
    sapply(t,f)
  })
}


############################################################################
################################PD functions################################
############################################################################

PD_model_chen=function(t,initial_values,PD_param_ind,
                       doses=c(10,25,50,100,100,100,100),t_doses=c(1,4,8,15,22,29,36),
                       t_inf=NULL,
                       PK_model,PK_param_ind) {
  
  n_admin=7
  
  C=PK_model(t,doses,t_doses,t_inf,PK_param_ind)
  
  with(as.list(c(PD_param_ind,initial_values)),{
    dAUC=R
    
    dR=Emax*C^H/(EC50^H+C^H)*(1-Imax*AUC/(IC50/K^(n_admin-1)+AUC))-kdeg*R
    
    list(c(dAUC,dR),"C"=C)
  }) 
}

PD_initial_values_chen=function(PD_param){
  AUC0=0
  R0=0
  return(c(AUC=AUC0,R=R0))
}

############################################################################
##############################PK/PD simulation##############################
############################################################################

simu_individual_param=function(pop_param=list(mu=c(k=0.1,V=2),w=matrix(c(0.1,0,0,0.1),2,2))){
  
  #Function description:
  ##Simulation of patient specific PK/PD parameters using the log-normal distribution
  
  #Parameters:
  ##pop_param: Population parameters as a list of the fixed effects and the variance/covariance matrix
  
  with(pop_param,{
    
    return(as.list(mu*exp(mvrnorm(n = 1,rep(0,length(mu)), w))))
    
  })
  
} 


PK_PD_simu=function(t=seq(0,50,by=1),
                    doses=c(10,25,50,100,100,100,100),
                    t_doses=c(1,4,8,15,22,29,36),
                    t_inf=NULL,
                    PK_model=c_bolus_1cpt_multiple,
                    PK_param_pop=list(mu=c(k=0.1,V=2),w=matrix(c(0.1,0,0,0.1),2,2)),
                    PD_model=PD_model_precursor,
                    PD_param_pop=list(
                      mu=c(Smax=1,IC50=500,k0=30,kp=0.1,ks=0,kout=2),
                      w=matrix(0,6,6)
                    ),
                    PD_initial_values=PD_initial_values_precursor,
                    method_ode="vode"){
  
  ##Simulation of the entire PK/PD profile (prediction from the PK/PD model)
  
  n_doses=ifelse(is.null(nrow(doses)),1,nrow(doses))
  
  doses=matrix(doses,nrow=n_doses)
  
  if(!is.null(t_inf)) {
    t_inf=matrix(t_inf,nrow=n_doses)
  }
  
  PK_param_ind=simu_individual_param(PK_param_pop)
  
  PD_param_ind=simu_individual_param(PD_param_pop)
  
  initial_values=PD_initial_values(PD_param_ind)
  
  C=matrix(NA,nrow=length(t),ncol=n_doses)
  R=matrix(NA,nrow=length(t),ncol=n_doses)
  
  for(i in 1:n_doses){
    res=ode(y=initial_values,times=t,func=PD_model,method=method_ode,
            parms=PD_param_ind,doses=doses[i,],t_doses=t_doses,
            t_inf=t_inf[i,],
            PK_model=PK_model,
            PK_param_ind=PK_param_ind)
    
    C[,i]=res[,"C"]
    R[,i]=res[,"R"]
    
  }
  
  return(list(t=t,C=C,R=R,param=unlist(c(PK_param_ind,PD_param_ind))))
  
}

PK_PD_simu_ind=function(t=seq(0,50,by=1),
                        doses=c(10,25,50,100,100,100,100),
                        t_doses=c(1,4,8,15,22,29,36),
                        t_inf=NULL,
                        PK_model=c_bolus_1cpt_multiple,
                        PK_param_ind=c(k=0.1,V=2),
                        PD_model=PD_model_precursor,
                        PD_param_ind=c(Smax=1,IC50=500,k0=30,kp=0.1,ks=0,kout=2),
                        PD_initial_values=PD_initial_values_precursor,
                        method_ode="vode"){
  
  ##Simulation of the entire PK/PD profile (prediction from the PK/PD model)
  
  initial_values=PD_initial_values(PD_param_ind)
  
  
  res=ode(y=initial_values,times=t,func=PD_model,method=method_ode,
          parms=PD_param_ind,doses=doses,t_doses=t_doses,
          t_inf=t_inf,
          PK_model=PK_model,
          PK_param_ind=PK_param_ind)
  
  C=res[,"C"]
  R=res[,"R"]
  
  
  
  return(list(t=t,C=C,R=R))
  
}


PK_PD_simu_Rmax=function(t=seq(0,50,by=0.01),
                         doses=c(10,25,50,100,100,100,100),
                         t_doses=c(1,4,8,15,22,29,36),
                         t_inf=NULL,
                         PK_model=c_bolus_1cpt_multiple,
                         PK_param_pop=list(mu=c(k=0.1,V=2),w=matrix(c(0.1,0,0,0.1),2,2)),
                         PD_model=PD_model_precursor,
                         PD_param_pop=list(
                           mu=c(Smax=1,IC50=500,k0=30,kp=0.1,ks=0,kout=2),
                           w=matrix(0,6,6)
                         ),
                         PD_initial_values=PD_initial_values_precursor,
                         method_ode="vode"){
  
  ind_profile=PK_PD_simu(t,doses,t_doses,t_inf,
                         PK_model,PK_param_pop,PD_model,PD_param_pop,
                         PD_initial_values,
                         method_ode=method_ode)$R
  
  n_doses=ifelse(is.null(nrow(doses)),1,nrow(doses))
  
  return(apply(matrix(ind_profile,ncol=n_doses),2,max))
}


extract_max=function(R,time,doses=c(1,0,3,rep(10,4)),t_doses=c(1,4,8,15,22,29,36)){
  
  time_cut=as.numeric(cut(time,c(t_doses,time[length(time)]),include.lowest=T,right = F))
  
  time_split=split(time,time_cut)
  
  R_split=lapply(time_split,FUN=function(x) R[time %in% x])
  
  R_max_local=sapply(R_split,max)
  
  t_max_local=sapply(1:length(t_doses),FUN=function(x) time_split[[x]][R_split[[x]]==R_max_local[x]])
  
  R_max_global=max(R_max_local)
  t_max_global=t_max_local[R_max_local==R_max_global]
  
  R_max_cum=cummax(R_max_local)
  
  R_max_local[doses==0]=NA
  t_max_local[doses==0]=NA
  R_max_cum[doses==0]=NA
  
  return(list(R_max_global=R_max_global,t_max_global=t_max_global,
              R_max_local=R_max_local,t_max_local=t_max_local,
              R_max_cum=R_max_cum))
}


extract_max_global=function(R,time,doses,t_doses){
  return(extract_max(R,time,doses,t_doses)$R_max_global)
}


extract_max_local=function(R,time,doses,t_doses){
  return(extract_max(R,time,doses,t_doses)$R_max_local)
}


PK_PD_simu_Rmax_local=function(t=seq(0,50,by=0.01),
                               doses=c(10,25,50,100,100,100,100),
                               t_doses=c(1,4,8,15,22,29,36),
                               t_inf=NULL,
                               PK_model=c_bolus_1cpt_multiple,
                               PK_param_pop=list(mu=c(k=0.1,V=2),w=matrix(c(0.1,0,0,0.1),2,2)),
                               PD_model=PD_model_precursor,
                               PD_param_pop=list(
                                 mu=c(Smax=1,IC50=500,k0=30,kp=0.1,ks=0,kout=2),
                                 w=matrix(0,6,6)
                               ),
                               PD_initial_values=PD_initial_values_precursor,
                               method_ode="vode"){
  
  ind_profile=PK_PD_simu(t,doses,t_doses,t_inf,PK_model,PK_param_pop,PD_model,PD_param_pop,
                         PD_initial_values,
                         method_ode=method_ode)$R
  
  
  n_doses=ifelse(is.null(nrow(doses)),1,nrow(doses))
  
  res=apply(matrix(ind_profile,ncol=n_doses),2,extract_max_local,t,doses=rep(1,length(t_doses)),t_doses)
  
  res[t(doses==0)]=NA
  
  return(res)
  
}

PK_PD_simu_Rmax_local_ind=function(t=seq(0,50,by=0.01),
                               doses=c(10,25,50,100,100,100,100),
                               t_doses=c(1,4,8,15,22,29,36),
                               t_inf=NULL,
                               PK_model=c_bolus_1cpt_multiple,
                               PK_param_ind=c(k=0.1,V=2),
                               PD_model=PD_model_precursor,
                               PD_param_ind=c(Smax=1,IC50=500,k0=30,kp=0.1,ks=0,kout=2),
                               PD_initial_values=PD_initial_values_precursor,
                               method_ode="vode"){
  
  ind_profile=PK_PD_simu_ind(t,doses,t_doses,t_inf,PK_model,PK_param_ind,PD_model,PD_param_ind,
                         PD_initial_values,
                         method_ode=method_ode)$R
  res=extract_max_local(ind_profile,t,doses=rep(1,length(t_doses)),t_doses)

  res[doses==0]=NA
  
  return(res)
  
}


