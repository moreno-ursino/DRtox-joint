source(paste0(path_simus,"/tools/pk_functions.R"))

library(dfcrm)
library(rstan)
library(parallel)
library(foreach)
library(doParallel)


################STAN MODELS

logit=function(x) {
  log(x/(1-x))
}

invlogit=function (x) {
  1/(1 + exp(-x))
}

Pi_logistic=function(beta0,beta1,pred){
  invlogit(beta0+beta1*pred)
}

pTox_indep=function(p1,p2){
  1-(1-p1)*(1-p2)
}

clayton_copula=function(p1,p2,gamma){
  p1+p2-pmax(p1^(-gamma)+p2^(-gamma)-1,rep(0,length(p1)))^(-1/gamma)
}


logistic_CRS_stan="
  functions{
    real invlogit(real x){
      return 1/(1 + exp(-x));
    }
  }
  data {
    int<lower=0> n_patients;         // number of patients
    int<lower=0,upper=1> y1[n_patients]; // CRS response
    vector[n_patients] logR;         // Log tranformation of the PD response
    vector[2] prior_beta0;  // param of the normal prior
    vector[2] prior_beta1;  // param of the gamma prior
  }
  parameters {
    real beta0;
    real<lower=0> beta1;
  }
  model {

    for (n in 1:n_patients) {
      target+= y1[n]*log(invlogit(beta0+beta1*logR[n]))+(1-y1[n])*log(1-invlogit(beta0+beta1*logR[n]));
    }

    beta0 ~ normal(prior_beta0[1],prior_beta0[2]);
    beta1 ~ gamma(prior_beta1[1],prior_beta1[2]);
  }"

sm_logistic_CRS=stan_model(model_code = logistic_CRS_stan,auto_write = T)

initfunction_CRS=function(chain_id = 1) {
  list(
    beta0=rnorm(1,0,1),
    beta1=rexp(1,1)
  )
} 

cum_oDLT_stan="
  functions{
    real invlogit(real x){
      return 1/(1 + exp(-x));
    }
  }
  data {
    int<lower=0> n_pat;         // number of patients
    vector[n_pat] last_admin;         // Patients admin
    int<lower=0,upper=1> y[n_pat]; // oDLT response
    vector[n_pat] logcumdose;         // Log tranformation of the cumulative dose
    vector[n_pat] logcumdose_1;         // Log tranformation of the cumulative dose without the last admin
    vector[2] prior_beta0;  // param of the normal prior
    vector[2] prior_beta1;  // param of the normal prior
  }
  parameters {
    real beta0;
    real beta1;
  }
  model {

    vector[n_pat] p_cum;

    vector[n_pat] p_cum_1;

    for (i in 1:n_pat) {
      p_cum[i]=invlogit(beta0+exp(beta1)*logcumdose[i]);

      if(last_admin[i]>1){
        p_cum_1[i]=invlogit(beta0+exp(beta1)*logcumdose_1[i]);
      } else{
        p_cum_1[i]=0;
      }
    }

    for (i in 1:n_pat) {
      target+= y[i]*log(p_cum[i]-p_cum_1[i])+(1-y[i])*log(1-p_cum[i]);
    }

    beta0 ~ normal(prior_beta0[1],prior_beta0[2]);
    beta1 ~ normal(prior_beta1[1],prior_beta1[2]);
  }"

sm_cum_oDLT=stan_model(model_code = cum_oDLT_stan,auto_write = T)

initfunction_oDLT=function(chain_id = 1) {
  list(
    beta0=rnorm(1,0,1),
    beta1=log(rexp(1,1))
  )
} 

clayton_stan="
  functions{
    real fun_clayton(real p1, real p2, real gamma){
      return fmax(p1^(-gamma)+p2^(-gamma)-1,0)^(-1/gamma);
    }

    real log_likelihood(real x1, real x2, real p1, real p2, real gamma){
      real p11;
      p11=fun_clayton(p1,p2,gamma);

      return x1*x2*log(p11)+x1*(1-x2)*log(p1-p11)+(1-x1)*x2*log(p2-p11)+(1-x1)*(1-x2)*log(1-p1-p2+p11);
    }
  }
  data {
    int<lower=0> K;         // number of dose-regimens
    int<lower=0> n;         // number of patients
    vector[K] p1;         // proba CRS
    vector[K] p2;         // proba oDLT

    int id_s[n];         // dose-regimens
    int x1[n];         // CRS data
    int x2[n];         // oDLT data

    real alpha;
    real beta;
  }
  parameters {
    real<lower=0> gamma;
  }
  model {

    for (i in 1:n) {
      target+= log_likelihood(x1[i],x2[i],p1[id_s[i]],p2[id_s[i]],gamma);
    }

    gamma ~ gamma(alpha,beta);

  }"

sm_clayton=stan_model(model_code = clayton_stan,auto_write = T)


initfunction_clayton=function(chain_id = 1) {
  gamma=rgamma(1,1,1)
  while(gamma>4){
    gamma=rgamma(1,1,1)
  }
  list(
    gamma=gamma
  )
} 


fun_mean_beta0_joint=function(wm,seq_ref){
  return(logit(1-sqrt(1-wm[seq_ref])))
}


optim_beta1_logistic=function(beta1,beta0,
                              pred,
                              tox_guess){
  sum(
    (Pi_logistic(beta0,beta1,pred)-tox_guess)^2
  )
  
}

t_function=function(t_doses,t_inf){
  t=NULL
  last_1=length(t_doses)-1
  
  t_break=10
  step1=0.5
  step2=0.1
  step3=24
  t_add=100
  
  if(last_1>0){
    t=c(
      t,
      sort(unique(as.numeric(sapply(1:last_1,function(x) 
        c(
          seq(t_doses[x],t_doses[x]+t_inf[x]/2,by=step1),
          seq(t_doses[x]+t_inf[x]/2,t_doses[x]+t_inf[x]/2+t_break,by=step2),
          seq(t_doses[x]+t_inf[x]/2+t_break,t_doses[x+1],by=step3)
        )
      ))))
    )
  }
  
  t=sort(unique(c(
    t,
    seq(t_doses[length(t_doses)],t_doses[length(t_doses)]+t_inf[length(t_doses)]/2,by=step1),
    seq(t_doses[length(t_doses)]+t_inf[length(t_doses)]/2,t_doses[length(t_doses)]+t_inf[length(t_doses)]/2+t_break,by=step2),
    seq(t_doses[length(t_doses)]+t_inf[length(t_doses)]/2+t_break,t_doses[length(t_doses)]+t_inf[length(t_doses)]/2+t_add,by=step3)
  )))
  
  
  return(t)
}

getNbTox=function(data,type_tox,n_doses){
  res=NULL
  
  if(type_tox=="CRS"){
    res=ifelse(
      rep(sum(data$CRS),n_doses)!=0,
      as.numeric(table(
        as.numeric(tapply(data$CRS,data$id,max)),
        factor(as.numeric(tapply(data$id_seq,data$id,max)),levels=1:n_doses)
      )[2,]),
      rep(0,n_doses))
  }
  
  if(type_tox=="oDLT"){
    res=ifelse(
      rep(sum(data$oDLT),n_doses)!=0,
      as.numeric(table(
        as.numeric(tapply(data$oDLT,data$id,max)),
        factor(as.numeric(tapply(data$id_seq,data$id,max)),levels=1:n_doses)
      )[2,]),
      rep(0,n_doses))
  }
  
  if(type_tox=="DLT"){
    res=ifelse(
      rep(sum(data$tox),n_doses)!=0,
      as.numeric(table(
        as.numeric(tapply(data$tox,data$id,max)),
        factor(as.numeric(tapply(data$id_seq,data$id,max)),levels=1:n_doses)
      )[2,]),
      rep(0,n_doses))
  }
  
  return(res)
}

return_PKPD_param=function(dir,
                           names_PK=c("Cl","V"),
                           names_PD=c("Emax","EC50","H","Imax","IC50","kdeg","K"),
                           fixed_param=c(EC50=10000,Imax=0.995,IC50=18200)){
  
  setwd(dir = dir)
  
  pop_param=read.table(file="results/populationParameters.txt",
                       sep=",",header=T)
  
  ind_param=read.table(file="results/IndividualParameters/estimatedIndividualParameters.txt",
                       sep=",",header=T)
  
  #############Population parameters#############
  
  ####Fixed effects vector####
  mu_PK=numeric(length(names_PK))
  names(mu_PK)=names_PK
  
  mu_PK_index=unlist(sapply(
    1:length(names_PK),
    function(x) grep(paste0(names_PK[x],"_pop"),pop_param$parameter)
  ))
  mu_PK_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[mu_PK_index])),
    function(x) grep(as.character(pop_param$parameter)[mu_PK_index[x]],paste0(names(mu_PK),"_pop"))
  ))
  
  mu_PK[mu_PK_index2]=pop_param$value[mu_PK_index]
  
  
  mu_PD=numeric(length(names_PD))
  names(mu_PD)=names_PD
  
  mu_PD_index=unlist(sapply(
    1:length(names_PD),
    function(x) grep(paste0(names_PD[x],"_pop"),pop_param$parameter)
  ))
  mu_PD_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[mu_PD_index])),
    function(x) grep(as.character(pop_param$parameter)[mu_PD_index[x]],paste0(names(mu_PD),"_pop"))
  ))
  
  mu_PD[mu_PD_index2]=pop_param$value[mu_PD_index]
  
  
  mu_PK_index_fixed=unlist(sapply(
    1:length(names_PK),
    function(x) grep(names_PK[x],names(fixed_param))
  ))
  
  if(length(mu_PK_index_fixed)>0){
    mu_PK_index_fixed2=unlist(sapply(
      1:length(fixed_param[mu_PK_index_fixed]),
      function(x) grep(names(fixed_param[mu_PK_index_fixed[x]]),names(mu_PK))
    ))
    
    mu_PK[mu_PK_index_fixed2]=fixed_param[mu_PK_index_fixed]
  }
  
  
  mu_PD_index_fixed=unlist(sapply(
    1:length(names_PD),
    function(x) grep(names_PD[x],names(fixed_param))
  ))
  
  if(length(mu_PD_index_fixed)>0){
    mu_PD_index_fixed2=unlist(sapply(
      1:length(fixed_param[mu_PD_index_fixed]),
      function(x) grep(names(fixed_param[mu_PD_index_fixed[x]]),names(mu_PD))
    ))
    
    mu_PD[mu_PD_index_fixed2]=fixed_param[mu_PD_index_fixed]
  }
  
  
  
  ####Variance covariance matrix####
  w_PK=matrix(0,length(mu_PK),length(mu_PK))
  rownames(w_PK)=colnames(w_PK)=names_PK
  
  w_PK_index=unlist(sapply(
    1:length(names_PK),
    function(x) grep(paste0("omega_",names_PK[x]),pop_param$parameter)
  ))
  w_PK_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[w_PK_index])),
    function(x) grep(as.character(pop_param$parameter)[w_PK_index[x]],paste0("omega_",names(mu_PK)))
  ))
  
  diag(w_PK)[w_PK_index2]=pop_param$value[w_PK_index]^2
  
  
  
  w_PD=matrix(0,length(mu_PD),length(mu_PD))
  rownames(w_PD)=colnames(w_PD)=names_PD
  
  w_PD_index=unlist(sapply(
    1:length(names_PD),
    function(x) grep(paste0("omega_",names_PD[x]),pop_param$parameter)
  ))
  w_PD_index2=unlist(sapply(
    1:length(as.character(pop_param$parameter[w_PD_index])),
    function(x) grep(as.character(pop_param$parameter)[w_PD_index[x]],paste0("omega_",names(mu_PD)))
  ))
  
  diag(w_PD)[w_PD_index2]=pop_param$value[w_PD_index]^2
  
  
  PK_param_pop=list(mu=mu_PK,
                    w=w_PK)
  PD_param_pop=list(mu=mu_PD,
                    w=w_PD)
  
  
  #############Individual parameters#############
  PK_param_ind=matrix(0,nrow=nrow(ind_param),ncol=length(names_PK)+1)
  
  colnames(PK_param_ind)=c("id",names_PK)
  PK_param_ind[,1]=ind_param$id
  
  ind_PK_index=unlist(sapply(
    1:length(colnames(PK_param_ind)),
    function(x) grep(paste0(colnames(PK_param_ind)[x],"_SAEM"),colnames(ind_param))
  ))
  
  ind_PK_index2=unlist(sapply(
    1:length(colnames(ind_param)[ind_PK_index]),
    function(x) grep(colnames(ind_param)[ind_PK_index[x]],paste0(colnames(PK_param_ind),"_SAEM"))
  ))
  
  PK_param_ind[,ind_PK_index2]=as.matrix(ind_param[,ind_PK_index])
  
  
  
  PD_param_ind=matrix(0,nrow=nrow(ind_param),ncol=length(names_PD)+1)
  
  colnames(PD_param_ind)=c("id",names_PD)
  PD_param_ind[,1]=ind_param$id
  
  ind_PD_index=unlist(sapply(
    1:length(colnames(PD_param_ind)),
    function(x) grep(paste0(colnames(PD_param_ind)[x],"_SAEM"),colnames(ind_param))
  ))
  
  ind_PD_index2=unlist(sapply(
    1:length(colnames(ind_param)[ind_PD_index]),
    function(x) grep(colnames(ind_param)[ind_PD_index[x]],paste0(colnames(PD_param_ind),"_SAEM"))
  ))
  
  PD_param_ind[,ind_PD_index2]=as.matrix(ind_param[,ind_PD_index])
  
  
  if(length(mu_PK_index_fixed)>0){
    
    PK_param_ind[,1+mu_PK_index_fixed2]=matrix(rep(fixed_param[mu_PK_index_fixed],each=nrow(PK_param_ind)),nrow=nrow(PK_param_ind))
    
  }
  
  if(length(mu_PD_index_fixed)>0){
    
    PD_param_ind[,1+mu_PD_index_fixed2]=matrix(rep(fixed_param[mu_PD_index_fixed],each=nrow(PD_param_ind)),nrow=nrow(PD_param_ind))
    
  }
  
  return(list(
    PK_param_pop=PK_param_pop,
    PD_param_pop=PD_param_pop,
    PK_param_ind=PK_param_ind,
    PD_param_ind=PD_param_ind
  ))
  
  
  
}

return_data=function(dir,
                     t,
                     DL,
                     t_doses,
                     t_inf,
                     PK_model,
                     PK_param_ind,
                     PD_model,
                     PD_param_ind,
                     PD_initial_values,
                     method_ode){
  
  data=read.table(file=paste0(dir,"/data.txt"),header=T)
  
  data$R=0
  
  
  for(i in 1:nrow(PK_param_ind)){
    
    Rmax_ind=PK_PD_simu_Rmax_local_ind(t=t,
                                       doses=DL[data[data$id==i,][1,"id_seq"],1:sum(data$id==i)],
                                       t_doses=t_doses[1:sum(data$id==i)],
                                       t_inf=t_inf[data[data$id==i,][1,"id_seq"],1:sum(data$id==i)],
                                       PK_model=PK_model,
                                       PK_param_ind=as.list(PK_param_ind[i,2:ncol(PK_param_ind)]),
                                       PD_model=PD_model,
                                       PD_param_ind=as.list(PD_param_ind[i,2:ncol(PD_param_ind)]),
                                       PD_initial_values=PD_initial_values,
                                       method_ode=method_ode)
    
    data$R[data$id==i]=as.numeric(Rmax_ind)
    
    
  }
  
  return(data=data)
  
  
}


estim_one_trial=function(
  itrial,
  sc,
  design,
  dir,
  names_PK,
  names_PD,
  fixed_param,
  t,
  DL,
  t_doses,
  t_inf,
  PK_model,
  PD_model,
  PD_initial_values,
  method_ode,
  seq_ref,
  seed,
  n_chains,
  n_warmup,
  n_PKPD,
  sd_beta0,
  alpha,
  sd_beta1,
  alpha_clayton,
  beta_clayton,
  wm,
  sim=T){
  
  DL_kg=DL/70
  
  ######Results definition
  seed_trial=seed+itrial
  n_doses=nrow(DL)
  
  n_predict=ifelse(!is.null(DL_predict),nrow(DL_predict),0)
  
  if(n_predict==0){
    DL_kg_predict=NULL
  }else{
    DL_kg_predict=DL_predict/70
  }
  
  
  ######Data definition
  PKPD_param=return_PKPD_param(dir,names_PK,names_PD,fixed_param)
  
  PK_param_pop=PKPD_param$PK_param_pop
  PD_param_pop=PKPD_param$PD_param_pop
  PK_param_ind=PKPD_param$PK_param_ind
  PD_param_ind=PKPD_param$PD_param_ind
  
  
  data=return_data(dir,t,DL,t_doses,t_inf,
                   PK_model,PK_param_ind,
                   PD_model,PD_param_ind,
                   PD_initial_values,
                   method_ode)
  
  t=t_function(t_doses,t_inf[1,])
  
  Rmax_pop=PK_PD_simu_Rmax(
    t,
    doses=DL,
    t_doses,
    t_inf,
    PK_model,
    PK_param_pop=list(
      mu=PK_param_pop$mu,
      w=matrix(0,length(PK_param_pop$mu),length(PK_param_pop$mu))
    ),
    PD_model,
    PD_param_pop=list(
      mu=PD_param_pop$mu,
      w=matrix(0,length(PD_param_pop$mu),length(PD_param_pop$mu))
    ),
    PD_initial_values,
    method_ode=method_ode
  )
  
  ###Dose/kg
  data$dose=data$dose/70
  data$cumdose=as.numeric(unlist(tapply(data$dose,data$id,cumsum)))
  data$cumdose_1=data$cumdose-data$dose
  
  data_aggreg=aggregate(data,by=list(data$id),max)[,-1]
  
  id_CRS=unique(data$id)[as.numeric(tapply(data$CRS,data$id,max))==1]
  
  data_noCRS=data[!data$id %in%id_CRS,]
  
  neighbor_ref=(seq_ref-1):(seq_ref+1)
  neighbor_ref=neighbor_ref[neighbor_ref %in% 1:n_doses]
  
  mean_beta0=fun_mean_beta0_joint(wm,seq_ref)
  
  optim_beta1_CRS=optim(
    par=c(0), fn=optim_beta1_logistic,
    pred=log(Rmax_pop[neighbor_ref]/Rmax_pop[seq_ref]),
    beta0=mean_beta0,
    tox_guess=1-sqrt(1-wm[neighbor_ref]),
    method="L-BFGS-B",
    lower = c(0), upper = c(Inf)
  )
  
  optim_beta1_oDLT=optim(
    par=c(0), fn=optim_beta1_logistic,
    pred=log(apply(DL_kg,1,sum)[neighbor_ref]/apply(DL_kg,1,sum)[seq_ref]),
    beta0=mean_beta0,
    tox_guess=1-sqrt(1-wm[neighbor_ref]),
    method="L-BFGS-B",
    lower = c(0), upper = c(Inf)
  )
  
  
  #####Initial values
  set.seed(seed_trial)
  init_CRS=lapply(1:n_chains, function(id) initfunction_CRS(chain_id = id))
  init_oDLT=lapply(1:n_chains, function(id) initfunction_oDLT(chain_id = id))
  init_clayton=lapply(1:n_chains, function(id) initfunction_clayton(chain_id = id))
  
  set.seed(seed_trial)
  
  #############CRS model
  sampl_CRS=sampling(sm_logistic_CRS,
                     data = list(
                       n_patients=length(unique(data$id)),
                       y1=as.numeric(tapply(data$CRS,data$id,max)),
                       logR=log(as.numeric(tapply(data$R,data$id,max))/Rmax_pop[seq_ref]),
                       prior_beta0=c(mean_beta0,sd_beta0),
                       prior_beta1=c(alpha,alpha/optim_beta1_CRS$par)
                     ),
                     iter=2*n_warmup, warmup=n_warmup,
                     chains = n_chains,control = list(adapt_delta = 0.99),
                     cores=1,seed=seed_trial,init=init_CRS)
  
  beta0=extract(sampl_CRS)$beta0
  beta1=extract(sampl_CRS)$beta1
  
  DL_tot=rbind(DL,DL_predict)
  DL_kg_tot=rbind(DL_kg,DL_kg_predict)
  t_inf_tot=rbind(t_inf,t_inf_predict)
  
  Rmax_function=function(i){
    cbind(
      rep(i,length(t_doses)),
      1:length(t_doses),
      PK_PD_simu_Rmax_local(t=t,
                            doses=DL_tot,
                            t_doses=t_doses,
                            t_inf=t_inf_tot,
                            PK_model=PK_model,
                            PK_param_pop=PK_param_pop,
                            PD_model=PD_model,
                            PD_param_pop=PD_param_pop,
                            PD_initial_values=PD_initial_values))
  }
  
  set.seed(seed_trial)
  Rmax_MC_local=do.call(rbind,lapply(1:n_PKPD,Rmax_function))
  if(n_predict==0){
    colnames(Rmax_MC_local)=c("id","admin",paste0("DL",1:n_doses))
  }else{
    colnames(Rmax_MC_local)=c("id","admin",paste0("DL",1:n_doses),paste0("DL_predict",1:n_predict))
  }
  
  
  Rmax_MC_global=aggregate(Rmax_MC_local[,grep("DL",colnames(Rmax_MC_local))],
                           by=list(Rmax_MC_local[,"id"]),
                           max,na.rm=T)[,-1]
  
  
  logR_MC=as.matrix(log(Rmax_MC_global/Rmax_pop[seq_ref]))
  
  p_estim_CRS=matrix(sapply(1:length(beta0),function(y){
    colMeans(matrix(sapply(1:nrow(logR_MC),function(x){
      Pi_logistic(beta0=beta0[y],beta1=beta1[y],
                  pred=logR_MC[x,])
    }),nrow=nrow(logR_MC),byrow=T))
  }),nrow=length(beta0),byrow = T)
  
  
  #############DLTo model
  sampl_oDLT=sampling(sm_cum_oDLT,
                      data = list(
                        n_pat=length(unique(data$id)),
                        last_admin=as.numeric(tapply(data$admin,data$id,max)),
                        y=as.numeric(tapply(data$oDLT,data$id,max)),
                        logcumdose=log(as.numeric(tapply(data$cumdose,data$id,max))/apply(DL_kg,1,sum)[seq_ref]),
                        logcumdose_1=log(as.numeric(tapply(data$cumdose_1,data$id,max))/apply(DL_kg,1,sum)[seq_ref]),
                        prior_beta0=c(mean_beta0,sd_beta0),
                        prior_beta1=c(log(optim_beta1_oDLT$par),sd_beta1)
                      ),
                      iter=2*n_warmup, warmup=n_warmup,
                      chains = n_chains,control = list(adapt_delta = 0.99),
                      cores=1,seed=seed_trial,init=init_oDLT)
  beta0=extract(sampl_oDLT)$beta0
  beta1=exp(extract(sampl_oDLT)$beta1)
  p_estim_oDLT=matrix(sapply(1:length(beta0),function(x){
    Pi_logistic(
      beta0=beta0[x],
      beta1=beta1[x],
      pred=log(apply(DL_kg_tot,1,sum)/apply(DL_kg,1,sum)[seq_ref])
    )
  }),nrow=length(beta0),byrow = T)
  
  sampl_oDLT_noCRS=sampling(sm_cum_oDLT,
                            data = list(
                              n_pat=length(unique(data_noCRS$id)),
                              last_admin=as.numeric(tapply(data_noCRS$admin,data_noCRS$id,max)),
                              y=as.numeric(tapply(data_noCRS$oDLT,data_noCRS$id,max)),
                              logcumdose=log(as.numeric(tapply(data_noCRS$cumdose,data_noCRS$id,max))/apply(DL_kg,1,sum)[seq_ref]),
                              logcumdose_1=log(as.numeric(tapply(data_noCRS$cumdose_1,data_noCRS$id,max))/apply(DL_kg,1,sum)[seq_ref]),
                              prior_beta0=c(mean_beta0,sd_beta0),
                              prior_beta1=c(log(optim_beta1_oDLT$par),sd_beta1)
                            ),
                            iter=2*n_warmup, warmup=n_warmup,
                            chains = n_chains,control = list(adapt_delta = 0.99),
                            cores=1,seed=seed_trial,init=init_oDLT)
  beta0=extract(sampl_oDLT_noCRS)$beta0
  beta1=exp(extract(sampl_oDLT_noCRS)$beta1)
  p_estim_oDLT_noCRS=matrix(sapply(1:length(beta0),function(x){
    Pi_logistic(
      beta0=beta0[x],
      beta1=beta1[x],
      pred=log(apply(DL_kg_tot,1,sum)/apply(DL_kg,1,sum)[seq_ref])
    )
  }),nrow=length(beta0),byrow = T)
  
  #############DLT model
  
  ####Independent
  p_estim_DLT_indep=matrix(sapply(1:nrow(p_estim_oDLT),function(x){
    pTox_indep(
      p1=p_estim_CRS[x,],
      p2=p_estim_oDLT[x,]
    )
  }),nrow=nrow(p_estim_oDLT),byrow = T)
  
  ####Conditional
  p_estim_DLT_cond=matrix(sapply(1:nrow(p_estim_oDLT_noCRS),function(x){
    pTox_indep(
      p1=p_estim_CRS[x,],
      p2=p_estim_oDLT_noCRS[x,]
    )
  }),nrow=nrow(p_estim_oDLT_noCRS),byrow = T)
  
  ####Copula
  data_aggreg$id_seq=factor(data_aggreg$id_seq,levels=1:n_doses)
  
  if(sum(data_aggreg$tox)){
    x=as.numeric(table(data_aggreg$id_seq,data_aggreg$tox)[,2])
  } else{
    x=rep(0,n_doses)
  }
  
  p1=colMeans(p_estim_CRS)[1:n_doses]
  p2=colMeans(p_estim_oDLT)[1:n_doses]
  
  p1_pred=colMeans(p_estim_CRS)
  p2_pred=colMeans(p_estim_oDLT)
  
  sampl_clayton=sampling(sm_clayton,
                         data = list(
                           K=n_doses,
                           n=nrow(data_aggreg),
                           p1=p1,
                           p2=p2,
                           id_s=as.numeric(as.character(data_aggreg$id_seq)),
                           x1=data_aggreg$CRS,
                           x2=data_aggreg$oDLT,
                           alpha=alpha_clayton,
                           beta=beta_clayton
                         ),
                         iter=2*n_warmup, warmup=n_warmup,
                         chains = n_chains,control = list(adapt_delta = 0.99),
                         cores=1,seed=seed_trial,init=init_clayton)
  
  gamma=extract(sampl_clayton)$gamma
  p_estim_DLT_copula=matrix(sapply(1:length(gamma),function(x){
    clayton_copula(
      p1=p1_pred,
      p2=p2_pred,
      gamma=gamma[x]
    )
  }),nrow=length(gamma),byrow = T)
  
  ###OUT
  
  if(sim){
    res=data.frame(matrix(NA,
                          nrow=4,
                          ncol=7*n_doses+11+3*n_predict))
    
    if(n_predict==0){
      colnames(res)=c("trial","Scenario","n","design",
                      "Method",
                      "Seed",
                      paste0("pTox_s",1:n_doses),
                      paste0("pTox_CRS_s",1:n_doses),
                      paste0("pTox_oDLT_s",1:n_doses),
                      "beta0_1","beta1_1","beta0_2","beta1_2","asso_param",
                      paste0("n_s",1:n_doses),
                      paste0("n_DLT_s",1:n_doses),
                      paste0("n_CRS_s",1:n_doses),
                      paste0("n_oDLT_s",1:n_doses))
    }else{
      colnames(res)=c("trial","Scenario","n","design",
                      "Method",
                      "Seed",
                      paste0("pTox_s",1:n_doses),
                      paste0("pTox_predict_s",1:n_predict),
                      paste0("pTox_CRS_s",1:n_doses),
                      paste0("pTox_predict_CRS_s",1:n_predict),
                      paste0("pTox_oDLT_s",1:n_doses),
                      paste0("pTox_predict_oDLT_s",1:n_predict),
                      "beta0_1","beta1_1","beta0_2","beta1_2","asso_param",
                      paste0("n_s",1:n_doses),
                      paste0("n_DLT_s",1:n_doses),
                      paste0("n_CRS_s",1:n_doses),
                      paste0("n_oDLT_s",1:n_doses))
    }
    
    res[1,]=c(itrial,sc,length(unique(data$id)),design,"DRtox_indep",seed_trial,
              colMeans(p_estim_DLT_indep),
              colMeans(p_estim_CRS),
              colMeans(p_estim_oDLT),
              mean(extract(sampl_CRS)$beta0),
              mean(extract(sampl_CRS)$beta1),
              mean(extract(sampl_oDLT)$beta0),
              mean(extract(sampl_oDLT)$beta1),
              NA,
              as.numeric(table(factor(data[data$admin==1,"id_seq"],levels=1:n_doses))),
              getNbTox(data=data,"DLT",n_doses),
              getNbTox(data=data,"CRS",n_doses),
              getNbTox(data=data,"oDLT",n_doses)
    )
    
    res[2,]=c(
      itrial,sc,length(unique(data$id)),design,"DRtox_copula",seed_trial,
      colMeans(p_estim_DLT_copula),
      colMeans(p_estim_CRS),
      colMeans(p_estim_oDLT),
      mean(extract(sampl_CRS)$beta0),
      mean(extract(sampl_CRS)$beta1),
      mean(extract(sampl_oDLT)$beta0),
      mean(extract(sampl_oDLT)$beta1),
      mean(extract(sampl_clayton)$gamma),
      as.numeric(table(factor(data[data$admin==1,"id_seq"],levels=1:n_doses))),
      getNbTox(data=data,"DLT",n_doses),
      getNbTox(data=data,"CRS",n_doses),
      getNbTox(data=data,"oDLT",n_doses)
    )
    
    res[3,]=c(itrial,sc,length(unique(data$id)),design,"DRtox_cond",seed_trial,
              colMeans(p_estim_DLT_cond),
              colMeans(p_estim_CRS),
              colMeans(p_estim_oDLT_noCRS),
              mean(extract(sampl_CRS)$beta0),
              mean(extract(sampl_CRS)$beta1),
              mean(extract(sampl_oDLT_noCRS)$beta0),
              mean(extract(sampl_oDLT_noCRS)$beta1),
              NA,
              as.numeric(table(factor(data[data$admin==1,"id_seq"],levels=1:n_doses))),
              getNbTox(data=data,"DLT",n_doses),
              getNbTox(data=data,"CRS",n_doses),
              getNbTox(data=data,"oDLT",n_doses)
    )
    
    
    res[4,]=c(
      itrial,sc,length(unique(data$id)),design,"crm",seed_trial,
      as.numeric(unlist(read.table(file=paste0(dir,"/pTox.txt")))),
      rep(NA,2*n_doses),
      rep(NA,3*n_predict),
      rep(NA,5),
      as.numeric(table(factor(data[data$admin==1,"id_seq"],levels=1:n_doses))),
      getNbTox(data=data,"DLT",n_doses),
      getNbTox(data=data,"CRS",n_doses),
      getNbTox(data=data,"oDLT",n_doses)
    )
    
  }else{
    res=list(
      sc=sc,
      itrial=itrial,
      seed_trial=seed_trial,
      wm=wm,
      seq_ref=seq_ref,
      n_chains=n_chains,
      n_warmup=n_warmup,
      n_PKPD=n_PKPD,
      weight_std=70,
      DL_kg=DL_kg,
      t_doses=t_doses,
      t_inf=t_inf,
      DL_kg_predict=DL_kg_predict,
      PK_param_pop=PK_param_pop,
      PD_param_pop=PD_param_pop,
      PK_param_ind=PK_param_ind,
      PD_param_ind=PD_param_ind,
      data=data,
      data_noCRS=data_noCRS,
      t=t,
      Rmax_pop=Rmax_pop,
      Rmax_MC_global=Rmax_MC_global,
      logR_MC=logR_MC,
      prior_beta0_CRS=c(mean_beta0,sd_beta0),
      prior_beta1_CRS=c(alpha,alpha/optim_beta1_CRS$par),
      init_CRS=init_CRS,
      prior_beta0_oDLT=c(mean_beta0,sd_beta0),
      prior_beta1_oDLT=c(log(optim_beta1_oDLT$par),sd_beta1),
      init_oDLT=init_oDLT,
      beta0_CRS=extract(sampl_CRS)$beta0,
      beta1_CRS=extract(sampl_CRS)$beta1,
      p_CRS=p_estim_CRS,
      beta0_oDLT=extract(sampl_oDLT)$beta0,
      beta1_oDLT=exp(extract(sampl_oDLT)$beta1),
      p_oDLT=p_estim_oDLT,
      beta0_oDLT_noCRS=extract(sampl_oDLT_noCRS)$beta0,
      beta1_oDLT_noCRS=exp(extract(sampl_oDLT_noCRS)$beta1),
      p_oDLT_noCRS=p_estim_oDLT_noCRS,
      init_copula=init_clayton,
      gamma=extract(sampl_clayton)$gamma,
      p_DLT_indep=p_estim_DLT_indep,
      p_DLT_cond=p_estim_DLT_cond,
      p_DLT_copula=p_estim_DLT_copula
    )
  }
  
  return(res)
}

estim_trials=function(
  n_trials,
  sc,
  design,
  names_PK,
  names_PD,
  fixed_param,
  DL,
  t_doses,
  t_inf,
  DL_predict=NULL,
  t_inf_predict=NULL,
  PK_model,
  PD_model,
  PD_initial_values,
  method_ode,
  seq_ref,
  seed,
  n_chains,
  n_warmup,
  n_PKPD,
  sd_beta0,
  alpha,
  sd_beta1,
  alpha_clayton,
  beta_clayton,
  wm,
  N_sim,
  no_cores,
  save_output,
  start=0,
  name_temp="temp"){
  
  
  setwd(dir_trials)
  
  vect_packages=c("deSolve","rstan","MASS","reshape2","dfcrm")
  
  N_pack=n_trials/N_sim
  
  dir.create(name_temp,showWarnings = FALSE)
  if(length(list.files(name_temp))>0){invisible(file.remove(paste(name_temp,"/",list.files(name_temp),sep="")))}
  t1 <- Sys.time()
  dump("t1",paste0(name_temp,"/Start.txt"))
  cat("Remaining",file=paste0(name_temp,"/Remain.txt"))
  
  ########################################
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  res_trials<- foreach(ip=1:N_pack,.combine=rbind,.inorder=FALSE,.packages = vect_packages,.export = ls(globalenv())) %dopar% { 
    res=NULL
    
    for(i in 1:N_sim){
      
      cat(ip,"/",N_pack,"-",i,"/",N_sim,"\n")
      
      res=rbind(
        res,
        estim_one_trial(
          itrial=(ip-1)*N_sim+i+start,
          sc,
          design,
          dir=paste0(dir_trials,"/trial",(ip-1)*N_sim+i+start),
          names_PK,
          names_PD,
          fixed_param,
          t,
          DL,
          t_doses,
          t_inf,
          PK_model,
          PD_model,
          PD_initial_values,
          method_ode,
          seq_ref,
          seed,
          n_chains,
          n_warmup,
          n_PKPD,
          sd_beta0,
          alpha,
          sd_beta1,
          alpha_clayton,
          beta_clayton,
          wm,
          sim=T)
      )
      
      
    } ## End loop i (simulations)
    
    setwd(dir_trials)
    
    sink(paste(tempfile(tmpdir=name_temp),"_",ip,".txt",sep=""))
    sink()
    
    t2 <- Sys.time()
    source(paste0(name_temp,"/Start.txt"))
    nfiles <- length(list.files(name_temp))-2
    dt <- as.numeric(difftime(t2,t1,units="mins"))/nfiles
    remain <- dt*(N_pack-nfiles)
    hrs <- remain %/% 60
    mins <- ceiling(remain %% 60)
    cat("Remaining ",hrs,"h",mins,"min",sep="",file=paste0(name_temp,"/Remain.txt"))
    
    res
    
  }
  
  stopCluster(cl)
  
  if(is.null(save_output)){
    return(res_trials)
  } else{
    
    setwd(dir_trials)
    
    if (!file.exists("results")){
      dir.create("results",showWarnings = FALSE)
    }
    
    save(res_trials,file=paste0("results/",save_output,".RData"))
    return("OK")
  }
  
}






