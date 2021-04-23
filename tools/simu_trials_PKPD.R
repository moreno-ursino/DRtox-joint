source(paste0(path_simus,"/tools/pk_functions.R"))
source(paste0(path_simus,"/tools/CRM_2_param/CRM_2p_general.R"))

library(dfcrm)
library(parallel)
library(foreach)
library(doParallel)

######PK/PD model

create_model_Chen = function(names_PD_param=names(mu_PD),fixed_PD_param=fixed_param,n_admin=7){
  PKPD_model = paste0(
    "DESCRIPTION: PK infusion + PD Chen

input =  {V, Cl ,", paste0(names_PD_param[-which(names_PD_param %in% names(fixed_PD_param))],collapse=", "),"}
    
EQUATION:
Cc = pkmodel(V, Cl)",
    paste0("\n",paste0(names(fixed_PD_param),"=",as.numeric(fixed_PD_param),collapse = "\n"),"\n")
    ,"AUC_0=0
R_0 = 0
ddt_AUC=R
ddt_R = Emax*Cc^H/(EC50^H+Cc^H)*(1-Imax*AUC/(IC50/K^(",n_admin,"-1)+AUC))-kdeg*R
    
OUTPUT:
output = {Cc, R}")
  return(PKPD_model)
}




create_mlxtran_Chen = function(path_data,path_result,path_model,
                               names_param,fixed_param,param_var,
                               init,
                               n_chains,n_smooth,n_explo){
  
  names_no_fixed=names_param[-which(names_param %in% names(fixed_param))]
  
  var_definition=rep("no-variability",length(names_no_fixed))
  var_definition[which(names_no_fixed %in% param_var)]=paste0("sd=omega_",names_no_fixed[which(names_no_fixed %in% param_var)])
  
  init_index=unlist(sapply(
    1:length(names_no_fixed),
    function(x) grep(names_no_fixed[x],names(init))
  ))
  
  mon_mlxtran = paste0(
    "; this script is generated automatically
    
    <DATAFILE>
    
    [FILEINFO]
    file = ",path_data,"
    delimiter = tab
    header = {ID, TIME, DV, DVID, AMT, TINF}
    
    [CONTENT]
    ID = {use=identifier}
    TIME = {use=time}
    DV = {use=observation, name={y_1, y_2}, yname={'1', '2'}, type={continuous, continuous}}
    DVID = {use=observationtype}
    AMT = {use=amount}
    TINF = {use=infusiontime}
    
    <MODEL>
    
    [INDIVIDUAL]
    input = {",paste0(c(
      paste0(names_no_fixed,"_pop"),
      paste0("omega_",param_var))
      ,collapse=", ")
    ,"}
    
    DEFINITION:
    ",paste0(paste0(names_no_fixed," ={distribution=logNormal, typical=",paste0(names_no_fixed,"_pop"),
                    ", ",var_definition,"}"),collapse="\n \t"),"
    
    [LONGITUDINAL]
    input = {b1, b2}
    
    file = ",path_model,"
    
    DEFINITION:
    y1 = {distribution=normal, prediction=Cc, errorModel=proportional(b1)}
    y2 = {distribution=normal, prediction=R, errorModel=proportional(b2)}
    
    <FIT>
    data = {y_1, y_2}
    model = {y1, y2}
    
    <PARAMETER>
    ",paste0(c(
      paste0(names_no_fixed,"_pop"," = {value = ",as.character(init[init_index]),", method=MLE}"),
      paste0("omega_",names_no_fixed[which(names_no_fixed %in% param_var)]," = {value = 1, method=MLE}")
    ),collapse="\n \t"),"
    b1 = {value=0.3, method=MLE}
    b2 = {value=0.3, method=MLE}

    <MONOLIX>
    
    [TASKS]
    populationParameters()
    individualParameters(method = conditionalMode)
    fim(method = StochasticApproximation)
    logLikelihood(run = false,method = ImportanceSampling)
    
    
    [SETTINGS]
    GLOBAL:
    exportpath = ",path_result,"
    nbchains = ",n_chains,"
    
    POPULATION:
    smoothingiterations = ",n_smooth,"
    exploratoryiterations = ",n_explo,"
    
    ")
  return(mon_mlxtran)
}


######Simu one trial
data_simu=function(id,id_seq,doses,t,t_doses,t_inf,
                   PK_model,PK_param_pop,
                   PD_model,PD_param_pop,
                   b_PK,b_PD,PD_initial_values,method_ode,
                   threshold,w_alpha,
                   pTox_oDLT_CRS_ind,pTox_oDLT_noCRS_ind){
  
  data_PKPD=NULL
  data_tox=NULL
  
  #Loop on the number of patients
  for(i in 1:nrow(doses)){

    data_true=PK_PD_simu(t=t,doses=doses[i,],
                         t_doses=t_doses,
                         t_inf=t_inf[i,],
                         PK_model=PK_model,
                         PK_param_pop=PK_param_pop,
                         PD_model=PD_model,
                         PD_param_pop=PD_param_pop,
                         PD_initial_values=PD_initial_values,
                         method_ode=method_ode)
    
    Rmax_true=extract_max_local(data_true$R,t,doses=rep(1,length(t_doses)),t_doses)
    Rmax_true[doses[i,]==0]=NA
    
    alpha=exp(rnorm(n = 1, sd=w_alpha))
    
    CRS=((alpha*Rmax_true)>=threshold)*1
    
    pTox_oDLT=array(0, dim = c(nrow(pTox_oDLT_CRS_ind),ncol(pTox_oDLT_CRS_ind),2))
    pTox_oDLT[,,1]=pTox_oDLT_noCRS_ind
    pTox_oDLT[,,2]=pTox_oDLT_CRS_ind
    
    oDLT=rbinom(ncol(pTox_oDLT),size=1,prob=pTox_oDLT[i,,(sum(CRS)>0)*1+1])
    
    tox=(CRS+oDLT>0)*1
    #######
    
    index_tox=ifelse(sum(tox,na.rm=T)>0,min(which(tox > 0)),length(Rmax_true))
    
    end_t=ifelse(index_tox<length(t_doses),t_doses[index_tox+1]-1,t[length(t)])
    
    data_patient_PK=data.frame(
      ID=rep(id[i],length(t[t<=end_t])),
      TIME=t[t<=end_t],
      DV=data_true$C[data_true$t<=end_t]*(1+b_PK*rnorm(length(t[t<=end_t]),0,1)),
      DVID=rep(1,length(t[t<=end_t])),
      AMT=NA,
      TINF=NA)
    
    data_patient_PD=data.frame(
      ID=rep(id[i],length(t[t<=end_t])),
      TIME=t[t<=end_t],
      DV=data_true$R[data_true$t<=end_t]*(1+b_PD*rnorm(length(t[t<=end_t]),0,1)),
      DVID=rep(2,length(t[t<=end_t])),
      AMT=NA,
      TINF=NA)
    
    data_patient_admin=data.frame(
      ID=rep(id[i],length(t_doses[t_doses<=end_t])),
      TIME=t_doses[t_doses<=end_t],
      DV=NA,
      DVID=rep(1,length(t_doses[t_doses<=end_t])),
      AMT=doses[i,1:sum(t_doses<=end_t)],
      TINF=t_inf[i,1:sum(t_doses<=end_t)]
    )
    
    data_patient=rbind(data_patient_admin,data_patient_PK,data_patient_PD)
    data_patient=data_patient[
      with(data_patient, order(TIME, DVID)),
      ]
    
    data_PKPD=rbind(data_PKPD,data_patient)
    
    data_tox=rbind(
      data_tox,
      data.frame(
        id=rep(id[i],index_tox),
        admin=1:index_tox,
        tox=tox[1:index_tox],
        CRS=CRS[1:index_tox],
        oDLT=oDLT[1:index_tox],
        dose=doses[i,1:index_tox],
        id_seq=rep(id_seq[i],index_tox)
      )
      
    )
  }
  
  return(list(data_PKPD=data_PKPD,data_tox=data_tox))
}

crm_empiric_trial=function(itrial,DL,t,t_doses,t_inf,
                           PK_model,PK_param_pop,
                           PD_model,PD_param_pop,
                           b_PK,b_PD,PD_initial_values,method_ode,
                           threshold,w_alpha,
                           pTox_oDLT_CRS_cond,
                           pTox_oDLT_noCRS_cond,
                           wm,
                           target,
                           coh_size,n_tot){
  
  set.seed(itrial)
  
  n_doses=nrow(DL)
  current_dose=1
  n_patients=0
  
  dose_donnees=c()
  y=c()
  
  data=NULL  
  data_PKPD=NULL
  
  
  while(n_patients+coh_size <= n_tot && current_dose > 0){
    
    if(is.null(data_PKPD)){
      id=1:coh_size
    } else{
      id=data_PKPD$ID[length(data_PKPD$ID)]+1:coh_size
    }
    
    simu_while=data_simu(id=id,id_seq=rep(current_dose,coh_size),
                         doses=DL[rep(current_dose,coh_size),],t=t,t_doses=t_doses,
                         t_inf=t_inf[rep(current_dose,coh_size),],
                         PK_model=PK_model,PK_param_pop=PK_param_pop,
                         PD_model=PD_model,PD_param_pop=PD_param_pop,
                         b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                         threshold=threshold,w_alpha=w_alpha,
                         pTox_oDLT_CRS_ind=pTox_oDLT_CRS_cond[rep(current_dose,coh_size),],
                         pTox_oDLT_noCRS_ind=pTox_oDLT_noCRS_cond[rep(current_dose,coh_size),])
    
    
    data_while=simu_while$data_PKPD
    
    data_tox=simu_while$data_tox
    
    tox_while_tot=tapply(data_tox$tox,data_tox$id,max)
    
    data_PKPD=rbind(data_PKPD,data_while)
    data=rbind(data,data_tox)
    
    n_patients=n_patients+coh_size
    dose_donnees=c(dose_donnees,rep(current_dose,coh_size))
    y=c(y,as.numeric(tox_while_tot))
    
    calcul_next_dose=crm(prior=wm,target=target,tox=y,level=dose_donnees,model="empiric")
    
    current_dose = min(calcul_next_dose$mtd, current_dose+1)
    
    pTox=calcul_next_dose$ptox
    
    
  }
  
  return(list(data_PKPD=data_PKPD,data=data,pTox=pTox))
  
}

crm_2p_trial=function(itrial,DL,t,t_doses,t_inf,
                      PK_model,PK_param_pop,
                      PD_model,PD_param_pop,
                      b_PK,b_PD,PD_initial_values,method_ode,
                      threshold,w_alpha,
                      pTox_oDLT_CRS_cond,
                      pTox_oDLT_noCRS_cond,
                      model,wm,
                      target,
                      coh_size,n_tot){
  
  set.seed(itrial)
  
  n_doses=nrow(DL)
  current_dose=1
  n_patients=0
  
  dose_donnees=c()
  y=c()
  
  data=NULL  
  data_PKPD=NULL
  
  
  while(n_patients+coh_size <= n_tot && current_dose > 0){
    
    if(is.null(data_PKPD)){
      id=1:coh_size
    } else{
      id=data_PKPD$ID[length(data_PKPD$ID)]+1:coh_size
    }
    
    simu_while=data_simu(id=id,id_seq=rep(current_dose,coh_size),
                         doses=DL[rep(current_dose,coh_size),],t=t,t_doses=t_doses,
                         t_inf=t_inf[rep(current_dose,coh_size),],
                         PK_model=PK_model,PK_param_pop=PK_param_pop,
                         PD_model=PD_model,PD_param_pop=PD_param_pop,
                         b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                         threshold=threshold,w_alpha=w_alpha,
                         pTox_oDLT_CRS_ind=pTox_oDLT_CRS_cond[rep(current_dose,coh_size),],
                         pTox_oDLT_noCRS_ind=pTox_oDLT_noCRS_cond[rep(current_dose,coh_size),])
    
    
    data_while=simu_while$data_PKPD
    
    data_tox=simu_while$data_tox

    tox_while_tot=tapply(data_tox$tox,data_tox$id,max)
    
    data_PKPD=rbind(data_PKPD,data_while)
    data=rbind(data,data_tox)
    
    n_patients=n_patients+coh_size
    dose_donnees=c(dose_donnees,rep(current_dose,coh_size))
    y=c(y,as.numeric(tox_while_tot))
    
    calcul_next_dose = calcul_next_dose_pTox(model=model, dose_donnees=dose_donnees, wm=wm, y=y, 
                                             interv=F, target=target,target_min=0.25,target_max=0.35, 
                                             dose_cur=current_dose, 
                                             saut=1,c_over=c_over, c_stop=c_stop,n_burn=2000, n_iter=4000, n_thin=1)
    
    current_dose = calcul_next_dose[[1]]
    
    pTox=calcul_next_dose[[2]]
    
    
  }
  
  return(list(data_PKPD=data_PKPD,data=data,pTox=pTox))
  
}


simu_one_trial_PKPD=function(itrial,DL,t,t_doses,t_inf,
                             PK_model,PK_param_pop,
                             PD_model,PD_param_pop,
                             fixed_param,param_var,init,
                             b_PK,b_PD,PD_initial_values,method_ode,
                             n_chains,n_smooth,n_explo,
                             create_model,
                             create_mlxtran,
                             threshold,w_alpha,
                             pTox_oDLT_CRS_cond,
                             pTox_oDLT_noCRS_cond,
                             design,
                             wm,
                             target,
                             coh_size,n_tot,
                             dir,
                             seed){
  
  set.seed(seed+itrial)
  
  dir_trial=paste0(dir,"/trial",itrial)
  
  dir.create(dir_trial)
  
  setwd(dir_trial)
  
  if(design=="empiric_crm"){
    
    trial_data=crm_empiric_trial(itrial=itrial,DL=DL,t=t,t_doses=t_doses,t_inf=t_inf,
                                 PK_model=PK_model,PK_param_pop=PK_param_pop,
                                 PD_model=PD_model,PD_param_pop=PD_param_pop,
                                 b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                                 threshold=threshold,w_alpha=w_alpha,
                                 pTox_oDLT_CRS_cond=pTox_oDLT_CRS_cond,
                                 pTox_oDLT_noCRS_cond=pTox_oDLT_noCRS_cond,
                                 wm=wm,
                                 target=target,
                                 coh_size=coh_size,n_tot=n_tot)
    
    data_PKPD=trial_data$data_PKPD
    data=trial_data$data
    pTox=trial_data$pTox
    
    write.table(pTox,
                file = "pTox.txt",sep="",row.names = FALSE, col.names = FALSE)
    
  }
  
  if(design=="crm"){
    trial_data=crm_2p_trial(itrial=itrial,DL=DL,t=t,t_doses=t_doses,t_inf=t_inf,
                            PK_model=PK_model,PK_param_pop=PK_param_pop,
                            PD_model=PD_model,PD_param_pop=PD_param_pop,
                            b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                            threshold=threshold,w_alpha=w_alpha,
                            pTox_oDLT_CRS_cond=pTox_oDLT_CRS_cond,
                            pTox_oDLT_noCRS_cond=pTox_oDLT_noCRS_cond,
                            model=logistic_crm,wm=wm,
                            target=target,
                            coh_size=coh_size,n_tot=n_tot)
    
    data_PKPD=trial_data$data_PKPD
    data=trial_data$data
    pTox=trial_data$pTox
    
    write.table(pTox,
                file = "pTox.txt",sep="",row.names = FALSE, col.names = FALSE)
    
  }
  
  
  
  write.table(data,
              file = "data.txt", 
              sep = "\t",
              row.names = FALSE, col.names = TRUE,na=".",quote=F)
  
  write.table(data_PKPD,
              file = "data_PKPD.txt", 
              sep = "\t",
              row.names = FALSE, col.names = TRUE,na=".",quote=F)
  
  
  model=create_model(names_PD_param=names(PD_param_pop$mu),fixed_PD_param=fixed_param,n_admin=length(t_doses))
  
  model_file="model.txt"
  
  sink(model_file)
  cat(model)
  sink()
  
  
  dir.create("results")
  
  mon_mlxtran = create_mlxtran(path_data=paste0("'data_PKPD.txt'"),
                               path_result=paste0("'results'"),
                               path_model=paste0("'model.txt'"),
                               names_param=c(names(PK_param_pop$mu),names(PD_param_pop$mu)),
                               fixed_param=fixed_param,
                               param_var=param_var,
                               init=init,
                               n_chains=n_chains,
                               n_smooth=n_smooth,
                               n_explo=n_explo)
  
  project="model_mlxtran.mlxtran"
  
  sink(project)
  cat(mon_mlxtran)
  sink()
  
  library(lixoftConnectors)
  initializeLixoftConnectors(software="monolix")
  
  lixoftConnectors::loadProject(projectFile = project)
  
  lixoftConnectors::runPopulationParameterEstimation()
  
  
}




simu_trials_PKPD=function(n_trials,DL,t,t_doses,t_inf,
                          PK_model,PK_param_pop,
                          PD_model,PD_param_pop,
                          fixed_param,param_var,init,
                          b_PK,b_PD,PD_initial_values,method_ode,
                          n_chains,n_smooth,n_explo,
                          create_model,
                          create_mlxtran,
                          threshold,w_alpha,
                          pTox_oDLT_CRS_cond,
                          pTox_oDLT_noCRS_cond,
                          design,
                          wm,
                          target,
                          coh_size,n_tot,
                          dir,
                          seed,
                          start=0){
  
  
  for(i in 1:n_trials){

    simu_one_trial_PKPD(itrial=i+start,DL=DL,t=t,t_doses=t_doses,t_inf=t_inf,
                        PK_model=PK_model,PK_param_pop=PK_param_pop,
                        PD_model=PD_model,PD_param_pop=PD_param_pop,
                        fixed_param=fixed_param,param_var=param_var,init=init,
                        b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                        n_chains=n_chains,n_smooth=n_smooth,n_explo=n_explo,
                        create_model=create_model,
                        create_mlxtran=create_mlxtran,
                        threshold=threshold,w_alpha=w_alpha,
                        pTox_oDLT_CRS_cond=pTox_oDLT_CRS_cond,
                        pTox_oDLT_noCRS_cond=pTox_oDLT_noCRS_cond,
                        design=design,
                        wm=wm,
                        target=target,
                        coh_size=coh_size,n_tot=n_tot,
                        dir=dir,seed=seed)
    
  }
  
}


simu_trials_PKPD_parallel=function(n_trials,DL,t,t_doses,t_inf,
                                   PK_model,PK_param_pop,
                                   PD_model,PD_param_pop,
                                   fixed_param,param_var,init,
                                   b_PK,b_PD,PD_initial_values,method_ode,
                                   n_chains,n_smooth,n_explo,
                                   create_model,
                                   create_mlxtran,
                                   threshold,w_alpha,
                                   pTox_oDLT_CRS_cond,
                                   pTox_oDLT_noCRS_cond,
                                   design,
                                   wm,
                                   target,
                                   coh_size,n_tot,
                                   dir,
                                   seed,
                                   no_cores,
                                   N_sim,
                                   start=0){
  
  N_pack=n_trials/N_sim
  
  cl=makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(ip=1:N_pack,.packages=c("deSolve","rstan","MASS","dfcrm","lixoftConnectors"),.export=ls(globalenv())) %dopar%{
    
    for(i in 1:N_sim){

      simu_one_trial_PKPD(itrial=(ip-1)*N_sim+i+start,DL=DL,t=t,t_doses=t_doses,t_inf=t_inf,
                          PK_model=PK_model,PK_param_pop=PK_param_pop,
                          PD_model=PD_model,PD_param_pop=PD_param_pop,
                          fixed_param=fixed_param,param_var=param_var,init=init,
                          b_PK=b_PK,b_PD=b_PD,PD_initial_values=PD_initial_values,method_ode=method_ode,
                          n_chains=n_chains,n_smooth=n_smooth,n_explo=n_explo,
                          create_model=create_model,
                          create_mlxtran=create_mlxtran,
                          threshold=threshold,w_alpha=w_alpha,
                          pTox_oDLT_CRS_cond=pTox_oDLT_CRS_cond,
                          pTox_oDLT_noCRS_cond=pTox_oDLT_noCRS_cond,
                          design=design,
                          wm=wm,
                          target=target,
                          coh_size=coh_size,n_tot=n_tot,
                          dir=dir,seed=seed)
    }
    
  }
  
  stopCluster(cl)
  
}




