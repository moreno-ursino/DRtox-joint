path_simus=""
path_simus="C:/Users/Emma/Dropbox/These/Redaction/joint_model_paper/Review1/code"

itrial=248
sc=2
design="crm"
seq_ref=4

dir_trials=paste0(path_simus,"/ex_one_trial")

load(paste0(dir_trials,"/trial_param.RData"))

#####

source(paste0(path_simus,"/tools/simu_trials_stat.R"))
source(paste0(path_simus,"/PK_PD_param.R"))

load(paste0(path_simus,"/scenario/scenario.Rdata"))

# n_chains=4
# n_warmup=2500
# n_PKPD=600
######Test
n_chains=4
n_warmup=100
n_PKPD=10

sd_beta0=2
alpha=5
sd_beta1=1

alpha_clayton=1
beta_clayton=1

DL_predict=matrix(data=c(
  5,10,rep(30,5),
  1,5,10,30,rep(60,3)
),ncol=7,byrow=T)

DL_predict=DL_predict*70

t_inf_predict=matrix(4,nrow=nrow(DL_predict),ncol=ncol(DL_predict))

names_PK=names(mu_PK)
names_PD=names(mu_PD)
DL=DL[scenario[[sc]]$id_DL,]
t_inf=t_inf[scenario[[sc]]$id_DL,]


res_one_trial=estim_one_trial(
  itrial=itrial,
  sc,
  design,
  dir=paste0(dir_trials,"/data"),
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
  sim=F)

##########

load(paste0(dir_trials,"/res_one_trial.RData"))

colMeans(res_one_trial$p_CRS)
colMeans(res_one_trial$p_oDLT)
colMeans(res_one_trial$p_oDLT_noCRS)
colMeans(res_one_trial$p_DLT_indep)
colMeans(res_one_trial$p_DLT_cond)
colMeans(res_one_trial$p_DLT_copula)

