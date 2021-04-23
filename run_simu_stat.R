path_simus=""
path_simus="C:/Users/Emma/Dropbox/These/Redaction/joint_model_paper/Review1/code"

sc=1
design="crm"
seq_ref=4

dir_trials=paste0(path_simus,"/",design,"_trials_sc",sc)

# N_sim=5
# no_cores=20
######Test
N_sim=1
no_cores=5

load(paste0(dir_trials,"/trial_param.RData"))
save_output=paste0(design,"_sc",index_scenario)
name_temp="temp"


#####

source(paste0(path_simus,"/tools/simu_trials_stat.R"))
source(paste0(path_simus,"/PK_PD_param.R"))

load(paste0(path_simus,"/scenario/scenario.Rdata"))

#n_trials=1000
######Test
n_trials=5

start=0

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
  1,5,10,30,45,rep(60,2),
  5,10,rep(30,5),
  1,5,10,30,rep(60,3),
  1,10,rep(45,5)
),ncol=7,byrow=T)

DL_predict=DL_predict*70

t_inf_predict=matrix(4,nrow=nrow(DL_predict),ncol=ncol(DL_predict))

names_PK=names(mu_PK)
names_PD=names(mu_PD)
DL=DL[scenario[[sc]]$id_DL,]
t_inf=t_inf[scenario[[sc]]$id_DL,]

estim_trials(
  n_trials,
  sc,
  design,
  names_PK,
  names_PD,
  fixed_param,
  DL,
  t_doses,
  t_inf,
  DL_predict,
  t_inf_predict,
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
  start,
  name_temp
)

