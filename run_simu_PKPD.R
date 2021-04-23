path_simus=""
path_simus="C:/Users/Emma/Dropbox/These/Redaction/joint_model_paper/Review1/code"

design="crm"
#design="empiric_crm"

#tot_index_scenario=1:6
#####Test
tot_index_scenario=1

# n_trials_tot=1000
# start=0
# n_trials=n_trials_tot-start
# no_cores=9
# N_sim=2
#####Test
n_trials_tot=10
start=0
n_trials=n_trials_tot-start
no_cores=5
N_sim=1

############Run

source(paste0(path_simus,"/tools/simu_trials_PKPD.R"))
source(paste0(path_simus,"/PK_PD_param.R"))

dir=path_simus

# n_chains=3
# n_smooth=200
# n_explo=500
##Test
n_chains=1
n_smooth=20
n_explo=50

target=0.3
wm=getprior(0.05,target,4,6) #According to sc1
coh_size=3
n_tot=30

seed=25011995

create_mlxtran=create_mlxtran_Chen
create_model=create_model_Chen

load(paste0(dir,"/scenario/scenario.Rdata"))

for(index_scenario in tot_index_scenario){
  
  dir.create(paste0(dir,"/",design,"_trials_sc",index_scenario))
  
  save(design,index_scenario,target,seed,wm, 
       file = paste0(dir,"/",design,"_trials_sc",index_scenario,"/trial_param.RData"))
  
  simu_trials_PKPD_parallel(n_trials,
                            DL=DL[scenario[[index_scenario]]$id_DL,],
                            t,t_doses,
                            t_inf=t_inf[scenario[[index_scenario]]$id_DL,],
                            PK_model,PK_param_pop,
                            PD_model,PD_param_pop,
                            fixed_param,param_var,init,
                            b_PK,b_PD,PD_initial_values,method_ode,
                            n_chains,n_smooth,n_explo,
                            create_model,
                            create_mlxtran,
                            threshold=scenario[[index_scenario]]$threshold,
                            w_alpha=scenario[[index_scenario]]$w_alpha,
                            pTox_oDLT_CRS_cond=scenario[[index_scenario]]$pTox_oDLT_CRS_cond,
                            pTox_oDLT_noCRS_cond=scenario[[index_scenario]]$pTox_oDLT_noCRS_cond,
                            design,
                            wm,
                            target,
                            coh_size,n_tot,
                            dir=paste0(dir,"/",design,"_trials_sc",index_scenario),
                            seed,
                            no_cores,
                            N_sim,
                            start)
}

