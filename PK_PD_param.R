#########################################################################################################
######################################Parameters of the PK/PD model######################################
#########################################################################################################

source(paste0(path_simus,"/tools/pk_functions.R"))

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

DL=DL_kg*70

t_inf=matrix(4,nrow=nrow(DL),ncol=ncol(DL))

t_doses_days=c(1,5,9,13,17,21,25)-1
t_doses=t_doses_days*24 ##In hours

PK_model=c_infusion_1cpt_multiple

mu_PK=c(Cl=1.36,V=3.40)
w_PK=matrix(0,length(mu_PK),length(mu_PK))
diag(w_PK)=c(0.419^2,0)
PK_param_pop=list(mu=mu_PK,w=w_PK)

PD_model=PD_model_chen

mu_PD=c(Emax=3.59*10^5,
        EC50=1*10^4,
        H=9.20*10^-1,
        Imax=0.995,
        IC50=1.82*10^4,
        kdeg=1.80*10^-1,
        K=2.83)


w_PD=matrix(0,length(mu_PD),length(mu_PD)) 
diag(w_PD)=c(
  0.14^2,
  0,
  0.03^2,
  0,
  0.12^2,
  0.13^2,
  0.36^2
)

PD_param_pop=list(mu=mu_PD,
                  w=w_PD)

PD_initial_values=PD_initial_values_chen

t=unique(sort(unlist(sapply(t_doses,function(x) c(x,x+1,x+2,x+5,x+7,x+24,x+48,x+72) ))))

method_ode="vode"


b_PK=0.1
b_PD=0.1


###PK/PD estimation with monolix
fixed_param=mu_PD[names(mu_PD) %in% c("EC50","Imax","IC50") ]
param_var=c("Cl","Emax","kdeg","K")
init=c(Cl="1",Emax="100000",H="1",K="2",V="4",kdeg="0.5")

