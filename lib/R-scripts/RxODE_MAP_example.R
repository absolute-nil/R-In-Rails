rm(list=ls())
#set working directory to current folder
curr.dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)

# Globally used objects and functions
# ------------------------------------------------------------------------------
# Load package libraries
library(RxODE)	#Differential equation solver
library(tidyverse)

# prior information from Dr.Amanda's previous work
# model file: M1_1bcomp_TPD18_mk_block_tol7.txt
POPpar <- c(8.21, 15.9, 0.673)
# Omega matrix (IIV)
OMEGA <- matrix(c(0.0466,0.0334,0.0301,
                      0.0334,0.0653,-0.0472,
                      0.0301,-0.0472,0.134),3,3)
# Sigma matrix (RUV)
SIGMA <- 0.229

#-----------------------------Define Functions--------------------------# ####
# 1.Cefalexin model containing differential equations 
#   for amount in each compartment (individual simulation)

mod <- RxODE({
  CL  = indCL; 
  V   = indV;
  MAT = indMAT;
  
  # Calculate rate-constants for differential equation solver
  KT  = 3/MAT;      # transit rate constant
  K   = CL/V;       # elimination rate constant
  
  #Transit compartment absorption model
  d/dt(depot) =-KT*depot;          # GUT depot compartment
  d/dt(T1)    = KT*(depot-T1);     # GUT transit compartment 1 amount
  d/dt(T2)    = KT*(T1-T2)   ;     # GUT transit compartment 2 amount
  d/dt(centr) = KT/V*T2 - K*centr; # central concentration
  d/dt(auc)   = centr;
})

# 2.Define objective function for MAP estimator
# Literature: [Wright, 2013; eq(1)]
obj.MAP <- function(log_ind_par,data,dose_data,weight,TV_par,var_omega,var_sigma) {
    
    # Define individual parameter values
    theta <- c(indCL  = exp(log_ind_par[1])*(weight/25)^0.75,      #Clearance (L/h)	
               indV   = exp(log_ind_par[2])*(weight/25),      #Volume (L) 
               indMAT = exp(log_ind_par[3]))     #Mean absorption time h
    
    #### Define event record
    ev <- et(time = dose_data$time,
             amt  = dose_data$dose,
             evid = 1, # dosing records
             cmt = "depot") %>%
      add.sampling(data$time)
    
    # Apply RxODE for simulation
    sim_data <- rxSolve(mod,theta,ev)
    
    # Individual prediction
    ipred <- sim_data$centr
    
    # log transformed for proportional residual error and exponential IIV
    dlog_par <- log_ind_par-log(TV_par)
    with(data, sum((log(ipred) - log(conc))^2)/var_sigma + as.numeric(t(dlog_par) %*% solve(var_omega) %*% dlog_par))
}

# 4. MAP estimator
MAP.cefa<- function(data, dose_data, weight) {
    # optimization with simplex method
    res <- optim(par=log(POPpar), fn=obj.MAP, data=data, dose_data=dose_data,
                 weight=weight, TV_par=POPpar, var_omega=OMEGA, 
                 var_sigma=SIGMA)
    
    # derive the individual PK parameters  
    log_ind_par <- res$par 
    CL  <- exp(log_ind_par[1])*(weight/25)^0.75
    V   <- exp(log_ind_par[2])*(weight/25)
    MAT <- exp(log_ind_par[3])
    
    # Collect the individual parameter values
    map_est <- c(indCL=CL, indV=V, indMAT=MAT)
    return(map_est)
}


#----------------------------MAP estimation-------------------------# ####

## user input
# weight 
weight<-10
# data information
data <- data.frame(time = c(22,27), conc = c(3,20))
# dose information
dose_data <- data.frame(time=c(0,12,24,36), dose=c(500,500,500,500))

## MAP estimation  
par_ind <- MAP.cefa(data, dose_data, weight)

#-------------Individual concentration profile simulation------------# ####

#### Define event record
ev <- et(time = dose_data$time,
         amt  = dose_data$dose,
         evid = 1, # dosing records
         cmt = "depot") %>%
  add.sampling(seq(from = 0, to = 48, by = 0.1))

# Apply RxODE for simulation
sim <- rxSolve(mod,params=par_ind,events=ev)

# calculate AUC from 24h to 48h  
AUC24_48 <- round(sim[sim$time==48,]$auc - sim[sim$time==24,]$auc,0)

## Output to screen
Y <- paste0("AUC from 24h to 48h: ", AUC24_48, " h*mg/L")
print(Y)