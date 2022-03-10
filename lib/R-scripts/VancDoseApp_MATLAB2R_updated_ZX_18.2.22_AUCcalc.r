## Simulation Code for Vancomycin Paed PK
# App for Amanda Gwee
# Simulation code: S Duffull 20 November 2018
# updated with AUC calculation: Derek Zhu 15 Feb 2022
#
# Purpose of the code: to predict the dose to achieve a target trough of 15
# mg/L (or as defined by the user)

remove(list=ls())

# load R packages
library(tidyverse) # for data visualisation and manipulation
library(RxODE) # for simulation

#-------------------------------------------------------------------
# Define model
#-------------------------------------------------------------------

# 2-comp oral absorption model + Emax model links Conc to Receptor occupancy

mod <- RxODE({
  CL = CL; 
  V1 = V1;
  Q  = Q;
  V2 = V2;
  C1 = centr/V1;
  C2 = peri/V2;
  d/dt(centr) = -CL*C1 - Q*C1 + Q*C2;
  d/dt(peri)  =          Q*C1 - Q*C2;
  d/dt(auc)   = C1;
})

## User input via APP

WT <- 2.5 # (kg)
CR <- 27 # (mcmol)
PMA <- 39 #(weeks)

## Variables that the user can change in settings of APP

Target_Trough <- 15 # [default=15] (mg/L)
Infusion_Duration <- 1 #[default=1] (hours)

## simulation environment variables (cannot be changed by user)
doseperkg <- 15 # mg/kg/dose
tk0 <- Infusion_Duration #IIVI
dd <- 2 # dose day (dd=2 means 24=48 hours) this is the day that observations are taken

## calculate dose interval
if (PMA < 29){
    di <- 24
} else if (PMA >= 29 && PMA < 36){
    di <- 12
} else if (PMA >= 36 && PMA < 44){
    di <- 8
} else if (PMA >= 44){
    di <- 6
}

## Get parameters
clmax <- 0.728
v1 <- 1.08
v2 <- 0.735
q <- 0.633
tm50 <- 47.7
hill <- 3.4
cl0 <- 0.012
wt_cl <- 0.75
wt_v1 <- 1.55
cr_cl <- 0.327

## generate patient
ID_MATMOD <- (PMA^hill)/(PMA^hill+tm50^hill)
ID_TVV1 <- v1*(WT/2.525)^wt_v1
ID_CLMAX <- ((clmax*(WT/2.525)^wt_cl-cl0)*ID_MATMOD)*((CR/27)^-1)^cr_cl
ID_TVCL <- cl0+ID_CLMAX
CL <- ID_TVCL
V1 <- ID_TVV1
V2 <- v2
Q <- q

## generate dosing schedules
dose <- doseperkg * WT

## reparameterise
K12 <- Q/V1
K21 <- Q/V2
K10 <- CL/V1
zzz <- K12+K21+K10
yyy <- K21*K10
alpha <- 0.5*(zzz+sqrt(zzz*zzz-4*yyy))
beta <- 0.5*(zzz-sqrt(zzz*zzz-4*yyy))
a <- dose*(K21-alpha)/(V1*(beta-alpha))
b <- dose*(K21-beta)/(V1*(alpha-beta))

## Calculate predicted trough concentration at first and final dose
Dose <- dose
trough <- di
dn <- dd*24/di

# first-dose
Temp_var <- trough < tk0
FP11 <- a/alpha * (1 - exp(-alpha * trough))*Temp_var
FP21 <- b/beta * (1 - exp(-beta * trough))*Temp_var
FP12 <- a/alpha * (1 - exp(-alpha * tk0)) * exp(-alpha * (trough-tk0))*(1-Temp_var)
FP22 <- b/beta * (1 - exp(-beta * tk0)) * exp(-beta * (trough-tk0))*(1-Temp_var)
FP1 <- FP11 + FP12
FP2 <- FP21 + FP22
F1D <- 1/tk0 * (FP1 + FP2)
C1D_temp <- F1D
First_Dose_Trough <- C1D_temp

# 48-hour
R1 <- (1-exp(-dn*alpha*di))/(1-exp(-alpha*di))
R2 <- (1-exp(-dn*beta*di))/(1-exp(-beta*di))
Temp_var <- trough < tk0
FP11 <- a/alpha * (1 - exp(-alpha * trough))*Temp_var
FP21 <- b/beta * (1 - exp(-beta * trough))*Temp_var
FP12 <- a/alpha * (1 - exp(-alpha * tk0)) * exp(-alpha * (trough-tk0))*R1*(1-Temp_var)
FP22 <- b/beta * (1 - exp(-beta * tk0)) * exp(-beta * (trough-tk0))*R2*(1-Temp_var)
FP1 <- FP11 + FP12
FP2 <- FP21 + FP22
FSSD <- 1/tk0 * (FP1 + FP2)
CSSD_temp <- FSSD
EndofDay2_Trough <- CSSD_temp

## What dose should I use?
factorDose <- Target_Trough / EndofDay2_Trough
dose <- doseperkg * WT * factorDose

## Calculate predicted trough concentration at first and final dose
DoseInterval <- di
RecommendedDose <- round(dose, 0)
RecDosePerKgPerDay <- round(dose/WT*24/di, 0)

## Output to screen
X <- paste0("Recommended Dose: ", RecommendedDose, " mg every ", DoseInterval, " hours (", RecDosePerKgPerDay, " mg/kg/day)")
print(X)

#----------------AUC Calculation---------------------#
#### Define fixed effect parameters
theta <- c(CL=CL, V1=V1, Q=Q, V2=V2) 

#### Define event record
ev <- et(amount.units="mg", time.units="hours") %>%
  add.dosing(dosing.to="centr",
             dose = RecommendedDose,
             rate = RecommendedDose, #1h infusion
             nbr.doses = dn, 
             dosing.interval= DoseInterval,
             start.time=0) %>%
  add.sampling(seq(from=24,to=48,by=0.1))

#### Perform simulation
sim  <- rxSolve(mod,theta,ev)

AUC24_48 <- round(sim[sim$time==48,]$auc - sim[sim$time==24,]$auc,0)

## Output to screen
Y <- paste0("AUC from 24h to 48h: ", AUC24_48, " h*mg/L")
print(Y)
Y %>% return(.)