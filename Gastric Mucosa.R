#####################################################################################
##########             Modelo CE progresión mucosa       #####################
#####################################################################################

# Based on the model developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)	
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Drug Policy Program, Center for Research and Teaching in Economics (CIDE) - CONACyT, 
#   Aguascalientes, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property 
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without written permission.
#####################################################################################

rm(list = ls())      # clear memory (removes all the variables from the workspace)

#### 01 Load packages ####
library(ggplot2)
devtools::install_github("DARTH-git/dampack")
library(dampack)

#### 02 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Prevention", "Prevention")  

## Number of strategies
n.str <- length(v.names.str)

## Markov model parameters
age     <- 20                                 # age at baseline
max.age <- 100                                # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles

v.n  <- c("Normalhpn", "Normalhpp", 
          "Gastritishpn", "Gastritishpp", 
          "Atrophyhpn", "Atrophyhpp", 
          "Intestinalhpn", "Intestinalhpp", 
          "Dysplasiahpn", "Dysplasiahpp", # Dyplasia state could be merged (hp - y hp +)
          "Dead")    # state names
n.s  <- length(v.n)                     # number of states

p.NnD <- 0.02  # probability to die when Normal mucosa with Helicobacter pilory negative

#Probability of transition without prevention
p.NnGn <- 0.05  # probability to become Gastritis hp(-) when Normal hp(-)
p.GnAn <- 0.05  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.AnIn <- 0.05  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.InDn <- 0.05  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.DnIn <- 0.05 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.InAn <- 0.05  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.AnGn <- 0.05  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.GnNn <- 0.05  # probability to become Normal hp(-) when Gastritis hp(-)
p.GnGp <- 0.05 # probability to become Gastritis hp(-) when Gastritis hp(+)
p.GpGn <- 0.05 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.AnAp <- 0.05 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.ApAn <- 0.05 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.InIp <- 0.05 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.IpIn <- 0.05 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.DnDp <- 0.05 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.DpDn <- 0.05 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.NpGp <- 0.05  # probability to become Gastritis hp(+) when Normal hp(+)
p.NnNp <- 0.05  # probability to become Normal hp(+) when Normal hp(-)
p.NpNn <- 0.05  # probability to become Normal hp(-) when Normal hp(+)
p.NpGp <- 0.05  # probability to become Gastritis hp(+) when Normanl hp(+)
p.GpAp <- 0.05  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.ApIp <- 0.05  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.IpDp <- 0.05  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.DpIp <- 0.05  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.IpAp <- 0.05  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.ApGp <- 0.05  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.GpNp <- 0.05  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 1 
t.effect.NnGn <- t.effect.GnAn <- 0.7 # This value should reflect treatment effectiveness in increase transition probability
t.ehpylori <- 0.8 # This value reflect probability of erradication after treatment
t.dhpylori <- 0.95 # This value reflect the probability of diagnostic strategy 1 to detect H. pylori cases (test sensitivity)
t.dehpylori <- 0.05 # This value reflect the probability of diagnostic strategy 1 to incorrectly clasify H. pylori - cases as + (test 1-specificity)

p.pNnGn <- p.NnGn * t.effect.NnGn  # probability to become Gastritis hp(-) when Normanl hp(-)
p.pGnAn <- p.GnAn * t.effect.GnAn  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.pAnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.pInDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.pDnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.pInAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.pAnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.pGnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.pGnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.pGpGn <- t.ehpylori # probability to become Gastritis hp(+) when Gastritis hp(-)
p.pAnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.pApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.pInIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.pIpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.pDnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.pDnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.pDpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.pDpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.pNpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.pNnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.pNpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.pNpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.pGpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.pApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.pIpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.pDpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.pIpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.pApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.pGpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

# Death rates (all equal. Same death riisk s assumed across states)
r.NnD    <- - log(1 - p.NnD) # rate of death in healthy
r.NpD   <- r.NnD  	 # rate of death in Normal helicobacter (+)
r.GnD   <- r.NnD  	 # rate of death in Gastritis helicobacter (-)
r.GpD   <- r.NnD  	 # rate of death in Gastritis helicobacter (+)
r.AnD   <- r.NnD  	 # rate of death in Atrophy helicobacter (-)
r.ApD   <- r.NnD  	 # rate of death in Atrophy helicobacter (+)
r.InD   <- r.NnD  	 # rate of death in Intestinal helicobacter (-)
r.IpD   <- r.NnD     # rate of death in Intestinal helicobacter (+)
r.DnD   <- r.NnD  	 # rate of death in Dysplasia helicobacter (-)
r.DpD   <- r.NnD     # rate of death in Dysplasia helicobacter (+)

p.NpD <- 1- exp(-r.NpD)  # probability to die when Normal helicobacter (+)
p.GnD <- 1- exp(-r.GnD)  # probability to die when Gastritis helicobacter (-)
p.GpD <- 1- exp(-r.GpD)  # probability to die when Gastritis helicobacter (+)
p.AnD <- 1- exp(-r.AnD)  # probability to die when Atrophy helicobacter (-)
p.ApD <- 1- exp(-r.ApD)  # probability to die when Atrophy helicobacter (+)
p.InD <- 1- exp(-r.InD)  # probability to die when Intestinal helicobacter (-) 
p.IpD <- 1- exp(-r.IpD)  # probability to die when Intestinal helicobacter (+) 
p.DnD <- 1- exp(-r.DnD)  # probability to die when Dysplasia helicobacter (-) 
p.DpD <- 1- exp(-r.DpD)  # probability to die when Dysplasia helicobacter (+) 

# Costs and utilities  
c.Nn  <- 400                     # cost of remaining one cycle Normal helicobacter (-)
c.Np  <- 400                     # cost of remaining one cycle Normal helicobacter (+)
c.Gn  <- 400                     # cost of remaining one cycle Gastritis helicobacter (-)
c.Gp  <- 400                     # cost of remaining one cycle Gastritis helicobacter (+)
c.An  <- 400                     # cost of remaining one cycle Atrophy helicobacter (-)
c.Ap  <- 400                     # cost of remaining one cycle Atrophy helicobacter (+)
c.In  <- 4000                    # cost of remaining one cycle Intestinal helicobacter (-) 
c.Ip  <- 8000                    # cost of remaining one cycle Intestinal helicobacter (+) 
c.Dn  <- 12000                   # cost of remaining one cycle Dysplasia helicobacter (-) 
c.Dp  <- 15000                   # cost of remaining one cycle Dysplasia helicobacter (+)
c.D  <- 0  # cost of remaining one cycle dead
c.diag <- 200 # cost of prevention strategy 
c.trt <- 300

# Single utility or hp+ / hp- utility? (pedir a Cristian)
u.Nn  <- 1                     # utility when Normal helicobacter (-)
u.Np  <- 1                     # utility when Normal helicobacter (+)
u.Gn  <- 1                     # utility when Gastritis helicobacter (-)
u.Gp  <- 0.7                   # utility when Gastritis helicobacter (+)
u.An  <- 0.8                   # utility when Atrophy helicobacter (-)
u.Ap  <- 0.5                   # utility when Atrophy helicobacter (+)
u.In  <- 0.4                   # utility when Intestinal helicobacter (-)
u.Ip  <- 0.3                   # utility when Intestinal helicobacter (+)
u.Dn  <- 0.2                   # utility when Dysplasia helicobacter (-)
u.Dp  <- 0.2                   # utility when Dysplasia helicobacter (+) 
u.D  <- 0                      # utility when dead

# Discounting factor
d.r  <- 0.03                    # equal discount of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


#### 03 Define and initialize matrices and vectors ####
#### 03.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_prev <- m.M_prev <- matrix(NA, 
                                  nrow = n.t + 1, ncol = n.s,
                                  dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_prev) # show first 6 rows of the matrix 

# The cohort starts as healthy (modify accordingly to prevalence of each state at t0)
m.M_no_prev[1, ] <- m.M_prev[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)                     # initialize first cycle of Markov trace

#### 03.2 Transition probability MATRIX ####
# create the transition probability matrix without prevention
m.P_noprev  <- matrix(0,
                      nrow = n.s,
                      ncol = n.s,
                      dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_noprev

# create the transition probability matrix witho prevention
m.P_prev  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_noprev


# fill in the transition probability matrix without prevention
### From Normal helicobacter (-)
m.P_noprev["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.NnNp - p.NnGn
m.P_noprev["Normalhpn", "Normalhpp"]    <- p.NnNp
m.P_noprev["Normalhpn", "Gastritishpn"]    <- p.NnGn
m.P_noprev["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_noprev["Normalhpp", "Normalhpp"]    <- 1 - p.NpGp - p.NpNn - p.NpD 
m.P_noprev["Normalhpp", "Gastritishpp"]    <- p.NpGp
m.P_noprev["Normalhpp", "Normalhpn"]    <- p.NpNn
m.P_noprev["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_noprev["Gastritishpn", "Gastritishpn"]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
  m.P_noprev["Gastritishpn", "Atrophyhpn"]    <- p.GnAn
m.P_noprev["Gastritishpn", "Gastritishpp"]    <- p.GnGp
m.P_noprev["Gastritishpn", "Normalhpn"]    <- p.GnNn
m.P_noprev["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_noprev["Gastritishpp", "Gastritishpp"]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD
m.P_noprev["Gastritishpp", "Atrophyhpp"]    <- p.GpAp
m.P_noprev["Gastritishpp", "Gastritishpn"]    <- p.GpGn
m.P_noprev["Gastritishpp", "Normalhpp"]    <- p.GpNp
m.P_noprev["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_noprev["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.AnIn - p.AnAp - p.AnGn -  p.AnD
m.P_noprev["Atrophyhpn", "Intestinalhpn"]    <- p.AnIn
m.P_noprev["Atrophyhpn", "Atrophyhpp"]    <- p.AnAp
m.P_noprev["Atrophyhpn", "Gastritishpn"]    <- p.AnGn
m.P_noprev["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (+) 
m.P_noprev["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD
m.P_noprev["Atrophyhpp", "Intestinalhpp"]    <- p.ApIp
m.P_noprev["Atrophyhpp", "Atrophyhpn"]    <- p.ApAn
m.P_noprev["Atrophyhpp", "Gastritishpp"]    <- p.ApGp
m.P_noprev["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_noprev["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD 
m.P_noprev["Intestinalhpn", "Dysplasiahpn"]    <- p.InDn
m.P_noprev["Intestinalhpn", "Intestinalhpp"]    <- p.InIp
m.P_noprev["Intestinalhpn", "Atrophyhpn"]    <- p.InAn
m.P_noprev["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_noprev["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD  
m.P_noprev["Intestinalhpp", "Dysplasiahpp"]    <- p.IpDp
m.P_noprev["Intestinalhpp", "Intestinalhpn"]    <- p.IpIn
m.P_noprev["Intestinalhpp", "Atrophyhpp"]    <- p.IpAp
m.P_noprev["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_noprev["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.DnDp - p.DnIn -  p.DnD
m.P_noprev["Dysplasiahpn", "Dysplasiahpp"]    <- p.DnDp
m.P_noprev["Dysplasiahpn", "Intestinalhpn"]    <- p.DnIn
m.P_noprev["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_noprev["Dysplasiahpp", "Dysplasiahpn"]    <- 1 - p.DpDn - p.DpIp - p.DpD
m.P_noprev["Dysplasiahpp", "Dysplasiahpp"]    <- p.DpDn
m.P_noprev["Dysplasiahpp", "Intestinalhpp"]    <- p.DpIp
m.P_noprev["Dysplasiahpp", "Dead"]    <- p.DpD

### From Dead
m.P_noprev["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_noprev)

# create transition probability matrix for prevention 
m.P_prev  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev

# fill in the transition probability matrix without prevention
### From Normal helicobacter (-)
m.P_prev["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.pNnNp - p.pNnGn
m.P_prev["Normalhpn", "Normalhpp"]    <- p.pNnNp
m.P_prev["Normalhpn", "Gastritishpn"]    <- p.pNnGn
m.P_prev["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev["Normalhpp", "Normalhpp"]    <- 1 - p.pNpGp - p.pNpNn - p.NpD 
m.P_prev["Normalhpp", "Gastritishpp"]    <- p.pNpGp
m.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn
m.P_prev["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev["Gastritishpn", "Gastritishpn"]    <- 1 - p.pGnAn - p.pGnGp - p.pGnNn - p.GnD
m.P_prev["Gastritishpn", "Atrophyhpn"]    <- p.pGnAn
m.P_prev["Gastritishpn", "Gastritishpp"]    <- p.pGnGp
m.P_prev["Gastritishpn", "Normalhpn"]    <- p.pGnNn
m.P_prev["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev["Gastritishpp", "Gastritishpp"]    <- 1 - p.pGpAp - p.pGpGn - p.pGpNp - p.GpD
m.P_prev["Gastritishpp", "Atrophyhpp"]    <- p.pGpAp
m.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
m.P_prev["Gastritishpp", "Normalhpp"]    <- p.pGpNp
m.P_prev["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.pAnIn - p.pAnAp - p.pAnGn -  p.AnD
m.P_prev["Atrophyhpn", "Intestinalhpn"]    <- p.pAnIn
m.P_prev["Atrophyhpn", "Atrophyhpp"]    <- p.pAnAp
m.P_prev["Atrophyhpn", "Gastritishpn"]    <- p.pAnGn
m.P_prev["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.pApIp - p.pApAn - p.pApGp - p.ApD
m.P_prev["Atrophyhpp", "Intestinalhpp"]    <- p.pApIp
m.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
m.P_prev["Atrophyhpp", "Gastritishpp"]    <- p.pApGp
m.P_prev["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.pInDn - p.pInIp - p.pInAn - p.InD 
m.P_prev["Intestinalhpn", "Dysplasiahpn"]    <- p.pInDn
m.P_prev["Intestinalhpn", "Intestinalhpp"]    <- p.pInIp
m.P_prev["Intestinalhpn", "Atrophyhpn"]    <- p.pInAn
m.P_prev["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.pIpDp - p.pIpIn - p.pIpAp - p.IpD  
m.P_prev["Intestinalhpp", "Dysplasiahpp"]    <- p.pIpDp
m.P_prev["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
m.P_prev["Intestinalhpp", "Atrophyhpp"]    <- p.pIpAp
m.P_prev["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.pDnDp - p.pDnIn -  p.DnD
m.P_prev["Dysplasiahpn", "Dysplasiahpp"]    <- p.pDnDp
m.P_prev["Dysplasiahpn", "Intestinalhpn"]    <- p.pDnIn
m.P_prev["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.pDpDn - p.pDpIp - p.DpD
m.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
m.P_prev["Dysplasiahpp", "Intestinalhpp"]    <- p.pDpIp
m.P_prev["Dysplasiahpp", "Dead"]    <- p.DpD

### From Dead
m.P_prev["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev)

#### 04 Run Markov model ####
v.trt <- rbinom(n.t,1,0.15) # Vector of treatment requires to be defined based on times of screening (3, 5 y 10 años)

for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_prev[t + 1, ] <- t(m.M_no_prev[t, ]) %*% m.P_noprev # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev[t + 1, ]    <- t(m.M_prev[t, ]) %*% if(v.trt[t]==1){m.P_prev} else {m.P_noprev} # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_prev)  # show the first 6 lines of the matrix

#### 05 Compute and Plot Epidemiological Outcomes ####
#### 05.1 Cohort trace #####
matplot(m.M_no_prev, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

# #### 05.2 Overall Survival (OS) #####
# v.os_no_prev <- 1 - m.M_no_prev[, "Dead"]       # calculate the overall survival (OS) probability for no prevention
# 
# plot(0:n.t, v.os_no_prev, type = 'l', 
#      ylim = c(0, 1),
#      ylab = "Survival probability",
#      xlab = "Cycle",
#      main = "Overall Survival")             # create a simple plot showing the OS
# grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

#### 05.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os_no_prev)                       # summing probablity of OS over time  (i.e. life expectancy)

#### 05.3 Disease prevalence ##### Revisar (para cada estado)
v.preva <- rowSums(m.M_no_prev[, c("Normalhpn", "Normalhpp", "Gastritishpn", "Gastritishpp", "Atrophyhpn", "Atrophyhpp", "Intestinalhpn", "Intestinalhpp", "Dysplasiahpn", "Dysplasiahpp")])/v.os_no_prev
plot(v.preva,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 05.4 Proportion of sick in Dysplasia Hp(+) #####
v.prop.Dpp <- m.M_no_prev[, "Dysplasiahpp"] / v.preva
plot(0:n.t, v.prop.Dpp,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Dysplasiahpp", 
     col = "black", type = "l")

#### 06 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_prev <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)

v.c_no_prev <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.D)
v.c_prev    <- c(c.Nn + c.diag, c.Np + c.diag, c.Gn + c.diag, c.Gp + c.diag, c.An + c.diag, 
                 c.Ap + c.diag, c.In + c.diag, c.Ip + c.diag, c.Dn + c.diag, c.Dp + c.diag, c.D) 
# Pendant to add treatment cost to positive diagnostic test (c.trt * p.pos). Pendant to define p.pos more accurately

#### 06.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_prev <- m.M_no_prev %*% v.u_no_prev
v.tu_prev <- m.M_prev %*% v.u_prev

v.tc_no_prev <- m.M_no_prev %*% v.c_no_prev
v.tc_prev    <- m.M_prev    %*% v.c_prev

#### 06.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_prev <- t(v.tu_no_prev) %*% v.dwe  
tu.d_prev    <- t(v.tu_prev)    %*% v.dwe

tc.d_no_prev <- t(v.tc_no_prev) %*% v.dwc
tc.d_prev    <- t(v.tc_prev)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_prev, tc.d_prev)
v.tu.d <- c(tu.d_no_prev, tu.d_prev)

# Matrix with discounted costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.tc.d,
                   Effect   = v.tu.d)
m.ce

#### 07 Compute ICERs of FONIS CaGa model #### Revisar función calculate_icers
m.cea <- calculate_icers(cost = m.ce$Cost,
                         effect = m.ce$Effect,
                         strategies = m.ce$Strategy)
m.cea

#### 08 Plot frontier of Sick-Sicker model ####
plot(m.cea, xlim = c(15.7, 16.5))
ggsave("figs/Markov-FONISCaGast-CEA-Frontier-gastricmucosa.png", width = 8, height = 6)