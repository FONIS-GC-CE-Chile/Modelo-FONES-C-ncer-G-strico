####################################################################################
##########             Modelo CE progresión mucosa       #####################
#####################################################################################

# Based on the model developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2,)	
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
library(dampack)

#### 02 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Prevention", "EDAUreasa", "UreaAire", "Serología", "Antígenofecal", "EDAbp", "SerologíaEDA", "UreaAirePSEDA", "AntígenohecesPSEDA" )

## Number of strategies
n.str <- length(v.names.str)

## Markov model parameters
age     <- 20                                 # age at baseline
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles

v.n  <- c("Normalhpn", "Normalhpp", "Gastritishpn", "Gastritishpp", "Atrophyhpn", "Atrophyhpp", "Intestinalhpn", "Intestinalhpp", "Dysplasiahpn", "Dysplasiahpp", "Dead")    # state names
n.s  <- length(v.n)                     # number of states

p.NnD <-  0.000616  # probability to die when healthy  ## Valor ficticio considerando 2018 muertes Chile 6,16 por 1000 hab

#Probability of transition without prevention

p.NnGn <- 0.0360531  # probability to become Gastritis hp(-) when Normal hp(-) --> Huang, 2020
p.GnAn <- 0.1386021  # probability to become Atrophy hp(-) when Gastritis hp(-) --> Huang, 2020
p.AnIn <- 0.2419695  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) --> Huang, 2020
p.InDn <- 0.1195749  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.DnIn <- 0.2573678  # probability to become Intestinal hp(-) when Dysplasia hp(-)--> Huang, 2020
p.InAn <- 0.0824981  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-) --> Huang, 2020
p.AnGn <- 0.1401598  # probability to become Gastritis hp(-) when Atrophy hp(-) --> Huang, 2020
p.GnNn <- 0.0043724  # probability to become Normal hp(-) when Gastritis hp(-) --> Huang, 2020
p.GnGp <- 0.5 # probability to become Gastritis hp(-) when Gastritis hp(+)
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
p.GpAp <- 5.2 * p.GnAn  # probability to become Atrophy hp(+) when Gastritis hp(+) --> Huang, 2020
p.ApIp <- 0.05  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.IpDp <- 0.05  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.DpIp <- 0  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.IpAp <- 0  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.ApGp <- 0.05  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.GpNp <- 0.05  # probability to become Normal hp (+) when Gastritis hp(+)


#Probability of transition with prevention strategy 1 (EDAUreasa)

p.pNnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.pGnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.pAnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.pInDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.pDnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.pInAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.pAnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.pGnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.pGnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.pGpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
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

#Probability of transition with prevention strategy 2 (UreaAire)

p.p2NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p2GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p2AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p2InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p2DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p2InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p2AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p2GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p2GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p2GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p2AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p2ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p2InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p2IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p2DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p2DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p2DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p2DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p2NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p2NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p2NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p2NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p2GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p2ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p2IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p2DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p2IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p2ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p2GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 3 (Serología)

p.p3NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p3GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p3AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p3InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p3DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p3InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p3AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p3GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p3GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p3GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p3AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p3ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p3InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p3IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p3DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p3DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p3DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p3DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p3NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p3NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p3NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p3NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p3GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p3ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p3IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p3DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p3IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p3ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p3GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 4 (Antígenofecal)

p.p4NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p4GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p4AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p4InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p4DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p4InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p4AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p4GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p4GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p4GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p4AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p4ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p4InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p4IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p4DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p4DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p4DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p4DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p4NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p4NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p4NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p4NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p4GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p4ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p4IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p4DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p4IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p4ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p4GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 5 (EDAbp)

p.p5NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p5GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p5AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p5InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p5DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p5InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p5AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p5GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p5GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p5GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p5AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p5ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p5InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p5IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p5DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p5DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p5DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p5DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p5NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p5NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p5NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p5NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p5GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p5ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p5IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p5DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p5IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p5ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p5GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 6 (SerologíaEDA)
p.p6NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p6GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p6AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p6InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p6DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p6InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p6AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p6GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p6GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p6GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p6AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p6ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p6InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p6IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p6DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p6DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p6DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p6DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p6NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p6NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p6NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p6NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p6GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p6ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p6IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p6DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p6IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p6ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p6GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 7 (UreaAirePSEDA) 

p.p7NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p7GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p7AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p7InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p7DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p7InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p7AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p7GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p7GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p7GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p7AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p7ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p7InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p7IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p7DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p7DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p7DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p7DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p7NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p7NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p7NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p7NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p7GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p7ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p7IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p7DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p7IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p7ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p7GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)

#Probability of transition with prevention strategy 8 (AntígenohecesPSEDA)

p.p8NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p8GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p8AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p8InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p8DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p8InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p8AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p8GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p8GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p8GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p8AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p8ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p8InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p8IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p8DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p8DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p8DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p8DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p8NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p8NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p8NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p8NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p8GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p8ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p8IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p8DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p8IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p8ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p8GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)


#Hazard ratio 

hr.Np <- 1 #hazard ratio of death in Normal helicobacter(+) vs healthy
hr.Gn <- 1 #hazard ratio of death in Gastritis helicobacter(-) vs healthy
hr.Gp <- 1 #hazard ratio of death in Gastritis helicobacter(+) vs healthy
hr.An <- 1 #hazard ratio of death in Atrophy helicobacter (-) vs healthy
hr.Ap <- 1 #hazard ratio of death in Atrophy helicobacter (+) vs healthy
hr.In <- 1 #hazard ratio of death in Intestional helicobacter (-) vs healthy
hr.Ip <- 1 #hazard ratio of death in Intestinal helicobacter (+) vs healthy
hr.Dn <- 1 #hazard ratio of death in Dysplasia helicobacter (-) vs healthy
hr.Dp <- 1 #hazard ratio of death in Dysplasia helicobacter (+) vs healthy

r.NnD    <- - log(1 - p.NnD) # rate of death in healthy
r.NpD   <- hr.Np * r.NnD  	 # rate of death in Normal helicobacter (+)
r.GnD   <- hr.Gn * r.NnD  	 # rate of death in Gastritis helicobacter (-)
r.GpD   <- hr.Gp * r.NnD  	 # rate of death in Gastritis helicobacter (+)
r.AnD   <- hr.An * r.NnD  	 # rate of death in Atrophy helicobacter (-)
r.ApD   <- hr.Ap * r.NnD  	 # rate of death in Atrophy helicobacter (+)
r.InD   <- hr.In * r.NnD  	 # rate of death in Intestinal helicobacter (-)
r.IpD   <- hr.Ip * r.NnD     # rate of death in Intestinal helicobacter (+)
r.DnD   <- hr.Dn * r.NnD  	 # rate of death in Dysplasia helicobacter (-)
r.DpD   <- hr.Dp * r.NnD     # rate of death in Dysplasia helicobacter (+)


p.NpD <- 1- exp(-r.NpD)  # probability to die when Normal helicobacter (+)
p.GnD <- 1- exp(-r.GnD) # probability to die when Gastritis helicobacter (-)
p.GpD <- 1- exp(-r.GpD) # probability to die when Gastritis helicobacter (+)
p.AnD <- 1- exp(-r.AnD) # probability to die when Atrophy helicobacter (-)
p.ApD <- 1- exp(-r.ApD) # probability to die when Atrophy helicobacter (+)
p.InD <- 1- exp(-r.InD)  # probability to die when Intestinal helicobacter (-) 
p.IpD <- 1- exp(-r.IpD)  # probability to die when Intestinal helicobacter (+) 
p.DnD <- 1- exp(-r.DnD)  # probability to die when Dysplasia helicobacter (-) 
p.DpD <- 1- exp(-r.DpD)  # probability to die when Dysplasia helicobacter (+) 


# Costs and utilities  
c.Nn  <- 0                     # cost of remaining one cycle Normal helicobacter (-)
c.Np  <- 0                     # cost of remaining one cycle Normal helicobacter (+)
c.Gn  <- 0                     # cost of remaining one cycle Gastritis helicobacter (-)
c.Gp  <- 0                     # cost of remaining one cycle Gastritis helicobacter (+)
c.An  <- 0                     # cost of remaining one cycle Atrophy helicobacter (-)
c.Ap  <- 0                     # cost of remaining one cycle Atrophy helicobacter (+)
c.In  <- 0                     # cost of remaining one cycle Intestinal helicobacter (-) 
c.Ip  <- 0                     # cost of remaining one cycle Intestinal helicobacter (+) 
c.Dn  <- 0                     # cost of remaining one cycle Dysplasia helicobacter (-) 
c.Dp  <- 0                     # cost of remaining one cycle Dysplasia helicobacter (+)
c.D  <- 0  # cost of remaining one cycle dead

#Prevention costs (UF)
c.p1 <- 1.93                  #Cost of Prevention Strategy 1 "EDAUreasa"  # Estudio de verificación
c.p2 <- 4.39                     #Cost of Prevention Strategy 2 "UreaAire" # Precio PUC convenio preferente  
c.p3 <- 2                  #Cost of Prevention Strategy 3 "Serología"  # Precio falso
c.p4 <- 1.72                  #Cost of Prevention Strategy 4 "Antígenofecal" # Precio PUC convenio preferente   
c.p5 <- 2.34                  #Cost of Prevention Strategy 5 "EDAbp" # Estudio de verificación
c.p6 <- 4                  #Cost of Prevention Strategy 6 "SerologíaEDA" # Precio falso  
c.p7 <- 6.73                  #Cost of Prevention Strategy 7 "UreaAirePSEDA"  # c.p2 + c.p5
c.p8 <- 4.06                  #Cost of Prevention Strategy 8 "AntígenohecesPSEDA" # c.p4 + c.p5


u.Nn  <- 1                     # utility when Normal helicobacter (-)
u.Np  <- 1                     # utility when Normal helicobacter (+)
u.Gn  <- 1                     # utility when Gastritis helicobacter (-)
u.Gp  <- 1                     # utility when Gastritis helicobacter (+)
u.An  <- 1                     # utility when Atrophy helicobacter (-)
u.Ap  <- 1                     # utility when Atrophy helicobacter (+)
u.In  <- 1                     # utility when Intestinal helicobacter (-)
u.Ip  <- 1                     # utility when Intestinal helicobacter (+)
u.Dn  <- 1                     # utility when Dysplasia helicobacter (-)
u.Dp  <- 1                     # utility when Dysplasia helicobacter (+) 
u.D  <- 0                       # utility when dead

# Discounting factor
d.r  <- 0.03                    # equal discount of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


#### 03 Define and initialize matrices and vectors ####
#### 03.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_prev <- m.M_prev <- m.M_prev2 <- m.M_prev3 <- m.M_prev4 <- m.M_prev5 <- m.M_prev6 <- m.M_prev7 <- m.M_prev8 <- matrix(NA, 
                                  nrow = n.t + 1, ncol = n.s,
                                  dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_prev) # show first 6 rows of the matrix 

# The cohort starts as healthy

m.M_no_prev[1, ] <- m.M_prev[1, ] <- m.M_prev2[1, ] <- m.M_prev3[1, ] <- m.M_prev4[1, ] <- m.M_prev5[1, ] <- m.M_prev6[1, ] <- m.M_prev7[1, ] <- m.M_prev8[1, ] <-  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)                     # initialize first cycle of Markov trace

#### 03.2 Transition probability MATRIX ####
# create the transition probability matrix without prevention
m.P_noprev  <- matrix(0,
                      nrow = n.s,
                      ncol = n.s,
                      dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_noprev


# create the transition probability matrix with prevention strategy 1 
m.P_prev  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev

# create the transition probability matrix with prevention strategy 2 
m.P_prev2  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev2

# create the transition probability matrix with prevention strategy 3 
m.P_prev3  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev3

# create the transition probability matrix with prevention strategy 4 
m.P_prev4  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev4

# create the transition probability matrix with prevention strategy 5 
m.P_prev5  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev5

# create the transition probability matrix with prevention strategy 6 
m.P_prev6  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev6

# create the transition probability matrix with prevention strategy 7 
m.P_prev7  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

m.P_prev7

# create the transition probability matrix with prevention strategy 8 
m.P_prev8  <- matrix(0,
                    nrow = n.s,
                    ncol = n.s,
                    dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev8

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

### From Atrophy Helicobacter (-) 
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

# fill in the transition probability matrix witho prevention strategy 1 
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




# fill in the transition probability matrix witho prevention strategy 2
### From Normal helicobacter (-)
m.P_prev2["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p2NnNp - p.p2NnGn
m.P_prev2["Normalhpn", "Normalhpp"]    <- p.p2NnNp
m.P_prev2["Normalhpn", "Gastritishpn"]    <- p.p2NnGn
m.P_prev2["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev2["Normalhpp", "Normalhpp"]    <- 1 - p.p2NpGp - p.p2NpNn - p.NpD 
m.P_prev2["Normalhpp", "Gastritishpp"]    <- p.p2NpGp
m.P_prev2["Normalhpp", "Normalhpn"]    <- p.p2NpNn
m.P_prev2["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev2["Gastritishpn", "Gastritishpn"]    <- 1 - p.p2GnAn - p.p2GnGp - p.p2GnNn - p.GnD
m.P_prev2["Gastritishpn", "Atrophyhpn"]    <- p.p2GnAn
m.P_prev2["Gastritishpn", "Gastritishpp"]    <- p.p2GnGp
m.P_prev2["Gastritishpn", "Normalhpn"]    <- p.p2GnNn
m.P_prev2["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev2["Gastritishpp", "Gastritishpp"]    <- 1 - p.p2GpAp - p.p2GpGn - p.p2GpNp - p.GpD
m.P_prev2["Gastritishpp", "Atrophyhpp"]    <- p.p2GpAp
m.P_prev2["Gastritishpp", "Gastritishpn"]    <- p.p2GpGn
m.P_prev2["Gastritishpp", "Normalhpp"]    <- p.p2GpNp
m.P_prev2["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev2["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p2AnIn - p.p2AnAp - p.p2AnGn -  p.AnD
m.P_prev2["Atrophyhpn", "Intestinalhpn"]    <- p.p2AnIn
m.P_prev2["Atrophyhpn", "Atrophyhpp"]    <- p.p2AnAp
m.P_prev2["Atrophyhpn", "Gastritishpn"]    <- p.p2AnGn
m.P_prev2["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev2["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p2ApIp - p.p2ApAn - p.p2ApGp - p.ApD
m.P_prev2["Atrophyhpp", "Intestinalhpp"]    <- p.p2ApIp
m.P_prev2["Atrophyhpp", "Atrophyhpn"]    <- p.p2ApAn
m.P_prev2["Atrophyhpp", "Gastritishpp"]    <- p.p2ApGp
m.P_prev2["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev2["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p2InDn - p.p2InIp - p.p2InAn - p.InD 
m.P_prev2["Intestinalhpn", "Dysplasiahpn"]    <- p.p2InDn
m.P_prev2["Intestinalhpn", "Intestinalhpp"]    <- p.p2InIp
m.P_prev2["Intestinalhpn", "Atrophyhpn"]    <- p.p2InAn
m.P_prev2["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev2["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p2IpDp - p.p2IpIn - p.p2IpAp - p.IpD  
m.P_prev2["Intestinalhpp", "Dysplasiahpp"]    <- p.p2IpDp
m.P_prev2["Intestinalhpp", "Intestinalhpn"]    <- p.p2IpIn
m.P_prev2["Intestinalhpp", "Atrophyhpp"]    <- p.p2IpAp
m.P_prev2["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev2["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p2DnDp - p.p2DnIn -  p.DnD
m.P_prev2["Dysplasiahpn", "Dysplasiahpp"]    <- p.p2DnDp
m.P_prev2["Dysplasiahpn", "Intestinalhpn"]    <- p.p2DnIn
m.P_prev2["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev2["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p2DpDn - p.p2DpIp - p.DpD
m.P_prev2["Dysplasiahpp", "Dysplasiahpn"]    <- p.p2DpDn
m.P_prev2["Dysplasiahpp", "Intestinalhpp"]    <- p.p2DpIp
m.P_prev2["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 3
### From Normal helicobacter (-)
m.P_prev3["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p3NnNp - p.p3NnGn
m.P_prev3["Normalhpn", "Normalhpp"]    <- p.p3NnNp
m.P_prev3["Normalhpn", "Gastritishpn"]    <- p.p3NnGn
m.P_prev3["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev3["Normalhpp", "Normalhpp"]    <- 1 - p.p3NpGp - p.p3NpNn - p.NpD 
m.P_prev3["Normalhpp", "Gastritishpp"]    <- p.p3NpGp
m.P_prev3["Normalhpp", "Normalhpn"]    <- p.p3NpNn
m.P_prev3["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev3["Gastritishpn", "Gastritishpn"]    <- 1 - p.p3GnAn - p.p3GnGp - p.p3GnNn - p.GnD
m.P_prev3["Gastritishpn", "Atrophyhpn"]    <- p.p3GnAn
m.P_prev3["Gastritishpn", "Gastritishpp"]    <- p.p3GnGp
m.P_prev3["Gastritishpn", "Normalhpn"]    <- p.p3GnNn
m.P_prev3["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev3["Gastritishpp", "Gastritishpp"]    <- 1 - p.p3GpAp - p.p3GpGn - p.p3GpNp - p.GpD
m.P_prev3["Gastritishpp", "Atrophyhpp"]    <- p.p3GpAp
m.P_prev3["Gastritishpp", "Gastritishpn"]    <- p.p3GpGn
m.P_prev3["Gastritishpp", "Normalhpp"]    <- p.p3GpNp
m.P_prev3["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev3["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p3AnIn - p.p3AnAp - p.p3AnGn -  p.AnD
m.P_prev3["Atrophyhpn", "Intestinalhpn"]    <- p.p3AnIn
m.P_prev3["Atrophyhpn", "Atrophyhpp"]    <- p.p3AnAp
m.P_prev3["Atrophyhpn", "Gastritishpn"]    <- p.p3AnGn
m.P_prev3["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev3["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p3ApIp - p.p3ApAn - p.p3ApGp - p.ApD
m.P_prev3["Atrophyhpp", "Intestinalhpp"]    <- p.p3ApIp
m.P_prev3["Atrophyhpp", "Atrophyhpn"]    <- p.p3ApAn
m.P_prev3["Atrophyhpp", "Gastritishpp"]    <- p.p3ApGp
m.P_prev3["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev3["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p3InDn - p.p3InIp - p.p3InAn - p.InD 
m.P_prev3["Intestinalhpn", "Dysplasiahpn"]    <- p.p3InDn
m.P_prev3["Intestinalhpn", "Intestinalhpp"]    <- p.p3InIp
m.P_prev3["Intestinalhpn", "Atrophyhpn"]    <- p.p3InAn
m.P_prev3["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev3["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p3IpDp - p.p3IpIn - p.p3IpAp - p.IpD  
m.P_prev3["Intestinalhpp", "Dysplasiahpp"]    <- p.p3IpDp
m.P_prev3["Intestinalhpp", "Intestinalhpn"]    <- p.p3IpIn
m.P_prev3["Intestinalhpp", "Atrophyhpp"]    <- p.p3IpAp
m.P_prev3["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev3["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p3DnDp - p.p3DnIn -  p.DnD
m.P_prev3["Dysplasiahpn", "Dysplasiahpp"]    <- p.p3DnDp
m.P_prev3["Dysplasiahpn", "Intestinalhpn"]    <- p.p3DnIn
m.P_prev3["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev3["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p3DpDn - p.p3DpIp - p.DpD
m.P_prev3["Dysplasiahpp", "Dysplasiahpn"]    <- p.p3DpDn
m.P_prev3["Dysplasiahpp", "Intestinalhpp"]    <- p.p3DpIp
m.P_prev3["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 4
### From Normal helicobacter (-)
m.P_prev4["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p4NnNp - p.p4NnGn
m.P_prev4["Normalhpn", "Normalhpp"]    <- p.p4NnNp
m.P_prev4["Normalhpn", "Gastritishpn"]    <- p.p4NnGn
m.P_prev4["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev4["Normalhpp", "Normalhpp"]    <- 1 - p.p4NpGp - p.p4NpNn - p.NpD 
m.P_prev4["Normalhpp", "Gastritishpp"]    <- p.p4NpGp
m.P_prev4["Normalhpp", "Normalhpn"]    <- p.p4NpNn
m.P_prev4["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev4["Gastritishpn", "Gastritishpn"]    <- 1 - p.p4GnAn - p.p4GnGp - p.p4GnNn - p.GnD
m.P_prev4["Gastritishpn", "Atrophyhpn"]    <- p.p4GnAn
m.P_prev4["Gastritishpn", "Gastritishpp"]    <- p.p4GnGp
m.P_prev4["Gastritishpn", "Normalhpn"]    <- p.p4GnNn
m.P_prev4["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev4["Gastritishpp", "Gastritishpp"]    <- 1 - p.p4GpAp - p.p4GpGn - p.p4GpNp - p.GpD
m.P_prev4["Gastritishpp", "Atrophyhpp"]    <- p.p4GpAp
m.P_prev4["Gastritishpp", "Gastritishpn"]    <- p.p4GpGn
m.P_prev4["Gastritishpp", "Normalhpp"]    <- p.p4GpNp
m.P_prev4["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev4["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p4AnIn - p.p4AnAp - p.p4AnGn -  p.AnD
m.P_prev4["Atrophyhpn", "Intestinalhpn"]    <- p.p4AnIn
m.P_prev4["Atrophyhpn", "Atrophyhpp"]    <- p.p4AnAp
m.P_prev4["Atrophyhpn", "Gastritishpn"]    <- p.p4AnGn
m.P_prev4["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev4["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p4ApIp - p.p4ApAn - p.p4ApGp - p.ApD
m.P_prev4["Atrophyhpp", "Intestinalhpp"]    <- p.p4ApIp
m.P_prev4["Atrophyhpp", "Atrophyhpn"]    <- p.p4ApAn
m.P_prev4["Atrophyhpp", "Gastritishpp"]    <- p.p4ApGp
m.P_prev4["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev4["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p4InDn - p.p4InIp - p.p4InAn - p.InD 
m.P_prev4["Intestinalhpn", "Dysplasiahpn"]    <- p.p4InDn
m.P_prev4["Intestinalhpn", "Intestinalhpp"]    <- p.p4InIp
m.P_prev4["Intestinalhpn", "Atrophyhpn"]    <- p.p4InAn
m.P_prev4["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev4["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p4IpDp - p.p4IpIn - p.p4IpAp - p.IpD  
m.P_prev4["Intestinalhpp", "Dysplasiahpp"]    <- p.p4IpDp
m.P_prev4["Intestinalhpp", "Intestinalhpn"]    <- p.p4IpIn
m.P_prev4["Intestinalhpp", "Atrophyhpp"]    <- p.p4IpAp
m.P_prev4["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev4["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p4DnDp - p.p4DnIn -  p.DnD
m.P_prev4["Dysplasiahpn", "Dysplasiahpp"]    <- p.p4DnDp
m.P_prev4["Dysplasiahpn", "Intestinalhpn"]    <- p.p4DnIn
m.P_prev4["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev4["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p4DpDn - p.p4DpIp - p.DpD
m.P_prev4["Dysplasiahpp", "Dysplasiahpn"]    <- p.p4DpDn
m.P_prev4["Dysplasiahpp", "Intestinalhpp"]    <- p.p4DpIp
m.P_prev4["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 5
### From Normal helicobacter (-)
m.P_prev5["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p5NnNp - p.p5NnGn
m.P_prev5["Normalhpn", "Normalhpp"]    <- p.p5NnNp
m.P_prev5["Normalhpn", "Gastritishpn"]    <- p.p5NnGn
m.P_prev5["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev5["Normalhpp", "Normalhpp"]    <- 1 - p.p5NpGp - p.p5NpNn - p.NpD 
m.P_prev5["Normalhpp", "Gastritishpp"]    <- p.p5NpGp
m.P_prev5["Normalhpp", "Normalhpn"]    <- p.p5NpNn
m.P_prev5["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev5["Gastritishpn", "Gastritishpn"]    <- 1 - p.p5GnAn - p.p5GnGp - p.p5GnNn - p.GnD
m.P_prev5["Gastritishpn", "Atrophyhpn"]    <- p.p5GnAn
m.P_prev5["Gastritishpn", "Gastritishpp"]    <- p.p5GnGp
m.P_prev5["Gastritishpn", "Normalhpn"]    <- p.p5GnNn
m.P_prev5["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev5["Gastritishpp", "Gastritishpp"]    <- 1 - p.p5GpAp - p.p5GpGn - p.p5GpNp - p.GpD
m.P_prev5["Gastritishpp", "Atrophyhpp"]    <- p.p5GpAp
m.P_prev5["Gastritishpp", "Gastritishpn"]    <- p.p5GpGn
m.P_prev5["Gastritishpp", "Normalhpp"]    <- p.p5GpNp
m.P_prev5["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev5["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p5AnIn - p.p5AnAp - p.p5AnGn -  p.AnD
m.P_prev5["Atrophyhpn", "Intestinalhpn"]    <- p.p5AnIn
m.P_prev5["Atrophyhpn", "Atrophyhpp"]    <- p.p5AnAp
m.P_prev5["Atrophyhpn", "Gastritishpn"]    <- p.p5AnGn
m.P_prev5["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev5["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p5ApIp - p.p5ApAn - p.p5ApGp - p.ApD
m.P_prev5["Atrophyhpp", "Intestinalhpp"]    <- p.p5ApIp
m.P_prev5["Atrophyhpp", "Atrophyhpn"]    <- p.p5ApAn
m.P_prev5["Atrophyhpp", "Gastritishpp"]    <- p.p5ApGp
m.P_prev5["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev5["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p5InDn - p.p5InIp - p.p5InAn - p.InD 
m.P_prev5["Intestinalhpn", "Dysplasiahpn"]    <- p.p5InDn
m.P_prev5["Intestinalhpn", "Intestinalhpp"]    <- p.p5InIp
m.P_prev5["Intestinalhpn", "Atrophyhpn"]    <- p.p5InAn
m.P_prev5["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev5["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p5IpDp - p.p5IpIn - p.p5IpAp - p.IpD  
m.P_prev5["Intestinalhpp", "Dysplasiahpp"]    <- p.p5IpDp
m.P_prev5["Intestinalhpp", "Intestinalhpn"]    <- p.p5IpIn
m.P_prev5["Intestinalhpp", "Atrophyhpp"]    <- p.p5IpAp
m.P_prev5["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev5["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p5DnDp - p.p5DnIn -  p.DnD
m.P_prev5["Dysplasiahpn", "Dysplasiahpp"]    <- p.p5DnDp
m.P_prev5["Dysplasiahpn", "Intestinalhpn"]    <- p.p5DnIn
m.P_prev5["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev5["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p5DpDn - p.p5DpIp - p.DpD
m.P_prev5["Dysplasiahpp", "Dysplasiahpn"]    <- p.p5DpDn
m.P_prev5["Dysplasiahpp", "Intestinalhpp"]    <- p.p5DpIp
m.P_prev5["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 6
### From Normal helicobacter (-)
m.P_prev6["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p6NnNp - p.p6NnGn
m.P_prev6["Normalhpn", "Normalhpp"]    <- p.p6NnNp
m.P_prev6["Normalhpn", "Gastritishpn"]    <- p.p6NnGn
m.P_prev6["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev6["Normalhpp", "Normalhpp"]    <- 1 - p.p6NpGp - p.p6NpNn - p.NpD 
m.P_prev6["Normalhpp", "Gastritishpp"]    <- p.p6NpGp
m.P_prev6["Normalhpp", "Normalhpn"]    <- p.p6NpNn
m.P_prev6["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev6["Gastritishpn", "Gastritishpn"]    <- 1 - p.p6GnAn - p.p6GnGp - p.p6GnNn - p.GnD
m.P_prev6["Gastritishpn", "Atrophyhpn"]    <- p.p6GnAn
m.P_prev6["Gastritishpn", "Gastritishpp"]    <- p.p6GnGp
m.P_prev6["Gastritishpn", "Normalhpn"]    <- p.p6GnNn
m.P_prev6["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev6["Gastritishpp", "Gastritishpp"]    <- 1 - p.p6GpAp - p.p6GpGn - p.p6GpNp - p.GpD
m.P_prev6["Gastritishpp", "Atrophyhpp"]    <- p.p6GpAp
m.P_prev6["Gastritishpp", "Gastritishpn"]    <- p.p6GpGn
m.P_prev6["Gastritishpp", "Normalhpp"]    <- p.p6GpNp
m.P_prev6["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev6["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p6AnIn - p.p6AnAp - p.p6AnGn -  p.AnD
m.P_prev6["Atrophyhpn", "Intestinalhpn"]    <- p.p6AnIn
m.P_prev6["Atrophyhpn", "Atrophyhpp"]    <- p.p6AnAp
m.P_prev6["Atrophyhpn", "Gastritishpn"]    <- p.p6AnGn
m.P_prev6["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev6["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p6ApIp - p.p6ApAn - p.p6ApGp - p.ApD
m.P_prev6["Atrophyhpp", "Intestinalhpp"]    <- p.p6ApIp
m.P_prev6["Atrophyhpp", "Atrophyhpn"]    <- p.p6ApAn
m.P_prev6["Atrophyhpp", "Gastritishpp"]    <- p.p6ApGp
m.P_prev6["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev6["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p6InDn - p.p6InIp - p.p6InAn - p.InD 
m.P_prev6["Intestinalhpn", "Dysplasiahpn"]    <- p.p6InDn
m.P_prev6["Intestinalhpn", "Intestinalhpp"]    <- p.p6InIp
m.P_prev6["Intestinalhpn", "Atrophyhpn"]    <- p.p6InAn
m.P_prev6["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev6["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p6IpDp - p.p6IpIn - p.p6IpAp - p.IpD  
m.P_prev6["Intestinalhpp", "Dysplasiahpp"]    <- p.p6IpDp
m.P_prev6["Intestinalhpp", "Intestinalhpn"]    <- p.p6IpIn
m.P_prev6["Intestinalhpp", "Atrophyhpp"]    <- p.p6IpAp
m.P_prev6["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev6["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p6DnDp - p.p6DnIn -  p.DnD
m.P_prev6["Dysplasiahpn", "Dysplasiahpp"]    <- p.p6DnDp
m.P_prev6["Dysplasiahpn", "Intestinalhpn"]    <- p.p6DnIn
m.P_prev6["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev6["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p6DpDn - p.p6DpIp - p.DpD
m.P_prev6["Dysplasiahpp", "Dysplasiahpn"]    <- p.p6DpDn
m.P_prev6["Dysplasiahpp", "Intestinalhpp"]    <- p.p6DpIp
m.P_prev6["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 7
### From Normal helicobacter (-)
m.P_prev7["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p7NnNp - p.p7NnGn
m.P_prev7["Normalhpn", "Normalhpp"]    <- p.p7NnNp
m.P_prev7["Normalhpn", "Gastritishpn"]    <- p.p7NnGn
m.P_prev7["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev7["Normalhpp", "Normalhpp"]    <- 1 - p.p7NpGp - p.p7NpNn - p.NpD 
m.P_prev7["Normalhpp", "Gastritishpp"]    <- p.p7NpGp
m.P_prev7["Normalhpp", "Normalhpn"]    <- p.p7NpNn
m.P_prev7["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev7["Gastritishpn", "Gastritishpn"]    <- 1 - p.p7GnAn - p.p7GnGp - p.p7GnNn - p.GnD
m.P_prev7["Gastritishpn", "Atrophyhpn"]    <- p.p7GnAn
m.P_prev7["Gastritishpn", "Gastritishpp"]    <- p.p7GnGp
m.P_prev7["Gastritishpn", "Normalhpn"]    <- p.p7GnNn
m.P_prev7["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev7["Gastritishpp", "Gastritishpp"]    <- 1 - p.p7GpAp - p.p7GpGn - p.p7GpNp - p.GpD
m.P_prev7["Gastritishpp", "Atrophyhpp"]    <- p.p7GpAp
m.P_prev7["Gastritishpp", "Gastritishpn"]    <- p.p7GpGn
m.P_prev7["Gastritishpp", "Normalhpp"]    <- p.p7GpNp
m.P_prev7["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev7["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p7AnIn - p.p7AnAp - p.p7AnGn -  p.AnD
m.P_prev7["Atrophyhpn", "Intestinalhpn"]    <- p.p7AnIn
m.P_prev7["Atrophyhpn", "Atrophyhpp"]    <- p.p7AnAp
m.P_prev7["Atrophyhpn", "Gastritishpn"]    <- p.p7AnGn
m.P_prev7["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev7["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p7ApIp - p.p7ApAn - p.p7ApGp - p.ApD
m.P_prev7["Atrophyhpp", "Intestinalhpp"]    <- p.p7ApIp
m.P_prev7["Atrophyhpp", "Atrophyhpn"]    <- p.p7ApAn
m.P_prev7["Atrophyhpp", "Gastritishpp"]    <- p.p7ApGp
m.P_prev7["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev7["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p7InDn - p.p7InIp - p.p7InAn - p.InD 
m.P_prev7["Intestinalhpn", "Dysplasiahpn"]    <- p.p7InDn
m.P_prev7["Intestinalhpn", "Intestinalhpp"]    <- p.p7InIp
m.P_prev7["Intestinalhpn", "Atrophyhpn"]    <- p.p7InAn
m.P_prev7["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev7["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p7IpDp - p.p7IpIn - p.p7IpAp - p.IpD  
m.P_prev7["Intestinalhpp", "Dysplasiahpp"]    <- p.p7IpDp
m.P_prev7["Intestinalhpp", "Intestinalhpn"]    <- p.p7IpIn
m.P_prev7["Intestinalhpp", "Atrophyhpp"]    <- p.p7IpAp
m.P_prev7["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev7["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p7DnDp - p.p7DnIn -  p.DnD
m.P_prev7["Dysplasiahpn", "Dysplasiahpp"]    <- p.p7DnDp
m.P_prev7["Dysplasiahpn", "Intestinalhpn"]    <- p.p7DnIn
m.P_prev7["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev7["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p7DpDn - p.p7DpIp - p.DpD
m.P_prev7["Dysplasiahpp", "Dysplasiahpn"]    <- p.p7DpDn
m.P_prev7["Dysplasiahpp", "Intestinalhpp"]    <- p.p7DpIp
m.P_prev7["Dysplasiahpp", "Dead"]    <- p.DpD

# fill in the transition probability matrix witho prevention strategy 8
### From Normal helicobacter (-)
m.P_prev8["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p8NnNp - p.p8NnGn
m.P_prev8["Normalhpn", "Normalhpp"]    <- p.p8NnNp
m.P_prev8["Normalhpn", "Gastritishpn"]    <- p.p8NnGn
m.P_prev8["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev8["Normalhpp", "Normalhpp"]    <- 1 - p.p8NpGp - p.p8NpNn - p.NpD 
m.P_prev8["Normalhpp", "Gastritishpp"]    <- p.p8NpGp
m.P_prev8["Normalhpp", "Normalhpn"]    <- p.p8NpNn
m.P_prev8["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev8["Gastritishpn", "Gastritishpn"]    <- 1 - p.p8GnAn - p.p8GnGp - p.p8GnNn - p.GnD
m.P_prev8["Gastritishpn", "Atrophyhpn"]    <- p.p8GnAn
m.P_prev8["Gastritishpn", "Gastritishpp"]    <- p.p8GnGp
m.P_prev8["Gastritishpn", "Normalhpn"]    <- p.p8GnNn
m.P_prev8["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev8["Gastritishpp", "Gastritishpp"]    <- 1 - p.p8GpAp - p.p8GpGn - p.p8GpNp - p.GpD
m.P_prev8["Gastritishpp", "Atrophyhpp"]    <- p.p8GpAp
m.P_prev8["Gastritishpp", "Gastritishpn"]    <- p.p8GpGn
m.P_prev8["Gastritishpp", "Normalhpp"]    <- p.p8GpNp
m.P_prev8["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev8["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p8AnIn - p.p8AnAp - p.p8AnGn -  p.AnD
m.P_prev8["Atrophyhpn", "Intestinalhpn"]    <- p.p8AnIn
m.P_prev8["Atrophyhpn", "Atrophyhpp"]    <- p.p8AnAp
m.P_prev8["Atrophyhpn", "Gastritishpn"]    <- p.p8AnGn
m.P_prev8["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev8["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p8ApIp - p.p8ApAn - p.p8ApGp - p.ApD
m.P_prev8["Atrophyhpp", "Intestinalhpp"]    <- p.p8ApIp
m.P_prev8["Atrophyhpp", "Atrophyhpn"]    <- p.p8ApAn
m.P_prev8["Atrophyhpp", "Gastritishpp"]    <- p.p8ApGp
m.P_prev8["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev8["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p8InDn - p.p8InIp - p.p8InAn - p.InD 
m.P_prev8["Intestinalhpn", "Dysplasiahpn"]    <- p.p8InDn
m.P_prev8["Intestinalhpn", "Intestinalhpp"]    <- p.p8InIp
m.P_prev8["Intestinalhpn", "Atrophyhpn"]    <- p.p8InAn
m.P_prev8["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev8["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p8IpDp - p.p8IpIn - p.p8IpAp - p.IpD  
m.P_prev8["Intestinalhpp", "Dysplasiahpp"]    <- p.p8IpDp
m.P_prev8["Intestinalhpp", "Intestinalhpn"]    <- p.p8IpIn
m.P_prev8["Intestinalhpp", "Atrophyhpp"]    <- p.p8IpAp
m.P_prev8["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev8["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p8DnDp - p.p8DnIn -  p.DnD
m.P_prev8["Dysplasiahpn", "Dysplasiahpp"]    <- p.p8DnDp
m.P_prev8["Dysplasiahpn", "Intestinalhpn"]    <- p.p8DnIn
m.P_prev8["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev8["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p8DpDn - p.p8DpIp - p.DpD
m.P_prev8["Dysplasiahpp", "Dysplasiahpn"]    <- p.p8DpDn
m.P_prev8["Dysplasiahpp", "Intestinalhpp"]    <- p.p8DpIp
m.P_prev8["Dysplasiahpp", "Dead"]    <- p.DpD


### From Dead
m.P_prev["Dead", "Dead"] <- 1
m.P_prev2["Dead", "Dead"] <- 1
m.P_prev3["Dead", "Dead"] <- 1
m.P_prev4["Dead", "Dead"] <- 1
m.P_prev5["Dead", "Dead"] <- 1
m.P_prev6["Dead", "Dead"] <- 1
m.P_prev7["Dead", "Dead"] <- 1
m.P_prev8["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev)
rowSums(m.P_prev2)
rowSums(m.P_prev3)
rowSums(m.P_prev4)
rowSums(m.P_prev5)
rowSums(m.P_prev6)
rowSums(m.P_prev7)
rowSums(m.P_prev8)


#### 04 Run Markov model ####
for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_prev[t + 1, ] <- t(m.M_no_prev[t, ]) %*% m.P_noprev # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev[t + 1, ]    <- t(m.M_prev[t, ])    %*% m.P_prev   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev2[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev2   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev3[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev3   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev4[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev4   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev5[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev5   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev6[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev6   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev7[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev7   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev8[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev8   # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_prev)  # show the first 6 lines of the matrix
head(m.M_prev)  # show the first 6 lines of the matrix
head(m.M_prev2)  # show the first 6 lines of the matrix
head(m.M_prev3)  # show the first 6 lines of the matrix
head(m.M_prev4)  # show the first 6 lines of the matrix
head(m.M_prev5)  # show the first 6 lines of the matrix
head(m.M_prev6)  # show the first 6 lines of the matrix
head(m.M_prev7)  # show the first 6 lines of the matrix
head(m.M_prev8)  # show the first 6 lines of the matrix

#### 05 Compute and Plot Epidemiological Outcomes ####
#### 05.1 Cohort trace #####
matplot(m.M_no_prev, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

#### 05.2 Overall Survival (OS) #####
v.os_no_prev <- 1 - m.M_no_prev[, "Dead"]       # calculate the overall survival (OS) probability for no prevention
v.os_prev <- 1 - m.M_prev[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 1
v.os_prev2 <- 1 - m.M_prev2[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 2
v.os_prev3 <- 1 - m.M_prev3[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 3 
v.os_prev4 <- 1 - m.M_prev4[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 4 
v.os_prev5 <- 1 - m.M_prev5[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 5
v.os_prev6 <- 1 - m.M_prev6[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 6
v.os_prev7 <- 1 - m.M_prev7[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 7
v.os_prev8 <- 1 - m.M_prev8[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 8 

# Plot OS no prevention

plot(0:n.t, v.os_no_prev, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 1 

plot(0:n.t, v.os_prev, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 1")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

#### 05.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os_no_prev)                       # summing probablity of OS over time  (i.e. life expectancy)

#### 05.3 H. Pylori prevalence without prevention#####
v.preva <- rowSums(m.M_prev[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_no_prev
plot(v.preva,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence")

#### 05.4 H. Pylori prevalence with prevention strategy 1#####
v.preva1 <- rowSums(m.M_prev[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev
plot(v.preva1,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence EDA-Ureasa")

#### 05.3 H. Pylori prevalence with prevention strategy 2 #####
v.preva2 <- rowSums(m.M_prev2[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev2
plot(v.preva2,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence UreaAire")

#### 05.3 H. Pylori prevalence with prevention strategy 3#####
v.preva3 <- rowSums(m.M_prev3[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev3
plot(v.preva3,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Serología")

#### 05.3 H. Pylori prevalence prevention strategy 4 #####
v.preva4 <- rowSums(m.M_prev4[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev4
plot(v.preva4,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Antígeno fecal")

#### 05.3 H. Pylori prevalence prevention strategy 5 #####
v.preva5 <- rowSums(m.M_prev5[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev5
plot(v.preva5,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence EDAbp")

#### 05.3 H. Pylori prevalence prevention strategy 6 #####
v.preva6 <- rowSums(m.M_prev6[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev6
plot(v.preva6,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence SerologíaEDA")

#### 05.3 H. Pylori prevalence prevention strategy 7 #####
v.preva7 <- rowSums(m.M_prev7[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev7
plot(v.preva7,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence UreaAirePSEDA")

#### 05.3 H. Pylori prevalence prevention strategy 8 #####
v.preva8 <- rowSums(m.M_prev8[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev8
plot(v.preva8,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence AntígenohecesPSEDA")


#### 05.4 Proportion of sick in Gastritis - No prevention
v.prop.G <- rowSums(m.M_no_prev[, c("Gastritishpp", "Gastritishpn")])/ v.os_no_prev
plot(0:n.t, v.prop.G,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis No prevention", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 1 
v.prop.G1 <- rowSums(m.M_prev[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev
plot(0:n.t, v.prop.G1,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis- EDA-Ureasa", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 2 
v.prop.G2 <- rowSums(m.M_prev2[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev2
plot(0:n.t, v.prop.G2,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - Urea Aire", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 3 
v.prop.G3 <- rowSums(m.M_prev3[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev3
plot(0:n.t, v.prop.G3,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis- Serología", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 4 
v.prop.G4 <- rowSums(m.M_prev4[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev4
plot(0:n.t, v.prop.G4,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - Antígeno fecal", 
     col = "black", type = "l")


#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 5 
v.prop.G5 <- rowSums(m.M_prev5[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev5
plot(0:n.t, v.prop.G5,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - EDAbd", 
     col = "black", type = "l")


#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 6 
v.prop.G6 <- rowSums(m.M_prev6[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev6
plot(0:n.t, v.prop.G6,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - SerologíaEDA", 
     col = "black", type = "l")


#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 7 
v.prop.G7 <- rowSums(m.M_prev7[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev7
plot(0:n.t, v.prop.G7,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - UreaAirePSEDA", 
     col = "black", type = "l")


#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 8 
v.prop.G8 <- rowSums(m.M_prev8[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev8
plot(0:n.t, v.prop.G8,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - AntígenohecesPSEDA", 
     col = "black", type = "l")


#### 06 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_prev <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev2    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev3    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev4    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev5    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev6    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev7    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)
v.u_prev8    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.D)

v.c_no_prev <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.D)
v.c_prev    <- c(c.Nn + c.p1, c.Np, c.Gn + c.p1, c.Gp, c.An + c.p1, c.Ap, c.In + c.p1, c.Ip, c.Dn + c.p1, c.Dp, c.D)
v.c_prev2    <- c(c.Nn + c.p2, c.Np, c.Gn + c.p2, c.Gp, c.An + c.p2, c.Ap, c.In + c.p2, c.Ip, c.Dn + c.p2, c.Dp, c.D)
v.c_prev3    <- c(c.Nn +  c.p3, c.Np, c.Gn +  c.p3, c.Gp, c.An +  c.p3, c.Ap, c.In +  c.p3, c.Ip, c.Dn +  c.p3, c.Dp, c.D)
v.c_prev4    <- c(c.Nn +  c.p4, c.Np, c.Gn +  c.p4, c.Gp, c.An +  c.p4, c.Ap, c.In +  c.p4, c.Ip, c.Dn +  c.p4, c.Dp, c.D)
v.c_prev5    <- c(c.Nn + c.p5, c.Np, c.Gn + c.p5, c.Gp, c.An + c.p5, c.Ap, c.In + c.p5, c.Ip, c.Dn + c.p5, c.Dp, c.D)
v.c_prev6    <- c(c.Nn +  c.p6, c.Np, c.Gn +  c.p6, c.Gp, c.An +  c.p6, c.Ap, c.In +  c.p6, c.Ip, c.Dn +  c.p6, c.Dp, c.D)
v.c_prev7    <- c(c.Nn +  c.p7, c.Np, c.Gn + c.p7, c.Gp, c.An + c.p7, c.Ap, c.In + c.p7, c.Ip, c.Dn + c.p7, c.Dp, c.D)
v.c_prev8    <- c(c.Nn +  c.p8, c.Np, c.Gn +  c.p8, c.Gp, c.An +  c.p8, c.Ap, c.In +  c.p8, c.Ip, c.Dn +  c.p8, c.Dp, c.D)

#### 06.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_prev <- m.M_no_prev %*% v.u_no_prev
v.tu_prev <- m.M_prev %*% v.u_prev
v.tu_prev2 <- m.M_prev2 %*% v.u_prev2
v.tu_prev3 <- m.M_prev3 %*% v.u_prev3
v.tu_prev4 <- m.M_prev4 %*% v.u_prev4
v.tu_prev5 <- m.M_prev5 %*% v.u_prev5
v.tu_prev6 <- m.M_prev6 %*% v.u_prev6
v.tu_prev7 <- m.M_prev7 %*% v.u_prev7
v.tu_prev8 <- m.M_prev8 %*% v.u_prev8

v.tc_no_prev <- m.M_no_prev %*% v.c_no_prev
v.tc_prev    <- m.M_prev    %*% v.c_prev
v.tc_prev2    <- m.M_prev2    %*% v.c_prev2
v.tc_prev3    <- m.M_prev3    %*% v.c_prev3
v.tc_prev4    <- m.M_prev4    %*% v.c_prev4
v.tc_prev5    <- m.M_prev5    %*% v.c_prev5
v.tc_prev6    <- m.M_prev6    %*% v.c_prev6
v.tc_prev7    <- m.M_prev7    %*% v.c_prev7
v.tc_prev8    <- m.M_prev8    %*% v.c_prev8

#### 06.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_prev <- t(v.tu_no_prev) %*% v.dwe  
tu.d_prev    <- t(v.tu_prev)    %*% v.dwe
tu.d_prev2    <- t(v.tu_prev2)    %*% v.dwe
tu.d_prev3    <- t(v.tu_prev3)    %*% v.dwe
tu.d_prev4    <- t(v.tu_prev4)    %*% v.dwe
tu.d_prev5    <- t(v.tu_prev5)    %*% v.dwe
tu.d_prev6    <- t(v.tu_prev6)    %*% v.dwe
tu.d_prev7    <- t(v.tu_prev7)    %*% v.dwe
tu.d_prev8    <- t(v.tu_prev8)    %*% v.dwe

tc.d_no_prev <- t(v.tc_no_prev) %*% v.dwc
tc.d_prev    <- t(v.tc_prev)    %*% v.dwc
tc.d_prev2    <- t(v.tc_prev2)    %*% v.dwc
tc.d_prev3    <- t(v.tc_prev3)    %*% v.dwc
tc.d_prev4    <- t(v.tc_prev4)    %*% v.dwc
tc.d_prev5    <- t(v.tc_prev5)    %*% v.dwc
tc.d_prev6    <- t(v.tc_prev6)    %*% v.dwc
tc.d_prev7    <- t(v.tc_prev7)    %*% v.dwc
tc.d_prev8    <- t(v.tc_prev8)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_prev, tc.d_prev, tc.d_prev2, tc.d_prev3, tc.d_prev4, tc.d_prev5, tc.d_prev6, tc.d_prev7, tc.d_prev8)
v.tu.d <- c(tu.d_no_prev, tu.d_prev, tu.d_prev2, tu.d_prev3, tu.d_prev4, tu.d_prev5, tu.d_prev6, tu.d_prev7, tu.d_prev8)

# Matrix with discounted costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.tc.d,
                   Effect   = v.tu.d)
m.ce

#### 07 Compute ICERs of FONIS CaGa model ####
m.cea <- calculate_icers(cost = m.ce$Cost,
                         effect = m.ce$Effect,
                         strategies = m.ce$Strategy)
m.cea

#### 08 Plot frontier of Sick-Sicker model ####
plot(m.cea, xlim = c(15.7, 16.5))
ggsave("figs/Markov-FONISCaGast-CEA-Frontier-gastricmucosa.png", width = 8, height = 6)