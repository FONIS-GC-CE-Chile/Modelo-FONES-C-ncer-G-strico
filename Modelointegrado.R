#####################################################################################
##########             Modelo CE integrado progresión mucosa - estadíos cáncer gástrico       #####################
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
library(dampack)

#### 02 Input Model Parameters ####

### Strategies explanation

## Prevention Strategies Helicobacter Pilory 
# Prevention Strategy 1 - EDAUreasa: EDA con prueba de UREASA
# Prevention Strategy 2 - UreaAire: Prueba de Urea en Aire espirado marcado con Carbono 13 
# Prevention Strategy 3 - Serología: Serología para H Pylori 
# Prevention Strategy 4 - Antígenofecal: Estudio de Antígeno en heces fecales 

#Preventions Strategies combined 
#Prevention Strategy 5 - EDAbp: EDA con biopsia para prueba de ureasa y estudio histopatológico
#Prevention Strategy 6 - SerologiaPSEDA: Serología para H. pylori con erradicación para casos + y PS con EDA para estudio de casos +
#Prevention Strategy 7 - UreaAirePSEDA: Urea aire espirado C13 con erradicación + PS con EDA para estudio histológico de casos +.
#Prevention Strategy 8 - AntígenofecalPSEDA: Antígeno fecal cone erradicación de casos +, + PS con EDA para estudio histológico de +. 

## Prevention Strategies Cancer 
#Prevention Strategy 9 - EDA: EDA para estudio histológico
#Prevention Strategy 10 - PSEDA: Pepsinógeno sérico con EDA para estudio histológico en casos positivos
#Prevention Strategy 11 - RxEDA: Radiografía de doble contraste baritado con EDA para estudio histológico de casos + 


## Strategy names
v.names.str <- c("No Prevention", 
                 "EDAUreasa", "UreaAire", "Serología", "Antígenofecal", 
                 "EDAbp", "SerologíaPSEDA", "UreaAirePSEDA", "AntígenofecalPSEDA", 
                 "EDA", "PSEDA", "RxEDA")

## Number of strategies
n.str <- length(v.names.str)


## Markov model parameters
age     <- 20                                 # age at baseline
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles

v.n  <- c("Normalhpn", "Normalhpp", 
          "Gastritishpn", "Gastritishpp", 
          "Atrophyhpn", "Atrophyhpp", 
          "Intestinalhpn", "Intestinalhpp", 
          "Dysplasiahpn", "Dysplasiahpp", 
          "1preclinical", "2preclinical", 
          "3preclinical", "4preclinical", 
          "1clinicala", "1clinicalb", "1clinicalc", 
          "2clinicala", "2clinicalb", "2clinicalc",
          "3clinicala", "3clinicalb", "3clinicalc", 
          "4clinicala", "4clinicalb", "4clinicalc", 
          "Dead")    # state names
n.s  <- length(v.n)                     # number of states

##DESCRIBIR ESTADOS

## Normalhpn = Normal mucosa Helicobacter Pilory Negative
## Normalhpp = Normal mucosa Helicobacter Pilory Positive 
## Gastritishpn = Gastritis Helicobacter Pilory Negative
## Gastritishpp = Gastritis Helicobacter Pilory Positive 
## Atrophyhpn = Atrophy Helicobacter Pilory Negative
## Atrophyhpp = Atrophy Gastritis Helicobacter Pilory Positive  
## Intestinalhpn= Intestinal Metaplasia Helicobacter Pilory Negative
## Intestinalhpp = Intestinal Metaplasia Helicobacter Pilory Positive  
## Dysplasiahpn= Dysplasia Helicobacter Pilory Negative
## Dysplasiahpp = Dysplasia Helicobacter Pilory Positive 
## 1 preclinical = Stage 1, preclinical 
## 2 preclinical = Stage 2, preclinical 
## 3 preclinical = Stage 3, preclinical 
## 4 preclinical = Stage 4, preclinical 
## 1 clinical a = Stage 1, 1st year after diagnosis
## 1 clinical b = Stage 1, 2 -5 years after diagnosis 
## 1 clinical c = Stage 1, 6 years after diagnosis and over 
## 2 clinical a = Stage 2, 1st year after diagnosis
## 2 clinical b = Stage 2, 2 -5 years after diagnosis 
## 2 clinical c = Stage 2, 6 years after diagnosis and over 
## 3 clinical a = Stage 3, 1st year after diagnosis
## 3 clinical b = Stage 3, 2 -5 years after diagnosis 
## 3 clinical c = Stage 3, 6 years after diagnosis and over 
## 4 clinical a = Stage 4, 1st year after diagnosis
## 4 clinical b = Stage 4, 2 -5 years after diagnosis 
## 4 clinical c = Stage 4, 6 years after diagnosis and over
## Death = Death 

## Probabilities of Transition 

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
p.Dp1p <- 0.05  # probability to become stage 1 preclinical when Dysplasia (+)
p.Dn1p <- 0.05  # probability to become stage 1 preclinical when Dysplasia (-)
p.1p1ca <- 0.05  # probability to become stage 1 clinical a when stage 1 preclinical
p.1ca1cb <- 0.05 # probability to become stage 1 clinical b when stage 1 clinical a 
p.1cb1cc <- 0.05 # probability to become stage 1 clinical c when stage 1 clinical b 
p.1p2p <- 0.05  # probability to become stage 2 preclinical when stage 1 preclinical
p.2p2ca <- 0.05  # probability to become stage 2 clinical a when stage 2 preclinical
p.2ca2cb <- 0.05 # probability to become stage 2 clinical b when stage 2 clinical a 
p.2cb2cc <- 0.05 # probability to become stage 2 clinical c when stage 2 clinical b 
p.2p3p <- 0.05  # probability to become stage 3 preclinical when stage 2 preclinical
p.3p3ca <- 0.05  # probability to become stage 3 clinical a when stage 3 preclinical
p.3ca3cb <- 0.05 # probability to become stage 3 clinical b when stage 3 clinical a 
p.3cb3cc <- 0.05 # probability to become stage 3 clinical c when stage 3 clinical b 
p.3p4p <- 0.05  # probability to become stage 4 preclinical when stage 3 preclinical
p.4p4ca <- 0.05  # probability to become stage 4 clinical a when stage 4 preclinical
p.4ca4cb <- 0.05 # probability to become stage 4 clinical b when stage 4 clinical a 
p.4cb4cc <- 0.05 # probability to become stage 1 clinical c when stage 4 clinical b 
p.1ccNn <- 0.08 # probability to become Normal hp (-) when stage 1 clinical c
p.2ccNn <- 0.05 # probability to become Normal hp (-) when stage 2 clinical c
p.3ccNn <- 0.04 # probability to become Normal hp (-) when stage 3 clinical c
p.4ccNn <- 0.03 # probability to become Normal hp (-) when stage 4 clinical c

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
p.pDn <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.pDp <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.pDn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.pDp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c

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
p.p2_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p2_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p2_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p2_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p2_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p2_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p2_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p2_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p2_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p2_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p2_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p2_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p2_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p2_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p2_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p2_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p2_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p2_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p2_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p2_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p2_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p2_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p2_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p2_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p2_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p2_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p2_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p2_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p2_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p2_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p2_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p2_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p2_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p2_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p2_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p2_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p2_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p2_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p2_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p2_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p3_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p3_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p3_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p3_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p3_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p3_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p3_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p3_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p3_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p3_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p3_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p3_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p3_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p3_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p3_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p3_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p3_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p3_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p3_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p3_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p3_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p3_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p3_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p3_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p3_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p3_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p3_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p3_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p3_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p3_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p3_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p3_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p3_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p3_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p3_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p3_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p4_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p4_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p4_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p4_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p4_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p4_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p4_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p4_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p4_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p4_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p4_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p4_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p4_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p4_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p4_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p4_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p4_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p4_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p4_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p4_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p4_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p4_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p4_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p4_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p4_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p4_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p4_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p4_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p4_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p4_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p4_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p4_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p4_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p4_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p4_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p4_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p5_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p5_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p5_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p5_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p5_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p5_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p5_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p5_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p5_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p5_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p5_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p5_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p5_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p5_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p5_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p5_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p5_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p5_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p5_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p5_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p5_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p5_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p5_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p5_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p5_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p5_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p5_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p5_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p5_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p5_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p5_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p5_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p5_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p5_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p5_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p5_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p5_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p5_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p5_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p5_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p6_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p6_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p6_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p6_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p6_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p6_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p6_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p6_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p6_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p6_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p6_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p6_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p6_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p6_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p6_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p6_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p6_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p6_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p6_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p6_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p6_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p6_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p6_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p6_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p6_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p6_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p6_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p6_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p6_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p6_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p6_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p6_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p6_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p6_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p6_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p6_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p6_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p6_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p6_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p6_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p7_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p7_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p7_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p7_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p7_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p7_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p7_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p7_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p7_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p7_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p7_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p7_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p7_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p7_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p7_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p7_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p7_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p7_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p7_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p7_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p7_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p7_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p7_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p7_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p7_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p7_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p7_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p7_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p7_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p7_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p7_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p7_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p7_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p7_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p7_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p7_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p7_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p7_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p7_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p7_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c


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
p.p8_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p8_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p8_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p8_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p8_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p8_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p8_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p8_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p8_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p8_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p8_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p8_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p8_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p8_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p8_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p8_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p8_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p8_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p8_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p8_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p8_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p8_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p8_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p8_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p8_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p8_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p8_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p8_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p8_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p8_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p8_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p8_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p8_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p8_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p8_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p8_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p8_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p8_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p8_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p8_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c

#Probabilities of transitions with Cancer prevention strategies 

# Probability of transition when prevention strategy 9 (EDA) 

p.p9NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p9GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p9AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p9InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p9DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p9InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p9AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p9GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p9GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p9GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p9AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p9ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p9InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p9IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p9DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p9DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p9DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p9DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p9NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p9NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p9NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p9NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p9GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p9ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p9IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p9DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p9IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p9ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p9GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)
p.p9_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p9_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p9_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p9_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p9_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p9_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p9_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p9_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p9_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p9_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p9_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p9_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p9_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p9_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p9_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p9_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p9_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p9_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p9_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p9_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) clinical c
p.p9_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) clinical c

# Probability of transition when prevention strategy 10 (PSEDA) 
p.p10NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p10GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p10AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p10InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p10DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p10InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p10AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p10GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p10GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p10GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p10AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p10ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p10InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p10IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p10DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p10DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p10DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p10DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p10NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p10NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p10NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p10NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p10GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p10ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p10IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p10DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p10IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p10ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p10GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)
p.p10_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p10_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p10_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p10_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p10_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p10_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p10_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p10_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p10_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p10_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p10_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p10_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p10_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p10_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p10_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p10_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p10_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p10_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p10_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) when stage 2 clinical c
p.p10_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) when stage 3 clinical c
p.p10_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) when stage 4 clinical c

# Probability of transition when prevention strategy 11 (RxEDA) 
p.p11NnGn <- p.NnGn * 0.5  # probability to become Gastritis hp(-) when Normanl hp(-)
p.p11GnAn <- p.GnAn * 0.5  # probability to become Atrophy hp(-) when Gastritis hp(-)
p.p11AnIn <- p.AnIn * 0.5  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.p11InDn <- p.InDn * 0.5  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.p11DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p11InAn <- p.InAn * 0.5  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-)
p.p11AnGn <- p.AnGn * 0.5  # probability to become Gastritis hp(-) when Atrophy hp(-)
p.p11GnNn <- p.GnNn * 0.5  # probability to become Normal hp(-) when Gastritis hp(-)
p.p11GnGp <- p.GnGp * 0.5  # probability to become Gastritis hp(-) when Gastritis hp(+)
p.p11GpGn <- p.GpGn * 0.5 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p11AnAp <- p.AnAp * 0.5 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.p11ApAn <- p.ApAn * 0.5 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p11InIp <- p.InIp * 0.5 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.p11IpIn <- p.IpIn * 0.5 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p11DnIn <- p.DnIn * 0.5 # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.p11DnDp <- p.DnDp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p11DpIp <- p.DpIp * 0.5 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.p11DpDn <- p.DpDn * 0.5 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p11NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normal hp(+)
p.p11NnNp <- p.NnNp * 0.5  # probability to become Normal hp(+) when Normal hp(-)
p.p11NpNn <- p.NpNn * 0.5  # probability to become Normal hp(-) when Normal hp(+)
p.p11NpGp <- p.NpGp * 0.5  # probability to become Gastritis hp(+) when Normanl hp(+)
p.p11GpAp <- p.GpAp * 0.5  # probability to become Atrophy hp(+) when Gastritis hp(+)
p.p11ApIp <- p.ApIp * 0.5  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.p11IpDp <- p.IpDp * 0.5  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)
p.p11DpIp <- p.DpIp * 0.5  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+)
p.p11IpAp <- p.IpAp * 0.5  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)
p.p11ApGp <-p.ApGp * 0.5  # probability to become Gastritis hp(+) when Atrophy hp(+) 
p.p11GpNp <- p.GpNp * 0.5  # probability to become Normal hp (+) when Gastritis hp(+)
p.p11_Dn1p <- 0.05 * p.Dn1p  # probability to become stage 1 preclinical when Dysplasia Hp (-)
p.p11_Dp1p <- 0.05 * p.Dp1p  # probability to become stage 1 preclinical when Dysplasia Hp (+)
p.p11_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p11_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p11_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p11_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p11_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p11_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p11_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p11_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p11_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p11_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p11_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p11_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p11_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p11_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p11_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p11_1ccNn <- 0.08 * p.1ccNn # probability to become Normal mucosa Hp (-) when stage 1 clinical c
p.p11_2ccNn <- 0.05 * p.2ccNn # probability to become Normal mucosa Hp (-) when stage 2 clinical c
p.p11_3ccNn <- 0.04 * p.3ccNn # probability to become Normal mucosa Hp (-) when stage 3 clinical c
p.p11_4ccNn <- 0.03 * p.4ccNn # probability to become Normal mucosa Hp (-) when stage 4 clinical c

#Hazard ratio 

hr.Np <- 1 #hazard ratio of death in Normal helicobacter(+) vs Normal mucosa Hp (-)
hr.Gn <- 1 #hazard ratio of death in Gastritis helicobacter(-) vs Normal mucosa Hp (-)
hr.Gp <- 1 #hazard ratio of death in Gastritis helicobacter(+) vs Normal mucosa Hp (-)
hr.An <- 1 #hazard ratio of death in Atrophy helicobacter (-) vs Normal mucosa Hp (-)
hr.Ap <- 1 #hazard ratio of death in Atrophy helicobacter (+) vs Normal mucosa Hp (-)
hr.In <- 1 #hazard ratio of death in Intestional helicobacter (-) vs Normal mucosa Hp (-)
hr.Ip <- 1 #hazard ratio of death in Intestinal helicobacter (+) vs Normal mucosa Hp (-)
hr.Dn <- 1 #hazard ratio of death in Dysplasia helicobacter (-) vs Normal mucosa Hp (-)
hr.Dp <- 1 #hazard ratio of death in Dysplasia helicobacter (+) vs Normal mucosa Hp (-)
hr.1p <- 1 #hazard ratio of death in stage 1 preclinical vs Normal mucosa Hp (-)
hr.2p <- 1 #hazard ratio of death in stage 2 preclinical vs Normal mucosa Hp (-)
hr.3p <- 1 #hazard ratio of death in stage 3 preclinical vs Normal mucosa Hp (-)
hr.4p <- 1 #hazard ratio of death in stage 4 preclinical vs Normal mucosa Hp (-)
hr.1ca <- 3 #hazard ratio of death in stage 1 clinical a vs Normal mucosa Hp (-)
hr.2ca <- 6  #hazard ratio of death in stage 2 clinical a vs Normal mucosa Hp (-)   
hr.3ca <- 9 #hazard ratio of death in stage 3 clinical a vs Normal mucosa Hp (-)
hr.4ca <- 10 #hazard ratio of death in stage 4 clinical a vs Normal mucosa Hp (-)
hr.1cb <- 3 #hazard ratio of death in stage 1 clinical b vs Normal mucosa Hp (-)
hr.2cb <- 6  #hazard ratio of death in stage 2 clinical b vs Normal mucosa Hp (-)   
hr.3cb <- 9 #hazard ratio of death in stage 3 clinical b vs Normal mucosa Hp (-)
hr.4cb <- 8 #hazard ratio of death in stage 4 clinical b vs Normal mucosa Hp (-)
hr.1cc <- 3 #hazard ratio of death in stage 1 clinical c vs Normal mucosa Hp (-)
hr.2cc <- 6  #hazard ratio of death in stage 2 clinical c vs Normal mucosa Hp (-)   
hr.3cc <- 9 #hazard ratio of death in stage 3 clinical c vs Normal mucosa Hp (-)
hr.4cc <- 8 #hazard ratio of death in stage 4 clinical c vs Normal mucosa Hp (-)

r.NnD    <- - log(1 - p.NnD) # rate of death in Normal mucosa Hp (-)
r.NpD   <- hr.Np * r.NnD  	 # rate of death in Normal helicobacter (+)
r.GnD   <- hr.Gn * r.NnD  	 # rate of death in Gastritis helicobacter (-)
r.GpD   <- hr.Gp * r.NnD  	 # rate of death in Gastritis helicobacter (+)
r.AnD   <- hr.An * r.NnD  	 # rate of death in Atrophy helicobacter (-)
r.ApD   <- hr.Ap * r.NnD  	 # rate of death in Atrophy helicobacter (+)
r.InD   <- hr.In * r.NnD  	 # rate of death in Intestinal helicobacter (-)
r.IpD   <- hr.Ip * r.NnD     # rate of death in Intestinal helicobacter (+)
r.DnD   <- hr.Dn * r.NnD  	 # rate of death in Dysplasia helicobacter (-)
r.DpD   <- hr.Dp * r.NnD     # rate of death in Dysplasia helicobacter (+)
r.1pD   <- hr.1p * r.NnD  	 # rate of death in 1 preclinical
r.2pD   <- hr.2p * r.NnD  	 # rate of death in 2 preclinical
r.3pD   <- hr.3p * r.NnD  	 # rate of death in 3 preclinical
r.4pD   <- hr.4p * r.NnD  	 # rate of death in 4 preclinical
r.1caD   <- hr.1ca * r.NnD  	 # rate of death in 1 clinical a
r.2caD   <- hr.2ca * r.NnD  	 # rate of death in 2 clinical a
r.3caD   <- hr.3ca * r.NnD  	 # rate of death in 3 clinical a
r.4caD   <- hr.4ca * r.NnD  	 # rate of death in 4 clinical a  
r.1cbD   <- hr.1cb * r.NnD  	 # rate of death in 1 clinical b
r.2cbD   <- hr.2cb * r.NnD  	 # rate of death in 2 clinical b
r.3cbD   <- hr.3cb * r.NnD  	 # rate of death in 3 clinical b
r.4cbD   <- hr.4cb * r.NnD  	 # rate of death in 4 clinical b  
r.1ccD   <- hr.1cc * r.NnD  	 # rate of death in 1 clinical c
r.2ccD   <- hr.2cc * r.NnD  	 # rate of death in 2 clinical c
r.3ccD   <- hr.3cc * r.NnD  	 # rate of death in 3 clinical c
r.4ccD   <- hr.4cc * r.NnD  	 # rate of death in 4 clinical c


p.NpD <- 1- exp(-r.NpD)  # probability to die when Normal helicobacter (+)
p.GnD <- 1- exp(-r.GnD) # probability to die when Gastritis helicobacter (-)
p.GpD <- 1- exp(-r.GpD) # probability to die when Gastritis helicobacter (+)
p.AnD <- 1- exp(-r.AnD) # probability to die when Atrophy helicobacter (-)
p.ApD <- 1- exp(-r.ApD) # probability to die when Atrophy helicobacter (+)
p.InD <- 1- exp(-r.InD)  # probability to die when Intestinal helicobacter (-) 
p.IpD <- 1- exp(-r.IpD)  # probability to die when Intestinal helicobacter (+) 
p.DnD <- 1- exp(-r.DnD)  # probability to die when Dysplasia helicobacter (-) 
p.DpD <- 1- exp(-r.DpD)  # probability to die when Dysplasia helicobacter (+) 
p.1pD <- 1- exp(-r.1pD)  # probability to die when stage 1 preclinical 
p.2pD <- 1- exp(-r.2pD) # probability to die when stage 2 preclinical
p.3pD <- 1- exp(-r.3pD) # probability to die when stage 3 preclinical
p.4pD <- 1- exp(-r.4pD) # probability to die when stage 4 preclinical
p.1caD <- 1- exp(-r.1caD)  # probability to die when stage 1 clinical a
p.2caD <- 1- exp(-r.2caD) # probability to die when stage 2 clinical a
p.3caD <- 1- exp(-r.3caD) # probability to die when stage 3 clinical a
p.4caD <- 1- exp(-r.4caD) # probability to die when stage 4 clinical a
p.1cbD <- 1- exp(-r.1cbD)  # probability to die when stage 1 clinical b 
p.2cbD <- 1- exp(-r.2cbD) # probability to die when stage 2 clinical b
p.3cbD <- 1- exp(-r.3cbD) # probability to die when stage 3 clinical b
p.4cbD <- 1- exp(-r.4cbD) # probability to die when stage 4 clinical b
p.1ccD <- 1- exp(-r.1ccD)  # probability to die when stage 1 clinical c
p.2ccD <- 1- exp(-r.2ccD) # probability to die when stage 2 clinical c
p.3ccD <- 1- exp(-r.3ccD) # probability to die when stage 3 clinical c
p.4ccD <- 1- exp(-r.4ccD) # probability to die when stage 4 clinical c


# Costs   
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
c.1p  <- 0                     # cost of remaining one cycle stage 1 preclinical
c.2p  <- 0                     # cost of remaining one cycle stage 2 preclinical
c.3p  <- 0                     # cost of remaining one cycle stage 3 preclinical
c.4p  <- 0                     # cost of remaining one cycle stage 4 preclinical
c.1ca  <- 0                     # cost of remaining one cycle stage 1 clinical a (treatment)
c.2ca  <- 0                     # cost of remaining one cycle stage 2 clinical a (treatment)
c.3ca  <- 0                     # cost of remaining one cycle stage 3 clinical a (treatment)
c.4ca  <- 0                     # cost of remaining one cycle stage 4 clinical a (treatment)
c.1cb  <- 0                     # cost of remaining one cycle stage 1 clinical b (treatment)
c.2cb  <- 0                     # cost of remaining one cycle stage 2 clinical b (treatment)
c.3cb  <- 0                     # cost of remaining one cycle stage 3 clinical b (treatment)
c.4cb  <- 0                     # cost of remaining one cycle stage 4 clinical b (treatment)
c.1cc  <- 0                     # cost of remaining one cycle stage 1 clinical c (treatment)
c.2cc  <- 0                     # cost of remaining one cycle stage 2 clinical c (treatment)
c.3cc  <- 0                     # cost of remaining one cycle stage 3 clinical c (treatment)
c.4cc  <- 0                     # cost of remaining one cycle stage 4 clinical c (treatment)

c.D  <- 0  # cost of remaining one cycle dead

#Prevention costs (UF)
c.EDAUreasa <- 1.93                  #Cost of Prevention "EDAUreasa"  # Estudio de verificación
c.UreaAire <- 4.39                     #Cost of Prevention "UreaAire" # Precio PUC convenio preferente  
c.EDABp <- 2.34                  #Cost of Prevention "EDAbp" # Estudio de verificación
c.EDA <- 2.34                   # cost of EDA # Estudio de verificación de costos 
c.PS <- 2                   # cost of PS # El costo de pepsinógeno sérico no está disponible
c.Rx <- 4                   # cost of Rx # El costo de la radiografía doble contrate no está disponible
c.UreaAire <- 4.39              #Cost of Prevention Strategy 2 "UreaAire" # Precio PUC convenio preferente  
c.Serología <- 2            # cost of serología # Precio falso
c.Antígenofecal <- 1.72                  #Cost of "Antígenofecal" # Precio PUC convenio preferente 

#Erradication costs (UF)
c.err <- 1   #Cost of H. pilory therapy 

#Utilities

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
u.1p  <- 1                     # utility when stage 1 preclinical
u.2p  <- 1                     # utility when stage 2 preclinical
u.3p  <- 0.7                     # utility when stage 3 preclinical
u.4p  <- 0.8                     # utility when stage 4 preclinical
u.1ca  <- 0.5                     # utility when stage 1 clinical a
u.2ca  <- 0.4                     # utility when stage 2 clinical a
u.3ca  <- 0.3                     # utility when stage 3 clinical a
u.4ca  <- 0.2                     # utility when stage 4 clinical a
u.1cb  <- 0.5                     # utility when stage 1 clinical b
u.2cb  <- 0.4                     # utility when stage 2 clinical b
u.3cb  <- 0.3                     # utility when stage 3 clinical b
u.4cb  <- 0.2                     # utility when stage 4 clinical b
u.1cc  <- 0.5                     # utility when stage 1 clinical c
u.2cc  <- 0.4                     # utility when stage 2 clinical c
u.3cc  <- 0.3                     # utility when stage 3 clinical c
u.4cc  <- 0.2                     # utility when stage 4 clinical c
u.D  <- 0                       # utility when dead


# Discounting factor
d.r  <- 0.03                    # equal discount of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


#### 03 Define and initialize matrices and vectors ####
#### 03.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_prev <- m.M_prev <- m.M_prev2 <- m.M_prev3 <- m.M_prev4 <- m.M_prev5 <- m.M_prev6 <- m.M_prev7 <- m.M_prev8 <- m.M_prev9 <- m.M_prev10 <- m.M_prev11 <- matrix(NA, 
                                                                                                                             nrow = n.t + 1, ncol = n.s,
                                                                                                                             dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_prev) # show first 6 rows of the matrix 

# The cohort starts as healthy

m.M_no_prev[1, ] <- m.M_prev[1, ] <- m.M_prev2[1, ] <- m.M_prev3[1, ] <- m.M_prev4[1, ] <- m.M_prev5[1, ] <- m.M_prev6[1, ] <- m.M_prev7[1, ] <- m.M_prev8[1, ] <- m.M_prev9[1, ] <- m.M_prev10[1, ] <- m.M_prev11[1, ]  <-  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0, 0, 0 , 0 , 0, 0)                     # initialize first cycle of Markov trace

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

# create the transition probability matrix with prevention strategy 1 cancer (EDA)
m.P_prev9  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev9

# create the transition probability matrix with prevention strategy 2 cancer (PSEDA)
m.P_prev10  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev10

# create the transition probability matrix with prevention strategy 3 cancer (RxEDA)
m.P_prev11  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prev11

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
m.P_noprev["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.DnDp - p.DnIn -  p.DnD - p.Dn1p
m.P_noprev["Dysplasiahpn", "Dysplasiahpp"]    <- p.DnDp
m.P_noprev["Dysplasiahpn", "Intestinalhpn"]    <- p.DnIn
m.P_noprev["Dysplasiahpn", "1preclinical"]     <- p.Dn1p
m.P_noprev["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_noprev["Dysplasiahpp", "Dysplasiahpn"]    <- 1 - p.DpDn - p.DpIp - p.DpD - p.Dp1p
m.P_noprev["Dysplasiahpp", "Dysplasiahpp"]    <- p.DpDn
m.P_noprev["Dysplasiahpp", "Intestinalhpp"]    <- p.DpIp
m.P_noprev["Dysplasiahpp", "1preclinical"]     <- p.Dp1p
m.P_noprev["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_noprev["1preclinical", "1preclinical"]    <- 1 - p.1p1ca - p.1p2p - p.1pD	
m.P_noprev["1preclinical", "1clinicala"]    <- p.1p1ca	
m.P_noprev["1preclinical", "2preclinical"]    <- p.1p2p	
m.P_noprev["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_noprev["2preclinical", "2preclinical"]    <- 1 - p.2p2ca - p.2p3p - p.2pD	
m.P_noprev["2preclinical", "2clinicala"]    <- p.2p2ca	
m.P_noprev["2preclinical", "3preclinical"]    <- p.2p3p	
m.P_noprev["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_noprev["3preclinical", "3preclinical"]    <- 1 - p.3p3ca - p.3p4p - p.3pD	
m.P_noprev["3preclinical", "3clinicala"]    <- p.3p3ca	
m.P_noprev["3preclinical", "4preclinical"]    <- p.3p4p	
m.P_noprev["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_noprev["4preclinical", "4preclinical"]    <- 1 - p.4p4ca - p.4pD	
m.P_noprev["4preclinical", "4clinicala"]    <- p.4p4ca	
m.P_noprev["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_noprev["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.1ca1cb	
m.P_noprev["1clinicala", "1clinicalb"] <- p.1ca1cb	
m.P_noprev["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_noprev["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.2ca2cb	
m.P_noprev["2clinicala", "2clinicalb"] <- p.2ca2cb	
m.P_noprev["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_noprev["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.3ca3cb 	
m.P_noprev["3clinicala", "3clinicalb"] <- p.3ca3cb	
m.P_noprev["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_noprev["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.4ca4cb	
m.P_noprev["4clinicala", "4clinicalb"] <- p.4ca4cb  	
m.P_noprev["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_noprev["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.1cb1cc	
m.P_noprev["1clinicalb", "1clinicalc"] <- p.1cb1cc	
m.P_noprev["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_noprev["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.2cb2cc	
m.P_noprev["2clinicalb", "2clinicalc"] <- p.2cb2cc	
m.P_noprev["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_noprev["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.3cb3cc 	
m.P_noprev["3clinicalb", "3clinicalc"] <- p.3cb3cc	
m.P_noprev["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_noprev["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.4cb4cc 	
m.P_noprev["4clinicalb", "4clinicalc"] <- p.4cb4cc	
m.P_noprev["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_noprev["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.1ccNn	
m.P_noprev["1clinicalc", "Normalhpn"] <- p.1ccNn	
m.P_noprev["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_noprev["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.2ccNn	
m.P_noprev["2clinicalc", "Normalhpn"] <- p.2ccNn	
m.P_noprev["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_noprev["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p1_3ccNn 	
m.P_noprev["3clinicalc", "Normalhpn"] <- p.p1_3ccNn
m.P_noprev["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_noprev["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p1_4ccNn 	
m.P_noprev["4clinicalc", "Normalhpn"] <- p.p1_4ccNn	
m.P_noprev["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_noprev["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_noprev)

# fill in the transition probability matrix with prevention strategy 1 

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
m.P_prev["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.pDnDp - p.pDnIn -  p.DnD - p.pDn1p
m.P_prev["Dysplasiahpn", "Dysplasiahpp"]    <- p.pDnDp
m.P_prev["Dysplasiahpn", "Intestinalhpn"]    <- p.pDnIn
m.P_prev["Dysplasiahpn", "1preclinical"]    <- p.pDn1p
m.P_prev["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.pDpDn - p.pDpIp - p.DpD - p.pDp1p
m.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
m.P_prev["Dysplasiahpp", "Intestinalhpp"]    <- p.pDpIp
m.P_prev["Dysplasiahpp", "1preclinical"]    <- p.pDp1p
m.P_prev["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev["1preclinical", "1preclinical"]    <- 1 - p.p1p1ca - p.p1p2p - p.1pD	
m.P_prev["1preclinical", "1clinicala"]    <- p.p1p1ca	
m.P_prev["1preclinical", "2preclinical"]    <- p.p1p2p	
m.P_prev["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev["2preclinical", "2preclinical"]    <- 1 - p.p2p2ca - p.p2p3p - p.2pD	
m.P_prev["2preclinical", "2clinicala"]    <- p.p2p2ca	
m.P_prev["2preclinical", "3preclinical"]    <- p.p2p3p	
m.P_prev["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev["3preclinical", "3preclinical"]    <- 1 - p.p3p3ca - p.p3p4p - p.p3pD	
m.P_prev["3preclinical", "3clinicala"]    <- p.p3p3ca	
m.P_prev["3preclinical", "4preclinical"]    <- p.p3p4p	
m.P_prev["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev["4preclinical", "4preclinical"]    <- 1 - p.p4p4ca - p.4pD	
m.P_prev["4preclinical", "4clinicala"]    <- p.p4p4ca	
m.P_prev["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p1ca1cb	
m.P_prev["1clinicala", "1clinicalb"] <- p.p1ca1cb	
m.P_prev["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p2ca2cb	
m.P_prev["2clinicala", "2clinicalb"] <- p.p2ca2cb	
m.P_prev["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p3ca3cb 	
m.P_prev["3clinicala", "3clinicalb"] <- p.p3ca3cb	
m.P_prev["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p4ca4cb	
m.P_prev["4clinicala", "4clinicalb"] <- p.p4ca4cb  	
m.P_prev["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p1cb1cc	
m.P_prev["1clinicalb", "1clinicalc"] <- p.p1cb1cc	
m.P_prev["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p2cb2cc	
m.P_prev["2clinicalb", "2clinicalc"] <- p.p2cb2cc	
m.P_prev["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p3cb3cc 	
m.P_prev["3clinicalb", "3clinicalc"] <- p.p3cb3cc	
m.P_prev["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p4cb4cc 	
m.P_prev["4clinicalb", "4clinicalc"] <- p.p4cb4cc	
m.P_prev["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p1ccNn	
m.P_prev["1clinicalc", "Normalhpn"] <- p.p1ccNn	
m.P_prev["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p2ccNn	
m.P_prev["2clinicalc", "Normalhpn"] <- p.p2ccNn	
m.P_prev["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p3ccNn 	
m.P_prev["3clinicalc", "Normalhpn"] <- p.p3ccNn
m.P_prev["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p4ccNn 	
m.P_prev["4clinicalc", "Normalhpn"] <- p.p4ccNn	
m.P_prev["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev)


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
m.P_prev2["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p2DnDp - p.p2DnIn -  p.DnD - p.Dn1p
m.P_prev2["Dysplasiahpn", "Dysplasiahpp"]    <- p.p2DnDp
m.P_prev2["Dysplasiahpn", "Intestinalhpn"]    <- p.p2DnIn
m.P_prev2["Dysplasiahpn", "1preclinical"]    <- p.Dn1p
m.P_prev2["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev2["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p2DpDn - p.p2DpIp - p.DpD - p.p2_Dp1p
m.P_prev2["Dysplasiahpp", "Dysplasiahpn"]    <- p.p2DpDn
m.P_prev2["Dysplasiahpp", "Intestinalhpp"]    <- p.p2DpIp
m.P_prev2["Dysplasiahpp", "1preclinical"]    <- p.p2_Dp1p
m.P_prev2["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev2["1preclinical", "1preclinical"]    <- 1 - p.p2_1p1ca - p.p2_1p2p - p.1pD	
m.P_prev2["1preclinical", "1clinicala"]    <- p.p2_1p1ca	
m.P_prev2["1preclinical", "2preclinical"]    <- p.p2_1p2p	
m.P_prev2["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev2["2preclinical", "2preclinical"]    <- 1 - p.p2_2p2ca - p.p2_2p3p - p.2pD	
m.P_prev2["2preclinical", "2clinicala"]    <- p.p2_2p2ca	
m.P_prev2["2preclinical", "3preclinical"]    <- p.p2_2p3p	
m.P_prev2["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev2["3preclinical", "3preclinical"]    <- 1 - p.p2_3p3ca - p.p2_3p4p - p.3pD	
m.P_prev2["3preclinical", "3clinicala"]    <- p.p2_3p3ca	
m.P_prev2["3preclinical", "4preclinical"]    <- p.p2_3p4p	
m.P_prev2["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev2["4preclinical", "4preclinical"]    <- 1 - p.p2_4p4ca - p.4pD	
m.P_prev2["4preclinical", "4clinicala"]    <- p.p2_4p4ca	
m.P_prev2["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev2["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p2_1ca1cb	
m.P_prev2["1clinicala", "1clinicalb"] <- p.p2_1ca1cb	
m.P_prev2["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev2["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p2_2ca2cb	
m.P_prev2["2clinicala", "2clinicalb"] <- p.p2_2ca2cb	
m.P_prev2["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev2["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p2_3ca3cb 	
m.P_prev2["3clinicala", "3clinicalb"] <- p.p2_3ca3cb	
m.P_prev2["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev2["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p2_4ca4cb	
m.P_prev2["4clinicala", "4clinicalb"] <- p.p2_4ca4cb  	
m.P_prev2["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev2["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p2_1cb1cc	
m.P_prev2["1clinicalb", "1clinicalc"] <- p.p2_1cb1cc	
m.P_prev2["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev2["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p2_2cb2cc	
m.P_prev2["2clinicalb", "2clinicalc"] <- p.p2_2cb2cc	
m.P_prev2["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev2["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p2_3cb3cc 	
m.P_prev2["3clinicalb", "3clinicalc"] <- p.p2_3cb3cc	
m.P_prev2["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev2["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p2_4cb4cc 	
m.P_prev2["4clinicalb", "4clinicalc"] <- p.p2_4cb4cc	
m.P_prev2["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev2["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p2_1ccNn	
m.P_prev2["1clinicalc", "Normalhpn"] <- p.p2_1ccNn	
m.P_prev2["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev2["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p2_2ccNn	
m.P_prev2["2clinicalc", "Normalhpn"] <- p.p2_2ccNn	
m.P_prev2["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev2["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p2_3ccNn 	
m.P_prev2["3clinicalc", "Normalhpn"] <- p.p2_3ccNn
m.P_prev2["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev2["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p2_4ccNn 	
m.P_prev2["4clinicalc", "Normalhpn"] <- p.p2_4ccNn	
m.P_prev2["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev2["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev2)


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
m.P_prev3["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p3DnDp - p.p3DnIn -  p.DnD -  p.p3_Dn1p
m.P_prev3["Dysplasiahpn", "Dysplasiahpp"]    <- p.p3DnDp
m.P_prev3["Dysplasiahpn", "Intestinalhpn"]    <- p.p3DnIn
m.P_prev3["Dysplasiahpn", "1preclinical"]    <- p.p3_Dn1p
m.P_prev3["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev3["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p3DpDn - p.p3DpIp - p.DpD - p.p3_Dp1p
m.P_prev3["Dysplasiahpp", "Dysplasiahpn"]    <- p.p3DpDn
m.P_prev3["Dysplasiahpp", "Intestinalhpp"]    <- p.p3DpIp
m.P_prev3["Dysplasiahpp", "1preclinical"]    <- p.p3_Dp1p
m.P_prev3["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev3["1preclinical", "1preclinical"]    <- 1 - p.p3_1p1ca - p.p3_1p2p - p.1pD	
m.P_prev3["1preclinical", "1clinicala"]    <- p.p3_1p1ca	
m.P_prev3["1preclinical", "2preclinical"]    <- p.p3_1p2p	
m.P_prev3["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev3["2preclinical", "2preclinical"]    <- 1 - p.p3_2p2ca - p.p3_2p3p - p.2pD	
m.P_prev3["2preclinical", "2clinicala"]    <- p.p3_2p2ca	
m.P_prev3["2preclinical", "3preclinical"]    <- p.p3_2p3p	
m.P_prev3["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev3["3preclinical", "3preclinical"]    <- 1 - p.p3_3p3ca - p.p3_3p4p - p.3pD	
m.P_prev3["3preclinical", "3clinicala"]    <- p.p3_3p3ca	
m.P_prev3["3preclinical", "4preclinical"]    <- p.p3_3p4p	
m.P_prev3["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev3["4preclinical", "4preclinical"]    <- 1 - p.p3_4p4ca - p.4pD	
m.P_prev3["4preclinical", "4clinicala"]    <- p.p3_4p4ca	
m.P_prev3["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev3["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p3_1ca1cb	
m.P_prev3["1clinicala", "1clinicalb"] <- p.p3_1ca1cb	
m.P_prev3["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev3["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p3_2ca2cb	
m.P_prev3["2clinicala", "2clinicalb"] <- p.p3_2ca2cb	
m.P_prev3["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev3["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p3_3ca3cb 	
m.P_prev3["3clinicala", "3clinicalb"] <- p.p3_3ca3cb	
m.P_prev3["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev3["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p3_4ca4cb	
m.P_prev3["4clinicala", "4clinicalb"] <- p.p3_4ca4cb  	
m.P_prev3["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev3["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p3_1cb1cc	
m.P_prev3["1clinicalb", "1clinicalc"] <- p.p3_1cb1cc	
m.P_prev3["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev3["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p3_2cb2cc	
m.P_prev3["2clinicalb", "2clinicalc"] <- p.p3_2cb2cc	
m.P_prev3["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev3["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p3_3cb3cc 	
m.P_prev3["3clinicalb", "3clinicalc"] <- p.p3_3cb3cc	
m.P_prev3["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev3["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p3_4cb4cc 	
m.P_prev3["4clinicalb", "4clinicalc"] <- p.p3_4cb4cc	
m.P_prev3["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev3["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p3_1ccNn	
m.P_prev3["1clinicalc", "Normalhpn"] <- p.p3_1ccNn	
m.P_prev3["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev3["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p3_2ccNn	
m.P_prev3["2clinicalc", "Normalhpn"] <- p.p3_2ccNn	
m.P_prev3["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev3["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p3_3ccNn 	
m.P_prev3["3clinicalc", "Normalhpn"] <- p.p3_3ccNn
m.P_prev3["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev3["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p3_4ccNn 	
m.P_prev3["4clinicalc", "Normalhpn"] <- p.p3_4ccNn	
m.P_prev3["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev3["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev3)


# fill in the transition probability matrix with prevention strategy 4

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
m.P_prev4["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p4DnDp - p.p4DnIn -  p.DnD - p.p4_Dn1p
m.P_prev4["Dysplasiahpn", "Dysplasiahpp"]    <- p.p4DnDp
m.P_prev4["Dysplasiahpn", "Intestinalhpn"]    <- p.p4DnIn
m.P_prev4["Dysplasiahpn", "1preclinical"]    <- p.p4_Dn1p
m.P_prev4["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev4["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p4DpDn - p.p4DpIp - p.DpD - p.p4_Dp1p
m.P_prev4["Dysplasiahpp", "Dysplasiahpn"]    <- p.p4DpDn
m.P_prev4["Dysplasiahpp", "Intestinalhpp"]    <- p.p4DpIp
m.P_prev4["Dysplasiahpp", "1preclinical"]    <- p.p4_Dp1p
m.P_prev4["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev4["1preclinical", "1preclinical"]    <- 1 - p.p4_1p1ca - p.p4_1p2p - p.1pD	
m.P_prev4["1preclinical", "1clinicala"]    <- p.p4_1p1ca	
m.P_prev4["1preclinical", "2preclinical"]    <- p.p4_1p2p	
m.P_prev4["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev4["2preclinical", "2preclinical"]    <- 1 - p.p4_2p2ca - p.p4_2p3p - p.2pD	
m.P_prev4["2preclinical", "2clinicala"]    <- p.p4_2p2ca	
m.P_prev4["2preclinical", "3preclinical"]    <- p.p4_2p3p	
m.P_prev4["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev4["3preclinical", "3preclinical"]    <- 1 - p.p4_3p3ca - p.p4_3p4p - p.3pD	
m.P_prev4["3preclinical", "3clinicala"]    <- p.p4_3p3ca	
m.P_prev4["3preclinical", "4preclinical"]    <- p.p4_3p4p	
m.P_prev4["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev4["4preclinical", "4preclinical"]    <- 1 - p.p4_4p4ca - p.4pD	
m.P_prev4["4preclinical", "4clinicala"]    <- p.p4_4p4ca	
m.P_prev4["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev4["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p4_1ca1cb	
m.P_prev4["1clinicala", "1clinicalb"] <- p.p4_1ca1cb	
m.P_prev4["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev4["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p4_2ca2cb	
m.P_prev4["2clinicala", "2clinicalb"] <- p.p4_2ca2cb	
m.P_prev4["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev4["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p4_3ca3cb 	
m.P_prev4["3clinicala", "3clinicalb"] <- p.p4_3ca3cb	
m.P_prev4["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev4["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p4_4ca4cb	
m.P_prev4["4clinicala", "4clinicalb"] <- p.p4_4ca4cb  	
m.P_prev4["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev4["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p4_1cb1cc	
m.P_prev4["1clinicalb", "1clinicalc"] <- p.p4_1cb1cc	
m.P_prev4["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev4["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p4_2cb2cc	
m.P_prev4["2clinicalb", "2clinicalc"] <- p.p4_2cb2cc	
m.P_prev4["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev4["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p4_3cb3cc 	
m.P_prev4["3clinicalb", "3clinicalc"] <- p.p4_3cb3cc	
m.P_prev4["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev4["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p4_4cb4cc 	
m.P_prev4["4clinicalb", "4clinicalc"] <- p.p4_4cb4cc	
m.P_prev4["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev4["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p4_1ccNn	
m.P_prev4["1clinicalc", "Normalhpn"] <- p.p4_1ccNn	
m.P_prev4["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev4["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p4_2ccNn	
m.P_prev4["2clinicalc", "Normalhpn"] <- p.p4_2ccNn	
m.P_prev4["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev4["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p4_3ccNn 	
m.P_prev4["3clinicalc", "Normalhpn"] <- p.p4_3ccNn
m.P_prev4["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev4["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p4_4ccNn 	
m.P_prev4["4clinicalc", "Normalhpn"] <- p.p4_4ccNn	
m.P_prev4["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev4["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev4)


# fill in the transition probability matrix witho prevention strategy 5
# Transition probabilities of the second part of the model (gastric cancer) are assumed as the ones with no prevention

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
m.P_prev5["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p5DnDp - p.p5DnIn -  p.DnD - p.p5_Dn1p
m.P_prev5["Dysplasiahpn", "Dysplasiahpp"]    <- p.p5DnDp
m.P_prev5["Dysplasiahpn", "Intestinalhpn"]    <- p.p5DnIn
m.P_prev5["Dysplasiahpn", "1preclinical"]    <- p.p5_Dn1p
m.P_prev5["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev5["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p5DpDn - p.p5DpIp - p.DpD - p.p5_Dp1p
m.P_prev5["Dysplasiahpp", "Dysplasiahpn"]    <- p.p5DpDn
m.P_prev5["Dysplasiahpp", "Intestinalhpp"]    <- p.p5DpIp
m.P_prev5["Dysplasiahpp", "1preclinical"]    <- p.p5_Dp1p
m.P_prev5["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev5["1preclinical", "1preclinical"]    <- 1 - p.p5_1p1ca - p.p5_1p2p - p.1pD	
m.P_prev5["1preclinical", "1clinicala"]    <- p.p5_1p1ca	
m.P_prev5["1preclinical", "2preclinical"]    <- p.p5_1p2p	
m.P_prev5["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev5["2preclinical", "2preclinical"]    <- 1 - p.p5_2p2ca - p.p5_2p3p - p.2pD	
m.P_prev5["2preclinical", "2clinicala"]    <- p.p5_2p2ca	
m.P_prev5["2preclinical", "3preclinical"]    <- p.p5_2p3p	
m.P_prev5["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev5["3preclinical", "3preclinical"]    <- 1 - p.p5_3p3ca - p.p5_3p4p - p.3pD	
m.P_prev5["3preclinical", "3clinicala"]    <- p.p5_3p3ca	
m.P_prev5["3preclinical", "4preclinical"]    <- p.p5_3p4p	
m.P_prev5["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev5["4preclinical", "4preclinical"]    <- 1 - p.p5_4p4ca - p.4pD	
m.P_prev5["4preclinical", "4clinicala"]    <- p.p5_4p4ca	
m.P_prev5["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev5["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p5_1ca1cb	
m.P_prev5["1clinicala", "1clinicalb"] <- p.p5_1ca1cb	
m.P_prev5["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev5["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p5_2ca2cb	
m.P_prev5["2clinicala", "2clinicalb"] <- p.p5_2ca2cb	
m.P_prev5["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev5["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p5_3ca3cb 	
m.P_prev5["3clinicala", "3clinicalb"] <- p.p5_3ca3cb	
m.P_prev5["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev5["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p5_4ca4cb	
m.P_prev5["4clinicala", "4clinicalb"] <- p.p5_4ca4cb  	
m.P_prev5["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev5["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p5_1cb1cc	
m.P_prev5["1clinicalb", "1clinicalc"] <- p.p5_1cb1cc	
m.P_prev5["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev5["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p5_2cb2cc	
m.P_prev5["2clinicalb", "2clinicalc"] <- p.p5_2cb2cc	
m.P_prev5["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev5["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p5_3cb3cc 	
m.P_prev5["3clinicalb", "3clinicalc"] <- p.p5_3cb3cc	
m.P_prev5["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev5["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p5_4cb4cc 	
m.P_prev5["4clinicalb", "4clinicalc"] <- p.p5_4cb4cc	
m.P_prev5["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev5["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p5_1ccNn	
m.P_prev5["1clinicalc", "Normalhpn"] <- p.p5_1ccNn	
m.P_prev5["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev5["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p5_2ccNn	
m.P_prev5["2clinicalc", "Normalhpn"] <- p.p5_2ccNn	
m.P_prev5["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev5["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p5_3ccNn 	
m.P_prev5["3clinicalc", "Normalhpn"] <- p.p5_3ccNn
m.P_prev5["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev5["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p5_4ccNn 	
m.P_prev5["4clinicalc", "Normalhpn"] <- p.p5_4ccNn	
m.P_prev5["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev5["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev5)


# fill in the transition probability matrix with prevention strategy 6
# Transition probabilities of the second part of the model (gastric cancer) are assumed as the ones with no prevention

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
m.P_prev6["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p6DnDp - p.p6DnIn -  p.DnD - p.p6_Dn1p
m.P_prev6["Dysplasiahpn", "Dysplasiahpp"]    <- p.p6DnDp
m.P_prev6["Dysplasiahpn", "Intestinalhpn"]    <- p.p6DnIn
m.P_prev6["Dysplasiahpn", "1preclinical"]    <- p.p6_Dn1p
m.P_prev6["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev6["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p6DpDn - p.p6DpIp - p.DpD - p.p6_Dp1p
m.P_prev6["Dysplasiahpp", "Dysplasiahpn"]    <- p.p6DpDn
m.P_prev6["Dysplasiahpp", "Intestinalhpp"]    <- p.p6DpIp
m.P_prev6["Dysplasiahpp", "1preclinical"]    <- p.p6_Dp1p
m.P_prev6["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev6["1preclinical", "1preclinical"]    <- 1 - p.p6_1p1ca - p.p6_1p2p - p.1pD	
m.P_prev6["1preclinical", "1clinicala"]    <- p.p6_1p1ca	
m.P_prev6["1preclinical", "2preclinical"]    <- p.p6_1p2p	
m.P_prev6["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev6["2preclinical", "2preclinical"]    <- 1 - p.p6_2p2ca - p.p6_2p3p - p.2pD	
m.P_prev6["2preclinical", "2clinicala"]    <- p.p6_2p2ca	
m.P_prev6["2preclinical", "3preclinical"]    <- p.p6_2p3p	
m.P_prev6["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev6["3preclinical", "3preclinical"]    <- 1 - p.p6_3p3ca - p.p6_3p4p - p.3pD	
m.P_prev6["3preclinical", "3clinicala"]    <- p.p6_3p3ca	
m.P_prev6["3preclinical", "4preclinical"]    <- p.p6_3p4p	
m.P_prev6["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev6["4preclinical", "4preclinical"]    <- 1 - p.p6_4p4ca - p.4pD	
m.P_prev6["4preclinical", "4clinicala"]    <- p.p6_4p4ca	
m.P_prev6["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev6["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p6_1ca1cb	
m.P_prev6["1clinicala", "1clinicalb"] <- p.p6_1ca1cb	
m.P_prev6["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev6["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p6_2ca2cb	
m.P_prev6["2clinicala", "2clinicalb"] <- p.p6_2ca2cb	
m.P_prev6["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev6["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p6_3ca3cb 	
m.P_prev6["3clinicala", "3clinicalb"] <- p.p6_3ca3cb	
m.P_prev6["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev6["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p6_4ca4cb	
m.P_prev6["4clinicala", "4clinicalb"] <- p.p6_4ca4cb  	
m.P_prev6["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev6["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p6_1cb1cc	
m.P_prev6["1clinicalb", "1clinicalc"] <- p.p6_1cb1cc	
m.P_prev6["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev6["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p6_2cb2cc	
m.P_prev6["2clinicalb", "2clinicalc"] <- p.p6_2cb2cc	
m.P_prev6["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev6["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p6_3cb3cc 	
m.P_prev6["3clinicalb", "3clinicalc"] <- p.p6_3cb3cc	
m.P_prev6["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev6["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p6_4cb4cc 	
m.P_prev6["4clinicalb", "4clinicalc"] <- p.p6_4cb4cc	
m.P_prev6["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev6["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p6_1ccNn	
m.P_prev6["1clinicalc", "Normalhpn"] <- p.p6_1ccNn	
m.P_prev6["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev6["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p6_2ccNn	
m.P_prev6["2clinicalc", "Normalhpn"] <- p.p6_2ccNn	
m.P_prev6["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev6["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p6_3ccNn 	
m.P_prev6["3clinicalc", "Normalhpn"] <- p.p6_3ccNn
m.P_prev6["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev6["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p6_4ccNn 	
m.P_prev6["4clinicalc", "Normalhpn"] <- p.p6_4ccNn	
m.P_prev6["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev6["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev6)


# fill in the transition probability matrix witho prevention strategy 7
# Transition probabilities of the second part of the model (gastric cancer) are assumed as the ones with no prevention

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
m.P_prev7["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p7DnDp - p.p7DnIn -  p.DnD -  p.p7_Dn1p
m.P_prev7["Dysplasiahpn", "Dysplasiahpp"]    <- p.p7DnDp
m.P_prev7["Dysplasiahpn", "Intestinalhpn"]    <- p.p7DnIn
m.P_prev7["Dysplasiahpn", "1preclinical"]    <- p.p7_Dn1p
m.P_prev7["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev7["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p7DpDn - p.p7DpIp - p.DpD - p.p7_Dp1p
m.P_prev7["Dysplasiahpp", "Dysplasiahpn"]    <- p.p7DpDn
m.P_prev7["Dysplasiahpp", "Intestinalhpp"]    <- p.p7DpIp
m.P_prev7["Dysplasiahpp", "1preclinical"]    <- p.p7_Dp1p
m.P_prev7["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev7["1preclinical", "1preclinical"]    <- 1 - p.p7_1p1ca - p.p7_1p2p - p.1pD	
m.P_prev7["1preclinical", "1clinicala"]    <- p.p7_1p1ca	
m.P_prev7["1preclinical", "2preclinical"]    <- p.p7_1p2p	
m.P_prev7["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev7["2preclinical", "2preclinical"]    <- 1 - p.p7_2p2ca - p.p7_2p3p - p.2pD	
m.P_prev7["2preclinical", "2clinicala"]    <- p.p7_2p2ca	
m.P_prev7["2preclinical", "3preclinical"]    <- p.p7_2p3p	
m.P_prev7["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev7["3preclinical", "3preclinical"]    <- 1 - p.p7_3p3ca - p.p7_3p4p - p.3pD	
m.P_prev7["3preclinical", "3clinicala"]    <- p.p7_3p3ca	
m.P_prev7["3preclinical", "4preclinical"]    <- p.p7_3p4p	
m.P_prev7["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev7["4preclinical", "4preclinical"]    <- 1 - p.p7_4p4ca - p.4pD	
m.P_prev7["4preclinical", "4clinicala"]    <- p.p7_4p4ca	
m.P_prev7["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev7["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p7_1ca1cb	
m.P_prev7["1clinicala", "1clinicalb"] <- p.p7_1ca1cb	
m.P_prev7["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev7["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p7_2ca2cb	
m.P_prev7["2clinicala", "2clinicalb"] <- p.p7_2ca2cb	
m.P_prev7["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev7["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p7_3ca3cb 	
m.P_prev7["3clinicala", "3clinicalb"] <- p.p7_3ca3cb	
m.P_prev7["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev7["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p7_4ca4cb	
m.P_prev7["4clinicala", "4clinicalb"] <- p.p7_4ca4cb  	
m.P_prev7["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev7["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p7_1cb1cc	
m.P_prev7["1clinicalb", "1clinicalc"] <- p.p7_1cb1cc	
m.P_prev7["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev7["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p7_2cb2cc	
m.P_prev7["2clinicalb", "2clinicalc"] <- p.p7_2cb2cc	
m.P_prev7["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev7["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p7_3cb3cc 	
m.P_prev7["3clinicalb", "3clinicalc"] <- p.p7_3cb3cc	
m.P_prev7["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev7["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p7_4cb4cc 	
m.P_prev7["4clinicalb", "4clinicalc"] <- p.p7_4cb4cc	
m.P_prev7["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev7["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p7_1ccNn	
m.P_prev7["1clinicalc", "Normalhpn"] <- p.p7_1ccNn	
m.P_prev7["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev7["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p7_2ccNn	
m.P_prev7["2clinicalc", "Normalhpn"] <- p.p7_2ccNn	
m.P_prev7["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev7["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p7_3ccNn 	
m.P_prev7["3clinicalc", "Normalhpn"] <- p.p7_3ccNn
m.P_prev7["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev7["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p7_4ccNn 	
m.P_prev7["4clinicalc", "Normalhpn"] <- p.p7_4ccNn	
m.P_prev7["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev7["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev7)


# fill in the transition probability matrix witho prevention strategy 8
# Transition probabilities of the second part of the model (gastric cancer) are assumed as the ones with no prevention


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
m.P_prev8["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p8DnDp - p.p8DnIn -  p.DnD - p.p8_Dn1p
m.P_prev8["Dysplasiahpn", "Dysplasiahpp"]    <- p.p8DnDp
m.P_prev8["Dysplasiahpn", "Intestinalhpn"]    <- p.p8DnIn
m.P_prev8["Dysplasiahpn", "1preclinical"]    <- p.p8_Dn1p
m.P_prev8["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev8["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p8DpDn - p.p8DpIp - p.DpD - p.p8_Dp1p
m.P_prev8["Dysplasiahpp", "Dysplasiahpn"]    <- p.p8DpDn
m.P_prev8["Dysplasiahpp", "Intestinalhpp"]    <- p.p8DpIp
m.P_prev8["Dysplasiahpp", "1preclinical"]    <- p.p8_Dp1p
m.P_prev8["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev8["1preclinical", "1preclinical"]    <- 1 - p.p8_1p1ca - p.p8_1p2p - p.1pD	
m.P_prev8["1preclinical", "1clinicala"]    <- p.p8_1p1ca	
m.P_prev8["1preclinical", "2preclinical"]    <- p.p8_1p2p	
m.P_prev8["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev8["2preclinical", "2preclinical"]    <- 1 - p.p8_2p2ca - p.p8_2p3p - p.2pD	
m.P_prev8["2preclinical", "2clinicala"]    <- p.p8_2p2ca	
m.P_prev8["2preclinical", "3preclinical"]    <- p.p8_2p3p	
m.P_prev8["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev8["3preclinical", "3preclinical"]    <- 1 - p.p8_3p3ca - p.p8_3p4p - p.3pD	
m.P_prev8["3preclinical", "3clinicala"]    <- p.p8_3p3ca	
m.P_prev8["3preclinical", "4preclinical"]    <- p.p8_3p4p	
m.P_prev8["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev8["4preclinical", "4preclinical"]    <- 1 - p.p8_4p4ca - p.4pD	
m.P_prev8["4preclinical", "4clinicala"]    <- p.p8_4p4ca	
m.P_prev8["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev8["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p8_1ca1cb	
m.P_prev8["1clinicala", "1clinicalb"] <- p.p8_1ca1cb	
m.P_prev8["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev8["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p8_2ca2cb	
m.P_prev8["2clinicala", "2clinicalb"] <- p.p8_2ca2cb	
m.P_prev8["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev8["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p8_3ca3cb 	
m.P_prev8["3clinicala", "3clinicalb"] <- p.p8_3ca3cb	
m.P_prev8["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev8["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p8_4ca4cb	
m.P_prev8["4clinicala", "4clinicalb"] <- p.p8_4ca4cb  	
m.P_prev8["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev8["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p8_1cb1cc	
m.P_prev8["1clinicalb", "1clinicalc"] <- p.p8_1cb1cc	
m.P_prev8["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev8["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p8_2cb2cc	
m.P_prev8["2clinicalb", "2clinicalc"] <- p.p8_2cb2cc	
m.P_prev8["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev8["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p8_3cb3cc 	
m.P_prev8["3clinicalb", "3clinicalc"] <- p.p8_3cb3cc	
m.P_prev8["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev8["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p8_4cb4cc 	
m.P_prev8["4clinicalb", "4clinicalc"] <- p.p8_4cb4cc	
m.P_prev8["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev8["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p8_1ccNn	
m.P_prev8["1clinicalc", "Normalhpn"] <- p.p8_1ccNn	
m.P_prev8["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev8["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p8_2ccNn	
m.P_prev8["2clinicalc", "Normalhpn"] <- p.p8_2ccNn	
m.P_prev8["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev8["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p8_3ccNn 	
m.P_prev8["3clinicalc", "Normalhpn"] <- p.p8_3ccNn
m.P_prev8["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev8["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p8_4ccNn 	
m.P_prev8["4clinicalc", "Normalhpn"] <- p.p8_4ccNn	
m.P_prev8["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
m.P_prev8["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev8)


## fill in the transition probability matrix with prevention strategy 9

### From Normal helicobacter (-)
m.P_prev9["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p9NnNp - p.p9NnGn
m.P_prev9["Normalhpn", "Normalhpp"]    <- p.p9NnNp
m.P_prev9["Normalhpn", "Gastritishpn"]    <- p.p9NnGn
m.P_prev9["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev9["Normalhpp", "Normalhpp"]    <- 1 - p.p9NpGp - p.p9NpNn - p.NpD 
m.P_prev9["Normalhpp", "Gastritishpp"]    <- p.p9NpGp
m.P_prev9["Normalhpp", "Normalhpn"]    <- p.p9NpNn
m.P_prev9["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev9["Gastritishpn", "Gastritishpn"]    <- 1 - p.p9GnAn - p.p9GnGp - p.p9GnNn - p.GnD
m.P_prev9["Gastritishpn", "Atrophyhpn"]    <- p.p9GnAn
m.P_prev9["Gastritishpn", "Gastritishpp"]    <- p.p9GnGp
m.P_prev9["Gastritishpn", "Normalhpn"]    <- p.p9GnNn
m.P_prev9["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev9["Gastritishpp", "Gastritishpp"]    <- 1 - p.p9GpAp - p.p9GpGn - p.p9GpNp - p.GpD
m.P_prev9["Gastritishpp", "Atrophyhpp"]    <- p.p9GpAp
m.P_prev9["Gastritishpp", "Gastritishpn"]    <- p.p9GpGn
m.P_prev9["Gastritishpp", "Normalhpp"]    <- p.p9GpNp
m.P_prev9["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev9["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p9AnIn - p.p9AnAp - p.p9AnGn -  p.AnD
m.P_prev9["Atrophyhpn", "Intestinalhpn"]    <- p.p9AnIn
m.P_prev9["Atrophyhpn", "Atrophyhpp"]    <- p.p9AnAp
m.P_prev9["Atrophyhpn", "Gastritishpn"]    <- p.p9AnGn
m.P_prev9["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev9["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p9ApIp - p.p9ApAn - p.p9ApGp - p.ApD
m.P_prev9["Atrophyhpp", "Intestinalhpp"]    <- p.p9ApIp
m.P_prev9["Atrophyhpp", "Atrophyhpn"]    <- p.p9ApAn
m.P_prev9["Atrophyhpp", "Gastritishpp"]    <- p.p9ApGp
m.P_prev9["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev9["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p9InDn - p.p9InIp - p.p9InAn - p.InD 
m.P_prev9["Intestinalhpn", "Dysplasiahpn"]    <- p.p9InDn
m.P_prev9["Intestinalhpn", "Intestinalhpp"]    <- p.p9InIp
m.P_prev9["Intestinalhpn", "Atrophyhpn"]    <- p.p9InAn
m.P_prev9["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev9["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p9IpDp - p.p9IpIn - p.p9IpAp - p.IpD  
m.P_prev9["Intestinalhpp", "Dysplasiahpp"]    <- p.p9IpDp
m.P_prev9["Intestinalhpp", "Intestinalhpn"]    <- p.p9IpIn
m.P_prev9["Intestinalhpp", "Atrophyhpp"]    <- p.p9IpAp
m.P_prev9["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev9["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p9DnDp - p.p9DnIn -  p.p9_Dn1p -p.DnD 
m.P_prev9["Dysplasiahpn", "Dysplasiahpp"]    <- p.p9DnDp
m.P_prev9["Dysplasiahpn", "Intestinalhpn"]    <- p.p9DnIn
m.P_prev9["Dysplasiahpn", "1preclinical"]     <- p.p9_Dn1p
m.P_prev9["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev9["Dysplasiahpp", "Dysplasiahpn"]    <- 1 - p.p9DpDn - p.p9DpIp - p.p9_Dp1p - p.DpD  
m.P_prev9["Dysplasiahpp", "Dysplasiahpp"]    <- p.p9DpDn
m.P_prev9["Dysplasiahpp", "Intestinalhpp"]    <- p.p9DpIp
m.P_prev9["Dysplasiahpp", "1preclinical"]     <- p.p9_Dp1p
m.P_prev9["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev9["1preclinical", "1preclinical"]    <- 1 - p.p9_1p1ca - p.p9_1p2p - p.1pD	
m.P_prev9["1preclinical", "1clinicala"]    <- p.p9_1p1ca	
m.P_prev9["1preclinical", "2preclinical"]    <- p.p9_1p2p	
m.P_prev9["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev9["2preclinical", "2preclinical"]    <- 1 - p.p9_2p2ca - p.p9_2p3p - p.2pD	
m.P_prev9["2preclinical", "2clinicala"]    <- p.p9_2p2ca	
m.P_prev9["2preclinical", "3preclinical"]    <- p.p9_2p3p	
m.P_prev9["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev9["3preclinical", "3preclinical"]    <- 1 - p.p9_3p3ca - p.p9_3p4p - p.3pD	
m.P_prev9["3preclinical", "3clinicala"]    <- p.p9_3p3ca	
m.P_prev9["3preclinical", "4preclinical"]    <- p.p9_3p4p	
m.P_prev9["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev9["4preclinical", "4preclinical"]    <- 1 - p.p9_4p4ca - p.4pD	
m.P_prev9["4preclinical", "4clinicala"]    <- p.p9_4p4ca	
m.P_prev9["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev9["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p9_1ca1cb	
m.P_prev9["1clinicala", "1clinicalb"] <- p.p9_1ca1cb	
m.P_prev9["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev9["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p9_2ca2cb	
m.P_prev9["2clinicala", "2clinicalb"] <- p.p9_2ca2cb	
m.P_prev9["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev9["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p9_3ca3cb 	
m.P_prev9["3clinicala", "3clinicalb"] <- p.p9_3ca3cb	
m.P_prev9["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev9["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p9_4ca4cb	
m.P_prev9["4clinicala", "4clinicalb"] <- p.p9_4ca4cb  	
m.P_prev9["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev9["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p9_1cb1cc	
m.P_prev9["1clinicalb", "1clinicalc"] <- p.p9_1cb1cc	
m.P_prev9["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev9["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p9_2cb2cc	
m.P_prev9["2clinicalb", "2clinicalc"] <- p.p9_2cb2cc	
m.P_prev9["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev9["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p9_3cb3cc 	
m.P_prev9["3clinicalb", "3clinicalc"] <- p.p9_3cb3cc	
m.P_prev9["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev9["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p9_4cb4cc 	
m.P_prev9["4clinicalb", "4clinicalc"] <- p.p9_4cb4cc	
m.P_prev9["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev9["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p9_1ccNn		
m.P_prev9["1clinicalc", "Normalhpn"] <- p.p9_1ccNn	
m.P_prev9["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev9["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p9_2ccNn	
m.P_prev9["2clinicalc", "Normalhpn"] <- p.p9_2ccNn	
m.P_prev9["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev9["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p9_3ccNn 	
m.P_prev9["3clinicalc", "Normalhpn"] <- p.p9_3ccNn	
m.P_prev9["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev9["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p9_4ccNn 	
m.P_prev9["4clinicalc", "Normalhpn"] <- p.p9_4ccNn	
m.P_prev9["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead	
m.P_prev9["Dead", "Dead"] <- 1	

# check rows add up to 1	
rowSums(m.P_prev9)	


## fill transitions prevention strategy 10

### From Normal helicobacter (-)
m.P_prev10["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p10NnNp - p.p10NnGn
m.P_prev10["Normalhpn", "Normalhpp"]    <- p.p10NnNp
m.P_prev10["Normalhpn", "Gastritishpn"]    <- p.p10NnGn
m.P_prev10["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev10["Normalhpp", "Normalhpp"]    <- 1 - p.p10NpGp - p.p10NpNn - p.NpD 
m.P_prev10["Normalhpp", "Gastritishpp"]    <- p.p10NpGp
m.P_prev10["Normalhpp", "Normalhpn"]    <- p.p10NpNn
m.P_prev10["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev10["Gastritishpn", "Gastritishpn"]    <- 1 - p.p10GnAn - p.p10GnGp - p.p10GnNn - p.GnD
m.P_prev10["Gastritishpn", "Atrophyhpn"]    <- p.p10GnAn
m.P_prev10["Gastritishpn", "Gastritishpp"]    <- p.p10GnGp
m.P_prev10["Gastritishpn", "Normalhpn"]    <- p.p10GnNn
m.P_prev10["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev10["Gastritishpp", "Gastritishpp"]    <- 1 - p.p10GpAp - p.p10GpGn - p.p10GpNp - p.GpD
m.P_prev10["Gastritishpp", "Atrophyhpp"]    <- p.p10GpAp
m.P_prev10["Gastritishpp", "Gastritishpn"]    <- p.p10GpGn
m.P_prev10["Gastritishpp", "Normalhpp"]    <- p.p10GpNp
m.P_prev10["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev10["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p10AnIn - p.p10AnAp - p.p10AnGn -  p.AnD
m.P_prev10["Atrophyhpn", "Intestinalhpn"]    <- p.p10AnIn
m.P_prev10["Atrophyhpn", "Atrophyhpp"]    <- p.p10AnAp
m.P_prev10["Atrophyhpn", "Gastritishpn"]    <- p.p10AnGn
m.P_prev10["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev10["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p10ApIp - p.p10ApAn - p.p10ApGp - p.ApD
m.P_prev10["Atrophyhpp", "Intestinalhpp"]    <- p.p10ApIp
m.P_prev10["Atrophyhpp", "Atrophyhpn"]    <- p.p10ApAn
m.P_prev10["Atrophyhpp", "Gastritishpp"]    <- p.p10ApGp
m.P_prev10["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev10["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p10InDn - p.p10InIp - p.p10InAn - p.InD 
m.P_prev10["Intestinalhpn", "Dysplasiahpn"]    <- p.p10InDn
m.P_prev10["Intestinalhpn", "Intestinalhpp"]    <- p.p10InIp
m.P_prev10["Intestinalhpn", "Atrophyhpn"]    <- p.p10InAn
m.P_prev10["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev10["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p10IpDp - p.p10IpIn - p.p10IpAp - p.IpD  
m.P_prev10["Intestinalhpp", "Dysplasiahpp"]    <- p.p10IpDp
m.P_prev10["Intestinalhpp", "Intestinalhpn"]    <- p.p10IpIn
m.P_prev10["Intestinalhpp", "Atrophyhpp"]    <- p.p10IpAp
m.P_prev10["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev10["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p10DnDp - p.p10DnIn -  p.p10_Dn1p -p.DnD 
m.P_prev10["Dysplasiahpn", "Dysplasiahpp"]    <- p.p10DnDp
m.P_prev10["Dysplasiahpn", "Intestinalhpn"]    <- p.p10DnIn
m.P_prev10["Dysplasiahpn", "1preclinical"]     <- p.p10_Dn1p
m.P_prev10["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev10["Dysplasiahpp", "Dysplasiahpn"]    <- 1 - p.p10DpDn - p.p10DpIp - p.p10_Dp1p - p.DpD  
m.P_prev10["Dysplasiahpp", "Dysplasiahpp"]    <- p.p10DpDn
m.P_prev10["Dysplasiahpp", "Intestinalhpp"]    <- p.p10DpIp
m.P_prev10["Dysplasiahpp", "1preclinical"]     <- p.p10_Dp1p
m.P_prev10["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev10["1preclinical", "1preclinical"]    <- 1 - p.p10_1p1ca - p.p10_1p2p - p.1pD	
m.P_prev10["1preclinical", "1clinicala"]    <- p.p10_1p1ca	
m.P_prev10["1preclinical", "2preclinical"]    <- p.p10_1p2p	
m.P_prev10["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev10["2preclinical", "2preclinical"]    <- 1 - p.p10_2p2ca - p.p10_2p3p - p.2pD	
m.P_prev10["2preclinical", "2clinicala"]    <- p.p10_2p2ca	
m.P_prev10["2preclinical", "3preclinical"]    <- p.p10_2p3p	
m.P_prev10["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev10["3preclinical", "3preclinical"]    <- 1 - p.p10_3p3ca - p.p10_3p4p - p.3pD	
m.P_prev10["3preclinical", "3clinicala"]    <- p.p10_3p3ca	
m.P_prev10["3preclinical", "4preclinical"]    <- p.p10_3p4p	
m.P_prev10["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev10["4preclinical", "4preclinical"]    <- 1 - p.p10_4p4ca - p.4pD	
m.P_prev10["4preclinical", "4clinicala"]    <- p.p10_4p4ca	
m.P_prev10["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev10["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p10_1ca1cb	
m.P_prev10["1clinicala", "1clinicalb"] <- p.p10_1ca1cb	
m.P_prev10["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev10["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p10_2ca2cb	
m.P_prev10["2clinicala", "2clinicalb"] <- p.p10_2ca2cb	
m.P_prev10["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev10["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p10_3ca3cb 	
m.P_prev10["3clinicala", "3clinicalb"] <- p.p10_3ca3cb	
m.P_prev10["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev10["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p10_4ca4cb	
m.P_prev10["4clinicala", "4clinicalb"] <- p.p10_4ca4cb  	
m.P_prev10["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev10["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p10_1cb1cc	
m.P_prev10["1clinicalb", "1clinicalc"] <- p.p10_1cb1cc	
m.P_prev10["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev10["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p10_2cb2cc	
m.P_prev10["2clinicalb", "2clinicalc"] <- p.p10_2cb2cc	
m.P_prev10["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev10["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p10_3cb3cc 	
m.P_prev10["3clinicalb", "3clinicalc"] <- p.p10_3cb3cc	
m.P_prev10["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev10["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p10_4cb4cc 	
m.P_prev10["4clinicalb", "4clinicalc"] <- p.p10_4cb4cc	
m.P_prev10["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev10["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p10_1ccNn		
m.P_prev10["1clinicalc", "Normalhpn"] <- p.p10_1ccNn	
m.P_prev10["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev10["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p10_2ccNn	
m.P_prev10["2clinicalc", "Normalhpn"] <- p.p10_2ccNn	
m.P_prev10["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev10["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p10_3ccNn 	
m.P_prev10["3clinicalc", "Normalhpn"] <- p.p10_3ccNn	
m.P_prev10["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev10["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p10_4ccNn 	
m.P_prev10["4clinicalc", "Normalhpn"] <- p.p10_4ccNn	
m.P_prev10["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead	
m.P_prev10["Dead", "Dead"] <- 1	

# check rows add up to 1	
rowSums(m.P_prev10)	

# fill transitions prevention strategy 11
### From Normal helicobacter (-)
m.P_prev11["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p11NnNp - p.p11NnGn
m.P_prev11["Normalhpn", "Normalhpp"]    <- p.p11NnNp
m.P_prev11["Normalhpn", "Gastritishpn"]    <- p.p11NnGn
m.P_prev11["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
m.P_prev11["Normalhpp", "Normalhpp"]    <- 1 - p.p11NpGp - p.p11NpNn - p.NpD 
m.P_prev11["Normalhpp", "Gastritishpp"]    <- p.p11NpGp
m.P_prev11["Normalhpp", "Normalhpn"]    <- p.p11NpNn
m.P_prev11["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
m.P_prev11["Gastritishpn", "Gastritishpn"]    <- 1 - p.p11GnAn - p.p11GnGp - p.p11GnNn - p.GnD
m.P_prev11["Gastritishpn", "Atrophyhpn"]    <- p.p11GnAn
m.P_prev11["Gastritishpn", "Gastritishpp"]    <- p.p11GnGp
m.P_prev11["Gastritishpn", "Normalhpn"]    <- p.p11GnNn
m.P_prev11["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
m.P_prev11["Gastritishpp", "Gastritishpp"]    <- 1 - p.p11GpAp - p.p11GpGn - p.p11GpNp - p.GpD
m.P_prev11["Gastritishpp", "Atrophyhpp"]    <- p.p11GpAp
m.P_prev11["Gastritishpp", "Gastritishpn"]    <- p.p11GpGn
m.P_prev11["Gastritishpp", "Normalhpp"]    <- p.p11GpNp
m.P_prev11["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
m.P_prev11["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p11AnIn - p.p11AnAp - p.p11AnGn -  p.AnD
m.P_prev11["Atrophyhpn", "Intestinalhpn"]    <- p.p11AnIn
m.P_prev11["Atrophyhpn", "Atrophyhpp"]    <- p.p11AnAp
m.P_prev11["Atrophyhpn", "Gastritishpn"]    <- p.p11AnGn
m.P_prev11["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
m.P_prev11["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p11ApIp - p.p11ApAn - p.p11ApGp - p.ApD
m.P_prev11["Atrophyhpp", "Intestinalhpp"]    <- p.p11ApIp
m.P_prev11["Atrophyhpp", "Atrophyhpn"]    <- p.p11ApAn
m.P_prev11["Atrophyhpp", "Gastritishpp"]    <- p.p11ApGp
m.P_prev11["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
m.P_prev11["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p11InDn - p.p11InIp - p.p11InAn - p.InD 
m.P_prev11["Intestinalhpn", "Dysplasiahpn"]    <- p.p11InDn
m.P_prev11["Intestinalhpn", "Intestinalhpp"]    <- p.p11InIp
m.P_prev11["Intestinalhpn", "Atrophyhpn"]    <- p.p11InAn
m.P_prev11["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
m.P_prev11["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p11IpDp - p.p11IpIn - p.p11IpAp - p.IpD  
m.P_prev11["Intestinalhpp", "Dysplasiahpp"]    <- p.p11IpDp
m.P_prev11["Intestinalhpp", "Intestinalhpn"]    <- p.p11IpIn
m.P_prev11["Intestinalhpp", "Atrophyhpp"]    <- p.p11IpAp
m.P_prev11["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
m.P_prev11["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p11DnDp - p.p11DnIn -  p.p11_Dn1p -p.DnD 
m.P_prev11["Dysplasiahpn", "Dysplasiahpp"]    <- p.p11DnDp
m.P_prev11["Dysplasiahpn", "Intestinalhpn"]    <- p.p11DnIn
m.P_prev11["Dysplasiahpn", "1preclinical"]     <- p.p11_Dn1p
m.P_prev11["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
m.P_prev11["Dysplasiahpp", "Dysplasiahpn"]    <- 1 - p.p11DpDn - p.p11DpIp - p.p11_Dp1p - p.DpD  
m.P_prev11["Dysplasiahpp", "Dysplasiahpp"]    <- p.p11DpDn
m.P_prev11["Dysplasiahpp", "Intestinalhpp"]    <- p.p11DpIp
m.P_prev11["Dysplasiahpp", "1preclinical"]     <- p.p11_Dp1p
m.P_prev11["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
m.P_prev11["1preclinical", "1preclinical"]    <- 1 - p.p11_1p1ca - p.p11_1p2p - p.1pD	
m.P_prev11["1preclinical", "1clinicala"]    <- p.p11_1p1ca	
m.P_prev11["1preclinical", "2preclinical"]    <- p.p11_1p2p	
m.P_prev11["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prev11["2preclinical", "2preclinical"]    <- 1 - p.p11_2p2ca - p.p11_2p3p - p.2pD	
m.P_prev11["2preclinical", "2clinicala"]    <- p.p11_2p2ca	
m.P_prev11["2preclinical", "3preclinical"]    <- p.p11_2p3p	
m.P_prev11["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prev11["3preclinical", "3preclinical"]    <- 1 - p.p11_3p3ca - p.p11_3p4p - p.3pD	
m.P_prev11["3preclinical", "3clinicala"]    <- p.p11_3p3ca	
m.P_prev11["3preclinical", "4preclinical"]    <- p.p11_3p4p	
m.P_prev11["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prev11["4preclinical", "4preclinical"]    <- 1 - p.p11_4p4ca - p.4pD	
m.P_prev11["4preclinical", "4clinicala"]    <- p.p11_4p4ca	
m.P_prev11["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prev11["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p11_1ca1cb	
m.P_prev11["1clinicala", "1clinicalb"] <- p.p11_1ca1cb	
m.P_prev11["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prev11["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p11_2ca2cb	
m.P_prev11["2clinicala", "2clinicalb"] <- p.p11_2ca2cb	
m.P_prev11["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
m.P_prev11["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p11_3ca3cb 	
m.P_prev11["3clinicala", "3clinicalb"] <- p.p11_3ca3cb	
m.P_prev11["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
m.P_prev11["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p11_4ca4cb	
m.P_prev11["4clinicala", "4clinicalb"] <- p.p11_4ca4cb  	
m.P_prev11["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prev11["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p11_1cb1cc	
m.P_prev11["1clinicalb", "1clinicalc"] <- p.p11_1cb1cc	
m.P_prev11["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prev11["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p11_2cb2cc	
m.P_prev11["2clinicalb", "2clinicalc"] <- p.p11_2cb2cc	
m.P_prev11["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prev11["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p11_3cb3cc 	
m.P_prev11["3clinicalb", "3clinicalc"] <- p.p11_3cb3cc	
m.P_prev11["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prev11["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p11_4cb4cc 	
m.P_prev11["4clinicalb", "4clinicalc"] <- p.p11_4cb4cc	
m.P_prev11["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prev11["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p11_1ccNn		
m.P_prev11["1clinicalc", "Normalhpn"] <- p.p11_1ccNn	
m.P_prev11["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prev11["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p11_2ccNn	
m.P_prev11["2clinicalc", "Normalhpn"] <- p.p11_2ccNn	
m.P_prev11["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prev11["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p11_3ccNn 	
m.P_prev11["3clinicalc", "Normalhpn"] <- p.p11_3ccNn	
m.P_prev11["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prev11["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p11_4ccNn 	
m.P_prev11["4clinicalc", "Normalhpn"] <- p.p11_4ccNn	
m.P_prev11["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead	
m.P_prev11["Dead", "Dead"] <- 1	

# check rows add up to 1	
rowSums(m.P_prev11)	


# double check rows add up to 1
rowSums(m.P_prev)
rowSums(m.P_prev2)
rowSums(m.P_prev3)
rowSums(m.P_prev4)
rowSums(m.P_prev5)
rowSums(m.P_prev6)
rowSums(m.P_prev7)
rowSums(m.P_prev8)
rowSums(m.P_prev9)
rowSums(m.P_prev10)
rowSums(m.P_prev11)



#### 04 Run Markov model ####
for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_prev[t + 1, ] <- t(m.M_no_prev[t, ]) %*% m.P_noprev # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev[t + 1, ]    <- t(m.M_prev[t, ])    %*% m.P_prev   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev2[t + 1, ]    <- t(m.M_prev2[t, ])    %*% m.P_prev2   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev3[t + 1, ]    <- t(m.M_prev3[t, ])    %*% m.P_prev3   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev4[t + 1, ]    <- t(m.M_prev4[t, ])    %*% m.P_prev4   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev5[t + 1, ]    <- t(m.M_prev5[t, ])    %*% m.P_prev5   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev6[t + 1, ]    <- t(m.M_prev6[t, ])    %*% m.P_prev6   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev7[t + 1, ]    <- t(m.M_prev7[t, ])    %*% m.P_prev7   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev8[t + 1, ]    <- t(m.M_prev8[t, ])    %*% m.P_prev8   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev9[t + 1, ]    <- t(m.M_prev9[t, ])    %*% m.P_prev9   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev10[t + 1, ]    <- t(m.M_prev10[t, ])    %*% m.P_prev10   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev11[t + 1, ]    <- t(m.M_prev11[t, ])    %*% m.P_prev11   # estimate the Markov trace for cycle the next cycle (t + 1)
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
head(m.M_prev9)  # show the first 6 lines of the matrix with prevention strategy 1
head(m.M_prev10)  # show the first 6 lines of the matrix with prevention strategy 2
head(m.M_prev11)  # show the first 6 lines of the matrix with prevention strategy 3

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
v.os_prev9 <- 1 - m.M_prev9[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 9
v.os_prev10 <- 1 - m.M_prev10[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 10
v.os_prev11 <- 1 - m.M_prev11[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 11



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

# Plot OS prevention strategy 2 

plot(0:n.t, v.os_prev2, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 2")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 3 

plot(0:n.t, v.os_prev3, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 3")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 4 

plot(0:n.t, v.os_prev4, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 4")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 5 

plot(0:n.t, v.os_prev5, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 5")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 6

plot(0:n.t, v.os_prev6, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 6")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 7 

plot(0:n.t, v.os_prev7, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 7")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 


# Plot OS prevention strategy 8 

plot(0:n.t, v.os_prev8, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 8")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

# Plot OS prevention strategy 9

plot(0:n.t, v.os_prev9, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 9")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

#Plot OS prevention strategy 10 

plot(0:n.t, v.os_prev10, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 10")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

#Plot OS prevention strategy 11

plot(0:n.t, v.os_prev11, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention Strategy 11")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid


#### 05.2.1 Life Expectancy (LE) #####
#Without prevention
v.le <- sum(v.os_no_prev)                       # summing probablity of OS over time  (i.e. life expectancy)
v.le1 <- sum(v.os_prev)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 1
v.le2 <- sum(v.os_prev2)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 2
v.le3 <- sum(v.os_prev3)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 3
v.le4 <- sum(v.os_prev4)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 4
v.le5 <- sum(v.os_prev5)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 5
v.le6 <- sum(v.os_prev6)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 6
v.le7 <- sum(v.os_prev7)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 7
v.le8 <- sum(v.os_prev8)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 8
v.le9 <- sum(v.os_prev9)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 9
v.le10 <- sum(v.os_prev10)                     # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 10
v.le11 <- sum(v.os_prev11)                     # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 11


#### 05.3 Disease Prevalence 
## H. Pylori prevalence without prevention#####
v.preva <- rowSums(m.M_prev[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_no_prev
plot(v.preva,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence")

#### 05.3 H. Pylori prevalence with prevention strategy 1#####
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

## 05.3 Cancer Stages with no No prevention
v.preva_no_prev <- rowSums(m.M_no_prev[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_no_prev
plot(v.preva_no_prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence with no prevention")

## 05.3 Cancer Stages with prevention
# May be necessary to add the ones for combined strategies (5-8)

##Prevention Strategy 9 
v.preva_prev9 <- rowSums(m.M_prev9[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prev9
plot(v.preva_prev9,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 10 
v.preva_prev10 <- rowSums(m.M_prev10[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prev10
plot(v.preva_prev10,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 11 
v.preva_prev11 <- rowSums(m.M_prev11[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prev11
plot(v.preva_prev11,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

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

#### 05.4 Proportion of sick in Stage 1 - No Prevention
v.prop.1c_noprev <- rowSums(m.M_no_prev[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_no_prev
plot(0:n.t, v.prop.1c_noprev,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 9 
v.prop.1c_prev9 <- rowSums(m.M_prev9[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_prev9
plot(0:n.t, v.prop.1c_prev9,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 10 
v.prop.1c_prev10 <- rowSums(m.M_prev10[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_prev10
plot(0:n.t, v.prop.1c_prev10,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 11 
v.prop.1c_prev11 <- rowSums(m.M_prev11[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_prev11
plot(0:n.t, v.prop.1c_prev11,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 06 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_prev <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev2    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev3    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev4    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev5    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev6    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev7    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev8    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev9    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev10    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prev11    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)

v.c_no_prev <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev    <- c(c.Nn + c.EDAUreasa, c.Np + c.EDAUreasa, c.Gn + c.EDAUreasa, c.Gp+ c.EDAUreasa, c.An + c.EDAUreasa, c.Ap + c.EDAUreasa, c.In + c.EDAUreasa, c.Ip + c.EDAUreasa, c.Dn + c.EDAUreasa, c.Dp + c.EDAUreasa, c.1p + c.EDAUreasa, c.2p + c.EDAUreasa, c.3p + c.EDAUreasa, c.4p + c.EDAUreasa, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev2    <- c(c.Nn + c.UreaAire, c.Np + c.UreaAire, c.Gn + c.UreaAire, c.Gp + c.UreaAire, c.An + c.UreaAire, c.Ap + c.UreaAire, c.In + c.UreaAire, c.Ip + c.UreaAire, c.Dn + c.UreaAire, c.Dp + c.UreaAire, c.1p + c.UreaAire, c.2p + c.UreaAire, c.3p + c.UreaAire, c.4p + c.UreaAire, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev3    <- c(c.Nn + c.Serología, c.Np + c.Serología, c.Gn +  c.Serología, c.Gp + c.Serología, c.An +  c.Serología, c.Ap + c.Serología, c.In + c.Serología, c.Ip + c.Serología, c.Dn + c.Serología, c.Dp + c.Serología, c.1p + c.Serología, c.2p + c.Serología, c.3p + c.Serología, c.4p + c.Serología, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev4    <- c(c.Nn +  c.Antígenofecal, c.Np + c.Antígenofecal, c.Gn + c.Antígenofecal, c.Gp + c.Antígenofecal, c.An + c.Antígenofecal, c.Ap + c.Antígenofecal, c.In +  c.Antígenofecal, c.Ip + c.Antígenofecal, c.Dn + c.Antígenofecal, c.Dp + c.Antígenofecal, c.1p + c.Antígenofecal, c.2p + c.Antígenofecal, c.3p + c.Antígenofecal, c.4p + c.Antígenofecal, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev5    <- c(c.Nn + c.EDABp, c.Np + c.EDABp, c.Gn + c.EDABp, c.Gp + c.EDABp, c.An + c.EDABp, c.Ap + c.EDABp, c.In + c.EDABp, c.Ip + c.EDABp, c.Dn + c.EDABp, c.Dp + c.EDABp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev6    <- c(c.Nn + c.Serología + c.EDA, c.Np + c.Serología + c.EDA, c.Gn +  c.Serología + c.EDA, c.Gp + c.Serología + c.EDA, c.An +  c.Serología + c.EDA, c.Ap + c.Serología + c.EDA, c.In + c.Serología + c.EDA, c.Ip + c.Serología + c.EDA, c.Dn +  c.Serología + c.EDA, c.Dp + c.Serología + c.EDA, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev7    <- c(c.Nn +  c.UreaAire + c.PS + c.EDA, c.Np + c.UreaAire + c.PS + c.EDA, c.Gn + c.UreaAire + c.PS + c.EDA, c.Gp + c.UreaAire + c.PS + c.EDA, c.An + c.UreaAire + c.PS + c.EDA, c.Ap + c.UreaAire + c.PS + c.EDA, c.In + c.UreaAire + c.PS + c.EDA, c.Ip + c.UreaAire + c.PS + c.EDA, c.Dn + c.UreaAire + c.PS + c.EDA, c.Dp + c.UreaAire + c.PS + c.EDA, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev8    <- c(c.Nn +  c.Antígenofecal + c.PS + c.EDA, c.Np + c.Antígenofecal + c.PS + c.EDA, c.Gn + c.Antígenofecal + c.PS + c.EDA, c.Gp + c.Antígenofecal + c.PS + c.EDA, c.An + c.Antígenofecal + c.PS + c.EDA, c.Ap + c.Antígenofecal + c.PS + c.EDA, c.In + c.Antígenofecal + c.PS + c.EDA, c.Ip + c.Antígenofecal + c.PS + c.EDA, c.Dn + c.Antígenofecal + c.PS + c.EDA, c.Dp + c.Antígenofecal + c.PS + c.EDA, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev9    <- c(c.Nn + c.EDA, c.Np + c.EDA, c.Gn + c.EDA, c.Gp + c.EDA, c.An + c.EDA, c.Ap + c.EDA, c.In + c.EDA, c.Ip + c.EDA, c.Dn + c.EDA, c.Dp + c.EDA, c.1p + c.EDA, c.2p + c.EDA, c.3p + c.EDA, c.4p + c.EDA, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev10    <- c(c.Nn + c.EDA + c.PS, c.Np + c.EDA + c.PS, c.Gn + c.EDA + c.PS, c.Gp + c.EDA + c.PS, c.An + c.EDA + c.PS, c.Ap + c.EDA + c.PS, c.In + c.EDA + c.PS, c.Ip + c.EDA + c.PS, c.Dn + c.EDA + c.PS, c.Dp + c.EDA + c.PS, c.1p + c.EDA + c.PS, c.2p + c.EDA + c.PS, c.3p + c.EDA + c.PS, c.4p + c.EDA + c.PS, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prev11    <- c(c.Nn + c.EDA + c.Rx, c.Np  + c.EDA + c.Rx, c.Gn  + c.EDA + c.Rx, c.Gp  + c.EDA + c.Rx, c.An  + c.EDA + c.Rx, c.Ap  + c.EDA + c.Rx, c.In  + c.EDA + c.Rx, c.Ip  + c.EDA + c.Rx, c.Dn  + c.EDA + c.Rx, c.Dp  + c.EDA + c.Rx, c.1p  + c.EDA + c.Rx, c.2p  + c.EDA + c.Rx, c.3p  + c.EDA + c.Rx, c.4p  + c.EDA + c.Rx, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)


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
v.tu_prev9 <- m.M_prev9 %*% v.u_prev9
v.tu_prev10 <- m.M_prev10 %*% v.u_prev10
v.tu_prev11 <- m.M_prev10 %*% v.u_prev11

v.tc_no_prev <- m.M_no_prev %*% v.c_no_prev
v.tc_prev    <- m.M_prev    %*% v.c_prev
v.tc_prev2    <- m.M_prev2    %*% v.c_prev2
v.tc_prev3    <- m.M_prev3    %*% v.c_prev3
v.tc_prev4    <- m.M_prev4    %*% v.c_prev4
v.tc_prev5    <- m.M_prev5    %*% v.c_prev5
v.tc_prev6    <- m.M_prev6    %*% v.c_prev6
v.tc_prev7    <- m.M_prev7    %*% v.c_prev7
v.tc_prev8    <- m.M_prev8    %*% v.c_prev8
v.tc_prev9    <- m.M_prev9    %*% v.c_prev9
v.tc_prev10    <- m.M_prev10    %*% v.c_prev10
v.tc_prev11    <- m.M_prev11    %*% v.c_prev11

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
tu.d_prev9    <- t(v.tu_prev9)    %*% v.dwe
tu.d_prev10    <- t(v.tu_prev10)    %*% v.dwe
tu.d_prev11    <- t(v.tu_prev11)    %*% v.dwe

tc.d_no_prev <- t(v.tc_no_prev) %*% v.dwc
tc.d_prev    <- t(v.tc_prev)    %*% v.dwc
tc.d_prev2    <- t(v.tc_prev2)    %*% v.dwc
tc.d_prev3    <- t(v.tc_prev3)    %*% v.dwc
tc.d_prev4    <- t(v.tc_prev4)    %*% v.dwc
tc.d_prev5    <- t(v.tc_prev5)    %*% v.dwc
tc.d_prev6    <- t(v.tc_prev6)    %*% v.dwc
tc.d_prev7    <- t(v.tc_prev7)    %*% v.dwc
tc.d_prev8    <- t(v.tc_prev8)    %*% v.dwc
tc.d_prev9    <- t(v.tc_prev9)    %*% v.dwc
tc.d_prev10    <- t(v.tc_prev10)    %*% v.dwc
tc.d_prev11    <- t(v.tc_prev11)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_prev, tc.d_prev, tc.d_prev2, tc.d_prev3, tc.d_prev4, tc.d_prev5, tc.d_prev6, tc.d_prev7, tc.d_prev8, tc.d_prev9, tc.d_prev10, tc.d_prev11)
v.tu.d <- c(tu.d_no_prev, tu.d_prev, tu.d_prev2, tu.d_prev3, tu.d_prev4, tu.d_prev5, tu.d_prev6, tu.d_prev7, tu.d_prev8, tu.d_prev9, tu.d_prev10, tu.d_prev11)

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