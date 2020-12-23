####################################################################################
##########             Modelo CE Cáncer Gástrico       #####################
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
library(dampack)
library(ggplot2)

#### 02 Input Model Parameters ####
## Strategy names
v.names.str <- c("No prevention", "EDA", "PSEDA", "RxEDA", "EDAbp", "SerologíaEDA", "UreaAirePSEDA", "AntígenohecesPSEDA")  

## Number of strategies
n.str <- length(v.names.str)
## Markov model parameters
age     <- 20                                 # age at baseline
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles

v.n  <- c("Healthy", "1preclinical", "2preclinical", "3preclinical", "4preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2clinicala", "2clinicalb", "2clinicalc", "3clinicala", "3clinicalb", "3clinicalc", "4clinicala", "4clinicalb", "4clinicalc", "Dead")    # state names

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

n.s  <- length(v.n)                     # number of states
n.t  <- 60                              # number of cycles

p.HD <- 0.000616  # probability to die when healthy  ## Valor ficticio considerando 2018 muertes Chile 6,16 por 1000 hab

#Probability of no prevention

p.H1p <- 0.05  # probability to become stage 1 preclinical when healthy
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
p.1ccH <- 0.08 # probability to become healthy when stage 1 clinical c
p.2ccH <- 0.05 # probability to become healthy when stage 2 clinical c
p.3ccH <- 0.04 # probability to become healthy when stage 3 clinical c
p.4ccH <- 0.03 # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 1 (EDA)

p.p1_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
p.p1_1p1ca <- 0.05 * p.1p1ca  # probability to become stage 1 clinical a when stage 1 preclinical
p.p1_1ca1cb <- 0.05 * p.1ca1cb # probability to become stage 1 clinical b when stage 1 clinical a 
p.p1_1cb1cc <- 0.05 * p.1cb1cc # probability to become stage 1 clinical c when stage 1 clinical b 
p.p1_1p2p <- 0.05 * p.1p2p # probability to become stage 2 preclinical when stage 1 preclinical
p.p1_2p2ca <- 0.05  * p.2p2ca # probability to become stage 2 clinical a when stage 2 preclinical
p.p1_2ca2cb <- 0.05 * p.2ca2cb # probability to become stage 2 clinical b when stage 2 clinical a 
p.p1_2cb2cc <- 0.05 * p.2cb2cc # probability to become stage 2 clinical c when stage 2 clinical b 
p.p1_2p3p <- 0.05 * p.2p3p # probability to become stage 3 preclinical when stage 2 preclinical
p.p1_3p3ca <- 0.05 * p.3p3ca  # probability to become stage 3 clinical a when stage 3 preclinical
p.p1_3ca3cb <- 0.05 * p.3ca3cb # probability to become stage 3 clinical b when stage 3 clinical a 
p.p1_3cb3cc <- 0.05 * p.3cb3cc # probability to become stage 3 clinical c when stage 3 clinical b 
p.p1_3p4p <- 0.05 * p.3p4p # probability to become stage 4 preclinical when stage 3 preclinical
p.p1_4p4ca <- 0.05 * p.4p4ca # probability to become stage 4 clinical a when stage 4 preclinical
p.p1_4ca4cb <- 0.05 * p.4ca4cb # probability to become stage 4 clinical b when stage 4 clinical a 
p.p1_4cb4cc <- 0.05 * p.4cb4cc # probability to become stage 1 clinical c when stage 4 clinical b 
p.p1_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p1_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p1_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p1_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 2 (PSEDA)
p.p2_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p2_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p2_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p2_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p2_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 3 (RxEDA)
p.p3_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p3_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p3_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p3_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p3_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c


# Probability of transition when prevention strategy 4 (EDAbp)
p.p4_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p4_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p4_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p4_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p4_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 5 "SerologíaEDA"
p.p5_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p5_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p5_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p5_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p5_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 6 "UreaAirePSEDA"

p.p6_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p6_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p6_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p6_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p6_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

# Probability of transition when prevention strategy 7 "AntígenohecesPSEDA"

p.p7_H1p <- 0.05 * p.H1p  # probability to become stage 1 preclinical when healthy
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
p.p7_1ccH <- 0.08 * p.1ccH # probability to become healthy when stage 1 clinical c
p.p7_2ccH <- 0.05 * p.2ccH # probability to become healthy when stage 2 clinical c
p.p7_3ccH <- 0.04 * p.3ccH # probability to become healthy when stage 3 clinical c
p.p7_4ccH <- 0.03 * p.4ccH # probability to become healthy when stage 4 clinical c

#Hazard ratio 

hr.1p <- 1 #hazard ratio of death in stage 1 preclinical vs healthy
hr.2p <- 1 #hazard ratio of death in stage 2 preclinical vs healthy
hr.3p <- 1 #hazard ratio of death in stage 3 preclinical vs healthy
hr.4p <- 1 #hazard ratio of death in stage 4 preclinical vs healthy
hr.1ca <- 3 #hazard ratio of death in stage 1 clinical a vs healthy
hr.2ca <- 6  #hazard ratio of death in stage 2 clinical a vs healthy   
hr.3ca <- 9 #hazard ratio of death in stage 3 clinical a vs healthy
hr.4ca <- 10 #hazard ratio of death in stage 4 clinical a vs healthy
hr.1cb <- 3 #hazard ratio of death in stage 1 clinical b vs healthy
hr.2cb <- 6  #hazard ratio of death in stage 2 clinical b vs healthy   
hr.3cb <- 9 #hazard ratio of death in stage 3 clinical b vs healthy
hr.4cb <- 8 #hazard ratio of death in stage 4 clinical b vs healthy
hr.1cc <- 3 #hazard ratio of death in stage 1 clinical c vs healthy
hr.2cc <- 6  #hazard ratio of death in stage 2 clinical c vs healthy   
hr.3cc <- 9 #hazard ratio of death in stage 3 clinical c vs healthy
hr.4cc <- 8 #hazard ratio of death in stage 4 clinical c vs healthy


r.HD    <- - log(1 - p.HD) # rate of death in healthy
r.1pD   <- hr.1p * r.HD  	 # rate of death in 1 preclinical
r.2pD   <- hr.2p * r.HD  	 # rate of death in 2 preclinical
r.3pD   <- hr.3p * r.HD  	 # rate of death in 3 preclinical
r.4pD   <- hr.4p * r.HD  	 # rate of death in 4 preclinical
r.1caD   <- hr.1ca * r.HD  	 # rate of death in 1 clinical a
r.2caD   <- hr.2ca * r.HD  	 # rate of death in 2 clinical a
r.3caD   <- hr.3ca * r.HD  	 # rate of death in 3 clinical a
r.4caD   <- hr.4ca * r.HD  	 # rate of death in 4 clinical a  
r.1cbD   <- hr.1cb * r.HD  	 # rate of death in 1 clinical b
r.2cbD   <- hr.2cb * r.HD  	 # rate of death in 2 clinical b
r.3cbD   <- hr.3cb * r.HD  	 # rate of death in 3 clinical b
r.4cbD   <- hr.4cb * r.HD  	 # rate of death in 4 clinical b  
r.1ccD   <- hr.1cc * r.HD  	 # rate of death in 1 clinical c
r.2ccD   <- hr.2cc * r.HD  	 # rate of death in 2 clinical c
r.3ccD   <- hr.3cc * r.HD  	 # rate of death in 3 clinical c
r.4ccD   <- hr.4cc * r.HD  	 # rate of death in 4 clinical c
    
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

# Treatment costs   
c.H  <- 0                     # cost of remaining one cycle healthy
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

# Prevention costs 
c.EDA <- 2.34                   # cost of EDA # Estudio de verificación de costos 
c.PS <- 2                   # cost of PS # El costo de pepsinógeno sérico no está disponible
c.Rx <- 4                   # cost of Rx # El costo de la radiografía doble contrate no está disponible
c.UreaAire <- 4.39              #Cost of Prevention Strategy 2 "UreaAire" # Precio PUC convenio preferente  
c.Serología <- 2            # cost of serología # Precio falso
c.Antígenofecal <- 1.72                  #Cost of "Antígenofecal" # Precio PUC convenio preferente 

# Utilities  

u.H  <- 1                     # utility when healthy 
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
m.M_no_prevca <- m.M_prevca1 <- m.M_prevca2 <- m.M_prevca3 <- m.M_prevca4 <- m.M_prevca5 <- m.M_prevca6 <- m.M_prevca7 <- matrix(NA, 
                                nrow = n.t + 1, ncol = n.s,
                                dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_prevca) # show first 6 rows of the matrix 

# The cohort starts as healthy

m.M_no_prevca[1, ] <- m.M_prevca1[1, ] <- m.M_prevca2[1, ] <- m.M_prevca3[1, ] <- m.M_prevca4[1, ] <- m.M_prevca5[1, ] <- m.M_prevca6[1, ] <- m.M_prevca7[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)                     # initialize first cycle of Markov trace

#### 03.2 Transition probability MATRIX ####
# create the transition probability matrix without prevention
m.P_noprevca  <- matrix(0,
               nrow = n.s,
               ncol = n.s,
               dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_noprevca

# create the transition probability matrix with prevention strategy 1 cancer (EDA)
m.P_prevca1  <- matrix(0,
                      nrow = n.s,
                      ncol = n.s,
                      dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca1

# create the transition probability matrix with prevention strategy 2 cancer (PSEDA)
m.P_prevca2  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca2

# create the transition probability matrix with prevention strategy 3 cancer (RxEDA)
m.P_prevca3  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca3

# create the transition probability matrix with prevention strategy 4 cancer (EDAbp)
m.P_prevca4  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca4

# create the transition probability matrix with prevention strategy 5 cancer (SerologíaEDA)
m.P_prevca5  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca5

# create the transition probability matrix with prevention strategy 6 cancer (UreaAirePSEDA)
m.P_prevca6  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca6

# create the transition probability matrix with prevention strategy 7 cancer (AntígenohecesPSEDA)
m.P_prevca7  <- matrix(0,
                       nrow = n.s,
                       ncol = n.s,
                       dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
m.P_prevca7

# fill in the transition probability matrix without prevention
### From Healthy
m.P_noprevca["Healthy", "Healthy"] <- 1 - p.HD - p.H1p
m.P_noprevca["Healthy", "1preclinical"]    <- p.H1p
m.P_noprevca["Healthy", "Dead"]    <- p.HD

### From 1 preclinical 
m.P_noprevca["1preclinical", "1preclinical"]    <- 1 - p.1p1ca - p.1p2p - p.1pD
m.P_noprevca["1preclinical", "1clinicala"]    <- p.1p1ca
m.P_noprevca["1preclinical", "2preclinical"]    <- p.1p2p
m.P_noprevca["1preclinical", "Dead"]    <- p.1pD

### From 2 preclinical
m.P_noprevca["2preclinical", "2preclinical"]    <- 1 - p.2p2ca - p.2p3p - p.2pD
m.P_noprevca["2preclinical", "2clinicala"]    <- p.2p2ca
m.P_noprevca["2preclinical", "3preclinical"]    <- p.2p3p
m.P_noprevca["2preclinical", "Dead"]    <- p.2pD

### From 3 preclinical
m.P_noprevca["3preclinical", "3preclinical"]    <- 1 - p.3p3ca - p.3p4p - p.3pD
m.P_noprevca["3preclinical", "3clinicala"]    <- p.3p3ca
m.P_noprevca["3preclinical", "4preclinical"]    <- p.3p4p
m.P_noprevca["3preclinical", "Dead"]    <- p.3pD

### From 4 preclinical
m.P_noprevca["4preclinical", "4preclinical"]    <- 1 - p.4p4ca - p.4pD
m.P_noprevca["4preclinical", "4clinicala"]    <- p.4p4ca
m.P_noprevca["4preclinical", "Dead"]    <- p.4pD

### From 1 clinical a
m.P_noprevca["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.1ca1cb
m.P_noprevca["1clinicala", "1clinicalb"] <- p.1ca1cb
m.P_noprevca["1clinicala", "Dead"]    <- p.1caD

### From 2 clinical a
m.P_noprevca["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.2ca2cb
m.P_noprevca["2clinicala", "2clinicalb"] <- p.2ca2cb
m.P_noprevca["2clinicala", "Dead"]    <- p.2caD


### From 3 clinical a
m.P_noprevca["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.3ca3cb 
m.P_noprevca["3clinicala", "3clinicalb"] <- p.3ca3cb
m.P_noprevca["3clinicala", "Dead"]    <- p.3caD


### From 4 clinical a 
m.P_noprevca["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.4ca4cb
m.P_noprevca["4clinicala", "4clinicalb"] <- p.4ca4cb  
m.P_noprevca["4clinicala", "Dead"]    <- p.4caD

### From 1 clinical b
m.P_noprevca["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.1cb1cc
m.P_noprevca["1clinicalb", "1clinicalc"] <- p.1cb1cc
m.P_noprevca["1clinicalb", "Dead"]    <- p.1cbD

### From 2 clinical b
m.P_noprevca["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.2cb2cc
m.P_noprevca["2clinicalb", "2clinicalc"] <- p.2cb2cc
m.P_noprevca["2clinicalb", "Dead"]    <- p.2cbD

### From 3 clinical b
m.P_noprevca["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.3cb3cc 
m.P_noprevca["3clinicalb", "3clinicalc"] <- p.3cb3cc
m.P_noprevca["3clinicalb", "Dead"]    <- p.3cbD

### From 4 clinical b
m.P_noprevca["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.4cb4cc 
m.P_noprevca["4clinicalb", "4clinicalc"] <- p.4cb4cc
m.P_noprevca["4clinicalb", "Dead"]    <- p.4cbD

### From 1 clinical c
m.P_noprevca["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.1ccH
m.P_noprevca["1clinicalc", "Healthy"] <- p.1ccH
m.P_noprevca["1clinicalc", "Dead"]    <- p.1ccD

### From 2 clinical c
m.P_noprevca["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.2ccH
m.P_noprevca["2clinicalc", "Healthy"] <- p.2ccH
m.P_noprevca["2clinicalc", "Dead"]    <- p.2ccD

### From 3 clinical c
m.P_noprevca["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.3ccH 
m.P_noprevca["3clinicalc", "Healthy"] <- p.3ccH
m.P_noprevca["3clinicalc", "Dead"]    <- p.3ccD

### From 4 clinical c
m.P_noprevca["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.4ccH 
m.P_noprevca["4clinicalc", "Healthy"] <- p.4ccH
m.P_noprevca["4clinicalc", "Dead"]    <- p.4ccD

### From Dead
m.P_noprevca["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_noprevca)

## fill in the transition probability matrix with prevention strategy 1

### From Healthy	
m.P_prevca1["Healthy", "Healthy"] <- 1 - p.HD - p.p1_H1p	
m.P_prevca1["Healthy", "1preclinical"]    <- p.p1_H1p	
m.P_prevca1["Healthy", "Dead"]    <- p.HD	

### From 1 preclinical 	
m.P_prevca1["1preclinical", "1preclinical"]    <- 1 - p.p1_1p1ca - p.p1_1p2p - p.1pD	
m.P_prevca1["1preclinical", "1clinicala"]    <- p.p1_1p1ca	
m.P_prevca1["1preclinical", "2preclinical"]    <- p.p1_1p2p	
m.P_prevca1["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
m.P_prevca1["2preclinical", "2preclinical"]    <- 1 - p.p1_2p2ca - p.p1_2p3p - p.2pD	
m.P_prevca1["2preclinical", "2clinicala"]    <- p.p1_2p2ca	
m.P_prevca1["2preclinical", "3preclinical"]    <- p.p1_2p3p	
m.P_prevca1["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
m.P_prevca1["3preclinical", "3preclinical"]    <- 1 - p.p1_3p3ca - p.p1_3p4p - p.3pD	
m.P_prevca1["3preclinical", "3clinicala"]    <- p.p1_3p3ca	
m.P_prevca1["3preclinical", "4preclinical"]    <- p.p1_3p4p	
m.P_prevca1["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
m.P_prevca1["4preclinical", "4preclinical"]    <- 1 - p.p1_4p4ca - p.4pD	
m.P_prevca1["4preclinical", "4clinicala"]    <- p.p1_4p4ca	
m.P_prevca1["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
m.P_prevca1["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p1_1ca1cb	
m.P_prevca1["1clinicala", "1clinicalb"] <- p.p1_1ca1cb	
m.P_prevca1["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
m.P_prevca1["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p1_2ca2cb	
m.P_prevca1["2clinicala", "2clinicalb"] <- p.p1_2ca2cb	
m.P_prevca1["2clinicala", "Dead"]    <- p.2caD	


### From 3 clinical a	
m.P_prevca1["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p1_3ca3cb 	
m.P_prevca1["3clinicala", "3clinicalb"] <- p.p1_3ca3cb	
m.P_prevca1["3clinicala", "Dead"]    <- p.3caD	


### From 4 clinical a 	
m.P_prevca1["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p1_4ca4cb	
m.P_prevca1["4clinicala", "4clinicalb"] <- p.p1_4ca4cb  	
m.P_prevca1["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
m.P_prevca1["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p1_1cb1cc	
m.P_prevca1["1clinicalb", "1clinicalc"] <- p.p1_1cb1cc	
m.P_prevca1["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
m.P_prevca1["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p1_2cb2cc	
m.P_prevca1["2clinicalb", "2clinicalc"] <- p.p1_2cb2cc	
m.P_prevca1["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
m.P_prevca1["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p1_3cb3cc 	
m.P_prevca1["3clinicalb", "3clinicalc"] <- p.p1_3cb3cc	
m.P_prevca1["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
m.P_prevca1["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p1_4cb4cc 	
m.P_prevca1["4clinicalb", "4clinicalc"] <- p.p1_4cb4cc	
m.P_prevca1["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
m.P_prevca1["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p1_1ccH	
m.P_prevca1["1clinicalc", "Healthy"] <- p.p1_1ccH	
m.P_prevca1["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
m.P_prevca1["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p1_2ccH	
m.P_prevca1["2clinicalc", "Healthy"] <- p.p1_2ccH	
m.P_prevca1["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
m.P_prevca1["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p1_3ccH 	
m.P_prevca1["3clinicalc", "Healthy"] <- p.p1_3ccH	
m.P_prevca1["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
m.P_prevca1["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p1_4ccH 	
m.P_prevca1["4clinicalc", "Healthy"] <- p.p1_4ccH	
m.P_prevca1["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead	
m.P_prevca1["Dead", "Dead"] <- 1	

# check rows add up to 1	
rowSums(m.P_prevca1)	

## fill transitions prevention strategy 2

### From 1 preclinical 
m.P_prevca2["1preclinical", "1preclinical"]    <- 1 - p.p2_1p1ca - p.p2_1p2p - p.1pD 
m.P_prevca2["1preclinical", "1clinicala"]    <- p.p2_1p1ca 
m.P_prevca2["1preclinical", "2preclinical"]    <- p.p2_1p2p 
m.P_prevca2["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca2["2preclinical", "2preclinical"]    <- 1 - p.p2_2p2ca - p.p2_2p3p - p.2pD 
m.P_prevca2["2preclinical", "2clinicala"]    <- p.p2_2p2ca 
m.P_prevca2["2preclinical", "3preclinical"]    <- p.p2_2p3p 
m.P_prevca2["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca2["3preclinical", "3preclinical"]    <- 1 - p.p2_3p3ca - p.p2_3p4p - p.3pD 
m.P_prevca2["3preclinical", "3clinicala"]    <- p.p2_3p3ca 
m.P_prevca2["3preclinical", "4preclinical"]    <- p.p2_3p4p 
m.P_prevca2["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca2["4preclinical", "4preclinical"]    <- 1 - p.p2_4p4ca - p.4pD 
m.P_prevca2["4preclinical", "4clinicala"]    <- p.p2_4p4ca 
m.P_prevca2["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca2["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p2_1ca1cb 
m.P_prevca2["1clinicala", "1clinicalb"] <- p.p2_1ca1cb 
m.P_prevca2["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca2["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p2_2ca2cb 
m.P_prevca2["2clinicala", "2clinicalb"] <- p.p2_2ca2cb 
m.P_prevca2["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca2["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p2_3ca3cb 
m.P_prevca2["3clinicala", "3clinicalb"] <- p.p2_3ca3cb 
m.P_prevca2["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca2["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p2_4ca4cb 
m.P_prevca2["4clinicala", "4clinicalb"] <- p.p2_4ca4cb  
m.P_prevca2["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca2["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p2_1cb1cc 
m.P_prevca2["1clinicalb", "1clinicalc"] <- p.p2_1cb1cc 
m.P_prevca2["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca2["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p2_2cb2cc 
m.P_prevca2["2clinicalb", "2clinicalc"] <- p.p2_2cb2cc 
m.P_prevca2["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca2["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p2_3cb3cc 
m.P_prevca2["3clinicalb", "3clinicalc"] <- p.p2_3cb3cc 
m.P_prevca2["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca2["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p2_4cb4cc 
m.P_prevca2["4clinicalb", "4clinicalc"] <- p.p2_4cb4cc 
m.P_prevca2["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca2["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p2_1ccH 
m.P_prevca2["1clinicalc", "Healthy"] <- p.p2_1ccH 
m.P_prevca2["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca2["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p2_2ccH 
m.P_prevca2["2clinicalc", "Healthy"] <- p.p2_2ccH 
m.P_prevca2["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca2["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p2_3ccH 
m.P_prevca2["3clinicalc", "Healthy"] <- p.p2_3ccH 
m.P_prevca2["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca2["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p2_4ccH 
m.P_prevca2["4clinicalc", "Healthy"] <- p.p2_4ccH 
m.P_prevca2["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca2["Dead", "Dead"] <- 1 

# check rows add up to 1 
rowSums(m.P_prevca2) 

# fill transitions prevention strategy 3

### From 1 preclinical 
m.P_prevca3["1preclinical", "1preclinical"]    <- 1 - p.p3_1p1ca - p.p3_1p2p - p.1pD 
m.P_prevca3["1preclinical", "1clinicala"]    <- p.p3_1p1ca 
m.P_prevca3["1preclinical", "2preclinical"]    <- p.p3_1p2p 
m.P_prevca3["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca3["2preclinical", "2preclinical"]    <- 1 - p.p3_2p2ca - p.p3_2p3p - p.2pD 
m.P_prevca3["2preclinical", "2clinicala"]    <- p.p3_2p2ca 
m.P_prevca3["2preclinical", "3preclinical"]    <- p.p3_2p3p 
m.P_prevca3["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca3["3preclinical", "3preclinical"]    <- 1 - p.p3_3p3ca - p.p3_3p4p - p.3pD 
m.P_prevca3["3preclinical", "3clinicala"]    <- p.p3_3p3ca 
m.P_prevca3["3preclinical", "4preclinical"]    <- p.p3_3p4p 
m.P_prevca3["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca3["4preclinical", "4preclinical"]    <- 1 - p.p3_4p4ca - p.4pD 
m.P_prevca3["4preclinical", "4clinicala"]    <- p.p3_4p4ca 
m.P_prevca3["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca3["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p3_1ca1cb 
m.P_prevca3["1clinicala", "1clinicalb"] <- p.p3_1ca1cb 
m.P_prevca3["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca3["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p3_2ca2cb 
m.P_prevca3["2clinicala", "2clinicalb"] <- p.p3_2ca2cb 
m.P_prevca3["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca3["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p3_3ca3cb 
m.P_prevca3["3clinicala", "3clinicalb"] <- p.p3_3ca3cb 
m.P_prevca3["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca3["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p3_4ca4cb 
m.P_prevca3["4clinicala", "4clinicalb"] <- p.p3_4ca4cb  
m.P_prevca3["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca3["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p3_1cb1cc 
m.P_prevca3["1clinicalb", "1clinicalc"] <- p.p3_1cb1cc 
m.P_prevca3["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca3["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p3_2cb2cc 
m.P_prevca3["2clinicalb", "2clinicalc"] <- p.p3_2cb2cc 
m.P_prevca3["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca3["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p3_3cb3cc 
m.P_prevca3["3clinicalb", "3clinicalc"] <- p.p3_3cb3cc 
m.P_prevca3["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca3["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p3_4cb4cc 
m.P_prevca3["4clinicalb", "4clinicalc"] <- p.p3_4cb4cc 
m.P_prevca3["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca3["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p3_1ccH 
m.P_prevca3["1clinicalc", "Healthy"] <- p.p3_1ccH 
m.P_prevca3["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca3["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p3_2ccH 
m.P_prevca3["2clinicalc", "Healthy"] <- p.p3_2ccH 
m.P_prevca3["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca3["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p3_3ccH 
m.P_prevca3["3clinicalc", "Healthy"] <- p.p3_3ccH 
m.P_prevca3["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca3["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p3_4ccH 
m.P_prevca3["4clinicalc", "Healthy"] <- p.p3_4ccH 
m.P_prevca3["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca3["Dead", "Dead"] <- 1 

# fill transitions prevention strategy 4 

### From 1 preclinical 
m.P_prevca4["1preclinical", "1preclinical"]    <- 1 - p.p4_1p1ca - p.p4_1p2p - p.1pD 
m.P_prevca4["1preclinical", "1clinicala"]    <- p.p4_1p1ca 
m.P_prevca4["1preclinical", "2preclinical"]    <- p.p4_1p2p 
m.P_prevca4["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca4["2preclinical", "2preclinical"]    <- 1 - p.p4_2p2ca - p.p4_2p3p - p.2pD 
m.P_prevca4["2preclinical", "2clinicala"]    <- p.p4_2p2ca 
m.P_prevca4["2preclinical", "3preclinical"]    <- p.p4_2p3p 
m.P_prevca4["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca4["3preclinical", "3preclinical"]    <- 1 - p.p4_3p3ca - p.p4_3p4p - p.3pD 
m.P_prevca4["3preclinical", "3clinicala"]    <- p.p4_3p3ca 
m.P_prevca4["3preclinical", "4preclinical"]    <- p.p4_3p4p 
m.P_prevca4["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca4["4preclinical", "4preclinical"]    <- 1 - p.p4_4p4ca - p.4pD 
m.P_prevca4["4preclinical", "4clinicala"]    <- p.p4_4p4ca 
m.P_prevca4["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca4["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p4_1ca1cb 
m.P_prevca4["1clinicala", "1clinicalb"] <- p.p4_1ca1cb 
m.P_prevca4["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca4["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p4_2ca2cb 
m.P_prevca4["2clinicala", "2clinicalb"] <- p.p4_2ca2cb 
m.P_prevca4["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca4["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p4_3ca3cb 
m.P_prevca4["3clinicala", "3clinicalb"] <- p.p4_3ca3cb 
m.P_prevca4["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca4["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p4_4ca4cb 
m.P_prevca4["4clinicala", "4clinicalb"] <- p.p4_4ca4cb  
m.P_prevca4["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca4["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p4_1cb1cc 
m.P_prevca4["1clinicalb", "1clinicalc"] <- p.p4_1cb1cc 
m.P_prevca4["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca4["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p4_2cb2cc 
m.P_prevca4["2clinicalb", "2clinicalc"] <- p.p4_2cb2cc 
m.P_prevca4["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca4["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p4_3cb3cc 
m.P_prevca4["3clinicalb", "3clinicalc"] <- p.p4_3cb3cc 
m.P_prevca4["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca4["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p4_4cb4cc 
m.P_prevca4["4clinicalb", "4clinicalc"] <- p.p4_4cb4cc 
m.P_prevca4["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca4["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p4_1ccH 
m.P_prevca4["1clinicalc", "Healthy"] <- p.p4_1ccH 
m.P_prevca4["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca4["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p4_2ccH 
m.P_prevca4["2clinicalc", "Healthy"] <- p.p4_2ccH 
m.P_prevca4["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca4["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p4_3ccH 
m.P_prevca4["3clinicalc", "Healthy"] <- p.p4_3ccH 
m.P_prevca4["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca4["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p4_4ccH 
m.P_prevca4["4clinicalc", "Healthy"] <- p.p4_4ccH 
m.P_prevca4["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca4["Dead", "Dead"] <- 1 

## fill transitions prevention strategy 5
### From 1 preclinical 
m.P_prevca5["1preclinical", "1preclinical"]    <- 1 - p.p5_1p1ca - p.p5_1p2p - p.1pD 
m.P_prevca5["1preclinical", "1clinicala"]    <- p.p5_1p1ca 
m.P_prevca5["1preclinical", "2preclinical"]    <- p.p5_1p2p 
m.P_prevca5["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca5["2preclinical", "2preclinical"]    <- 1 - p.p5_2p2ca - p.p5_2p3p - p.2pD 
m.P_prevca5["2preclinical", "2clinicala"]    <- p.p5_2p2ca 
m.P_prevca5["2preclinical", "3preclinical"]    <- p.p5_2p3p 
m.P_prevca5["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca5["3preclinical", "3preclinical"]    <- 1 - p.p5_3p3ca - p.p5_3p4p - p.3pD 
m.P_prevca5["3preclinical", "3clinicala"]    <- p.p5_3p3ca 
m.P_prevca5["3preclinical", "4preclinical"]    <- p.p5_3p4p 
m.P_prevca5["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca5["4preclinical", "4preclinical"]    <- 1 - p.p5_4p4ca - p.4pD 
m.P_prevca5["4preclinical", "4clinicala"]    <- p.p5_4p4ca 
m.P_prevca5["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca5["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p5_1ca1cb 
m.P_prevca5["1clinicala", "1clinicalb"] <- p.p5_1ca1cb 
m.P_prevca5["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca5["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p5_2ca2cb 
m.P_prevca5["2clinicala", "2clinicalb"] <- p.p5_2ca2cb 
m.P_prevca5["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca5["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p5_3ca3cb 
m.P_prevca5["3clinicala", "3clinicalb"] <- p.p5_3ca3cb 
m.P_prevca5["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca5["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p5_4ca4cb 
m.P_prevca5["4clinicala", "4clinicalb"] <- p.p5_4ca4cb  
m.P_prevca5["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca5["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p5_1cb1cc 
m.P_prevca5["1clinicalb", "1clinicalc"] <- p.p5_1cb1cc 
m.P_prevca5["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca5["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p5_2cb2cc 
m.P_prevca5["2clinicalb", "2clinicalc"] <- p.p5_2cb2cc 
m.P_prevca5["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca5["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p5_3cb3cc 
m.P_prevca5["3clinicalb", "3clinicalc"] <- p.p5_3cb3cc 
m.P_prevca5["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca5["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p5_4cb4cc 
m.P_prevca5["4clinicalb", "4clinicalc"] <- p.p5_4cb4cc 
m.P_prevca5["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca5["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p5_1ccH 
m.P_prevca5["1clinicalc", "Healthy"] <- p.p5_1ccH 
m.P_prevca5["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca5["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p5_2ccH 
m.P_prevca5["2clinicalc", "Healthy"] <- p.p5_2ccH 
m.P_prevca5["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca5["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p5_3ccH 
m.P_prevca5["3clinicalc", "Healthy"] <- p.p5_3ccH 
m.P_prevca5["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca5["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p5_4ccH 
m.P_prevca5["4clinicalc", "Healthy"] <- p.p5_4ccH 
m.P_prevca5["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca5["Dead", "Dead"] <- 1 

## fill transition probabilities prevention strategy 6 

### From 1 preclinical 
m.P_prevca6["1preclinical", "1preclinical"]    <- 1 - p.p6_1p1ca - p.p6_1p2p - p.1pD 
m.P_prevca6["1preclinical", "1clinicala"]    <- p.p6_1p1ca 
m.P_prevca6["1preclinical", "2preclinical"]    <- p.p6_1p2p 
m.P_prevca6["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca6["2preclinical", "2preclinical"]    <- 1 - p.p6_2p2ca - p.p6_2p3p - p.2pD 
m.P_prevca6["2preclinical", "2clinicala"]    <- p.p6_2p2ca 
m.P_prevca6["2preclinical", "3preclinical"]    <- p.p6_2p3p 
m.P_prevca6["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca6["3preclinical", "3preclinical"]    <- 1 - p.p6_3p3ca - p.p6_3p4p - p.3pD 
m.P_prevca6["3preclinical", "3clinicala"]    <- p.p6_3p3ca 
m.P_prevca6["3preclinical", "4preclinical"]    <- p.p6_3p4p 
m.P_prevca6["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca6["4preclinical", "4preclinical"]    <- 1 - p.p6_4p4ca - p.4pD 
m.P_prevca6["4preclinical", "4clinicala"]    <- p.p6_4p4ca 
m.P_prevca6["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca6["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p6_1ca1cb 
m.P_prevca6["1clinicala", "1clinicalb"] <- p.p6_1ca1cb 
m.P_prevca6["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca6["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p6_2ca2cb 
m.P_prevca6["2clinicala", "2clinicalb"] <- p.p6_2ca2cb 
m.P_prevca6["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca6["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p6_3ca3cb 
m.P_prevca6["3clinicala", "3clinicalb"] <- p.p6_3ca3cb 
m.P_prevca6["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca6["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p6_4ca4cb 
m.P_prevca6["4clinicala", "4clinicalb"] <- p.p6_4ca4cb  
m.P_prevca6["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca6["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p6_1cb1cc 
m.P_prevca6["1clinicalb", "1clinicalc"] <- p.p6_1cb1cc 
m.P_prevca6["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca6["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p6_2cb2cc 
m.P_prevca6["2clinicalb", "2clinicalc"] <- p.p6_2cb2cc 
m.P_prevca6["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca6["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p6_3cb3cc 
m.P_prevca6["3clinicalb", "3clinicalc"] <- p.p6_3cb3cc 
m.P_prevca6["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca6["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p6_4cb4cc 
m.P_prevca6["4clinicalb", "4clinicalc"] <- p.p6_4cb4cc 
m.P_prevca6["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca6["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p6_1ccH 
m.P_prevca6["1clinicalc", "Healthy"] <- p.p6_1ccH 
m.P_prevca6["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca6["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p6_2ccH 
m.P_prevca6["2clinicalc", "Healthy"] <- p.p6_2ccH 
m.P_prevca6["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca6["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p6_3ccH 
m.P_prevca6["3clinicalc", "Healthy"] <- p.p6_3ccH 
m.P_prevca6["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca6["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p6_4ccH 
m.P_prevca6["4clinicalc", "Healthy"] <- p.p6_4ccH 
m.P_prevca6["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca6["Dead", "Dead"] <- 1 

## fill the transition probability 7

### From 1 preclinical 
m.P_prevca7["1preclinical", "1preclinical"]    <- 1 - p.p7_1p1ca - p.p7_1p2p - p.1pD 
m.P_prevca7["1preclinical", "1clinicala"]    <- p.p7_1p1ca 
m.P_prevca7["1preclinical", "2preclinical"]    <- p.p7_1p2p 
m.P_prevca7["1preclinical", "Dead"]    <- p.1pD 

### From 2 preclinical 
m.P_prevca7["2preclinical", "2preclinical"]    <- 1 - p.p7_2p2ca - p.p7_2p3p - p.2pD 
m.P_prevca7["2preclinical", "2clinicala"]    <- p.p7_2p2ca 
m.P_prevca7["2preclinical", "3preclinical"]    <- p.p7_2p3p 
m.P_prevca7["2preclinical", "Dead"]    <- p.2pD 

### From 3 preclinical 
m.P_prevca7["3preclinical", "3preclinical"]    <- 1 - p.p7_3p3ca - p.p7_3p4p - p.3pD 
m.P_prevca7["3preclinical", "3clinicala"]    <- p.p7_3p3ca 
m.P_prevca7["3preclinical", "4preclinical"]    <- p.p7_3p4p 
m.P_prevca7["3preclinical", "Dead"]    <- p.3pD 

### From 4 preclinical 
m.P_prevca7["4preclinical", "4preclinical"]    <- 1 - p.p7_4p4ca - p.4pD 
m.P_prevca7["4preclinical", "4clinicala"]    <- p.p7_4p4ca 
m.P_prevca7["4preclinical", "Dead"]    <- p.4pD 

### From 1 clinical a 
m.P_prevca7["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p7_1ca1cb 
m.P_prevca7["1clinicala", "1clinicalb"] <- p.p7_1ca1cb 
m.P_prevca7["1clinicala", "Dead"]    <- p.1caD 

### From 2 clinical a 
m.P_prevca7["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p7_2ca2cb 
m.P_prevca7["2clinicala", "2clinicalb"] <- p.p7_2ca2cb 
m.P_prevca7["2clinicala", "Dead"]    <- p.2caD 


### From 3 clinical a 
m.P_prevca7["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p7_3ca3cb 
m.P_prevca7["3clinicala", "3clinicalb"] <- p.p7_3ca3cb 
m.P_prevca7["3clinicala", "Dead"]    <- p.3caD 


### From 4 clinical a 
m.P_prevca7["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p7_4ca4cb 
m.P_prevca7["4clinicala", "4clinicalb"] <- p.p7_4ca4cb  
m.P_prevca7["4clinicala", "Dead"]    <- p.4caD 

### From 1 clinical b 
m.P_prevca7["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p7_1cb1cc 
m.P_prevca7["1clinicalb", "1clinicalc"] <- p.p7_1cb1cc 
m.P_prevca7["1clinicalb", "Dead"]    <- p.1cbD 

### From 2 clinical b 
m.P_prevca7["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p7_2cb2cc 
m.P_prevca7["2clinicalb", "2clinicalc"] <- p.p7_2cb2cc 
m.P_prevca7["2clinicalb", "Dead"]    <- p.2cbD 

### From 3 clinical b 
m.P_prevca7["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p7_3cb3cc 
m.P_prevca7["3clinicalb", "3clinicalc"] <- p.p7_3cb3cc 
m.P_prevca7["3clinicalb", "Dead"]    <- p.3cbD 

### From 4 clinical b 
m.P_prevca7["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p7_4cb4cc 
m.P_prevca7["4clinicalb", "4clinicalc"] <- p.p7_4cb4cc 
m.P_prevca7["4clinicalb", "Dead"]    <- p.4cbD 

### From 1 clinical c 
m.P_prevca7["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p7_1ccH 
m.P_prevca7["1clinicalc", "Healthy"] <- p.p7_1ccH 
m.P_prevca7["1clinicalc", "Dead"]    <- p.1ccD 

### From 2 clinical c 
m.P_prevca7["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p7_2ccH 
m.P_prevca7["2clinicalc", "Healthy"] <- p.p7_2ccH 
m.P_prevca7["2clinicalc", "Dead"]    <- p.2ccD 

### From 3 clinical c 
m.P_prevca7["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p7_3ccH 
m.P_prevca7["3clinicalc", "Healthy"] <- p.p7_3ccH 
m.P_prevca7["3clinicalc", "Dead"]    <- p.3ccD 

### From 4 clinical c 
m.P_prevca7["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p7_4ccH 
m.P_prevca7["4clinicalc", "Healthy"] <- p.p7_4ccH 
m.P_prevca7["4clinicalc", "Dead"]    <- p.4ccD 

### From Dead 
m.P_prevca7["Dead", "Dead"] <- 1 


#### 04 Run Markov model ####
for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_prevca[t + 1, ] <- t(m.M_no_prevca[t, ]) %*% m.P_noprevca # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca1[t + 1, ]    <- t(m.M_prevca1[t, ])    %*% m.P_prevca1   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca2[t + 1, ]    <- t(m.M_prevca2[t, ])    %*% m.P_prevca2   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca3[t + 1, ]    <- t(m.M_prevca3[t, ])    %*% m.P_prevca3   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca4[t + 1, ]    <- t(m.M_prevca4[t, ])    %*% m.P_prevca4   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca5[t + 1, ]    <- t(m.M_prevca5[t, ])    %*% m.P_prevca5   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca6[t + 1, ]    <- t(m.M_prevca6[t, ])    %*% m.P_prevca6   # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prevca7[t + 1, ]    <- t(m.M_prevca7[t, ])    %*% m.P_prevca7   # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_prevca)  # show the first 6 lines of the matrix with no prevention
head(m.M_prevca1)  # show the first 6 lines of the matrix with prevention strategy 1
head(m.M_prevca2)  # show the first 6 lines of the matrix with prevention strategy 2
head(m.M_prevca3)  # show the first 6 lines of the matrix with prevention strategy 3
head(m.M_prevca4)  # show the first 6 lines of the matrix with prevention strategy 4
head(m.M_prevca5)  # show the first 6 lines of the matrix with prevention strategy 5
head(m.M_prevca6)  # show the first 6 lines of the matrix with prevention strategy 6
head(m.M_prevca7)  # show the first 6 lines of the matrix with prevention strategy 7


#### 05 Compute and Plot Epidemiological Outcomes ####
#### 05.1 Cohort trace #####
matplot(m.M_no_prevca, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

#### 05.2 Overall Survival (OS) #####

## No prevention 

v.os_no_prevca <- 1 - m.M_no_prevca[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_no_prevca, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

## With prevention 
v.os_prevca1 <- 1 - m.M_prevca1[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca1, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 1")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca2 <- 1 - m.M_prevca2[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca2, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 2")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca3 <- 1 - m.M_prevca3[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca3, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 3")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca4 <- 1 - m.M_prevca4[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca4, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 4")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca5 <- 1 - m.M_prevca5[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca5, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 5")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca6 <- 1 - m.M_prevca6[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca6, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 6")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

v.os_prevca7 <- 1 - m.M_prevca7[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_prevca7, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival Prevention 7")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid

#### 05.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os_no_prevca)                       # summing probablity of OS over time  (i.e. life expectancy) No prevention
v.le <- sum(v.os_prevca1)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 1
v.le <- sum(v.os_prevca2)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 2
v.le <- sum(v.os_prevca3)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 3
v.le <- sum(v.os_prevca4)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 4
v.le <- sum(v.os_prevca5)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 5
v.le <- sum(v.os_prevca6)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 6
v.le <- sum(v.os_prevca7)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 7

#### 05.3 Disease prevalence #####

## No prevention
v.preva_no_prevca <- rowSums(m.M_no_prevca[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_no_prevca
plot(v.preva_no_prevca,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

## Prevention Strategy 1 
v.preva_prevca1 <- rowSums(m.M_prevca1[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca1
plot(v.preva_prevca1,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca2 <- rowSums(m.M_prevca2[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca2
plot(v.preva_prevca2,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca3 <- rowSums(m.M_prevca3[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca3
plot(v.preva_prevca3,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca4<- rowSums(m.M_prevca4[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca4
plot(v.preva_prevca4,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca5 <- rowSums(m.M_prevca5[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca5
plot(v.preva_prevca5,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca6 <- rowSums(m.M_prevca6[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca6
plot(v.preva_prevca6,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

v.preva_prevca7 <- rowSums(m.M_prevca7[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "3preclinical",  "3clinicala", "3clinicalb", "3clinicalc","4preclinical", "4clinicala", "4clinicalb" , "4clinicalc" )])/v.os_prevca7
plot(v.preva_prevca7,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 05.4 Proportion of sick in stage 1 state #####
v.prop.1c_noprev <- rowSums(m.M_no_prevca[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_no_prevca
plot(0:n.t, v.prop.1c_noprev,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")


#### 06 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_prevca <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca1    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca2    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca3    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca4    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca5    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca6    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)
v.u_prevca7    <- c(u.H, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.2ca, u.2cb, u.2cc, u.3ca, u.3cb, u.3cc, u.4ca, u.4cb, u.4cc, u.D)

v.c_no_prevca <- c(c.H, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca1    <- c(c.H + c.EDA, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca2    <- c(c.H + c.EDA + c.PS, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca3    <- c(c.H + c.EDA + c.Rx, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca4    <- c(c.H + c.EDA , c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca5    <- c(c.H + c.EDA + c.Serología, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca6    <- c(c.H + c.EDA + c.UreaAire, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)
v.c_prevca7    <- c(c.H + c.EDA + c.PS + c.Antígenofecal, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.2ca, c.2cb, c.2cc, c.3ca, c.3cb, c.3cc, c.4ca, c.4cb, c.4cc, c.D)

#### 06.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_prevca <- m.M_no_prevca %*% v.u_no_prevca
v.tu_prevca1 <- m.M_prevca1 %*% v.u_prevca1
v.tu_prevca2 <- m.M_prevca2 %*% v.u_prevca2
v.tu_prevca3 <- m.M_prevca3 %*% v.u_prevca3
v.tu_prevca4 <- m.M_prevca4 %*% v.u_prevca4
v.tu_prevca5 <- m.M_prevca5 %*% v.u_prevca5
v.tu_prevca6 <- m.M_prevca6 %*% v.u_prevca6
v.tu_prevca7 <- m.M_prevca7 %*% v.u_prevca7


v.tc_no_prevca <- m.M_no_prevca %*% v.c_no_prevca
v.tc_prevca1    <- m.M_prevca1    %*% v.c_prevca1
v.tc_prevca2    <- m.M_prevca2    %*% v.c_prevca2
v.tc_prevca3    <- m.M_prevca3    %*% v.c_prevca3
v.tc_prevca4    <- m.M_prevca4    %*% v.c_prevca4
v.tc_prevca5    <- m.M_prevca5    %*% v.c_prevca5
v.tc_prevca6    <- m.M_prevca6    %*% v.c_prevca6
v.tc_prevca7    <- m.M_prevca7    %*% v.c_prevca7

#### 06.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_prevca <- t(v.tu_no_prevca) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
tu.d_prevca1    <- t(v.tu_prevca1)    %*% v.dwe
tu.d_prevca2    <- t(v.tu_prevca2)    %*% v.dwe
tu.d_prevca3    <- t(v.tu_prevca3)    %*% v.dwe
tu.d_prevca4    <- t(v.tu_prevca4)    %*% v.dwe
tu.d_prevca5    <- t(v.tu_prevca5)    %*% v.dwe
tu.d_prevca6    <- t(v.tu_prevca6)    %*% v.dwe
tu.d_prevca7    <- t(v.tu_prevca7)    %*% v.dwe

tc.d_no_prevca <- t(v.tc_no_prevca) %*% v.dwc
tc.d_prevca1    <- t(v.tc_prevca1)    %*% v.dwc
tc.d_prevca2    <- t(v.tc_prevca2)    %*% v.dwc
tc.d_prevca3    <- t(v.tc_prevca3)    %*% v.dwc
tc.d_prevca4    <- t(v.tc_prevca4)    %*% v.dwc
tc.d_prevca5    <- t(v.tc_prevca5)    %*% v.dwc
tc.d_prevca6    <- t(v.tc_prevca6)    %*% v.dwc
tc.d_prevca7    <- t(v.tc_prevca7)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_prevca, tc.d_prevca1, tc.d_prevca2, tc.d_prevca3, tc.d_prevca4, tc.d_prevca5, tc.d_prevca6, tc.d_prevca7)
v.tu.d <- c(tu.d_no_prevca, tu.d_prevca1, tu.d_prevca2, tu.d_prevca3, tu.d_prevca4, tu.d_prevca5, tu.d_prevca6, tu.d_prevca7)

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

#### 08 Plot frontier of FONIS CaGa model ####
plot(m.cea, xlim = c(15.7, 16.5))
ggsave("figs/Markov-FONISCaGast-CEA-Frontier.png", width = 8, height = 6)

