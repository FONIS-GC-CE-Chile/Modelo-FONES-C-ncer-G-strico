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

setwd("/Users/franc/Desktop/ESP/FONIScancer ")

#### 01 Load packages ####
library(dampack)
library(reshape2)
library(truncnorm)
library(ggplot2)
library(dplyr)
library(readxl)
library(wesanderson)


#### 02 Input Model Parameters ####


  ##    02.1 General Parameters ------------------------------------------------------

### Strategies 
# Prevention Strategy 1 - UreaAire: Prueba de Urea en Aire espirado marcado con Carbono 13 a los 30 años 1 vez 
# Prevention Strategy 2 - Antígenofecal Estudio de serología + antígeno en heces fecales a los 30 años 1 vez  
#Prevention Strategy 3 - EDA: EDA para estudio histológico a los 45 años 
#Prevention Strategy 4 - PSEDA: Pepsinógeno sérico + EDA para estudio histológico en casos positivos a los 40 años
#Prevention Strategy 5 - Serología: Serología + antígeno en deposiciones

## Strategy names
v.names.str <- c("No Prevention", "UreaAire", "Antígenofecal", "EDA", "PSEDA", "Serología" )

## Number of strategies
n.str <- length(v.names.str)

## Markov model parameters
age     <- 30                                 # age at baseline
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles
cycle_length <- 1 

v.n  <- c("Normalhpn", "Normalhpp", 
          "Gastritishpn", "Gastritishpp", 
          "Atrophyhpn", "Atrophyhpp", 
          "Intestinalhpn", "Intestinalhpp", 
          "Dysplasiahpn", "Dysplasiahpp", 
          "1preclinical", "2preclinical", 
          "3preclinical", "4preclinical", 
          "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", 
          "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf",
          "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", 
          "4clinicala", "4clinicalb", "4clinicalc", "4clinicald", "4clinicale", "4clinicalf",
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
## 1 clinical b = Stage 1, 2 years after diagnosis 
## 1 clinical c = Stage 1, 3 years after diagnosis 
## 1 clinical d = Stage 1, 4 years after diagnosis 
## 1 clinical e = Stage 1, 5 years after diagnosis 
## 1 clinical f = Stage 1, more than 6 years after diagnosis 
## 2 clinical a = Stage 2, 1st year after diagnosis
## 2 clinical b = Stage 2, 2 years after diagnosis 
## 2 clinical c = Stage 2, 3 years after diagnosis 
## 2 clinical d = Stage 2, 4 years after diagnosis 
## 2 clinical e = Stage 2, 5 years after diagnosis 
## 2 clinical f = Stage 2, more than 6 years after diagnosis 
## 3 clinical a = Stage 3, 1st year after diagnosis
## 3 clinical b = Stage 3, 2 years after diagnosis 
## 3 clinical c = Stage 3, 3 years after diagnosis 
## 3 clinical d = Stage 3, 4 years after diagnosis 
## 3 clinical e = Stage 3, 5 years after diagnosis 
## 3 clinical f = Stage 3, more than 6 years after diagnosis 
## 4 clinical a = Stage 4, 1st year after diagnosis
## 4 clinical b = Stage 4, 2 years after diagnosis 
## 4 clinical c = Stage 4, 3 years after diagnosis 
## 4 clinical d = Stage 4, 4 years after diagnosis 
## 4 clinical e = Stage 4, 5 years after diagnosis 
## 4 clinical f = Stage 4, more than 6 years after diagnosis 

## Death = Death 

  ##    02.2 Probabilidades de muerte ------------------------------------------------


## Load life table mortality rates from `.csv` format
LifeTable_Chile <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile.csv")

## Load hazards from the relative survival analysis 
SR1a <- read_excel("SRfinal.xlsx", sheet = "1a", range = "A1:B72")
SR1b <- read_excel("SRfinal.xlsx", sheet = "1b", range = "A1:B72")
SR1c <- read_excel("SRfinal.xlsx", sheet = "1c", range = "A1:B72")
SR1d <- read_excel("SRfinal.xlsx", sheet = "1d", range = "A1:B72")
SR1e <- read_excel("SRfinal.xlsx", sheet = "1e", range = "A1:B72")
SR1f <- read_excel("SRfinal.xlsx", sheet = "1f", range = "A1:B72")

SR2a <- read_excel("SRfinal.xlsx", sheet = "2a", range = "A1:B72")
SR2b <- read_excel("SRfinal.xlsx", sheet = "2b", range = "A1:B72")
SR2c <- read_excel("SRfinal.xlsx", sheet = "2c", range = "A1:B72")
SR2d <- read_excel("SRfinal.xlsx", sheet = "2d", range = "A1:B72")
SR2e <- read_excel("SRfinal.xlsx", sheet = "2e", range = "A1:B72")
SR2f <- read_excel("SRfinal.xlsx", sheet = "2f", range = "A1:B72")

SR3a <- read_excel("SRfinal.xlsx", sheet = "3a", range = "A1:B72")
SR3b <- read_excel("SRfinal.xlsx", sheet = "3b", range = "A1:B72")
SR3c <- read_excel("SRfinal.xlsx", sheet = "3c", range = "A1:B72")
SR3d <- read_excel("SRfinal.xlsx", sheet = "3d", range = "A1:B72")
SR3e <- read_excel("SRfinal.xlsx", sheet = "3e", range = "A1:B72")
SR3f <- read_excel("SRfinal.xlsx", sheet = "3f", range = "A1:B72")

SR4a <- read_excel("SRfinal.xlsx", sheet = "4a", range = "A1:B72")
SR4b <- read_excel("SRfinal.xlsx", sheet = "4b", range = "A1:B72")
SR4c <- read_excel("SRfinal.xlsx", sheet = "4c", range = "A1:B72")
SR4d <- read_excel("SRfinal.xlsx", sheet = "4d", range = "A1:B72")
SR4e <- read_excel("SRfinal.xlsx", sheet = "4e", range = "A1:B72")
SR4f <- read_excel("SRfinal.xlsx", sheet = "4f", range = "A1:B72")

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age <- SR1a %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age <- SR1b %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age <- SR1c %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age <- SR1d %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age <- SR1e %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age <- SR1f %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age <- SR2a %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age <- SR2b %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age <- SR2c %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age <- SR2d %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age <- SR2e %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age <- SR2f %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age <- SR3a %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age <- SR3b %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age <- SR3c %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age <- SR3d %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age <- SR3e %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age <- SR3f %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age <- SR4a %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age <- SR4b %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age <- SR4c %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age <- SR4d %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age <- SR4e %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age <- SR4f %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

## Extract age-specific relativa survival 
v_r_mort_by_age <- LifeTable_Chile %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()


# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage <- rep(v_r_mort_by_age, each = 1/cycle_length)
v_r_1aDage <- rep(v_r_sr1a_by_age, each = 1/cycle_length)
v_r_1bDage <- rep(v_r_sr1b_by_age, each = 1/cycle_length)
v_r_1cDage <- rep(v_r_sr1c_by_age, each = 1/cycle_length)
v_r_1dDage <- rep(v_r_sr1d_by_age, each = 1/cycle_length)
v_r_1eDage <- rep(v_r_sr1e_by_age, each = 1/cycle_length)
v_r_1fDage <- rep(v_r_sr1f_by_age, each = 1/cycle_length)
v_r_2aDage <- rep(v_r_sr2a_by_age, each = 1/cycle_length)
v_r_2bDage <- rep(v_r_sr2b_by_age, each = 1/cycle_length)
v_r_2cDage <- rep(v_r_sr2c_by_age, each = 1/cycle_length)
v_r_2dDage <- rep(v_r_sr2d_by_age, each = 1/cycle_length)
v_r_2eDage <- rep(v_r_sr2e_by_age, each = 1/cycle_length)
v_r_2fDage <- rep(v_r_sr2f_by_age, each = 1/cycle_length)
v_r_3aDage <- rep(v_r_sr3a_by_age, each = 1/cycle_length)
v_r_3bDage <- rep(v_r_sr3b_by_age, each = 1/cycle_length)
v_r_3cDage <- rep(v_r_sr3c_by_age, each = 1/cycle_length)
v_r_3dDage <- rep(v_r_sr3d_by_age, each = 1/cycle_length)
v_r_3eDage <- rep(v_r_sr3e_by_age, each = 1/cycle_length)
v_r_3fDage <- rep(v_r_sr3f_by_age, each = 1/cycle_length)
v_r_4aDage <- rep(v_r_sr4a_by_age, each = 1/cycle_length)
v_r_4bDage <- rep(v_r_sr4b_by_age, each = 1/cycle_length)
v_r_4cDage <- rep(v_r_sr4c_by_age, each = 1/cycle_length)
v_r_4dDage <- rep(v_r_sr4d_by_age, each = 1/cycle_length)
v_r_4eDage <- rep(v_r_sr4e_by_age, each = 1/cycle_length)
v_r_4fDage <- rep(v_r_sr4f_by_age, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage <- 1 - exp(-v_r_HDage * cycle_length)
v_p_1aDage <- v_p_HDage + (1 - exp(-v_r_1aDage * cycle_length))
v_p_1bDage <- v_p_HDage + (1 - exp(-v_r_1bDage * cycle_length))
v_p_1cDage <- v_p_HDage + (1 - exp(-v_r_1cDage * cycle_length))
v_p_1dDage <- v_p_HDage + (1 - exp(-v_r_1dDage * cycle_length))
v_p_1eDage <- v_p_HDage + (1 - exp(-v_r_1eDage * cycle_length))
v_p_1fDage <- v_p_HDage + (1 - exp(-v_r_1fDage * cycle_length))
v_p_2aDage <- v_p_HDage + (1 - exp(-v_r_2aDage * cycle_length))
v_p_2bDage <- v_p_HDage + (1 - exp(-v_r_2bDage * cycle_length))
v_p_2cDage <- v_p_HDage + (1 - exp(-v_r_2cDage * cycle_length))
v_p_2dDage <- v_p_HDage + (1 - exp(-v_r_2dDage * cycle_length))
v_p_2eDage <- v_p_HDage + (1 - exp(-v_r_2eDage * cycle_length))
v_p_2fDage <- v_p_HDage + (1 - exp(-v_r_2fDage * cycle_length))
v_p_3aDage <- v_p_HDage + (1 - exp(-v_r_3aDage * cycle_length))
v_p_3bDage <- v_p_HDage + (1 - exp(-v_r_3bDage * cycle_length))
v_p_3cDage <- v_p_HDage + (1 - exp(-v_r_3cDage * cycle_length))
v_p_3dDage <- v_p_HDage + (1 - exp(-v_r_3dDage * cycle_length))
v_p_3eDage <- v_p_HDage + (1 - exp(-v_r_3eDage * cycle_length))
v_p_3fDage <- v_p_HDage + (1 - exp(-v_r_3fDage * cycle_length))
v_p_4aDage <- v_p_HDage + (1 - exp(-v_r_4aDage * cycle_length))
v_p_4bDage <- v_p_HDage + (1 - exp(-v_r_4bDage * cycle_length))
v_p_4cDage <- v_p_HDage + (1 - exp(-v_r_4cDage * cycle_length))
v_p_4dDage <- v_p_HDage + (1 - exp(-v_r_4dDage * cycle_length))
v_p_4eDage <- v_p_HDage + (1 - exp(-v_r_4eDage * cycle_length))
v_p_4fDage <- v_p_HDage + (1 - exp(-v_r_4fDage * cycle_length))

## Probabilidades de muerte específica por edad 

## Load hazards from the relative survival analysis EDADES 31 -100
SR1a_31 <- read_excel("SR31.xlsx", sheet = "1a", range = "A1:B71")
SR1b_31 <- read_excel("SR31.xlsx", sheet = "1b", range = "A1:B71")
SR1c_31 <- read_excel("SR31.xlsx", sheet = "1c", range = "A1:B71")
SR1d_31 <- read_excel("SR31.xlsx", sheet = "1d", range = "A1:B71")
SR1e_31 <- read_excel("SR31.xlsx", sheet = "1e", range = "A1:B71")
SR1f_31 <- read_excel("SR31.xlsx", sheet = "1f", range = "A1:B71")

SR2a_31 <- read_excel("SR31.xlsx", sheet = "2a", range = "A1:B71")
SR2b_31 <- read_excel("SR31.xlsx", sheet = "2b", range = "A1:B71")
SR2c_31 <- read_excel("SR31.xlsx", sheet = "2c", range = "A1:B71")
SR2d_31 <- read_excel("SR31.xlsx", sheet = "2d", range = "A1:B71")
SR2e_31 <- read_excel("SR31.xlsx", sheet = "2e", range = "A1:B71")
SR2f_31 <- read_excel("SR31.xlsx", sheet = "2f", range = "A1:B71")

SR3a_31 <- read_excel("SR31.xlsx", sheet = "3a", range = "A1:B71")
SR3b_31 <- read_excel("SR31.xlsx", sheet = "3b", range = "A1:B71")
SR3c_31 <- read_excel("SR31.xlsx", sheet = "3c", range = "A1:B71")
SR3d_31 <- read_excel("SR31.xlsx", sheet = "3d", range = "A1:B71")
SR3e_31 <- read_excel("SR31.xlsx", sheet = "3e", range = "A1:B71")
SR3f_31 <- read_excel("SR31.xlsx", sheet = "3f", range = "A1:B71")

SR4a_31 <- read_excel("SR31.xlsx", sheet = "4a", range = "A1:B71")
SR4b_31 <- read_excel("SR31.xlsx", sheet = "4b", range = "A1:B71")
SR4c_31 <- read_excel("SR31.xlsx", sheet = "4c", range = "A1:B71")
SR4d_31 <- read_excel("SR31.xlsx", sheet = "4d", range = "A1:B71")
SR4e_31 <- read_excel("SR31.xlsx", sheet = "4e", range = "A1:B71")
SR4f_31 <- read_excel("SR31.xlsx", sheet = "4f", range = "A1:B71")

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age31 <- SR1a_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age31 <- SR1b_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age31 <- SR1c_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age31 <- SR1d_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age31 <- SR1e_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age31 <- SR1f_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age31 <- SR2a_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age31 <- SR2b_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age31 <- SR2c_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age31 <- SR2d_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age31 <- SR2e_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age31 <- SR2f_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age31 <- SR3a_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age31 <- SR3b_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age31 <- SR3c_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age31 <- SR3d_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age31 <- SR3e_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age31 <- SR3f_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age31 <- SR4a_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age31 <- SR4b_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age31 <- SR4c_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age31 <- SR4d_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age31 <- SR4e_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age31 <- SR4f_31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

##load table 31 
LifeTable_Chile31 <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile31.csv")

## Extract age-specific relativa survival 
v_r_mort_by_age31 <- LifeTable_Chile31 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()

# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage31 <- rep(v_r_mort_by_age31, each = 1/cycle_length)
v_r_1aDage31 <- rep(v_r_sr1a_by_age31, each = 1/cycle_length)
v_r_1bDage31 <- rep(v_r_sr1b_by_age31, each = 1/cycle_length)
v_r_1cDage31 <- rep(v_r_sr1c_by_age31, each = 1/cycle_length)
v_r_1dDage31 <- rep(v_r_sr1d_by_age31, each = 1/cycle_length)
v_r_1eDage31 <- rep(v_r_sr1e_by_age31, each = 1/cycle_length)
v_r_1fDage31 <- rep(v_r_sr1f_by_age31, each = 1/cycle_length)
v_r_2aDage31 <- rep(v_r_sr2a_by_age31, each = 1/cycle_length)
v_r_2bDage31 <- rep(v_r_sr2b_by_age31, each = 1/cycle_length)
v_r_2cDage31 <- rep(v_r_sr2c_by_age31, each = 1/cycle_length)
v_r_2dDage31 <- rep(v_r_sr2d_by_age31, each = 1/cycle_length)
v_r_2eDage31 <- rep(v_r_sr2e_by_age31, each = 1/cycle_length)
v_r_2fDage31 <- rep(v_r_sr2f_by_age31, each = 1/cycle_length)
v_r_3aDage31 <- rep(v_r_sr3a_by_age31, each = 1/cycle_length)
v_r_3bDage31 <- rep(v_r_sr3b_by_age31, each = 1/cycle_length)
v_r_3cDage31 <- rep(v_r_sr3c_by_age31, each = 1/cycle_length)
v_r_3dDage31 <- rep(v_r_sr3d_by_age31, each = 1/cycle_length)
v_r_3eDage31 <- rep(v_r_sr3e_by_age31, each = 1/cycle_length)
v_r_3fDage31 <- rep(v_r_sr3f_by_age31, each = 1/cycle_length)
v_r_4aDage31 <- rep(v_r_sr4a_by_age31, each = 1/cycle_length)
v_r_4bDage31 <- rep(v_r_sr4b_by_age31, each = 1/cycle_length)
v_r_4cDage31 <- rep(v_r_sr4c_by_age31, each = 1/cycle_length)
v_r_4dDage31 <- rep(v_r_sr4d_by_age31, each = 1/cycle_length)
v_r_4eDage31 <- rep(v_r_sr4e_by_age31, each = 1/cycle_length)
v_r_4fDage31 <- rep(v_r_sr4f_by_age31, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage31 <- 1 - exp(-v_r_HDage31 * cycle_length)
v_p_1aDage31 <- v_p_HDage31 + (1 - exp(-v_r_1aDage31 * cycle_length))
v_p_1bDage31 <- v_p_HDage31 + (1 - exp(-v_r_1bDage31 * cycle_length))
v_p_1cDage31 <- v_p_HDage31 + (1 - exp(-v_r_1cDage31 * cycle_length))
v_p_1dDage31 <- v_p_HDage31 + (1 - exp(-v_r_1dDage31 * cycle_length))
v_p_1eDage31 <- v_p_HDage31 + (1 - exp(-v_r_1eDage31 * cycle_length))
v_p_1fDage31 <- v_p_HDage31 + (1 - exp(-v_r_1fDage31 * cycle_length))
v_p_2aDage31 <- v_p_HDage31 + (1 - exp(-v_r_2aDage31 * cycle_length))
v_p_2bDage31 <- v_p_HDage31 + (1 - exp(-v_r_2bDage31 * cycle_length))
v_p_2cDage31 <- v_p_HDage31 + (1 - exp(-v_r_2cDage31 * cycle_length))
v_p_2dDage31 <- v_p_HDage31 + (1 - exp(-v_r_2dDage31 * cycle_length))
v_p_2eDage31 <- v_p_HDage31 + (1 - exp(-v_r_2eDage31 * cycle_length))
v_p_2fDage31 <- v_p_HDage31 + (1 - exp(-v_r_2fDage31 * cycle_length))
v_p_3aDage31 <- v_p_HDage31 + (1 - exp(-v_r_3aDage31 * cycle_length))
v_p_3bDage31 <- v_p_HDage31 + (1 - exp(-v_r_3bDage31 * cycle_length))
v_p_3cDage31 <- v_p_HDage31 + (1 - exp(-v_r_3cDage31 * cycle_length))
v_p_3dDage31 <- v_p_HDage31 + (1 - exp(-v_r_3dDage31 * cycle_length))
v_p_3eDage31 <- v_p_HDage31 + (1 - exp(-v_r_3eDage31 * cycle_length))
v_p_3fDage31 <- v_p_HDage31 + (1 - exp(-v_r_3fDage31 * cycle_length))
v_p_4aDage31 <- v_p_HDage31 + (1 - exp(-v_r_4aDage31 * cycle_length))
v_p_4bDage31 <- v_p_HDage31 + (1 - exp(-v_r_4bDage31 * cycle_length))
v_p_4cDage31 <- v_p_HDage31 + (1 - exp(-v_r_4cDage31 * cycle_length))
v_p_4dDage31 <- v_p_HDage31 + (1 - exp(-v_r_4dDage31 * cycle_length))
v_p_4eDage31 <- v_p_HDage31 + (1 - exp(-v_r_4eDage31 * cycle_length))
v_p_4fDage31 <- v_p_HDage31 + (1 - exp(-v_r_4fDage31 * cycle_length))

### Para 30 - 44 años 
## Load hazards from the relative survival analysis EDADES 30 -44
SR1a_3044 <- read_excel("SRfinal.xlsx", sheet = "1a", range = "A1:B16")
SR1b_3044 <- read_excel("SRfinal.xlsx", sheet = "1b", range = "A1:B16")
SR1c_3044 <- read_excel("SRfinal.xlsx", sheet = "1c", range = "A1:B16")
SR1d_3044 <- read_excel("SRfinal.xlsx", sheet = "1d", range = "A1:B16")
SR1e_3044 <- read_excel("SRfinal.xlsx", sheet = "1e", range = "A1:B16")
SR1f_3044 <- read_excel("SRfinal.xlsx", sheet = "1f", range = "A1:B16")

SR2a_3044 <- read_excel("SRfinal.xlsx", sheet = "2a", range = "A1:B16")
SR2b_3044 <- read_excel("SRfinal.xlsx", sheet = "2b", range = "A1:B16")
SR2c_3044 <- read_excel("SRfinal.xlsx", sheet = "2c", range = "A1:B16")
SR2d_3044 <- read_excel("SRfinal.xlsx", sheet = "2d", range = "A1:B16")
SR2e_3044 <- read_excel("SRfinal.xlsx", sheet = "2e", range = "A1:B16")
SR2f_3044 <- read_excel("SRfinal.xlsx", sheet = "2f", range = "A1:B16")

SR3a_3044 <- read_excel("SRfinal.xlsx", sheet = "3a", range = "A1:B16")
SR3b_3044 <- read_excel("SRfinal.xlsx", sheet = "3b", range = "A1:B16")
SR3c_3044 <- read_excel("SRfinal.xlsx", sheet = "3c", range = "A1:B16")
SR3d_3044 <- read_excel("SRfinal.xlsx", sheet = "3d", range = "A1:B16")
SR3e_3044 <- read_excel("SRfinal.xlsx", sheet = "3e", range = "A1:B16")
SR3f_3044 <- read_excel("SRfinal.xlsx", sheet = "3f", range = "A1:B16")

SR4a_3044 <- read_excel("SRfinal.xlsx", sheet = "4a", range = "A1:B16")
SR4b_3044 <- read_excel("SRfinal.xlsx", sheet = "4b", range = "A1:B16")
SR4c_3044 <- read_excel("SRfinal.xlsx", sheet = "4c", range = "A1:B16")
SR4d_3044 <- read_excel("SRfinal.xlsx", sheet = "4d", range = "A1:B16")
SR4e_3044 <- read_excel("SRfinal.xlsx", sheet = "4e", range = "A1:B16")
SR4f_3044 <- read_excel("SRfinal.xlsx", sheet = "4f", range = "A1:B16")

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age3044 <- SR1a_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age3044 <- SR1b_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age3044 <- SR1c_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age3044 <- SR1d_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age3044 <- SR1e_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age3044 <- SR1f_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age3044 <- SR2a_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age3044 <- SR2b_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age3044 <- SR2c_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age3044 <- SR2d_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age3044 <- SR2e_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age3044 <- SR2f_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age3044 <- SR3a_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age3044 <- SR3b_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age3044 <- SR3c_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age3044 <- SR3d_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age3044 <- SR3e_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age3044 <- SR3f_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age3044 <- SR4a_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age3044 <- SR4b_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age3044 <- SR4c_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age3044 <- SR4d_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age3044 <- SR4e_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age3044 <- SR4f_3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

##load table 30-44 
LifeTable_Chile3044 <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile3044.csv")

## Extract age-specific relativa survival 
v_r_mort_by_age3044 <- LifeTable_Chile3044 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()

# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage3044 <- rep(v_r_mort_by_age3044, each = 1/cycle_length)
v_r_1aDage3044 <- rep(v_r_sr1a_by_age3044, each = 1/cycle_length)
v_r_1bDage3044 <- rep(v_r_sr1b_by_age3044, each = 1/cycle_length)
v_r_1cDage3044 <- rep(v_r_sr1c_by_age3044, each = 1/cycle_length)
v_r_1dDage3044 <- rep(v_r_sr1d_by_age3044, each = 1/cycle_length)
v_r_1eDage3044 <- rep(v_r_sr1e_by_age3044, each = 1/cycle_length)
v_r_1fDage3044 <- rep(v_r_sr1f_by_age3044, each = 1/cycle_length)
v_r_2aDage3044 <- rep(v_r_sr2a_by_age3044, each = 1/cycle_length)
v_r_2bDage3044 <- rep(v_r_sr2b_by_age3044, each = 1/cycle_length)
v_r_2cDage3044 <- rep(v_r_sr2c_by_age3044, each = 1/cycle_length)
v_r_2dDage3044 <- rep(v_r_sr2d_by_age3044, each = 1/cycle_length)
v_r_2eDage3044 <- rep(v_r_sr2e_by_age3044, each = 1/cycle_length)
v_r_2fDage3044 <- rep(v_r_sr2f_by_age3044, each = 1/cycle_length)
v_r_3aDage3044 <- rep(v_r_sr3a_by_age3044, each = 1/cycle_length)
v_r_3bDage3044 <- rep(v_r_sr3b_by_age3044, each = 1/cycle_length)
v_r_3cDage3044 <- rep(v_r_sr3c_by_age3044, each = 1/cycle_length)
v_r_3dDage3044 <- rep(v_r_sr3d_by_age3044, each = 1/cycle_length)
v_r_3eDage3044 <- rep(v_r_sr3e_by_age3044, each = 1/cycle_length)
v_r_3fDage3044 <- rep(v_r_sr3f_by_age3044, each = 1/cycle_length)
v_r_4aDage3044 <- rep(v_r_sr4a_by_age3044, each = 1/cycle_length)
v_r_4bDage3044 <- rep(v_r_sr4b_by_age3044, each = 1/cycle_length)
v_r_4cDage3044 <- rep(v_r_sr4c_by_age3044, each = 1/cycle_length)
v_r_4dDage3044 <- rep(v_r_sr4d_by_age3044, each = 1/cycle_length)
v_r_4eDage3044 <- rep(v_r_sr4e_by_age3044, each = 1/cycle_length)
v_r_4fDage3044 <- rep(v_r_sr4f_by_age3044, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage3044 <- 1 - exp(-v_r_HDage3044 * cycle_length)
v_p_1aDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1aDage3044 * cycle_length))
v_p_1bDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1bDage3044 * cycle_length))
v_p_1cDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1cDage3044 * cycle_length))
v_p_1dDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1dDage3044 * cycle_length))
v_p_1eDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1eDage3044 * cycle_length))
v_p_1fDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_1fDage3044 * cycle_length))
v_p_2aDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2aDage3044 * cycle_length))
v_p_2bDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2bDage3044 * cycle_length))
v_p_2cDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2cDage3044 * cycle_length))
v_p_2dDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2dDage3044 * cycle_length))
v_p_2eDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2eDage3044 * cycle_length))
v_p_2fDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_2fDage3044 * cycle_length))
v_p_3aDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3aDage3044 * cycle_length))
v_p_3bDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3bDage3044 * cycle_length))
v_p_3cDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3cDage3044 * cycle_length))
v_p_3dDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3dDage3044 * cycle_length))
v_p_3eDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3eDage3044 * cycle_length))
v_p_3fDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_3fDage3044 * cycle_length))
v_p_4aDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4aDage3044 * cycle_length))
v_p_4bDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4bDage3044 * cycle_length))
v_p_4cDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4cDage3044 * cycle_length))
v_p_4dDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4dDage3044 * cycle_length))
v_p_4eDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4eDage3044 * cycle_length))
v_p_4fDage3044 <- v_p_HDage3044 + (1 - exp(-v_r_4fDage3044 * cycle_length))



## Para 46 - 100 
## Load hazards from the relative survival analysis EDADES 46 -100

v1 <- c("Age", "H")

SR1a_46 <- read_excel("SRfinal.xlsx", sheet = "1a", range = "A18:B72", col_names = v1)
SR1b_46 <- read_excel("SRfinal.xlsx", sheet = "1b", range = "A18:B72", col_names = v1)
SR1c_46 <- read_excel("SRfinal.xlsx", sheet = "1c", range = "A18:B72", col_names = v1)
SR1d_46 <- read_excel("SRfinal.xlsx", sheet = "1d", range = "A18:B72", col_names = v1)
SR1e_46 <- read_excel("SRfinal.xlsx", sheet = "1e", range = "A18:B72", col_names = v1)
SR1f_46 <- read_excel("SRfinal.xlsx", sheet = "1f", range = "A18:B72", col_names = v1)

SR2a_46 <- read_excel("SRfinal.xlsx", sheet = "2a", range = "A18:B72", col_names = v1)
SR2b_46 <- read_excel("SRfinal.xlsx", sheet = "2b", range = "A18:B72", col_names = v1)
SR2c_46 <- read_excel("SRfinal.xlsx", sheet = "2c", range = "A18:B72", col_names = v1)
SR2d_46 <- read_excel("SRfinal.xlsx", sheet = "2d", range = "A18:B72", col_names = v1)
SR2e_46 <- read_excel("SRfinal.xlsx", sheet = "2e", range = "A18:B72", col_names = v1)
SR2f_46 <- read_excel("SRfinal.xlsx", sheet = "2f", range = "A18:B72", col_names = v1)

SR3a_46 <- read_excel("SRfinal.xlsx", sheet = "3a", range = "A18:B72", col_names = v1)
SR3b_46 <- read_excel("SRfinal.xlsx", sheet = "3b", range = "A18:B72", col_names = v1)
SR3c_46 <- read_excel("SRfinal.xlsx", sheet = "3c", range = "A18:B72", col_names = v1)
SR3d_46 <- read_excel("SRfinal.xlsx", sheet = "3d", range = "A18:B72", col_names = v1)
SR3e_46 <- read_excel("SRfinal.xlsx", sheet = "3e", range = "A18:B72", col_names = v1)
SR3f_46 <- read_excel("SRfinal.xlsx", sheet = "3f", range = "A18:B72", col_names = v1)
SR4a_46 <- read_excel("SRfinal.xlsx", sheet = "4a", range = "A18:B72", col_names = v1)
SR4b_46 <- read_excel("SRfinal.xlsx", sheet = "4b", range = "A18:B72", col_names = v1)
SR4c_46 <- read_excel("SRfinal.xlsx", sheet = "4c", range = "A18:B72", col_names = v1)
SR4d_46 <- read_excel("SRfinal.xlsx", sheet = "4d", range = "A18:B72", col_names = v1)
SR4e_46 <- read_excel("SRfinal.xlsx", sheet = "4e", range = "A18:B72", col_names = v1)
SR4f_46 <- read_excel("SRfinal.xlsx", sheet = "4f", range = "A18:B72", col_names = v1)

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age46 <- SR1a_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age46 <- SR1b_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age46 <- SR1c_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age46 <- SR1d_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age46 <- SR1e_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age46 <- SR1f_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age46 <- SR2a_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age46 <- SR2b_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age46 <- SR2c_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age46 <- SR2d_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age46 <- SR2e_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age46 <- SR2f_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age46 <- SR3a_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age46 <- SR3b_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age46 <- SR3c_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age46 <- SR3d_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age46 <- SR3e_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age46 <- SR3f_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age46 <- SR4a_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age46 <- SR4b_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age46 <- SR4c_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age46 <- SR4d_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age46 <- SR4e_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age46 <- SR4f_46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

##load table 31 
LifeTable_Chile46 <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile46.csv")

## Extract age-specific relativa survival 
v_r_mort_by_age46 <- LifeTable_Chile46 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()

# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage46 <- rep(v_r_mort_by_age46, each = 1/cycle_length)
v_r_1aDage46 <- rep(v_r_sr1a_by_age46, each = 1/cycle_length)
v_r_1bDage46 <- rep(v_r_sr1b_by_age46, each = 1/cycle_length)
v_r_1cDage46 <- rep(v_r_sr1c_by_age46, each = 1/cycle_length)
v_r_1dDage46 <- rep(v_r_sr1d_by_age46, each = 1/cycle_length)
v_r_1eDage46 <- rep(v_r_sr1e_by_age46, each = 1/cycle_length)
v_r_1fDage46 <- rep(v_r_sr1f_by_age46, each = 1/cycle_length)
v_r_2aDage46 <- rep(v_r_sr2a_by_age46, each = 1/cycle_length)
v_r_2bDage46 <- rep(v_r_sr2b_by_age46, each = 1/cycle_length)
v_r_2cDage46 <- rep(v_r_sr2c_by_age46, each = 1/cycle_length)
v_r_2dDage46 <- rep(v_r_sr2d_by_age46, each = 1/cycle_length)
v_r_2eDage46 <- rep(v_r_sr2e_by_age46, each = 1/cycle_length)
v_r_2fDage46 <- rep(v_r_sr2f_by_age46, each = 1/cycle_length)
v_r_3aDage46 <- rep(v_r_sr3a_by_age46, each = 1/cycle_length)
v_r_3bDage46 <- rep(v_r_sr3b_by_age46, each = 1/cycle_length)
v_r_3cDage46 <- rep(v_r_sr3c_by_age46, each = 1/cycle_length)
v_r_3dDage46 <- rep(v_r_sr3d_by_age46, each = 1/cycle_length)
v_r_3eDage46 <- rep(v_r_sr3e_by_age46, each = 1/cycle_length)
v_r_3fDage46 <- rep(v_r_sr3f_by_age46, each = 1/cycle_length)
v_r_4aDage46 <- rep(v_r_sr4a_by_age46, each = 1/cycle_length)
v_r_4bDage46 <- rep(v_r_sr4b_by_age46, each = 1/cycle_length)
v_r_4cDage46 <- rep(v_r_sr4c_by_age46, each = 1/cycle_length)
v_r_4dDage46 <- rep(v_r_sr4d_by_age46, each = 1/cycle_length)
v_r_4eDage46 <- rep(v_r_sr4e_by_age46, each = 1/cycle_length)
v_r_4fDage46 <- rep(v_r_sr4f_by_age46, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage46 <- 1 - exp(-v_r_HDage46 * cycle_length)
v_p_1aDage46 <- v_p_HDage46 + (1 - exp(-v_r_1aDage46 * cycle_length))
v_p_1bDage46 <- v_p_HDage46 + (1 - exp(-v_r_1bDage46 * cycle_length))
v_p_1cDage46 <- v_p_HDage46 + (1 - exp(-v_r_1cDage46 * cycle_length))
v_p_1dDage46 <- v_p_HDage46 + (1 - exp(-v_r_1dDage46 * cycle_length))
v_p_1eDage46 <- v_p_HDage46 + (1 - exp(-v_r_1eDage46 * cycle_length))
v_p_1fDage46 <- v_p_HDage46 + (1 - exp(-v_r_1fDage46 * cycle_length))
v_p_2aDage46 <- v_p_HDage46 + (1 - exp(-v_r_2aDage46 * cycle_length))
v_p_2bDage46 <- v_p_HDage46 + (1 - exp(-v_r_2bDage46 * cycle_length))
v_p_2cDage46 <- v_p_HDage46 + (1 - exp(-v_r_2cDage46 * cycle_length))
v_p_2dDage46 <- v_p_HDage46 + (1 - exp(-v_r_2dDage46 * cycle_length))
v_p_2eDage46 <- v_p_HDage46 + (1 - exp(-v_r_2eDage46 * cycle_length))
v_p_2fDage46 <- v_p_HDage46 + (1 - exp(-v_r_2fDage46 * cycle_length))
v_p_3aDage46 <- v_p_HDage46 + (1 - exp(-v_r_3aDage46 * cycle_length))
v_p_3bDage46 <- v_p_HDage46 + (1 - exp(-v_r_3bDage46 * cycle_length))
v_p_3cDage46 <- v_p_HDage46 + (1 - exp(-v_r_3cDage46 * cycle_length))
v_p_3dDage46 <- v_p_HDage46 + (1 - exp(-v_r_3dDage46 * cycle_length))
v_p_3eDage46 <- v_p_HDage46 + (1 - exp(-v_r_3eDage46 * cycle_length))
v_p_3fDage46 <- v_p_HDage46 + (1 - exp(-v_r_3fDage46 * cycle_length))
v_p_4aDage46 <- v_p_HDage46 + (1 - exp(-v_r_4aDage46 * cycle_length))
v_p_4bDage46 <- v_p_HDage46 + (1 - exp(-v_r_4bDage46 * cycle_length))
v_p_4cDage46 <- v_p_HDage46 + (1 - exp(-v_r_4cDage46 * cycle_length))
v_p_4dDage46 <- v_p_HDage46 + (1 - exp(-v_r_4dDage46 * cycle_length))
v_p_4eDage46 <- v_p_HDage46 + (1 - exp(-v_r_4eDage46 * cycle_length))
v_p_4fDage46 <- v_p_HDage46 + (1 - exp(-v_r_4fDage46 * cycle_length))


### Para 30 - 39 años 
## Load hazards from the relative survival analysis EDADES 30 -44
SR1a_3039 <- read_excel("SRfinal.xlsx", sheet = "1a", range = "A1:B11")
SR1b_3039 <- read_excel("SRfinal.xlsx", sheet = "1b", range = "A1:B11")
SR1c_3039 <- read_excel("SRfinal.xlsx", sheet = "1c", range = "A1:B11")
SR1d_3039 <- read_excel("SRfinal.xlsx", sheet = "1d", range = "A1:B11")
SR1e_3039 <- read_excel("SRfinal.xlsx", sheet = "1e", range = "A1:B11")
SR1f_3039 <- read_excel("SRfinal.xlsx", sheet = "1f", range = "A1:B11")

SR2a_3039 <- read_excel("SRfinal.xlsx", sheet = "2a", range = "A1:B11")
SR2b_3039 <- read_excel("SRfinal.xlsx", sheet = "2b", range = "A1:B11")
SR2c_3039 <- read_excel("SRfinal.xlsx", sheet = "2c", range = "A1:B11")
SR2d_3039 <- read_excel("SRfinal.xlsx", sheet = "2d", range = "A1:B11")
SR2e_3039 <- read_excel("SRfinal.xlsx", sheet = "2e", range = "A1:B11")
SR2f_3039 <- read_excel("SRfinal.xlsx", sheet = "2f", range = "A1:B11")

SR3a_3039 <- read_excel("SRfinal.xlsx", sheet = "3a", range = "A1:B11")
SR3b_3039 <- read_excel("SRfinal.xlsx", sheet = "3b", range = "A1:B11")
SR3c_3039 <- read_excel("SRfinal.xlsx", sheet = "3c", range = "A1:B11")
SR3d_3039 <- read_excel("SRfinal.xlsx", sheet = "3d", range = "A1:B11")
SR3e_3039 <- read_excel("SRfinal.xlsx", sheet = "3e", range = "A1:B11")
SR3f_3039 <- read_excel("SRfinal.xlsx", sheet = "3f", range = "A1:B11")

SR4a_3039 <- read_excel("SRfinal.xlsx", sheet = "4a", range = "A1:B11")
SR4b_3039 <- read_excel("SRfinal.xlsx", sheet = "4b", range = "A1:B11")
SR4c_3039 <- read_excel("SRfinal.xlsx", sheet = "4c", range = "A1:B11")
SR4d_3039 <- read_excel("SRfinal.xlsx", sheet = "4d", range = "A1:B11")
SR4e_3039 <- read_excel("SRfinal.xlsx", sheet = "4e", range = "A1:B11")
SR4f_3039 <- read_excel("SRfinal.xlsx", sheet = "4f", range = "A1:B11")

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age3039 <- SR1a_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age3039 <- SR1b_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age3039 <- SR1c_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age3039 <- SR1d_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age3039 <- SR1e_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age3039 <- SR1f_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age3039 <- SR2a_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age3039 <- SR2b_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age3039 <- SR2c_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age3039 <- SR2d_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age3039 <- SR2e_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age3039 <- SR2f_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age3039 <- SR3a_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age3039 <- SR3b_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age3039 <- SR3c_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age3039 <- SR3d_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age3039 <- SR3e_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age3039 <- SR3f_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age3039 <- SR4a_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age3039 <- SR4b_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age3039 <- SR4c_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age3039 <- SR4d_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age3039 <- SR4e_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age3039 <- SR4f_3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

##load table 30-39 
LifeTable_Chile3039 <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile3039.csv")

## Extract age-specific relativa survival 
v_r_mort_by_age3039 <- LifeTable_Chile3039 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()

# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage3039 <- rep(v_r_mort_by_age3039, each = 1/cycle_length)
v_r_1aDage3039 <- rep(v_r_sr1a_by_age3039, each = 1/cycle_length)
v_r_1bDage3039 <- rep(v_r_sr1b_by_age3039, each = 1/cycle_length)
v_r_1cDage3039 <- rep(v_r_sr1c_by_age3039, each = 1/cycle_length)
v_r_1dDage3039 <- rep(v_r_sr1d_by_age3039, each = 1/cycle_length)
v_r_1eDage3039 <- rep(v_r_sr1e_by_age3039, each = 1/cycle_length)
v_r_1fDage3039 <- rep(v_r_sr1f_by_age3039, each = 1/cycle_length)
v_r_2aDage3039 <- rep(v_r_sr2a_by_age3039, each = 1/cycle_length)
v_r_2bDage3039 <- rep(v_r_sr2b_by_age3039, each = 1/cycle_length)
v_r_2cDage3039 <- rep(v_r_sr2c_by_age3039, each = 1/cycle_length)
v_r_2dDage3039 <- rep(v_r_sr2d_by_age3039, each = 1/cycle_length)
v_r_2eDage3039 <- rep(v_r_sr2e_by_age3039, each = 1/cycle_length)
v_r_2fDage3039 <- rep(v_r_sr2f_by_age3039, each = 1/cycle_length)
v_r_3aDage3039 <- rep(v_r_sr3a_by_age3039, each = 1/cycle_length)
v_r_3bDage3039 <- rep(v_r_sr3b_by_age3039, each = 1/cycle_length)
v_r_3cDage3039 <- rep(v_r_sr3c_by_age3039, each = 1/cycle_length)
v_r_3dDage3039 <- rep(v_r_sr3d_by_age3039, each = 1/cycle_length)
v_r_3eDage3039 <- rep(v_r_sr3e_by_age3039, each = 1/cycle_length)
v_r_3fDage3039 <- rep(v_r_sr3f_by_age3039, each = 1/cycle_length)
v_r_4aDage3039 <- rep(v_r_sr4a_by_age3039, each = 1/cycle_length)
v_r_4bDage3039 <- rep(v_r_sr4b_by_age3039, each = 1/cycle_length)
v_r_4cDage3039 <- rep(v_r_sr4c_by_age3039, each = 1/cycle_length)
v_r_4dDage3039 <- rep(v_r_sr4d_by_age3039, each = 1/cycle_length)
v_r_4eDage3039 <- rep(v_r_sr4e_by_age3039, each = 1/cycle_length)
v_r_4fDage3039 <- rep(v_r_sr4f_by_age3039, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage3039 <- 1 - exp(-v_r_HDage3039 * cycle_length)
v_p_1aDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1aDage3039 * cycle_length))
v_p_1bDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1bDage3039 * cycle_length))
v_p_1cDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1cDage3039 * cycle_length))
v_p_1dDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1dDage3039 * cycle_length))
v_p_1eDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1eDage3039 * cycle_length))
v_p_1fDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_1fDage3039 * cycle_length))
v_p_2aDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2aDage3039 * cycle_length))
v_p_2bDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2bDage3039 * cycle_length))
v_p_2cDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2cDage3039 * cycle_length))
v_p_2dDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2dDage3039 * cycle_length))
v_p_2eDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2eDage3039 * cycle_length))
v_p_2fDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_2fDage3039 * cycle_length))
v_p_3aDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3aDage3039 * cycle_length))
v_p_3bDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3bDage3039 * cycle_length))
v_p_3cDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3cDage3039 * cycle_length))
v_p_3dDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3dDage3039 * cycle_length))
v_p_3eDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3eDage3039 * cycle_length))
v_p_3fDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_3fDage3039 * cycle_length))
v_p_4aDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4aDage3039 * cycle_length))
v_p_4bDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4bDage3039 * cycle_length))
v_p_4cDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4cDage3039 * cycle_length))
v_p_4dDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4dDage3039 * cycle_length))
v_p_4eDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4eDage3039 * cycle_length))
v_p_4fDage3039 <- v_p_HDage3039 + (1 - exp(-v_r_4fDage3039 * cycle_length))


## Para 41 - 100 
## Load hazards from the relative survival analysis EDADES 41 -100

SR1a_41 <- read_excel("SRfinal.xlsx", sheet = "1a", range = "A13:B72", col_names = v1)
SR1b_41 <- read_excel("SRfinal.xlsx", sheet = "1b", range = "A13:B72", col_names = v1)
SR1c_41 <- read_excel("SRfinal.xlsx", sheet = "1c", range = "A13:B72", col_names = v1)
SR1d_41 <- read_excel("SRfinal.xlsx", sheet = "1d", range = "A13:B72", col_names = v1)
SR1e_41 <- read_excel("SRfinal.xlsx", sheet = "1e", range = "A13:B72", col_names = v1)
SR1f_41 <- read_excel("SRfinal.xlsx", sheet = "1f", range = "A13:B72", col_names = v1)

SR2a_41 <- read_excel("SRfinal.xlsx", sheet = "2a", range = "A13:B72", col_names = v1)
SR2b_41 <- read_excel("SRfinal.xlsx", sheet = "2b", range = "A13:B72", col_names = v1)
SR2c_41 <- read_excel("SRfinal.xlsx", sheet = "2c", range = "A13:B72", col_names = v1)
SR2d_41 <- read_excel("SRfinal.xlsx", sheet = "2d", range = "A13:B72", col_names = v1)
SR2e_41 <- read_excel("SRfinal.xlsx", sheet = "2e", range = "A13:B72", col_names = v1)
SR2f_41 <- read_excel("SRfinal.xlsx", sheet = "2f", range = "A13:B72", col_names = v1)

SR3a_41 <- read_excel("SRfinal.xlsx", sheet = "3a", range = "A13:B72", col_names = v1)
SR3b_41 <- read_excel("SRfinal.xlsx", sheet = "3b", range = "A13:B72", col_names = v1)
SR3c_41 <- read_excel("SRfinal.xlsx", sheet = "3c", range = "A13:B72", col_names = v1)
SR3d_41 <- read_excel("SRfinal.xlsx", sheet = "3d", range = "A13:B72", col_names = v1)
SR3e_41 <- read_excel("SRfinal.xlsx", sheet = "3e", range = "A13:B72", col_names = v1)
SR3f_41 <- read_excel("SRfinal.xlsx", sheet = "3f", range = "A13:B72", col_names = v1)
SR4a_41 <- read_excel("SRfinal.xlsx", sheet = "4a", range = "A13:B72", col_names = v1)
SR4b_41 <- read_excel("SRfinal.xlsx", sheet = "4b", range = "A13:B72", col_names = v1)
SR4c_41 <- read_excel("SRfinal.xlsx", sheet = "4c", range = "A13:B72", col_names = v1)
SR4d_41 <- read_excel("SRfinal.xlsx", sheet = "4d", range = "A13:B72", col_names = v1)
SR4e_41 <- read_excel("SRfinal.xlsx", sheet = "4e", range = "A13:B72", col_names = v1)
SR4f_41 <- read_excel("SRfinal.xlsx", sheet = "4f", range = "A13:B72", col_names = v1)

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_sr1a_by_age41 <- SR1a_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1b_by_age41 <- SR1b_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1c_by_age41 <- SR1c_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1d_by_age41 <- SR1d_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1e_by_age41 <- SR1e_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr1f_by_age41 <- SR1f_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2a_by_age41 <- SR2a_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2b_by_age41 <- SR2b_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2c_by_age41 <- SR2c_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2d_by_age41 <- SR2d_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2e_by_age41 <- SR2e_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr2f_by_age41 <- SR2f_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3a_by_age41 <- SR3a_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3b_by_age41 <- SR3b_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3c_by_age41 <- SR3c_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3d_by_age41 <- SR3d_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3e_by_age41 <- SR3e_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr3f_by_age41 <- SR3f_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4a_by_age41 <- SR4a_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4b_by_age41 <- SR4b_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4c_by_age41 <- SR4c_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4d_by_age41 <- SR4d_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4e_by_age41 <- SR4e_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

v_r_sr4f_by_age41 <- SR4f_41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(H) %>%
  as.matrix()

##load table 41 
LifeTable_Chile41 <- read.csv2("~/Desktop/ESP/FONIScancer /LifeTable_Chile41.csv")

## Extract age-specific relativa survival 
v_r_mort_by_age41 <- LifeTable_Chile41 %>%
  dplyr::filter(Age >= age & Age < max.age) %>%
  dplyr::select(TotPop) %>%
  as.matrix()

# Age-specific mortality rate in the Healthy state (background mortality)
# for all cycles
v_r_HDage41 <- rep(v_r_mort_by_age41, each = 1/cycle_length)
v_r_1aDage41 <- rep(v_r_sr1a_by_age41, each = 1/cycle_length)
v_r_1bDage41 <- rep(v_r_sr1b_by_age41, each = 1/cycle_length)
v_r_1cDage41 <- rep(v_r_sr1c_by_age41, each = 1/cycle_length)
v_r_1dDage41 <- rep(v_r_sr1d_by_age41, each = 1/cycle_length)
v_r_1eDage41 <- rep(v_r_sr1e_by_age41, each = 1/cycle_length)
v_r_1fDage41 <- rep(v_r_sr1f_by_age41, each = 1/cycle_length)
v_r_2aDage41 <- rep(v_r_sr2a_by_age41, each = 1/cycle_length)
v_r_2bDage41 <- rep(v_r_sr2b_by_age41, each = 1/cycle_length)
v_r_2cDage41 <- rep(v_r_sr2c_by_age41, each = 1/cycle_length)
v_r_2dDage41 <- rep(v_r_sr2d_by_age41, each = 1/cycle_length)
v_r_2eDage41 <- rep(v_r_sr2e_by_age41, each = 1/cycle_length)
v_r_2fDage41 <- rep(v_r_sr2f_by_age41, each = 1/cycle_length)
v_r_3aDage41 <- rep(v_r_sr3a_by_age41, each = 1/cycle_length)
v_r_3bDage41 <- rep(v_r_sr3b_by_age41, each = 1/cycle_length)
v_r_3cDage41 <- rep(v_r_sr3c_by_age41, each = 1/cycle_length)
v_r_3dDage41 <- rep(v_r_sr3d_by_age41, each = 1/cycle_length)
v_r_3eDage41 <- rep(v_r_sr3e_by_age41, each = 1/cycle_length)
v_r_3fDage41 <- rep(v_r_sr3f_by_age41, each = 1/cycle_length)
v_r_4aDage41 <- rep(v_r_sr4a_by_age41, each = 1/cycle_length)
v_r_4bDage41 <- rep(v_r_sr4b_by_age41, each = 1/cycle_length)
v_r_4cDage41 <- rep(v_r_sr4c_by_age41, each = 1/cycle_length)
v_r_4dDage41 <- rep(v_r_sr4d_by_age41, each = 1/cycle_length)
v_r_4eDage41 <- rep(v_r_sr4e_by_age41, each = 1/cycle_length)
v_r_4fDage41 <- rep(v_r_sr4f_by_age41, each = 1/cycle_length)

# Transform to age-specific background mortality risk for all cycles adjusting
# by cycle length
v_p_HDage41 <- 1 - exp(-v_r_HDage41 * cycle_length)
v_p_1aDage41 <- v_p_HDage41 + (1 - exp(-v_r_1aDage41 * cycle_length))
v_p_1bDage41 <- v_p_HDage41 + (1 - exp(-v_r_1bDage41 * cycle_length))
v_p_1cDage41 <- v_p_HDage41 + (1 - exp(-v_r_1cDage41 * cycle_length))
v_p_1dDage41 <- v_p_HDage41 + (1 - exp(-v_r_1dDage41 * cycle_length))
v_p_1eDage41 <- v_p_HDage41 + (1 - exp(-v_r_1eDage41 * cycle_length))
v_p_1fDage41 <- v_p_HDage41 + (1 - exp(-v_r_1fDage41 * cycle_length))
v_p_2aDage41 <- v_p_HDage41 + (1 - exp(-v_r_2aDage41 * cycle_length))
v_p_2bDage41 <- v_p_HDage41 + (1 - exp(-v_r_2bDage41 * cycle_length))
v_p_2cDage41 <- v_p_HDage41 + (1 - exp(-v_r_2cDage41 * cycle_length))
v_p_2dDage41 <- v_p_HDage41 + (1 - exp(-v_r_2dDage41 * cycle_length))
v_p_2eDage41 <- v_p_HDage41 + (1 - exp(-v_r_2eDage41 * cycle_length))
v_p_2fDage41 <- v_p_HDage41 + (1 - exp(-v_r_2fDage41 * cycle_length))
v_p_3aDage41 <- v_p_HDage41 + (1 - exp(-v_r_3aDage41 * cycle_length))
v_p_3bDage41 <- v_p_HDage41 + (1 - exp(-v_r_3bDage41 * cycle_length))
v_p_3cDage41 <- v_p_HDage41 + (1 - exp(-v_r_3cDage41 * cycle_length))
v_p_3dDage41 <- v_p_HDage41 + (1 - exp(-v_r_3dDage41 * cycle_length))
v_p_3eDage41 <- v_p_HDage41 + (1 - exp(-v_r_3eDage41 * cycle_length))
v_p_3fDage41 <- v_p_HDage41 + (1 - exp(-v_r_3fDage41 * cycle_length))
v_p_4aDage41 <- v_p_HDage41 + (1 - exp(-v_r_4aDage41 * cycle_length))
v_p_4bDage41 <- v_p_HDage41 + (1 - exp(-v_r_4bDage41 * cycle_length))
v_p_4cDage41 <- v_p_HDage41 + (1 - exp(-v_r_4cDage41 * cycle_length))
v_p_4dDage41 <- v_p_HDage41 + (1 - exp(-v_r_4dDage41 * cycle_length))
v_p_4eDage41 <- v_p_HDage41 + (1 - exp(-v_r_4eDage41 * cycle_length))
v_p_4fDage41 <- v_p_HDage41 + (1 - exp(-v_r_4fDage41 * cycle_length))


  ##    02.3 Probabilidades de transición --------------------------------------------

## Probabilities of Transition 
#Probability of transition without prevention --> Valores Huang

#p.NnGn <- 0.0360531  # probability to become Gastritis hp(-) when Normal hp(-) --> Huang, 2020
#p.GnAn <- 0.07787315  # probability to become Atrophy hp(-) when Gastritis hp(-) --> Promedio extraído de Huang, 2020 / Yeh 0.001 
#p.AnIn <- 0.12133965  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) --> Huang, 2020
#p.InDn <- 0.06380055  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
#p.DnIn <- 0.1544834  # probability to become Intestinal hp(-) when Dysplasia hp(-)--> Promedio Huang, 2020 / 0.0081 Yeh 2016 monthly
#p.InAn <- 0.04664585  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-) --> Promedio Huang, 2020 / 0,0029 Yeh 2016 monthly
#p.AnGn <- 0.0930757  # probability to become Gastritis hp(-) when Atrophy hp(-) --> promedio Huang, 2020 / 0,0086 Yeh 2016 monthly
#p.GnNn <- 0.0043724  # probability to become Normal hp(-) when Gastritis hp(-) --> Huang, 2020
#p.GnGp <- 0.0 # probability to become Gastritis hp(+) when Gastritis hp(-)
#p.GpGn <- 0.01 # probability to become Gastritis hp(-) when Gastritis hp(+)
#p.AnAp <- 0.0 # probability to become Atrophy hp(+) when Atrophy hp(-)
#p.ApAn <- 0.01 # probability to become Atrophy hp(-) when Atrophy hp(+)
#p.InIp <- 0.0 # probability to become Intestinal hp(+) when Intestinal hp(-)
#p.IpIn <- 0.01 # probability to become Intestinal hp(-) when Intestinal hp(+)
#p.DnDp <- 0.0 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
#p.DpDn <- 0.01 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
#p.NpGp <- 0.0  # probability to become Gastritis hp(+) when Normal hp(+)
#p.NnNp <- 0.0  # probability to become Normal hp(+) when Normal hp(-)
#p.NpNn <- 0.01  # probability to become Normal hp(-) when Normal hp(+)
#p.NpGp <- 0.0  # probability to become Gastritis hp(+) when Normanl hp(+)
#p.GpAp <- 5 * p.GnAn  # probability to become Atrophy hp(+) when Gastritis hp(+) --> Abdullahi 2011, metanalisis utilizado por Huang 
#p.ApIp <- p.AnIn  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) --> FALTA 
#p.IpDp <- p.InDn  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+) --> FALTA 
#p.DpIp <- 0  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+) --> Se decidió 
#p.IpAp <- 0  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)  --> Se decidió 
#p.ApGp <- 0.0  # probability to become Gastritis hp(+) when Atrophy hp(+) --> CONSULTAR, acá asumo que los positivos no regresan 
#p.GpNp <- 0.0  # probability to become Normal hp (+) when Gastritis hp(+) --> CONSULTAR, acá asumo que los positivos no regresan
#p.Dn1p <- 0.034925  # probability to become stage 1 preclinical when Dysplasia (-) --> Huang 2020, promedio 
#p.Dp1p <- p.Dn1p  # probability to become stage 1 preclinical when Dysplasia (+) --> Es distinta esta probabilidad con HP +? 
#p.1p1ca <- 0.12312  # probability to become stage 1 clinical a when stage 1 preclinical --> Yeh 2016, probabilidad mensual * 12
#p.1ca1cb <- 1 # probability to become stage 1 clinical b when stage 1 clinical a 
#p.1cb1cc <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
#p.1cc1cd <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
#p.1cd1ce <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
#p.1ce1cf <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
#p.1cfNn <- 0.5
#p.1p2p <- 0.19176665  # probability to become stage 2 preclinical when stage 1 preclinical --> Huang 2020, promedio
#p.2p2ca <- 0.48  # probability to become stage 2 clinical a when stage 2 preclinical --> Yeh 2016, probabilidad mensual * 12
#p.2ca2cb <- 1 # probability to become stage 2 clinical b when stage 2 clinical a 
#p.2cb2cc <- 0.2 # probability to become stage 2 clinical c when stage 2 clinical b 
#p.2cc2cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.2cd2ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.2ce2cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.2cfNn <- 0.4 
#p.2p3p <- 0.19176665  # probability to become stage 3 preclinical when stage 2 preclinical --> Huang 2020, promedio
#p.3p3ca <- 0.6  # probability to become stage 3 clinical a when stage 3 preclinical --> INVENTO
#p.3ca3cb <- 1 # probability to become stage 3 clinical b when stage 3 clinical a 
#p.3cb3cc <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.3cc3cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.3cd3ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.3ce3cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.3cfNn <- 0.3
#p.3p4p <- 0.3934693  # probability to become stage 4 preclinical when stage 3 preclinical --> Huang 2020
#p.4p4ca <- 0.7  # probability to become stage 4 clinical a when stage 4 preclinical --> INVENTO
#p.4ca4cb <- 1 # probability to become stage 4 clinical b when stage 4 clinical a 
#p.4cb4cc <- 1 # probability to become stage 4 clinical c when stage 4 clinical b 
#p.4cc4cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.4cd4ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.4ce4cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
#p.4cfNn <- 0.2

### Probabilidades de transición - Valores Arnoldo Riquelme 

p.NnGn <- 0.002008048339 # probability to become Gastritis hp(-) when Normal hp(-) 
p.GnAn <- 0.001002006021  # probability to become Atrophy hp(-) when Gastritis hp(-) 
p.AnIn <- 0.002008048339  # probability to become Intestinal Metaplasia hp (-) when Atrophy hp(-) 
p.InDn <- 0.006073323852  # probability to become Dysplasia hp(-) when Intestinal Metaplasia hp(-)
p.DnIn <- 0.4124841223  # probability to become Intestinal hp(-) when Dysplasia hp(-)
p.InAn <- 0.06885008491  # probability to become Atrophy hp(-) when Intestinal Metaplasia hp (-) 
p.AnGn <- 0.1294494367  # probability to become Gastritis hp(-) when Atrophy hp(-) 
p.GnNn <- 0.04364750021  # probability to become Normal hp(-) when Gastritis hp(-) 
p.GnGp <- 0.0 # probability to become Gastritis hp(+) when Gastritis hp(-)
p.GpGn <- 0.0 # probability to become Gastritis hp(-) when Gastritis hp(+)
p.AnAp <- 0.0 # probability to become Atrophy hp(+) when Atrophy hp(-)
p.ApAn <- 0.01 # probability to become Atrophy hp(-) when Atrophy hp(+)
p.InIp <- 0.0 # probability to become Intestinal hp(+) when Intestinal hp(-)
p.IpIn <- 0.01 # probability to become Intestinal hp(-) when Intestinal hp(+)
p.DnDp <- 0.0 # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.DpDn <- 0.01 # probability to become Dysplasia hp(-) when Dysplasia hp(+)
p.NpGp <- 0.4507197283  # probability to become Gastritis hp(+) when Normal hp(+)
p.NnNp <- 0.0  # probability to become Normal hp(+) when Normal hp(-)
p.NpNn <- 0.0  # probability to become Normal hp(-) when Normal hp(+)
p.NpGp <- 0.0  # probability to become Gastritis hp(+) when Normal hp(+)
p.GpAp <- 0.01020621831  # probability to become Atrophy hp(+) when Gastritis hp(+)  
p.ApIp <- 0.01440931621  # probability to become Intestinal Metaplasia hp (+) when Atrophy hp(+) 
p.IpDp <- 0.01440931621  # probability to become Dysplasia hp(+) when Intestinal Metaplasia hp(+)  
p.DpIp <- 0.3157445711  # probability to become Intestinal Metaplasia hp(+) when Dysplasia hp(+) 
p.IpAp <- 0.02085163764  # probability to become Atrophy hp(+) when Intestinal Metaplasia hp (+)   
p.ApGp <- 0.02085163764  # probability to become Gastritis hp(+) when Atrophy hp(+)  
p.GpNp <- 0.02085163764  # probability to become Normal hp (+) when Gastritis hp(+) 
p.Dn1p <- 0.01440931621  # probability to become stage 1 preclinical when Dysplasia (-) --> Huang 2020, promedio 
p.Dp1p <- 0.031981215  # probability to become stage 1 preclinical when Dysplasia (+) --> Es distinta esta probabilidad con HP +? 
p.1p1ca <- 0.1164045522  # probability to become stage 1 clinical a when stage 1 preclinical --> Yeh 2016, probabilidad mensual * 12
p.1ca1cb <- 1 # probability to become stage 1 clinical b when stage 1 clinical a 
p.1cb1cc <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
p.1cc1cd <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
p.1cd1ce <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
p.1ce1cf <- 1 # probability to become stage 1 clinical c when stage 1 clinical b 
p.1cfNn <- 1
p.1p2p <- 0.1229316  # probability to become stage 2 preclinical when stage 1 preclinical --> Huang 2020, límite inferior
p.2p2ca <- 0.3872902427  # probability to become stage 2 clinical a when stage 2 preclinical --> Yeh 2016, probabilidad mensual * 12
p.2ca2cb <- 1 # probability to become stage 2 clinical b when stage 2 clinical a 
p.2cb2cc <- 0.2 # probability to become stage 2 clinical c when stage 2 clinical b 
p.2cc2cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.2cd2ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.2ce2cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.2cfNn <- 1 
p.2p3p <- 0.1229316  # probability to become stage 3 preclinical when stage 2 preclinical --> Huang 2020, límite inferior
p.3p3ca <- 0.3872902427  # probability to become stage 3 clinical a when stage 3 preclinical --> Yeh, 2016 
p.3ca3cb <- 1 # probability to become stage 3 clinical b when stage 3 clinical a 
p.3cb3cc <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.3cc3cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.3cd3ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.3ce3cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.3cfNn <- 1
p.3p4p <- 0.3934693  # probability to become stage 4 preclinical when stage 3 preclinical --> Huang 2020
p.4p4ca <- 0.9688267694  # probability to become stage 4 clinical a when stage 4 preclinical --> Yeh, 2016
p.4ca4cb <- 1 # probability to become stage 4 clinical b when stage 4 clinical a 
p.4cb4cc <- 1 # probability to become stage 4 clinical c when stage 4 clinical b 
p.4cc4cd <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.4cd4ce <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.4ce4cf <- 1 # probability to become stage 3 clinical c when stage 3 clinical b 
p.4cfNn <- 1


#Probability of transition with prevention strategy 1 (UreaAire)

p.pGpGn <- 0.938   # probability to become Gastritis hp(-) when Gastritis hp(+)
p.pApAn <- 0.938   # probability to become Atrophy hp(+) when Atrophy hp(-)
p.pIpIn <- 0.938   # probability to become Intestinal hp(+) when Intestinal hp(-)
p.pDpDn <- 0.938   # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.pNpNn <- 0.938   # probability to become Normal hp(-) when Normal hp(+)

#Probability of transition with prevention strategy 2 (Antígeno Fecal)

p.p2GpGn <- 0.9204   # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p2ApAn <- 0.9204   # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p2IpIn <- 0.9204   # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p2DpDn <- 0.9204   # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p2NpNn <- 0.9204      # probability to become Normal hp(-) when Normal hp(+) AUMENTA CON LA INTERVENCIÓN

#Probability of transition with prevention strategy 3 (EDA)

p.p3GpGn <- 0.938   # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p3ApAn <- 0.938   # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p3IpIn <- 0.938   # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p3DpDn <- 0.938   # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p3NpNn <- 0.938      # probability to become Normal hp(-) when Normal hp(+) 
p.p3_1p1ca <- 0.99  # probability to become stage 1 clinical a when stage 1 preclinical 
p.p3_2p2ca <- 0.99 # probability to become stage 2 clinical a when stage 2 preclinical 
p.p3_3p3ca <- 0.99  # probability to become stage 3 clinical a when stage 3 preclinical  
p.p3_4p4ca <- 0.99 # probability to become stage 4 clinical a when stage 4 preclinical  

#Probability of transition with prevention strategy 4 (PSEDA)
#PENDIENTE 
p.p4GpGn <- 0.01*0.95   # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p4ApAn <- 0.06*0.95   # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p4IpIn <- 0.06*0.95   # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p4DpDn <- 0.07*0.95   # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p4NpNn <- 0.01*0.95      # probability to become Normal hp(-) when Normal hp(+) 
p.p4_1p1ca <- 0.7  # probability to become stage 1 clinical a when stage 1 preclinical
p.p4_2p2ca <- 0.7 # probability to become stage 2 clinical a when stage 2 preclinical
p.p4_3p3ca <- 0.7  # probability to become stage 3 clinical a when stage 3 preclinical
p.p4_4p4ca <- 0.7 # probability to become stage 4 clinical a when stage 4 preclinical

#Probability of transition with prevention strategy 5 (Serología)

p.p5GpGn <- 0.9468   # probability to become Gastritis hp(+) when Gastritis hp(-)
p.p5ApAn <- 0.9468   # probability to become Atrophy hp(+) when Atrophy hp(-)
p.p5IpIn <- 0.9468   # probability to become Intestinal hp(+) when Intestinal hp(-)
p.p5DpDn <- 0.9468   # probability to become Dysplasia hp(+) when Dysplasia hp(-)
p.p5NpNn <- 0.9468      # probability to become Normal hp(-) when Normal hp(+) AUMENTA CON LA INTERVENCIÓN


  ##    02.4 Probabilities of dying --------------------------------------------------


p.NnD <- v_p_HDage  # probability to die when healthy 
p.NpD <- v_p_HDage  # probability to die when Normal helicobacter (+)
p.GnD <- v_p_HDage # probability to die when Gastritis helicobacter (-)
p.GpD <- v_p_HDage # probability to die when Gastritis helicobacter (+)
p.AnD <- v_p_HDage # probability to die when Atrophy helicobacter (-)
p.ApD <- v_p_HDage # probability to die when Atrophy helicobacter (+)
p.InD <- v_p_HDage  # probability to die when Intestinal helicobacter (-) 
p.IpD <- v_p_HDage  # probability to die when Intestinal helicobacter (+) 
p.DnD <- v_p_HDage  # probability to die when Dysplasia helicobacter (-) 
p.DpD <- v_p_HDage  # probability to die when Dysplasia helicobacter (+) 
p.1pD <- v_p_1aDage  # probability to die when stage 1 preclinical 
p.2pD <- v_p_2aDage # probability to die when stage 2 preclinical
p.3pD <- v_p_3aDage # probability to die when stage 3 preclinical
p.4pD <- v_p_4aDage # probability to die when stage 4 preclinical
p.1caD <- v_p_1aDage  # probability to die when stage 1 clinical a
p.2caD <- v_p_2aDage # probability to die when stage 2 clinical a
p.3caD <- v_p_3aDage # probability to die when stage 3 clinical a
p.4caD <- v_p_4aDage # probability to die when stage 4 clinical a
p.1cbD <- v_p_1bDage # probability to die when stage 1 clinical b 
p.2cbD <- v_p_2bDage # probability to die when stage 2 clinical b
p.3cbD <- v_p_3bDage # probability to die when stage 3 clinical b
p.4cbD <- v_p_4bDage # probability to die when stage 4 clinical b
p.1ccD <- v_p_1cDage  # probability to die when stage 1 clinical c
p.2ccD <- v_p_2cDage # probability to die when stage 2 clinical c
p.3ccD <- v_p_3cDage # probability to die when stage 3 clinical c
p.4ccD <- v_p_4cDage # probability to die when stage 4 clinical c
p.1cdD <- v_p_1dDage  # probability to die when stage 1 clinical c
p.2cdD <- v_p_2dDage # probability to die when stage 2 clinical c
p.3cdD <- v_p_3dDage # probability to die when stage 3 clinical c
p.4cdD <- v_p_4dDage # probability to die when stage 4 clinical c
p.1ceD <- v_p_1eDage  # probability to die when stage 1 clinical c
p.2ceD <- v_p_2eDage # probability to die when stage 2 clinical c
p.3ceD <- v_p_3eDage # probability to die when stage 3 clinical c
p.4ceD <- v_p_4eDage # probability to die when stage 4 clinical c
p.1cfD <- v_p_1fDage  # probability to die when stage 1 clinical c
p.2cfD <- v_p_2fDage # probability to die when stage 2 clinical c
p.3cfD <- v_p_3fDage # probability to die when stage 3 clinical c
p.4cfD <- v_p_4fDage # probability to die when stage 4 clinical c

## Probability to die at 30 

p.NnD30 <- 0.0007297762825  # probability to die when healthy at 30 
p.NpD30 <- 0.0007297762825  # probability to die when Normal helicobacter (+)
p.GnD30 <- 0.0007297762825 # probability to die when Gastritis helicobacter (-)
p.GpD30 <- 0.0007297762825 # probability to die when Gastritis helicobacter (+)
p.AnD30 <- 0.0007297762825 # probability to die when Atrophy helicobacter (-)
p.ApD30 <- 0.0007297762825 # probability to die when Atrophy helicobacter (+)
p.InD30 <- 0.0007297762825  # probability to die when Intestinal helicobacter (-) 
p.IpD30 <- 0.0007297762825  # probability to die when Intestinal helicobacter (+) 
p.DnD30 <- 0.0007297762825  # probability to die when Dysplasia helicobacter (-) 
p.DpD30 <- 0.0007297762825  # probability to die when Dysplasia helicobacter (+) 
p.1pD30 <- 0.012214253 + 0.0007297762825  # probability to die when stage 1 preclinical 
p.2pD30 <- 0.08030691 + 0.0007297762825 # probability to die when stage 2 preclinical
p.3pD30 <- 0.174171454 + 0.0007297762825  # probability to die when stage 3 preclinical
p.4pD30 <- 0.228382927 + 0.0007297762825 # probability to die when stage 4 preclinical
p.1caD30 <- 0.012214253 + 0.0007297762825  # probability to die when stage 1 clinical a
p.2caD30 <- 0.08030691 + 0.0007297762825 # probability to die when stage 2 clinical a
p.3caD30 <- 0.174171454 + 0.0007297762825 # probability to die when stage 3 clinical a
p.4caD30 <- 0.228382927 + 0.0007297762825 # probability to die when stage 4 clinical a
p.1cbD30 <- 0.004444335 + 0.0007297762825 # probability to die when stage 1 clinical b 
p.2cbD30 <- 0.114891303 + 0.0007297762825 # probability to die when stage 2 clinical b
p.3cbD30 <- 0.258072382 + 0.0007297762825 # probability to die when stage 3 clinical b
p.4cbD30 <- 0.652838594 + 0.0007297762825 # probability to die when stage 4 clinical b
p.1ccD30 <- 0.003510454 + 0.0007297762825  # probability to die when stage 1 clinical c
p.2ccD30 <- 0.091858077 + 0.0007297762825 # probability to die when stage 2 clinical c
p.3ccD30 <- 0.209957968 + 0.0007297762825 # probability to die when stage 3 clinical c
p.4ccD30 <- 0.566240775 + 0.0007297762825 # probability to die when stage 4 clinical c
p.1cdD30 <- 0.002796383 + 0.0007297762825
p.1ceD30 <- 0.002018909 + 0.0007297762825
p.1cfD30 <- 0.006118674 + 0.0007297762825
p.2cdD30 <- 0.073857548 + 0.0007297762825
p.2ceD30 <- 0.053868251 + 0.0007297762825
p.2cfD30 <- 0.040946243 + 0.0007297762825
p.3cdD30 <- 0.171106376 + 0.0007297762825
p.3ceD30 <- 0.126664317 + 0.0007297762825
p.3cfD30 <- 0.091146106 + 0.0007297762825
p.4cdD30 <- 0.566240775 + 0.0007297762825
p.4ceD30 <- 0.566240775 + 0.0007297762825
p.4cfD30 <- 0.566240775 + 0.0007297762825

## Probability to die at 45 

p.NnD45 <- 0.001915352  # probability to die when healthy at 30 
p.NpD45 <- 0.001915352  # probability to die when Normal helicobacter (+)
p.GnD45 <- 0.001915352 # probability to die when Gastritis helicobacter (-)
p.GpD45 <- 0.001915352 # probability to die when Gastritis helicobacter (+)
p.AnD45 <- 0.001915352 # probability to die when Atrophy helicobacter (-)
p.ApD45 <- 0.001915352 # probability to die when Atrophy helicobacter (+)
p.InD45 <- 0.001915352  # probability to die when Intestinal helicobacter (-) 
p.IpD45 <- 0.001915352  # probability to die when Intestinal helicobacter (+) 
p.DnD45 <- 0.001915352  # probability to die when Dysplasia helicobacter (-) 
p.DpD45 <- 0.001915352  # probability to die when Dysplasia helicobacter (+) 
p.1pD45 <- 0.002686375 + 0.001915352   # probability to die when stage 1 preclinical 
p.2pD45 <- 0.071054215 + 0.001915352  # probability to die when stage 2 preclinical
p.3pD45 <- 0.164956367 + 0.001915352   # probability to die when stage 3 preclinical
p.4pD45 <- 0.472139366 + 0.001915352  # probability to die when stage 4 preclinical
p.1caD45 <- 0.002686375 + 0.001915352   # probability to die when stage 1 clinical a
p.2caD45 <- 0.071054215 + 0.001915352  # probability to die when stage 2 clinical a
p.3caD45 <- 0.164956367 + 0.001915352  # probability to die when stage 3 clinical a
p.4caD45 <- 0.472139366 + 0.001915352  # probability to die when stage 4 clinical a
p.1cbD45 <- 0.004444335 + 0.001915352  # probability to die when stage 1 clinical b 
p.2cbD45 <- 0.114891303 + 0.001915352  # probability to die when stage 2 clinical b
p.3cbD45 <- 0.258072382 + 0.001915352  # probability to die when stage 3 clinical b
p.4cbD45 <- 0.652838594 + 0.001915352  # probability to die when stage 4 clinical b
p.1ccD45 <- 0.003510454 + 0.001915352   # probability to die when stage 1 clinical c
p.2ccD45 <- 0.091858077 + 0.001915352  # probability to die when stage 2 clinical c
p.3ccD45 <- 0.209957968 + 0.001915352  # probability to die when stage 3 clinical c
p.4ccD45 <- 0.566240775 + 0.001915352  # probability to die when stage 4 clinical c
p.1cdD45 <- 0.002796383 + 0.001915352 
p.1ceD45 <- 0.002018909 + 0.001915352 
p.1cfD45 <- 0.006118674 + 0.001915352 
p.2cdD45 <- 0.073857548 + 0.001915352 
p.2ceD45 <- 0.053868251 + 0.001915352 
p.2cfD45 <- 0.040946243 + 0.001915352 
p.3cdD45 <- 0.171106376 + 0.001915352 
p.3ceD45 <- 0.126664317 + 0.001915352 
p.3cfD45 <- 0.091146106 + 0.001915352 
p.4cdD45 <- 0.566240775 + 0.001915352 
p.4ceD45 <- 0.566240775 + 0.001915352 
p.4cfD45 <- 0.566240775 + 0.001915352 

## Probability to die at 40 

p.NnD40 <- 0.001259204  # probability to die when healthy at 30 
p.NpD40 <- 0.001259204  # probability to die when Normal helicobacter (+)
p.GnD40 <- 0.001259204 # probability to die when Gastritis helicobacter (-)
p.GpD40 <- 0.001259204 # probability to die when Gastritis helicobacter (+)
p.AnD40 <- 0.001259204 # probability to die when Atrophy helicobacter (-)
p.ApD40 <- 0.001259204 # probability to die when Atrophy helicobacter (+)
p.InD40 <- 0.001259204  # probability to die when Intestinal helicobacter (-) 
p.IpD40 <- 0.001259204  # probability to die when Intestinal helicobacter (+) 
p.DnD40 <- 0.001259204  # probability to die when Dysplasia helicobacter (-) 
p.DpD40 <- 0.001259204  # probability to die when Dysplasia helicobacter (+) 
p.1pD40 <- 0.002686375 + 0.001259204   # probability to die when stage 1 preclinical 
p.2pD40 <- 0.071054215 + 0.001259204  # probability to die when stage 2 preclinical
p.3pD40 <- 0.164956367 + 0.001259204   # probability to die when stage 3 preclinical
p.4pD40 <- 0.472139366 + 0.001259204  # probability to die when stage 4 preclinical
p.1caD40 <- 0.002686375 + 0.001259204   # probability to die when stage 1 clinical a
p.2caD40 <- 0.071054215 + 0.001259204  # probability to die when stage 2 clinical a
p.3caD40 <- 0.164956367 + 0.001259204  # probability to die when stage 3 clinical a
p.4caD40 <- 0.472139366 + 0.001259204  # probability to die when stage 4 clinical a
p.1cbD40 <- 0.004444335 + 0.001259204  # probability to die when stage 1 clinical b 
p.2cbD40 <- 0.114891303 + 0.001259204  # probability to die when stage 2 clinical b
p.3cbD40 <- 0.258072382 + 0.001259204  # probability to die when stage 3 clinical b
p.4cbD40 <- 0.652838594 + 0.001259204  # probability to die when stage 4 clinical b
p.1ccD40 <- 0.003510454 + 0.001259204   # probability to die when stage 1 clinical c
p.2ccD40 <- 0.091858077 + 0.001259204  # probability to die when stage 2 clinical c
p.3ccD40 <- 0.209957968 + 0.001259204  # probability to die when stage 3 clinical c
p.4ccD40 <- 0.566240775 + 0.001259204  # probability to die when stage 4 clinical c
p.1cdD40 <- 0.002796383 + 0.001259204 
p.1ceD40 <- 0.002018909 + 0.001259204 
p.1cfD40 <- 0.006118674 + 0.001259204 
p.2cdD40 <- 0.073857548 + 0.001259204 
p.2ceD40 <- 0.053868251 + 0.001259204 
p.2cfD40 <- 0.040946243 + 0.001259204 
p.3cdD40 <- 0.171106376 + 0.001259204 
p.3ceD40 <- 0.126664317 + 0.001259204 
p.3cfD40 <- 0.091146106 + 0.001259204 
p.4cdD40 <- 0.566240775 + 0.001259204 
p.4ceD40 <- 0.566240775 + 0.001259204 
p.4cfD40 <- 0.566240775 + 0.001259204 


# Probability to die from 31 to 100
p.NnD31 <- v_p_HDage31  # probability to die when healthy 
p.NpD31 <- v_p_HDage31  # probability to die when Normal helicobacter (+)
p.GnD31 <- v_p_HDage31 # probability to die when Gastritis helicobacter (-)
p.GpD31 <- v_p_HDage31 # probability to die when Gastritis helicobacter (+)
p.AnD31 <- v_p_HDage31 # probability to die when Atrophy helicobacter (-)
p.ApD31 <- v_p_HDage31 # probability to die when Atrophy helicobacter (+)
p.InD31 <- v_p_HDage31  # probability to die when Intestinal helicobacter (-) 
p.IpD31 <- v_p_HDage31  # probability to die when Intestinal helicobacter (+) 
p.DnD31 <- v_p_HDage31  # probability to die when Dysplasia helicobacter (-) 
p.DpD31 <- v_p_HDage31  # probability to die when Dysplasia helicobacter (+) 
p.1pD31 <- v_p_1aDage31  # probability to die when stage 1 preclinical 
p.2pD31 <- v_p_2aDage31 # probability to die when stage 2 preclinical
p.3pD31 <- v_p_3aDage31 # probability to die when stage 3 preclinical
p.4pD31 <- v_p_4aDage31 # probability to die when stage 4 preclinical
p.1caD31 <- v_p_1aDage31  # probability to die when stage 1 clinical a
p.2caD31 <- v_p_2aDage31 # probability to die when stage 2 clinical a
p.3caD31 <- v_p_3aDage31 # probability to die when stage 3 clinical a
p.4caD31 <- v_p_4aDage31 # probability to die when stage 4 clinical a
p.1cbD31 <- v_p_1bDage31 # probability to die when stage 1 clinical b 
p.2cbD31 <- v_p_2bDage31 # probability to die when stage 2 clinical b
p.3cbD31 <- v_p_3bDage31 # probability to die when stage 3 clinical b
p.4cbD31 <- v_p_4bDage31 # probability to die when stage 4 clinical b
p.1ccD31 <- v_p_1cDage31  # probability to die when stage 1 clinical c
p.2ccD31 <- v_p_2cDage31 # probability to die when stage 2 clinical c
p.3ccD31 <- v_p_3cDage31 # probability to die when stage 3 clinical c
p.4ccD31 <- v_p_4cDage31 # probability to die when stage 4 clinical c
p.1cdD31 <- v_p_1dDage31  # probability to die when stage 1 clinical c
p.2cdD31 <- v_p_2dDage31 # probability to die when stage 2 clinical c
p.3cdD31 <- v_p_3dDage31 # probability to die when stage 3 clinical c
p.4cdD31 <- v_p_4dDage31 # probability to die when stage 4 clinical c
p.1ceD31 <- v_p_1eDage31  # probability to die when stage 1 clinical c
p.2ceD31 <- v_p_2eDage31 # probability to die when stage 2 clinical c
p.3ceD31 <- v_p_3eDage31 # probability to die when stage 3 clinical c
p.4ceD31 <- v_p_4eDage31 # probability to die when stage 4 clinical c
p.1cfD31 <- v_p_1fDage31  # probability to die when stage 1 clinical c
p.2cfD31 <- v_p_2fDage31 # probability to die when stage 2 clinical c
p.3cfD31 <- v_p_3fDage31 # probability to die when stage 3 clinical c
p.4cfD31 <- v_p_4fDage31 # probability to die when stage 4 clinical c

#Probability to die from 30 to 44 

p.NnD3044 <- v_p_HDage3044  # probability to die when healthy 
p.NpD3044 <- v_p_HDage3044  # probability to die when Normal helicobacter (+)
p.GnD3044 <- v_p_HDage3044 # probability to die when Gastritis helicobacter (-)
p.GpD3044 <- v_p_HDage3044 # probability to die when Gastritis helicobacter (+)
p.AnD3044 <- v_p_HDage3044 # probability to die when Atrophy helicobacter (-)
p.ApD3044 <- v_p_HDage3044 # probability to die when Atrophy helicobacter (+)
p.InD3044 <- v_p_HDage3044  # probability to die when Intestinal helicobacter (-) 
p.IpD3044 <- v_p_HDage3044  # probability to die when Intestinal helicobacter (+) 
p.DnD3044 <- v_p_HDage3044  # probability to die when Dysplasia helicobacter (-) 
p.DpD3044 <- v_p_HDage3044  # probability to die when Dysplasia helicobacter (+) 
p.1pD3044 <- v_p_1aDage3044  # probability to die when stage 1 preclinical 
p.2pD3044 <- v_p_2aDage3044 # probability to die when stage 2 preclinical
p.3pD3044 <- v_p_3aDage3044 # probability to die when stage 3 preclinical
p.4pD3044 <- v_p_4aDage3044 # probability to die when stage 4 preclinical
p.1caD3044 <- v_p_1aDage3044  # probability to die when stage 1 clinical a
p.2caD3044 <- v_p_2aDage3044 # probability to die when stage 2 clinical a
p.3caD3044 <- v_p_3aDage3044 # probability to die when stage 3 clinical a
p.4caD3044 <- v_p_4aDage3044 # probability to die when stage 4 clinical a
p.1cbD3044 <- v_p_1bDage3044 # probability to die when stage 1 clinical b 
p.2cbD3044 <- v_p_2bDage3044 # probability to die when stage 2 clinical b
p.3cbD3044 <- v_p_3bDage3044 # probability to die when stage 3 clinical b
p.4cbD3044 <- v_p_4bDage3044 # probability to die when stage 4 clinical b
p.1ccD3044 <- v_p_1cDage3044  # probability to die when stage 1 clinical c
p.2ccD3044 <- v_p_2cDage3044 # probability to die when stage 2 clinical c
p.3ccD3044 <- v_p_3cDage3044 # probability to die when stage 3 clinical c
p.4ccD3044 <- v_p_4cDage3044 # probability to die when stage 4 clinical c
p.1cdD3044 <- v_p_1dDage3044  # probability to die when stage 1 clinical c
p.2cdD3044 <- v_p_2dDage3044 # probability to die when stage 2 clinical c
p.3cdD3044 <- v_p_3dDage3044 # probability to die when stage 3 clinical c
p.4cdD3044 <- v_p_4dDage3044 # probability to die when stage 4 clinical c
p.1ceD3044 <- v_p_1eDage3044  # probability to die when stage 1 clinical c
p.2ceD3044 <- v_p_2eDage3044 # probability to die when stage 2 clinical c
p.3ceD3044 <- v_p_3eDage3044 # probability to die when stage 3 clinical c
p.4ceD3044 <- v_p_4eDage3044 # probability to die when stage 4 clinical c
p.1cfD3044 <- v_p_1fDage3044  # probability to die when stage 1 clinical c
p.2cfD3044 <- v_p_2fDage3044 # probability to die when stage 2 clinical c
p.3cfD3044 <- v_p_3fDage3044 # probability to die when stage 3 clinical c
p.4cfD3044 <- v_p_4fDage3044 # probability to die when stage 4 clinical c

# Probability to die from 46 to 100
p.NnD46 <- v_p_HDage46  # probability to die when healthy 
p.NpD46 <- v_p_HDage46  # probability to die when Normal helicobacter (+)
p.GnD46 <- v_p_HDage46 # probability to die when Gastritis helicobacter (-)
p.GpD46 <- v_p_HDage46 # probability to die when Gastritis helicobacter (+)
p.AnD46 <- v_p_HDage46 # probability to die when Atrophy helicobacter (-)
p.ApD46 <- v_p_HDage46 # probability to die when Atrophy helicobacter (+)
p.InD46 <- v_p_HDage46  # probability to die when Intestinal helicobacter (-) 
p.IpD46 <- v_p_HDage46  # probability to die when Intestinal helicobacter (+) 
p.DnD46 <- v_p_HDage46  # probability to die when Dysplasia helicobacter (-) 
p.DpD46 <- v_p_HDage46  # probability to die when Dysplasia helicobacter (+) 
p.1pD46 <- v_p_1aDage46  # probability to die when stage 1 preclinical 
p.2pD46 <- v_p_2aDage46 # probability to die when stage 2 preclinical
p.3pD46 <- v_p_3aDage46 # probability to die when stage 3 preclinical
p.4pD46 <- v_p_4aDage46 # probability to die when stage 4 preclinical
p.1caD46 <- v_p_1aDage46  # probability to die when stage 1 clinical a
p.2caD46 <- v_p_2aDage46 # probability to die when stage 2 clinical a
p.3caD46 <- v_p_3aDage46 # probability to die when stage 3 clinical a
p.4caD46 <- v_p_4aDage46 # probability to die when stage 4 clinical a
p.1cbD46 <- v_p_1bDage46 # probability to die when stage 1 clinical b 
p.2cbD46 <- v_p_2bDage46 # probability to die when stage 2 clinical b
p.3cbD46 <- v_p_3bDage46 # probability to die when stage 3 clinical b
p.4cbD46 <- v_p_4bDage46 # probability to die when stage 4 clinical b
p.1ccD46 <- v_p_1cDage46  # probability to die when stage 1 clinical c
p.2ccD46 <- v_p_2cDage46 # probability to die when stage 2 clinical c
p.3ccD46 <- v_p_3cDage46 # probability to die when stage 3 clinical c
p.4ccD46 <- v_p_4cDage46 # probability to die when stage 4 clinical c
p.1cdD46 <- v_p_1dDage46  # probability to die when stage 1 clinical c
p.2cdD46 <- v_p_2dDage46 # probability to die when stage 2 clinical c
p.3cdD46 <- v_p_3dDage46 # probability to die when stage 3 clinical c
p.4cdD46 <- v_p_4dDage46 # probability to die when stage 4 clinical c
p.1ceD46 <- v_p_1eDage46  # probability to die when stage 1 clinical c
p.2ceD46 <- v_p_2eDage46 # probability to die when stage 2 clinical c
p.3ceD46 <- v_p_3eDage46 # probability to die when stage 3 clinical c
p.4ceD46 <- v_p_4eDage46 # probability to die when stage 4 clinical c
p.1cfD46 <- v_p_1fDage46  # probability to die when stage 1 clinical c
p.2cfD46 <- v_p_2fDage46 # probability to die when stage 2 clinical c
p.3cfD46 <- v_p_3fDage46 # probability to die when stage 3 clinical c
p.4cfD46 <- v_p_4fDage46 # probability to die when stage 4 clinical c


#Probability to die from 30 to 39 

p.NnD3039 <- v_p_HDage3039  # probability to die when healthy 
p.NpD3039 <- v_p_HDage3039  # probability to die when Normal helicobacter (+)
p.GnD3039 <- v_p_HDage3039 # probability to die when Gastritis helicobacter (-)
p.GpD3039 <- v_p_HDage3039 # probability to die when Gastritis helicobacter (+)
p.AnD3039 <- v_p_HDage3039 # probability to die when Atrophy helicobacter (-)
p.ApD3039 <- v_p_HDage3039 # probability to die when Atrophy helicobacter (+)
p.InD3039 <- v_p_HDage3039  # probability to die when Intestinal helicobacter (-) 
p.IpD3039 <- v_p_HDage3039  # probability to die when Intestinal helicobacter (+) 
p.DnD3039 <- v_p_HDage3039  # probability to die when Dysplasia helicobacter (-) 
p.DpD3039 <- v_p_HDage3039  # probability to die when Dysplasia helicobacter (+) 
p.1pD3039 <- v_p_1aDage3039  # probability to die when stage 1 preclinical 
p.2pD3039 <- v_p_2aDage3039 # probability to die when stage 2 preclinical
p.3pD3039 <- v_p_3aDage3039 # probability to die when stage 3 preclinical
p.4pD3039 <- v_p_4aDage3039 # probability to die when stage 4 preclinical
p.1caD3039 <- v_p_1aDage3039  # probability to die when stage 1 clinical a
p.2caD3039 <- v_p_2aDage3039 # probability to die when stage 2 clinical a
p.3caD3039 <- v_p_3aDage3039 # probability to die when stage 3 clinical a
p.4caD3039 <- v_p_4aDage3039 # probability to die when stage 4 clinical a
p.1cbD3039 <- v_p_1bDage3039 # probability to die when stage 1 clinical b 
p.2cbD3039 <- v_p_2bDage3039 # probability to die when stage 2 clinical b
p.3cbD3039 <- v_p_3bDage3039 # probability to die when stage 3 clinical b
p.4cbD3039 <- v_p_4bDage3039 # probability to die when stage 4 clinical b
p.1ccD3039 <- v_p_1cDage3039  # probability to die when stage 1 clinical c
p.2ccD3039 <- v_p_2cDage3039 # probability to die when stage 2 clinical c
p.3ccD3039 <- v_p_3cDage3039 # probability to die when stage 3 clinical c
p.4ccD3039 <- v_p_4cDage3039 # probability to die when stage 4 clinical c
p.1cdD3039 <- v_p_1dDage3039  # probability to die when stage 1 clinical c
p.2cdD3039 <- v_p_2dDage3039 # probability to die when stage 2 clinical c
p.3cdD3039 <- v_p_3dDage3039 # probability to die when stage 3 clinical c
p.4cdD3039 <- v_p_4dDage3039 # probability to die when stage 4 clinical c
p.1ceD3039 <- v_p_1eDage3039  # probability to die when stage 1 clinical c
p.2ceD3039 <- v_p_2eDage3039 # probability to die when stage 2 clinical c
p.3ceD3039 <- v_p_3eDage3039 # probability to die when stage 3 clinical c
p.4ceD3039 <- v_p_4eDage3039 # probability to die when stage 4 clinical c
p.1cfD3039 <- v_p_1fDage3039  # probability to die when stage 1 clinical c
p.2cfD3039 <- v_p_2fDage3039 # probability to die when stage 2 clinical c
p.3cfD3039 <- v_p_3fDage3039 # probability to die when stage 3 clinical c
p.4cfD3039 <- v_p_4fDage3039 # probability to die when stage 4 clinical c

# Probability to die from 41 to 100
p.NnD41 <- v_p_HDage41  # probability to die when healthy 
p.NpD41 <- v_p_HDage41  # probability to die when Normal helicobacter (+)
p.GnD41 <- v_p_HDage41 # probability to die when Gastritis helicobacter (-)
p.GpD41 <- v_p_HDage41 # probability to die when Gastritis helicobacter (+)
p.AnD41 <- v_p_HDage41 # probability to die when Atrophy helicobacter (-)
p.ApD41 <- v_p_HDage41 # probability to die when Atrophy helicobacter (+)
p.InD41 <- v_p_HDage41  # probability to die when Intestinal helicobacter (-) 
p.IpD41 <- v_p_HDage41  # probability to die when Intestinal helicobacter (+) 
p.DnD41 <- v_p_HDage41  # probability to die when Dysplasia helicobacter (-) 
p.DpD41 <- v_p_HDage41  # probability to die when Dysplasia helicobacter (+) 
p.1pD41 <- v_p_1aDage41  # probability to die when stage 1 preclinical 
p.2pD41 <- v_p_2aDage41 # probability to die when stage 2 preclinical
p.3pD41 <- v_p_3aDage41 # probability to die when stage 3 preclinical
p.4pD41 <- v_p_4aDage41 # probability to die when stage 4 preclinical
p.1caD41 <- v_p_1aDage41  # probability to die when stage 1 clinical a
p.2caD41 <- v_p_2aDage41 # probability to die when stage 2 clinical a
p.3caD41 <- v_p_3aDage41 # probability to die when stage 3 clinical a
p.4caD41 <- v_p_4aDage41 # probability to die when stage 4 clinical a
p.1cbD41 <- v_p_1bDage41 # probability to die when stage 1 clinical b 
p.2cbD41 <- v_p_2bDage41 # probability to die when stage 2 clinical b
p.3cbD41 <- v_p_3bDage41 # probability to die when stage 3 clinical b
p.4cbD41 <- v_p_4bDage41 # probability to die when stage 4 clinical b
p.1ccD41 <- v_p_1cDage41  # probability to die when stage 1 clinical c
p.2ccD41 <- v_p_2cDage41 # probability to die when stage 2 clinical c
p.3ccD41 <- v_p_3cDage41 # probability to die when stage 3 clinical c
p.4ccD41 <- v_p_4cDage41 # probability to die when stage 4 clinical c
p.1cdD41 <- v_p_1dDage41  # probability to die when stage 1 clinical c
p.2cdD41 <- v_p_2dDage41 # probability to die when stage 2 clinical c
p.3cdD41 <- v_p_3dDage41 # probability to die when stage 3 clinical c
p.4cdD41 <- v_p_4dDage41 # probability to die when stage 4 clinical c
p.1ceD41 <- v_p_1eDage41  # probability to die when stage 1 clinical c
p.2ceD41 <- v_p_2eDage41 # probability to die when stage 2 clinical c
p.3ceD41 <- v_p_3eDage41 # probability to die when stage 3 clinical c
p.4ceD41 <- v_p_4eDage41 # probability to die when stage 4 clinical c
p.1cfD41 <- v_p_1fDage41  # probability to die when stage 1 clinical c
p.2cfD41 <- v_p_2fDage41 # probability to die when stage 2 clinical c
p.3cfD41 <- v_p_3fDage41 # probability to die when stage 3 clinical c
p.4cfD41 <- v_p_4fDage41 # probability to die when stage 4 clinical c



  ##    02.5 Costs -------------------------------------------------------------------


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
c.1ca  <- 66.64033534                     # cost of remaining one cycle stage 1 clinical a (treatment)
c.2ca  <- 124.1583862                    # cost of remaining one cycle stage 2 clinical a (treatment)
c.3ca  <- 326.1405023                     # cost of remaining one cycle stage 3 clinical a (treatment)
c.4ca  <- 197.5762362                     # cost of remaining one cycle stage 4 clinical a (treatment)
c.1cb  <- 2.003344331                    # cost of remaining one cycle stage 1 clinical b (treatment)
c.2cb  <- 8.761306628                     # cost of remaining one cycle stage 2 clinical b (treatment)
c.3cb  <- 8.130164988                     # cost of remaining one cycle stage 3 clinical b (treatment)
c.4cb  <- 8.130164988                     # cost of remaining one cycle stage 4 clinical b (treatment)
c.1cc  <- 2.003344331                     # cost of remaining one cycle stage 1 clinical c (treatment)
c.2cc  <- 8.020245578                     # cost of remaining one cycle stage 2 clinical c (treatment)
c.3cc  <- 7.389103938                     # cost of remaining one cycle stage 3 clinical c (treatment)
c.4cc  <- 7.389103938
c.1cd  <- 1.653509662                     # cost of remaining one cycle stage 1 clinical c (treatment)
c.2cd  <- 5.476832307                     # cost of remaining one cycle stage 2 clinical c (treatment)
c.3cd  <- 4.214549026                     # cost of remaining one cycle stage 3 clinical c (treatment)
c.4cd  <- 4.214549026
c.1ce  <- 1.653509662                     # cost of remaining one cycle stage 1 clinical c (treatment)
c.2ce  <- 5.476832307                    # cost of remaining one cycle stage 2 clinical c (treatment)
c.3ce  <- 4.214549026                     # cost of remaining one cycle stage 3 clinical c (treatment)
c.4ce  <- 4.214549026
c.1cf  <- 1.653509662                     # cost of remaining one cycle stage 1 clinical c (treatment)
c.2cf  <- 2.279706887                     # cost of remaining one cycle stage 2 clinical c (treatment)
c.3cf  <- 2.279706887                     # cost of remaining one cycle stage 3 clinical c (treatment)
c.4cf  <- 2.279706887

c.D  <- 0  # cost of remaining one cycle dead

#Prevention costs (UF)
c.UreaAire <- 2.86291                        #Costeo  
c.EDABp <- 1.96434                           #Costeo 
c.PSEda <- 0.8887518696                                # Costeo 
c.Antigenofecal <- 0.15167                  #Costeo  
c.Serología <- 0.395619914                       #Costeo 

#Erradication costs (UF)
c.errp.1 <- 0.2812988       # Costeo erradicación para estados positivos estrategia 1 
c.errp.2 <- 0.27537672      # Costeo erradicación para estados positivos estrategia 2
c.errp.3 <- 0.2812988       # Costeo erradicación para estados positivos estrategia 3
c.errp.4 <- 0.01706795866
c.errp.5 <- 0.2643616512    # Costeo erradicación para estados positivos estrategia 5 
c.errn.1 <- 0.0148052       # Costeo erradicación para estados negativos estrategia 1 
c.errn.2 <- 0.0736374       # Costeo erradicación para estados negativos estrategia 2 
c.errn.5 <- 0.067559328     # Costeo erradicación para estados negativos estrategia 5 
  

  ##    02.6 Utilities ---------------------------------------------------------------


u.Nn  <- 0.95                     # utility when Normal helicobacter (-)
u.Np  <- 0.95                     # utility when Normal helicobacter (+)
u.Gn  <- 0.95                     # utility when Gastritis helicobacter (-)
u.Gp  <- 0.95                     # utility when Gastritis helicobacter (+)
u.An  <- 0.95                     # utility when Atrophy helicobacter (-)
u.Ap  <- 0.95                     # utility when Atrophy helicobacter (+)
u.In  <- 0.95                     # utility when Intestinal helicobacter (-)
u.Ip  <- 0.95                     # utility when Intestinal helicobacter (+)
u.Dn  <- 0.95                     # utility when Dysplasia helicobacter (-)
u.Dp  <- 0.95                     # utility when Dysplasia helicobacter (+) 
u.1p  <- 0.8070                     # utility when stage 1 preclinical
u.2p  <- 0.7238                     # utility when stage 2 preclinical
u.3p  <- 0.6477                     # utility when stage 3 preclinical
u.4p  <- 0.5345                       # utility when stage 4 preclinical
u.1ca  <- 0.8070                     # utility when stage 1 clinical a
u.2ca  <- 0.7238                     # utility when stage 2 clinical a
u.3ca  <- 0.6477                     # utility when stage 3 clinical a
u.4ca  <- 0.5345                       # utility when stage 4 clinical a
u.1cb  <- 0.8070                     # utility when stage 1 clinical b
u.2cb  <- 0.7238                      # utility when stage 2 clinical b
u.3cb  <- 0.6477                     # utility when stage 3 clinical b
u.4cb  <- 0.5345                       # utility when stage 4 clinical b
u.1cc  <- 0.8070                     # utility when stage 1 clinical c
u.2cc  <- 0.7238                     # utility when stage 2 clinical c
u.3cc  <- 0.6477                     # utility when stage 3 clinical c
u.4cc  <- 0.5345                       # utility when stage 4 clinical c
u.1cd  <- 0.8070                     # utility when stage 1 clinical c
u.2cd  <- 0.7238                    # utility when stage 2 clinical c
u.3cd  <- 0.6477                     # utility when stage 3 clinical c
u.4cd  <- 0.5345                       # utility when stage 4 clinical c
u.1ce  <- 0.8070                     # utility when stage 1 clinical c
u.2ce  <- 0.7238                     # utility when stage 2 clinical c
u.3ce  <- 0.6477                     # utility when stage 3 clinical c
u.4ce  <- 0.5345                       # utility when stage 4 clinical c
u.1cf  <- 0.8070                     # utility when stage 1 clinical c
u.2cf  <- 0.7238                     # utility when stage 2 clinical c
u.3cf  <- 0.6477                     # utility when stage 3 clinical c
u.4cf  <- 0.5345                       # utility when stage 4 clinical c

u.D  <- 0                       # utility when dead


  ##    02.7 Discounting factor  -----------------------------------------------------

d.r  <- 0.03                    # equal discount of costs and QALYs by 3%
v.dwc <- 1/((1 + d.r) ^ (0:n.t)) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1/((1 + d.r) ^ (0:n.t)) # calculate discount weights for effectiveness for each cycle based on discount rate d.r


#### 03 Define and initialize matrices and vectors ####
  ##    03.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_prev <- m.M_prev <- m.M_prev2 <- m.M_prev3 <- m.M_prev4 <- m.M_prev5  <- matrix(NA, 
                                                    nrow = n.t + 1, ncol = n.s,
                                            dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_prev) # show first 6 rows of the matrix 

#  

m.M_no_prev[1, ] <- m.M_prev[1, ] <- m.M_prev2[1, ] <- m.M_prev3[1, ] <- m.M_prev4[1, ] <- m.M_prev5[1, ]  <-  c(0.044, 0.039 , 0.176 , 0.702 , 0, 0.0195, 0, 0.0195, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0, 0, 0 , 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)                     # initialize first cycle of Markov trace



  ##    03.2 Transition probability MATRIX ####
# create the transition probability matrix without prevention
#m.P_noprev  <- matrix(0,
#nrow = n.s,
#ncol = n.s,
#dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
#m.P_noprev


# create the transition probability matrix with prevention strategy 1 
#m.P_prev  <- matrix(0,
#nrow = n.s,
#ncol = n.s,
# dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix
#m.P_prev

# create the transition probability matrix with prevention strategy 2 
#m.P_prev2  <- matrix(0,
#nrow = n.s,
#ncol = n.s,
#dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

#m.P_prev2

# create the transition probability matrix with prevention strategy 3 
#m.P_prev3  <- matrix(0,
#nrow = n.s,
#ncol = n.s,
#dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

#m.P_prev3

# create the transition probability matrix with prevention strategy 4 
#m.P_prev4  <- matrix(0,
# nrow = n.s,
#ncol = n.s,
#dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

#m.P_prev4

# create the transition probability matrix with prevention strategy 5 
#m.P_prev5  <- matrix(0,
#nrow = n.s,
#ncol = n.s,
# dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

#m.P_prev5

#create the transition probability array with no prevention 

a.P_noprev <- array(0, 
                  dim = c(n.s, n.s, n.t), 
                  dimnames = list( from  = v.n, 
                                   to    = v.n, 
                                   cycle = 1:n.t))

# fill in the transition probability array without prevention
### From Normal helicobacter (-)
a.P_noprev["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_noprev["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_noprev["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_noprev["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
a.P_noprev["Normalhpp", "Normalhpp", ]    <- 1 - p.NpGp - p.NpNn - p.NpD 
a.P_noprev["Normalhpp", "Gastritishpp", ]    <- p.NpGp
a.P_noprev["Normalhpp", "Normalhpn", ]    <- p.NpNn
a.P_noprev["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_noprev["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_noprev["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_noprev["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_noprev["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_noprev["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_noprev["Gastritishpp", "Gastritishpp", ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD
a.P_noprev["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
a.P_noprev["Gastritishpp", "Gastritishpn", ]    <- p.GpGn
a.P_noprev["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_noprev["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_noprev["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn - p.AnD
a.P_noprev["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_noprev["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_noprev["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_noprev["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (-) 
a.P_noprev["Atrophyhpp", "Atrophyhpp", ]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD
a.P_noprev["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
a.P_noprev["Atrophyhpp", "Atrophyhpn", ]    <- p.ApAn
a.P_noprev["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_noprev["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_noprev["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD
a.P_noprev["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_noprev["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_noprev["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_noprev["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_noprev["Intestinalhpp", "Intestinalhpp", ]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD  
a.P_noprev["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
a.P_noprev["Intestinalhpp", "Intestinalhpn", ]    <- p.IpIn
a.P_noprev["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_noprev["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_noprev["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn - p.Dn1p - p.DnD
a.P_noprev["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_noprev["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_noprev["Dysplasiahpn", "1preclinical", ]     <- p.Dn1p
a.P_noprev["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_noprev["Dysplasiahpp", "Dysplasiahpp", ]    <- 1 - p.DpDn - p.DpIp - p.Dp1p - p.DpD
a.P_noprev["Dysplasiahpp", "Dysplasiahpn", ]    <- p.DpDn
a.P_noprev["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_noprev["Dysplasiahpp", "1preclinical", ]     <- p.Dp1p
a.P_noprev["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_noprev["1preclinical", "1preclinical",]    <- 1 - p.1p1ca - p.1p2p - p.1pD	
a.P_noprev["1preclinical", "1clinicala",]    <- p.1p1ca	
a.P_noprev["1preclinical", "2preclinical",]    <- p.1p2p	
a.P_noprev["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_noprev["2preclinical", "2preclinical",]    <- 1 - p.2p2ca - p.2p3p - p.2pD	
a.P_noprev["2preclinical", "2clinicala", ]    <- p.2p2ca	
a.P_noprev["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_noprev["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_noprev["3preclinical", "3preclinical", ]    <- 1 - p.3p3ca - p.3p4p - p.3pD	
a.P_noprev["3preclinical", "3clinicala", ]    <- p.3p3ca	
a.P_noprev["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_noprev["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_noprev["4preclinical", "4preclinical", ]    <- 1 - p.4p4ca - p.4pD	
a.P_noprev["4preclinical", "4clinicala", ]    <- p.4p4ca	
a.P_noprev["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_noprev["1clinicala", "1clinicala", ]    <- 0	
a.P_noprev["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_noprev["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_noprev["2clinicala", "2clinicala", ]    <- 0	
a.P_noprev["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_noprev["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_noprev["3clinicala", "3clinicala", ]    <- 0 	
a.P_noprev["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_noprev["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_noprev["4clinicala", "4clinicala", ]    <- 0	
a.P_noprev["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_noprev["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_noprev["1clinicalb", "1clinicalb", ]    <- 0	
a.P_noprev["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_noprev["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_noprev["2clinicalb", "2clinicalb", ]    <- 0	
a.P_noprev["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_noprev["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_noprev["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_noprev["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_noprev["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_noprev["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_noprev["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_noprev["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_noprev["1clinicalc", "1clinicalc", ]    <- 0
a.P_noprev["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_noprev["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_noprev["1clinicald", "1clinicald", ]    <- 0
a.P_noprev["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_noprev["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_noprev["1clinicale", "1clinicale", ]    <- 0
a.P_noprev["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_noprev["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_noprev["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_noprev["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_noprev["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_noprev["2clinicalc", "2clinicalc", ]    <- 0
a.P_noprev["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_noprev["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_noprev["2clinicald", "2clinicald", ]    <- 0
a.P_noprev["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_noprev["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_noprev["2clinicale", "2clinicale", ]    <- 0
a.P_noprev["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_noprev["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_noprev["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_noprev["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_noprev["2clinicalf", "Dead", ]    <- p.2cfD	

### From 3 clinical c	
a.P_noprev["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_noprev["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_noprev["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_noprev["3clinicald", "3clinicald", ]    <- 0
a.P_noprev["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_noprev["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_noprev["3clinicale", "3clinicale", ]    <- 0
a.P_noprev["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_noprev["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_noprev["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_noprev["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_noprev["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_noprev["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_noprev["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_noprev["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_noprev["4clinicald", "4clinicald", ]    <- 0
a.P_noprev["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_noprev["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_noprev["4clinicale", "4clinicale", ]    <- 0
a.P_noprev["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_noprev["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_noprev["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_noprev["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_noprev["4clinicalf", "Dead", ]    <- p.4cfD	

### From Dead
a.P_noprev["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_noprev)

# fill in the transition probability matrix with prevention strategy 1 

a.P_prev <- array(0, 
                  dim = c(n.s, n.s, n.t), 
                  dimnames = list( from  = v.n, 
                                   to    = v.n, 
                                   cycle = 1:n.t))


a.P_prev["Gastritishpp", "Gastritishpn", 1] <- p.pGpGn
a.P_prev["Atrophyhpp", "Atrophyhpn", 1] <- p.pApAn
a.P_prev["Intestinalhpp", "Intestinalhpn", 1] <- p.pIpIn
a.P_prev["Dysplasiahpp", "Dysplasiahpn", 1] <- p.pDpDn
a.P_prev["Gastritishpp", "Gastritishpn", 2:n.t] <- p.GpGn
a.P_prev["Atrophyhpp", "Atrophyhpn", 2:n.t] <- p.ApAn
a.P_prev["Intestinalhpp", "Intestinalhpn", 2:n.t] <- p.IpIn
a.P_prev["Dysplasiahpp", "Dysplasiahpn", 2:n.t] <- p.DpDn

### From Normal helicobacter (-)
a.P_prev["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_prev["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_prev["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_prev["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
largosecuencia <- length(p.NpD)
a.P_prev["Normalhpp", "Normalhpn", 1] <- p.pNpNn
a.P_prev["Normalhpp", "Normalhpn", 2:n.t] <- p.NpNn
a.P_prev["Normalhpp", "Normalhpp", 1]  <- 1 - p.NpD30 - p.NpGp - p.pNpNn
a.P_prev["Normalhpp", "Normalhpp", 2:n.t]    <- 1 - p.NpGp - p.NpNn - p.NpD31 
a.P_prev["Normalhpp", "Gastritishpp", ]    <- p.NpGp
##a.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn --> Probabilidad que cambia, está arriba 
a.P_prev["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_prev["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_prev["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_prev["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_prev["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_prev["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_prev["Gastritishpp", "Gastritishpp", 1 ]    <- 1 - p.GpAp - p.pGpGn - p.GpNp - p.GpD30
a.P_prev["Gastritishpp", "Gastritishpp", 2:n.t ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD31
a.P_prev["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
##a.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
a.P_prev["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_prev["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_prev["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn -  p.AnD
a.P_prev["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_prev["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_prev["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_prev["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (+) 
a.P_prev["Atrophyhpp", "Atrophyhpp", 1]    <- 1 - p.ApIp - p.pApAn - p.ApGp - p.ApD30
a.P_prev["Atrophyhpp", "Atrophyhpp", 2:n.t]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD31
a.P_prev["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
##a.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
a.P_prev["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_prev["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_prev["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD 
a.P_prev["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_prev["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_prev["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_prev["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_prev["Intestinalhpp", "Intestinalhpp", 1]    <- 1 - p.IpDp - p.pIpIn - p.IpAp - p.IpD30
a.P_prev["Intestinalhpp", "Intestinalhpp", 2:n.t]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD31  
a.P_prev["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
#a.P_prev["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
a.P_prev["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_prev["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_prev["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn -  p.DnD - p.Dn1p
a.P_prev["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_prev["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_prev["Dysplasiahpn", "1preclinical", ]    <- p.Dn1p
a.P_prev["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_prev["Dysplasiahpp", "Dysplasiahpp", 1]    <- 1 - p.pDpDn - p.DpIp - p.DpD30 - p.Dp1p
a.P_prev["Dysplasiahpp", "Dysplasiahpp", 2:n.t]    <- 1 - p.DpDn - p.DpIp - p.DpD31 - p.Dp1p
##a.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
a.P_prev["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_prev["Dysplasiahpp", "1preclinical", ]    <- p.Dp1p
a.P_prev["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_prev["1preclinical", "1preclinical",]    <- 1 - p.1p1ca - p.1p2p - p.1pD	
a.P_prev["1preclinical", "1clinicala",]    <- p.1p1ca	
a.P_prev["1preclinical", "2preclinical",]    <- p.1p2p	
a.P_prev["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_prev["2preclinical", "2preclinical",]    <- 1 - p.2p2ca - p.2p3p - p.2pD	
a.P_prev["2preclinical", "2clinicala", ]    <- p.2p2ca	
a.P_prev["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_prev["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_prev["3preclinical", "3preclinical", ]    <- 1 - p.3p3ca - p.3p4p - p.3pD	
a.P_prev["3preclinical", "3clinicala", ]    <- p.3p3ca	
a.P_prev["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_prev["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_prev["4preclinical", "4preclinical", ]    <- 1 - p.4p4ca - p.4pD	
a.P_prev["4preclinical", "4clinicala", ]    <- p.4p4ca	
a.P_prev["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_prev["1clinicala", "1clinicala", ]    <- 0	
a.P_prev["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_prev["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_prev["2clinicala", "2clinicala", ]    <- 0	
a.P_prev["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_prev["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_prev["3clinicala", "3clinicala", ]    <- 0 	
a.P_prev["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_prev["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_prev["4clinicala", "4clinicala", ]    <- 0	
a.P_prev["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_prev["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_prev["1clinicalb", "1clinicalb", ]    <- 0	
a.P_prev["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_prev["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_prev["2clinicalb", "2clinicalb", ]    <- 0	
a.P_prev["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_prev["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_prev["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_prev["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_prev["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_prev["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_prev["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_prev["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_prev["1clinicalc", "1clinicalc", ]    <- 0
a.P_prev["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_prev["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_prev["1clinicald", "1clinicald", ]    <- 0
a.P_prev["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_prev["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_prev["1clinicale", "1clinicale", ]    <- 0
a.P_prev["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_prev["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_prev["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_prev["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_prev["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_prev["2clinicalc", "2clinicalc", ]    <- 0
a.P_prev["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_prev["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_prev["2clinicald", "2clinicald", ]    <- 0
a.P_prev["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_prev["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_prev["2clinicale", "2clinicale", ]    <- 0
a.P_prev["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_prev["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_prev["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_prev["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_prev["2clinicalf", "Dead", ]    <- p.2cfD	

### From 3 clinical c	
a.P_prev["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_prev["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_prev["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_prev["3clinicald", "3clinicald", ]    <- 0
a.P_prev["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_prev["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_prev["3clinicale", "3clinicale", ]    <- 0
a.P_prev["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_prev["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_prev["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_prev["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_prev["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_prev["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_prev["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_prev["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_prev["4clinicald", "4clinicald", ]    <- 0
a.P_prev["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_prev["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_prev["4clinicale", "4clinicale", ]    <- 0
a.P_prev["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_prev["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_prev["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_prev["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_prev["4clinicalf", "Dead", ]    <- p.4cfD	

### From Dead
a.P_prev["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_prev)
###


### From Normal helicobacter (-)
#m.P_prev["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.pNnNp - p.pNnGn
#m.P_prev["Normalhpn", "Normalhpp"]    <- p.pNnNp
#m.P_prev["Normalhpn", "Gastritishpn"]    <- p.pNnGn
#m.P_prev["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
#m.P_prev["Normalhpp", "Normalhpp"]    <- 1 - p.pNpGp - p.pNpNn - p.NpD 
#m.P_prev["Normalhpp", "Gastritishpp"]    <- p.pNpGp
#m.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn
#m.P_prev["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
#m.P_prev["Gastritishpn", "Gastritishpn"]    <- 1 - p.pGnAn - p.pGnGp - p.pGnNn - p.GnD
#m.P_prev["Gastritishpn", "Atrophyhpn"]    <- p.pGnAn
#m.P_prev["Gastritishpn", "Gastritishpp"]    <- p.pGnGp
#m.P_prev["Gastritishpn", "Normalhpn"]    <- p.pGnNn
#m.P_prev["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
#m.P_prev["Gastritishpp", "Gastritishpp"]    <- 1 - p.pGpAp - p.pGpGn - p.pGpNp - p.GpD
#m.P_prev["Gastritishpp", "Atrophyhpp"]    <- p.pGpAp
#m.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
#m.P_prev["Gastritishpp", "Normalhpp"]    <- p.pGpNp
#m.P_prev["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
#m.P_prev["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.pAnIn - p.pAnAp - p.pAnGn -  p.AnD
#m.P_prev["Atrophyhpn", "Intestinalhpn"]    <- p.pAnIn
#m.P_prev["Atrophyhpn", "Atrophyhpp"]    <- p.pAnAp
#m.P_prev["Atrophyhpn", "Gastritishpn"]    <- p.pAnGn
#m.P_prev["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
#m.P_prev["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.pApIp - p.pApAn - p.pApGp - p.ApD
#m.P_prev["Atrophyhpp", "Intestinalhpp"]    <- p.pApIp
#m.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
#m.P_prev["Atrophyhpp", "Gastritishpp"]    <- p.pApGp
#m.P_prev["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
#m.P_prev["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.pInDn - p.pInIp - p.pInAn - p.InD 
#m.P_prev["Intestinalhpn", "Dysplasiahpn"]    <- p.pInDn
#m.P_prev["Intestinalhpn", "Intestinalhpp"]    <- p.pInIp
#m.P_prev["Intestinalhpn", "Atrophyhpn"]    <- p.pInAn
#m.P_prev["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
#m.P_prev["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.pIpDp - p.pIpIn - p.pIpAp - p.IpD  
#m.P_prev["Intestinalhpp", "Dysplasiahpp"]    <- p.pIpDp
#m.P_prev["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
#m.P_prev["Intestinalhpp", "Atrophyhpp"]    <- p.pIpAp
#m.P_prev["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
#m.P_prev["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.pDnDp - p.pDnIn -  p.DnD - p.pDn1p
#m.P_prev["Dysplasiahpn", "Dysplasiahpp"]    <- p.pDnDp
#m.P_prev["Dysplasiahpn", "Intestinalhpn"]    <- p.pDnIn
#m.P_prev["Dysplasiahpn", "1preclinical"]    <- p.pDn1p
#m.P_prev["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
#m.P_prev["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.pDpDn - p.pDpIp - p.DpD - p.pDp1p
#m.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
#m.P_prev["Dysplasiahpp", "Intestinalhpp"]    <- p.pDpIp
#m.P_prev["Dysplasiahpp", "1preclinical"]    <- p.pDp1p
#m.P_prev["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
#m.P_prev["1preclinical", "1preclinical"]    <- 1 - p.p1p1ca - p.p1p2p - p.1pD	
#m.P_prev["1preclinical", "1clinicala"]    <- p.p1p1ca	
#m.P_prev["1preclinical", "2preclinical"]    <- p.p1p2p	
#m.P_prev["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
#m.P_prev["2preclinical", "2preclinical"]    <- 1 - p.p2p2ca - p.p2p3p - p.2pD	
#m.P_prev["2preclinical", "2clinicala"]    <- p.p2p2ca	
#m.P_prev["2preclinical", "3preclinical"]    <- p.p2p3p	
#m.P_prev["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
#m.P_prev["3preclinical", "3preclinical"]    <- 1 - p.p3p3ca - p.p3p4p - p.3pD	
#m.P_prev["3preclinical", "3clinicala"]    <- p.p3p3ca	
#m.P_prev["3preclinical", "4preclinical"]    <- p.p3p4p	
#m.P_prev["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
#m.P_prev["4preclinical", "4preclinical"]    <- 1 - p.p4p4ca - p.4pD	
#m.P_prev["4preclinical", "4clinicala"]    <- p.p4p4ca	
#m.P_prev["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
#m.P_prev["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p1ca1cb	
#m.P_prev["1clinicala", "1clinicalb"] <- p.p1ca1cb	
#m.P_prev["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
#m.P_prev["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p2ca2cb	
#m.P_prev["2clinicala", "2clinicalb"] <- p.p2ca2cb	
#m.P_prev["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
#m.P_prev["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p3ca3cb 	
#m.P_prev["3clinicala", "3clinicalb"] <- p.p3ca3cb	
#m.P_prev["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
#m.P_prev["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p4ca4cb	
#m.P_prev["4clinicala", "4clinicalb"] <- p.p4ca4cb  	
#m.P_prev["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
#m.P_prev["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p1cb1cc	
#m.P_prev["1clinicalb", "1clinicalc"] <- p.p1cb1cc	
#m.P_prev["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
#m.P_prev["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p2cb2cc	
#m.P_prev["2clinicalb", "2clinicalc"] <- p.p2cb2cc	
#m.P_prev["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
#m.P_prev["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p3cb3cc 	
#m.P_prev["3clinicalb", "3clinicalc"] <- p.p3cb3cc	
#m.P_prev["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
#m.P_prev["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p4cb4cc 	
#m.P_prev["4clinicalb", "4clinicalc"] <- p.p4cb4cc	
#m.P_prev["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
#m.P_prev["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p1ccNn	
#m.P_prev["1clinicalc", "Normalhpn"] <- p.p1ccNn	
#m.P_prev["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
#m.P_prev["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p2ccNn	
#m.P_prev["2clinicalc", "Normalhpn"] <- p.p2ccNn	
#m.P_prev["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
#m.P_prev["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p3ccNn 	
#m.P_prev["3clinicalc", "Normalhpn"] <- p.p3ccNn
#m.P_prev["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
#m.P_prev["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p4ccNn 	
#m.P_prev["4clinicalc", "Normalhpn"] <- p.p4ccNn	
#m.P_prev["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
#m.P_prev["Dead", "Dead"] <- 1

# check rows add up to 1
#rowSums(m.P_prev)


# fill in the transition probability matrix with prevention strategy 2

# fill in the transition probability array with prevention strategy 2 

a.P_prev2 <- array(0, 
                  dim = c(n.s, n.s, n.t), 
                  dimnames = list( from  = v.n, 
                                   to    = v.n, 
                                   cycle = 1:n.t))

a.P_prev2["Normalhpp", "Normalhpn", 1] <- p.p2NpNn
a.P_prev2["Gastritishpp", "Gastritishpn", 1] <- p.p2GpGn
a.P_prev2["Atrophyhpp", "Atrophyhpn", 1] <- p.p2ApAn
a.P_prev2["Intestinalhpp", "Intestinalhpn", 1] <- p.p2IpIn
a.P_prev2["Dysplasiahpp", "Dysplasiahpn", 1] <- p.p2DpDn
a.P_prev2["Normalhpp", "Normalhpn", 2:n.t] <- p.NpNn
a.P_prev2["Gastritishpp", "Gastritishpn", 2:n.t] <- p.GpGn
a.P_prev2["Atrophyhpp", "Atrophyhpn", 2:n.t] <- p.ApAn
a.P_prev2["Intestinalhpp", "Intestinalhpn", 2:n.t] <- p.IpIn
a.P_prev2["Dysplasiahpp", "Dysplasiahpn", 2:n.t] <- p.DpDn

### From Normal helicobacter (-)
a.P_prev2["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_prev2["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_prev2["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_prev2["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
a.P_prev2["Normalhpp", "Normalhpp", 1 ]    <- 1 - p.NpGp - p.p2NpNn - p.NpD30
a.P_prev2["Normalhpp", "Normalhpp", 2:n.t]    <- 1 - p.NpGp - p.NpNn - p.NpD31 
a.P_prev2["Normalhpp", "Gastritishpp", ]    <- p.NpGp
##a.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn --> Probabilidad que cambia, está arriba 
a.P_prev2["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_prev2["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_prev2["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_prev2["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_prev2["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_prev2["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_prev2["Gastritishpp", "Gastritishpp", 1 ]    <- 1 - p.GpAp - p.p2GpGn - p.GpNp - p.GpD30
a.P_prev2["Gastritishpp", "Gastritishpp", 2:n.t ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD31
a.P_prev2["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
##a.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
a.P_prev2["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_prev2["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_prev2["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn -  p.AnD
a.P_prev2["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_prev2["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_prev2["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_prev2["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (+) 
a.P_prev2["Atrophyhpp", "Atrophyhpp", 1]    <- 1 - p.ApIp - p.p2ApAn - p.ApGp - p.ApD30
a.P_prev2["Atrophyhpp", "Atrophyhpp", 2:n.t]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD31
a.P_prev2["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
##a.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
a.P_prev2["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_prev2["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_prev2["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD 
a.P_prev2["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_prev2["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_prev2["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_prev2["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_prev2["Intestinalhpp", "Intestinalhpp", 1]    <- 1 - p.IpDp - p.p2IpIn - p.IpAp - p.IpD30
a.P_prev2["Intestinalhpp", "Intestinalhpp", 2:n.t]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD31  
a.P_prev2["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
#a.P_prev["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
a.P_prev2["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_prev2["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_prev2["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn -  p.DnD - p.Dn1p
a.P_prev2["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_prev2["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_prev2["Dysplasiahpn", "1preclinical", ]    <- p.Dn1p
a.P_prev2["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_prev2["Dysplasiahpp", "Dysplasiahpp", 1]    <- 1 - p.p2DpDn - p.DpIp - p.DpD30 - p.Dp1p
a.P_prev2["Dysplasiahpp", "Dysplasiahpp", 2:n.t]    <- 1 - p.DpDn - p.DpIp - p.DpD31 - p.Dp1p
##a.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
a.P_prev2["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_prev2["Dysplasiahpp", "1preclinical", ]    <- p.Dp1p
a.P_prev2["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_prev2["1preclinical", "1preclinical",]    <- 1 - p.1p1ca - p.1p2p - p.1pD	
a.P_prev2["1preclinical", "1clinicala",]    <- p.1p1ca	
a.P_prev2["1preclinical", "2preclinical",]    <- p.1p2p	
a.P_prev2["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_prev2["2preclinical", "2preclinical",]    <- 1 - p.2p2ca - p.2p3p - p.2pD	
a.P_prev2["2preclinical", "2clinicala", ]    <- p.2p2ca	
a.P_prev2["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_prev2["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_prev2["3preclinical", "3preclinical", ]    <- 1 - p.3p3ca - p.3p4p - p.3pD	
a.P_prev2["3preclinical", "3clinicala", ]    <- p.3p3ca	
a.P_prev2["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_prev2["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_prev2["4preclinical", "4preclinical", ]    <- 1 - p.4p4ca - p.4pD	
a.P_prev2["4preclinical", "4clinicala", ]    <- p.4p4ca	
a.P_prev2["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_prev2["1clinicala", "1clinicala", ]    <- 0	
a.P_prev2["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_prev2["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_prev2["2clinicala", "2clinicala", ]    <- 0	
a.P_prev2["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_prev2["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_prev2["3clinicala", "3clinicala", ]    <- 0 	
a.P_prev2["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_prev2["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_prev2["4clinicala", "4clinicala", ]    <- 0	
a.P_prev2["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_prev2["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_prev2["1clinicalb", "1clinicalb", ]    <- 0	
a.P_prev2["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_prev2["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_prev2["2clinicalb", "2clinicalb", ]    <- 0	
a.P_prev2["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_prev2["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_prev2["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_prev2["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_prev2["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_prev2["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_prev2["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_prev2["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_prev2["1clinicalc", "1clinicalc", ]    <- 0
a.P_prev2["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_prev2["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_prev2["1clinicald", "1clinicald", ]    <- 0
a.P_prev2["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_prev2["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_prev2["1clinicale", "1clinicale", ]    <- 0
a.P_prev2["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_prev2["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_prev2["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_prev2["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_prev2["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_prev2["2clinicalc", "2clinicalc", ]    <- 0
a.P_prev2["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_prev2["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_prev2["2clinicald", "2clinicald", ]    <- 0
a.P_prev2["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_prev2["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_prev2["2clinicale", "2clinicale", ]    <- 0
a.P_prev2["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_prev2["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_prev2["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_prev2["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_prev2["2clinicalf", "Dead", ]    <- p.2cfD	

### From 3 clinical c	
a.P_prev2["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_prev2["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_prev2["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_prev2["3clinicald", "3clinicald", ]    <- 0
a.P_prev2["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_prev2["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_prev2["3clinicale", "3clinicale", ]    <- 0
a.P_prev2["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_prev2["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_prev2["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_prev2["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_prev2["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_prev2["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_prev2["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_prev2["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_prev2["4clinicald", "4clinicald", ]    <- 0
a.P_prev2["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_prev2["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_prev2["4clinicale", "4clinicale", ]    <- 0
a.P_prev2["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_prev2["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_prev2["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_prev2["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_prev2["4clinicalf", "Dead", ]    <- p.4cfD	

### From Dead
a.P_prev2["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_prev2)


### Matriz de transición que ya no uso porque se cambió por el arreglo 


### From Normal helicobacter (-)
#m.P_prev2["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p2NnNp - p.p2NnGn
#m.P_prev2["Normalhpn", "Normalhpp"]    <- p.p2NnNp
#m.P_prev2["Normalhpn", "Gastritishpn"]    <- p.p2NnGn
#m.P_prev2["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
#m.P_prev2["Normalhpp", "Normalhpp"]    <- 1 - p.p2NpGp - p.p2NpNn - p.NpD 
#m.P_prev2["Normalhpp", "Gastritishpp"]    <- p.p2NpGp
#m.P_prev2["Normalhpp", "Normalhpn"]    <- p.p2NpNn
#m.P_prev2["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
#m.P_prev2["Gastritishpn", "Gastritishpn"]    <- 1 - p.p2GnAn - p.p2GnGp - p.p2GnNn - p.GnD
#m.P_prev2["Gastritishpn", "Atrophyhpn"]    <- p.p2GnAn
#m.P_prev2["Gastritishpn", "Gastritishpp"]    <- p.p2GnGp
#m.P_prev2["Gastritishpn", "Normalhpn"]    <- p.p2GnNn
#m.P_prev2["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
#m.P_prev2["Gastritishpp", "Gastritishpp"]    <- 1 - p.p2GpAp - p.p2GpGn - p.p2GpNp - p.GpD
#m.P_prev2["Gastritishpp", "Atrophyhpp"]    <- p.p2GpAp
#m.P_prev2["Gastritishpp", "Gastritishpn"]    <- p.p2GpGn
#m.P_prev2["Gastritishpp", "Normalhpp"]    <- p.p2GpNp
#m.P_prev2["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
#m.P_prev2["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p2AnIn - p.p2AnAp - p.p2AnGn -  p.AnD
#m.P_prev2["Atrophyhpn", "Intestinalhpn"]    <- p.p2AnIn
#m.P_prev2["Atrophyhpn", "Atrophyhpp"]    <- p.p2AnAp
#m.P_prev2["Atrophyhpn", "Gastritishpn"]    <- p.p2AnGn
#m.P_prev2["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
#m.P_prev2["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p2ApIp - p.p2ApAn - p.p2ApGp - p.ApD
#m.P_prev2["Atrophyhpp", "Intestinalhpp"]    <- p.p2ApIp
#m.P_prev2["Atrophyhpp", "Atrophyhpn"]    <- p.p2ApAn
#m.P_prev2["Atrophyhpp", "Gastritishpp"]    <- p.p2ApGp
#m.P_prev2["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
#m.P_prev2["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p2InDn - p.p2InIp - p.p2InAn - p.InD 
#m.P_prev2["Intestinalhpn", "Dysplasiahpn"]    <- p.p2InDn
#m.P_prev2["Intestinalhpn", "Intestinalhpp"]    <- p.p2InIp
#m.P_prev2["Intestinalhpn", "Atrophyhpn"]    <- p.p2InAn
#m.P_prev2["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
#m.P_prev2["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p2IpDp - p.p2IpIn - p.p2IpAp - p.IpD  
#m.P_prev2["Intestinalhpp", "Dysplasiahpp"]    <- p.p2IpDp
#m.P_prev2["Intestinalhpp", "Intestinalhpn"]    <- p.p2IpIn
#m.P_prev2["Intestinalhpp", "Atrophyhpp"]    <- p.p2IpAp
#m.P_prev2["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
#m.P_prev2["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p2DnDp - p.p2DnIn -  p.DnD - p.Dn1p
#m.P_prev2["Dysplasiahpn", "Dysplasiahpp"]    <- p.p2DnDp
#m.P_prev2["Dysplasiahpn", "Intestinalhpn"]    <- p.p2DnIn
#m.P_prev2["Dysplasiahpn", "1preclinical"]    <- p.Dn1p
#m.P_prev2["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
#m.P_prev2["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p2DpDn - p.p2DpIp - p.DpD - p.p2_Dp1p
#m.P_prev2["Dysplasiahpp", "Dysplasiahpn"]    <- p.p2DpDn
#m.P_prev2["Dysplasiahpp", "Intestinalhpp"]    <- p.p2DpIp
#m.P_prev2["Dysplasiahpp", "1preclinical"]    <- p.p2_Dp1p
#m.P_prev2["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
#m.P_prev2["1preclinical", "1preclinical"]    <- 1 - p.p2_1p1ca - p.p2_1p2p - p.1pD	
#m.P_prev2["1preclinical", "1clinicala"]    <- p.p2_1p1ca	
#m.P_prev2["1preclinical", "2preclinical"]    <- p.p2_1p2p	
#m.P_prev2["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
#m.P_prev2["2preclinical", "2preclinical"]    <- 1 - p.p2_2p2ca - p.p2_2p3p - p.2pD	
#m.P_prev2["2preclinical", "2clinicala"]    <- p.p2_2p2ca	
#m.P_prev2["2preclinical", "3preclinical"]    <- p.p2_2p3p	
#m.P_prev2["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
#m.P_prev2["3preclinical", "3preclinical"]    <- 1 - p.p2_3p3ca - p.p2_3p4p - p.3pD	
#m.P_prev2["3preclinical", "3clinicala"]    <- p.p2_3p3ca	
#m.P_prev2["3preclinical", "4preclinical"]    <- p.p2_3p4p	
#m.P_prev2["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
#m.P_prev2["4preclinical", "4preclinical"]    <- 1 - p.p2_4p4ca - p.4pD	
#m.P_prev2["4preclinical", "4clinicala"]    <- p.p2_4p4ca	
#m.P_prev2["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
#m.P_prev2["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p2_1ca1cb	
#m.P_prev2["1clinicala", "1clinicalb"] <- p.p2_1ca1cb	
#m.P_prev2["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
#m.P_prev2["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p2_2ca2cb	
#m.P_prev2["2clinicala", "2clinicalb"] <- p.p2_2ca2cb	
#m.P_prev2["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
#m.P_prev2["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p2_3ca3cb 	
#m.P_prev2["3clinicala", "3clinicalb"] <- p.p2_3ca3cb	
#m.P_prev2["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
#m.P_prev2["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p2_4ca4cb	
#m.P_prev2["4clinicala", "4clinicalb"] <- p.p2_4ca4cb  	
#m.P_prev2["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
#m.P_prev2["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p2_1cb1cc	
#m.P_prev2["1clinicalb", "1clinicalc"] <- p.p2_1cb1cc	
#m.P_prev2["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
#m.P_prev2["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p2_2cb2cc	
#m.P_prev2["2clinicalb", "2clinicalc"] <- p.p2_2cb2cc	
#m.P_prev2["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
#m.P_prev2["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p2_3cb3cc 	
#m.P_prev2["3clinicalb", "3clinicalc"] <- p.p2_3cb3cc	
#m.P_prev2["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
#m.P_prev2["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p2_4cb4cc 	
#m.P_prev2["4clinicalb", "4clinicalc"] <- p.p2_4cb4cc	
#m.P_prev2["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
#m.P_prev2["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p2_1ccNn	
#m.P_prev2["1clinicalc", "Normalhpn"] <- p.p2_1ccNn	
#m.P_prev2["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
#m.P_prev2["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p2_2ccNn	
#m.P_prev2["2clinicalc", "Normalhpn"] <- p.p2_2ccNn	
#m.P_prev2["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
#m.P_prev2["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p2_3ccNn 	
#m.P_prev2["3clinicalc", "Normalhpn"] <- p.p2_3ccNn
#m.P_prev2["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
#m.P_prev2["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p2_4ccNn 	
#m.P_prev2["4clinicalc", "Normalhpn"] <- p.p2_4ccNn	
#m.P_prev2["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
#m.P_prev2["Dead", "Dead"] <- 1

# check rows add up to 1
#rowSums(m.P_prev2)


# fill in the transition probability matrix with prevention strategy 3
a.P_prev3 <- array(0, 
                   dim = c(n.s, n.s, n.t), 
                   dimnames = list( from  = v.n, 
                                    to    = v.n, 
                                    cycle = 1:n.t))

a.P_prev3["Normalhpp", "Normalhpn", 16] <- p.p3NpNn
a.P_prev3["Gastritishpp", "Gastritishpn", 16] <- p.p3GpGn
a.P_prev3["Atrophyhpp", "Atrophyhpn", 16] <- p.p3ApAn
a.P_prev3["Intestinalhpp", "Intestinalhpn", 16] <- p.p3IpIn
a.P_prev3["Dysplasiahpp", "Dysplasiahpn", 16] <- p.p3DpDn
a.P_prev3["Normalhpp", "Normalhpn", 17:n.t] <- p.NpNn
a.P_prev3["Gastritishpp", "Gastritishpn", 17:n.t] <- p.GpGn
a.P_prev3["Atrophyhpp", "Atrophyhpn", 17:n.t] <- p.ApAn
a.P_prev3["Intestinalhpp", "Intestinalhpn", 17:n.t] <- p.IpIn
a.P_prev3["Dysplasiahpp", "Dysplasiahpn", 17:n.t] <- p.DpDn
a.P_prev3["Normalhpp", "Normalhpn", 1:15] <- p.NpNn
a.P_prev3["Gastritishpp", "Gastritishpn", 1:15] <- p.GpGn
a.P_prev3["Atrophyhpp", "Atrophyhpn", 1:15] <- p.ApAn
a.P_prev3["Intestinalhpp", "Intestinalhpn", 1:15] <- p.IpIn
a.P_prev3["Dysplasiahpp", "Dysplasiahpn", 1:15] <- p.DpDn

### From Normal helicobacter (-)
a.P_prev3["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_prev3["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_prev3["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_prev3["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
a.P_prev3["Normalhpp", "Normalhpp", 1:15]    <- 1 - p.NpGp - p.NpNn - p.NpD3044 
a.P_prev3["Normalhpp", "Normalhpp", 16 ]    <- 1 - p.NpGp - p.p3NpNn - p.NpD45
a.P_prev3["Normalhpp", "Normalhpp", 17:n.t]    <- 1 - p.NpGp - p.NpNn - p.NpD46 
a.P_prev3["Normalhpp", "Gastritishpp", ]    <- p.NpGp
##a.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn --> Probabilidad que cambia, está arriba 
a.P_prev3["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_prev3["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_prev3["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_prev3["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_prev3["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_prev3["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_prev3["Gastritishpp", "Gastritishpp", 1:15 ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD3044
a.P_prev3["Gastritishpp", "Gastritishpp", 16 ]    <- 1 - p.GpAp - p.p3GpGn - p.GpNp - p.GpD45
a.P_prev3["Gastritishpp", "Gastritishpp", 17:n.t ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD46
a.P_prev3["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
##a.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
a.P_prev3["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_prev3["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_prev3["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn -  p.AnD
a.P_prev3["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_prev3["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_prev3["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_prev3["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (+) 
a.P_prev3["Atrophyhpp", "Atrophyhpp", 1:15]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD3044
a.P_prev3["Atrophyhpp", "Atrophyhpp", 16]    <- 1 - p.ApIp - p.p3ApAn - p.ApGp - p.ApD45
a.P_prev3["Atrophyhpp", "Atrophyhpp", 17:n.t]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD46
a.P_prev3["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
##a.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
a.P_prev3["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_prev3["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_prev3["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD 
a.P_prev3["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_prev3["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_prev3["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_prev3["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_prev3["Intestinalhpp", "Intestinalhpp", 1:15]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD3044  
a.P_prev3["Intestinalhpp", "Intestinalhpp", 16]    <- 1 - p.IpDp - p.p3IpIn - p.IpAp - p.IpD45
a.P_prev3["Intestinalhpp", "Intestinalhpp", 17:n.t]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD46  
a.P_prev3["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
#a.P_preV["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
a.P_prev3["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_prev3["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_prev3["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn -  p.DnD - p.Dn1p
a.P_prev3["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_prev3["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_prev3["Dysplasiahpn", "1preclinical", ]    <- p.Dn1p
a.P_prev3["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_prev3["Dysplasiahpp", "Dysplasiahpp", 1:15]    <- 1 - p.DpDn - p.DpIp - p.DpD3044 - p.Dp1p
a.P_prev3["Dysplasiahpp", "Dysplasiahpp", 16]    <- 1 - p.p3DpDn - p.DpIp - p.DpD45 - p.Dp1p
a.P_prev3["Dysplasiahpp", "Dysplasiahpp", 17:n.t]    <- 1 - p.DpDn - p.DpIp - p.DpD46 - p.Dp1p
##a.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
a.P_prev3["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_prev3["Dysplasiahpp", "1preclinical", ]    <- p.Dp1p
a.P_prev3["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_prev3["1preclinical", "1preclinical", 1:15]    <- 1 - p.1p1ca - p.1p2p - p.1pD3044	
a.P_prev3["1preclinical", "1preclinical", 16]    <- 1 - p.p3_1p1ca - p.1p2p - p.1pD45
a.P_prev3["1preclinical", "1preclinical", 17:n.t]    <- 1 - p.1p1ca - p.1p2p - p.1pD46
a.P_prev3["1preclinical", "1clinicala", 1:15 ]    <- p.1p1ca	
a.P_prev3["1preclinical", "1clinicala", 16 ]    <- p.p3_1p1ca
a.P_prev3["1preclinical", "1clinicala", 17:n.t ]    <- p.1p1ca	
a.P_prev3["1preclinical", "2preclinical", ]    <- p.1p2p	
a.P_prev3["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_prev3["2preclinical", "2preclinical", 1:15]    <- 1 - p.2p2ca - p.2p3p - p.2pD3044	
a.P_prev3["2preclinical", "2preclinical", 16]    <- 1 - p.p3_2p2ca - p.2p3p - p.2pD45
a.P_prev3["2preclinical", "2preclinical", 17:n.t]    <- 1 - p.2p2ca - p.2p3p - p.2pD46
a.P_prev3["2preclinical", "2clinicala", 1:15]    <- p.2p2ca	
a.P_prev3["2preclinical", "2clinicala", 16]    <- p.p3_2p2ca	
a.P_prev3["2preclinical", "2clinicala", 17:n.t]    <- p.2p2ca	
a.P_prev3["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_prev3["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_prev3["3preclinical", "3preclinical", 1:15]    <- 1 - p.3p3ca - p.3p4p - p.3pD3044	
a.P_prev3["3preclinical", "3preclinical", 16]    <- 1 - p.p3_3p3ca - p.3p4p - p.3pD45
a.P_prev3["3preclinical", "3preclinical", 17:n.t]    <- 1 - p.3p3ca - p.3p4p - p.3pD46
a.P_prev3["3preclinical", "3clinicala", 1:15 ]    <- p.3p3ca	
a.P_prev3["3preclinical", "3clinicala", 16]    <- p.p3_3p3ca	
a.P_prev3["3preclinical", "3clinicala", 17:n.t]    <- p.3p3ca	
a.P_prev3["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_prev3["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_prev3["4preclinical", "4preclinical", 1:15]    <- 1 - p.4p4ca - p.4pD3044	
a.P_prev3["4preclinical", "4preclinical", 16]    <- 1 - p.p3_4p4ca - p.4pD45	
a.P_prev3["4preclinical", "4preclinical", 17:n.t]    <- 1 - p.4p4ca - p.4pD46	
a.P_prev3["4preclinical", "4clinicala", 1:15]    <- p.4p4ca	
a.P_prev3["4preclinical", "4clinicala", 16]    <- p.p3_4p4ca	
a.P_prev3["4preclinical", "4clinicala", 17:n.t]    <- p.4p4ca	
a.P_prev3["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_prev3["1clinicala", "1clinicala", ]    <- 0	
a.P_prev3["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_prev3["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_prev3["2clinicala", "2clinicala", ]    <- 0	
a.P_prev3["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_prev3["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_prev3["3clinicala", "3clinicala", ]    <- 0 	
a.P_prev3["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_prev3["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_prev3["4clinicala", "4clinicala", ]    <- 0	
a.P_prev3["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_prev3["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_prev3["1clinicalb", "1clinicalb", ]    <- 0	
a.P_prev3["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_prev3["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_prev3["2clinicalb", "2clinicalb", ]    <- 0	
a.P_prev3["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_prev3["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_prev3["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_prev3["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_prev3["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_prev3["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_prev3["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_prev3["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_prev3["1clinicalc", "1clinicalc", ]    <- 0
a.P_prev3["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_prev3["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_prev3["1clinicald", "1clinicald", ]    <- 0
a.P_prev3["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_prev3["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_prev3["1clinicale", "1clinicale", ]    <- 0
a.P_prev3["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_prev3["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_prev3["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_prev3["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_prev3["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_prev3["2clinicalc", "2clinicalc", ]    <- 0
a.P_prev3["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_prev3["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_prev3["2clinicald", "2clinicald", ]    <- 0
a.P_prev3["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_prev3["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_prev3["2clinicale", "2clinicale", ]    <- 0
a.P_prev3["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_prev3["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_prev3["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_prev3["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_prev3["2clinicalf", "Dead", ]    <- p.2cfD	

### From 3 clinical c	
a.P_prev3["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_prev3["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_prev3["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_prev3["3clinicald", "3clinicald", ]    <- 0
a.P_prev3["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_prev3["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_prev3["3clinicale", "3clinicale", ]    <- 0
a.P_prev3["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_prev3["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_prev3["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_prev3["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_prev3["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_prev3["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_prev3["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_prev3["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_prev3["4clinicald", "4clinicald", ]    <- 0
a.P_prev3["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_prev3["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_prev3["4clinicale", "4clinicale", ]    <- 0
a.P_prev3["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_prev3["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_prev3["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_prev3["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_prev3["4clinicalf", "Dead", ]    <- p.4cfD	

### From Dead
a.P_prev3["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_prev3)
###


### From Normal helicobacter (-)
#m.P_prev3["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p3NnNp - p.p3NnGn
#m.P_prev3["Normalhpn", "Normalhpp"]    <- p.p3NnNp
#m.P_prev3["Normalhpn", "Gastritishpn"]    <- p.p3NnGn
#m.P_prev3["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
#m.P_prev3["Normalhpp", "Normalhpp"]    <- 1 - p.p3NpGp - p.p3NpNn - p.NpD 
#m.P_prev3["Normalhpp", "Gastritishpp"]    <- p.p3NpGp
#m.P_prev3["Normalhpp", "Normalhpn"]    <- p.p3NpNn
#m.P_prev3["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
#m.P_prev3["Gastritishpn", "Gastritishpn"]    <- 1 - p.p3GnAn - p.p3GnGp - p.p3GnNn - p.GnD
#m.P_prev3["Gastritishpn", "Atrophyhpn"]    <- p.p3GnAn
#m.P_prev3["Gastritishpn", "Gastritishpp"]    <- p.p3GnGp
#m.P_prev3["Gastritishpn", "Normalhpn"]    <- p.p3GnNn
#m.P_prev3["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
#m.P_prev3["Gastritishpp", "Gastritishpp"]    <- 1 - p.p3GpAp - p.p3GpGn - p.p3GpNp - p.GpD
#m.P_prev3["Gastritishpp", "Atrophyhpp"]    <- p.p3GpAp
#m.P_prev3["Gastritishpp", "Gastritishpn"]    <- p.p3GpGn
#m.P_prev3["Gastritishpp", "Normalhpp"]    <- p.p3GpNp
#m.P_prev3["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
#m.P_prev3["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p3AnIn - p.p3AnAp - p.p3AnGn -  p.AnD
#m.P_prev3["Atrophyhpn", "Intestinalhpn"]    <- p.p3AnIn
#m.P_prev3["Atrophyhpn", "Atrophyhpp"]    <- p.p3AnAp
#m.P_prev3["Atrophyhpn", "Gastritishpn"]    <- p.p3AnGn
#m.P_prev3["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
#m.P_prev3["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p3ApIp - p.p3ApAn - p.p3ApGp - p.ApD
#m.P_prev3["Atrophyhpp", "Intestinalhpp"]    <- p.p3ApIp
#m.P_prev3["Atrophyhpp", "Atrophyhpn"]    <- p.p3ApAn
#m.P_prev3["Atrophyhpp", "Gastritishpp"]    <- p.p3ApGp
#m.P_prev3["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
#m.P_prev3["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p3InDn - p.p3InIp - p.p3InAn - p.InD 
#m.P_prev3["Intestinalhpn", "Dysplasiahpn"]    <- p.p3InDn
#m.P_prev3["Intestinalhpn", "Intestinalhpp"]    <- p.p3InIp
#m.P_prev3["Intestinalhpn", "Atrophyhpn"]    <- p.p3InAn
#m.P_prev3["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
#m.P_prev3["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p3IpDp - p.p3IpIn - p.p3IpAp - p.IpD  
#m.P_prev3["Intestinalhpp", "Dysplasiahpp"]    <- p.p3IpDp
#m.P_prev3["Intestinalhpp", "Intestinalhpn"]    <- p.p3IpIn
#m.P_prev3["Intestinalhpp", "Atrophyhpp"]    <- p.p3IpAp
#m.P_prev3["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
#m.P_prev3["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p3DnDp - p.p3DnIn -  p.DnD -  p.p3_Dn1p
#m.P_prev3["Dysplasiahpn", "Dysplasiahpp"]    <- p.p3DnDp
#m.P_prev3["Dysplasiahpn", "Intestinalhpn"]    <- p.p3DnIn
#m.P_prev3["Dysplasiahpn", "1preclinical"]    <- p.p3_Dn1p
#m.P_prev3["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
#m.P_prev3["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p3DpDn - p.p3DpIp - p.DpD - p.p3_Dp1p
#m.P_prev3["Dysplasiahpp", "Dysplasiahpn"]    <- p.p3DpDn
#m.P_prev3["Dysplasiahpp", "Intestinalhpp"]    <- p.p3DpIp
#m.P_prev3["Dysplasiahpp", "1preclinical"]    <- p.p3_Dp1p
#m.P_prev3["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
#m.P_prev3["1preclinical", "1preclinical"]    <- 1 - p.p3_1p1ca - p.p3_1p2p - p.1pD	
#m.P_prev3["1preclinical", "1clinicala"]    <- p.p3_1p1ca	
#m.P_prev3["1preclinical", "2preclinical"]    <- p.p3_1p2p	
#m.P_prev3["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
#m.P_prev3["2preclinical", "2preclinical"]    <- 1 - p.p3_2p2ca - p.p3_2p3p - p.2pD	
#m.P_prev3["2preclinical", "2clinicala"]    <- p.p3_2p2ca	
#m.P_prev3["2preclinical", "3preclinical"]    <- p.p3_2p3p	
#m.P_prev3["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
#m.P_prev3["3preclinical", "3preclinical"]    <- 1 - p.p3_3p3ca - p.p3_3p4p - p.3pD	
#m.P_prev3["3preclinical", "3clinicala"]    <- p.p3_3p3ca	
#m.P_prev3["3preclinical", "4preclinical"]    <- p.p3_3p4p	
#m.P_prev3["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
#m.P_prev3["4preclinical", "4preclinical"]    <- 1 - p.p3_4p4ca - p.4pD	
#m.P_prev3["4preclinical", "4clinicala"]    <- p.p3_4p4ca	
#m.P_prev3["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
#m.P_prev3["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p3_1ca1cb	
#m.P_prev3["1clinicala", "1clinicalb"] <- p.p3_1ca1cb	
#m.P_prev3["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
#m.P_prev3["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p3_2ca2cb	
#m.P_prev3["2clinicala", "2clinicalb"] <- p.p3_2ca2cb	
#m.P_prev3["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
#m.P_prev3["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p3_3ca3cb 	
#m.P_prev3["3clinicala", "3clinicalb"] <- p.p3_3ca3cb	
#m.P_prev3["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
#m.P_prev3["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p3_4ca4cb	
#m.P_prev3["4clinicala", "4clinicalb"] <- p.p3_4ca4cb  	
#m.P_prev3["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
#m.P_prev3["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p3_1cb1cc	
#m.P_prev3["1clinicalb", "1clinicalc"] <- p.p3_1cb1cc	
#m.P_prev3["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
#m.P_prev3["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p3_2cb2cc	
#m.P_prev3["2clinicalb", "2clinicalc"] <- p.p3_2cb2cc	
#m.P_prev3["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
#m.P_prev3["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p3_3cb3cc 	
#m.P_prev3["3clinicalb", "3clinicalc"] <- p.p3_3cb3cc	
#m.P_prev3["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
#m.P_prev3["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p3_4cb4cc 	
#m.P_prev3["4clinicalb", "4clinicalc"] <- p.p3_4cb4cc	
#m.P_prev3["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
#m.P_prev3["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p3_1ccNn	
#m.P_prev3["1clinicalc", "Normalhpn"] <- p.p3_1ccNn	
#m.P_prev3["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
#m.P_prev3["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p3_2ccNn	
#m.P_prev3["2clinicalc", "Normalhpn"] <- p.p3_2ccNn	
#m.P_prev3["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
#m.P_prev3["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p3_3ccNn 	
#m.P_prev3["3clinicalc", "Normalhpn"] <- p.p3_3ccNn
#m.P_prev3["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
#m.P_prev3["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p3_4ccNn 	
#m.P_prev3["4clinicalc", "Normalhpn"] <- p.p3_4ccNn	
#m.P_prev3["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
#m.P_prev3["Dead", "Dead"] <- 1

# check rows add up to 1
#rowSums(m.P_prev3)


# fill in the transition probability array with prevention strategy 4

a.P_prev4 <- array(0, 
                    dim = c(n.s, n.s, n.t), 
                    dimnames = list( from  = v.n, 
                                     to    = v.n, 
                                     cycle = 1:n.t))

# fill in the transition probability array without prevention
### From Normal helicobacter (-)
a.P_prev4["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_prev4["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_prev4["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_prev4["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
a.P_prev4["Normalhpp", "Normalhpp", 1:10 ]    <- 1 - p.NpGp - p.NpNn - p.NpD3039
a.P_prev4["Normalhpp", "Normalhpp", 11]    <- 1 - p.NpGp - p.p4NpNn - p.NpD40 
a.P_prev4["Normalhpp", "Normalhpp", 12:n.t]    <- 1 - p.NpGp - p.NpNn - p.NpD41 
a.P_prev4["Normalhpp", "Normalhpn", 1:10]    <- p.NpNn
a.P_prev4["Normalhpp", "Normalhpn", 11]    <- p.p4NpNn
a.P_prev4["Normalhpp", "Normalhpn", 12:n.t]    <- p.NpNn
a.P_prev4["Normalhpp", "Gastritishpp", ]    <- p.NpGp
a.P_prev4["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_prev4["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_prev4["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_prev4["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_prev4["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_prev4["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_prev4["Gastritishpp", "Gastritishpp", 1:10]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD3039
a.P_prev4["Gastritishpp", "Gastritishpp", 11]    <- 1 - p.GpAp - p.p4GpGn - p.GpNp - p.GpD40
a.P_prev4["Gastritishpp", "Gastritishpp", 12:n.t]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD41
a.P_prev4["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
a.P_prev4["Gastritishpp", "Gastritishpn", 1:10]    <- p.GpGn
a.P_prev4["Gastritishpp", "Gastritishpn", 11]    <- p.p4GpGn
a.P_prev4["Gastritishpp", "Gastritishpn", 12:n.t]    <- p.GpGn
a.P_prev4["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_prev4["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_prev4["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn - p.AnD
a.P_prev4["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_prev4["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_prev4["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_prev4["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (+) 
a.P_prev4["Atrophyhpp", "Atrophyhpp", 1:10]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD3039
a.P_prev4["Atrophyhpp", "Atrophyhpp", 11]    <- 1 - p.ApIp - p.p4ApAn - p.ApGp - p.ApD40
a.P_prev4["Atrophyhpp", "Atrophyhpp", 12:n.t]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD41
a.P_prev4["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
a.P_prev4["Atrophyhpp", "Atrophyhpn", 1:10]    <- p.ApAn
a.P_prev4["Atrophyhpp", "Atrophyhpn", 11]    <- p.p4ApAn
a.P_prev4["Atrophyhpp", "Atrophyhpn", 12:n.t]    <- p.ApAn
a.P_prev4["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_prev4["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_prev4["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD
a.P_prev4["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_prev4["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_prev4["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_prev4["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_prev4["Intestinalhpp", "Intestinalhpp", 1:10]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD3039  
a.P_prev4["Intestinalhpp", "Intestinalhpp", 11]    <- 1 - p.IpDp - p.p4IpIn - p.IpAp - p.IpD40 
a.P_prev4["Intestinalhpp", "Intestinalhpp", 12:n.t]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD41 
a.P_prev4["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
a.P_prev4["Intestinalhpp", "Intestinalhpn", 1:10]    <- p.IpIn
a.P_prev4["Intestinalhpp", "Intestinalhpn", 11]    <- p.p4IpIn
a.P_prev4["Intestinalhpp", "Intestinalhpn", 12:n.t]    <- p.IpIn
a.P_prev4["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_prev4["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_prev4["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn - p.Dn1p - p.DnD
a.P_prev4["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_prev4["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_prev4["Dysplasiahpn", "1preclinical", ]     <- p.Dn1p
a.P_prev4["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_prev4["Dysplasiahpp", "Dysplasiahpp", 1:10]    <- 1 - p.DpDn - p.DpIp - p.Dp1p - p.DpD3039
a.P_prev4["Dysplasiahpp", "Dysplasiahpp", 11]    <- 1 - p.p4DpDn - p.DpIp - p.Dp1p - p.DpD40
a.P_prev4["Dysplasiahpp", "Dysplasiahpp", 12:n.t]    <- 1 - p.DpDn - p.DpIp - p.Dp1p - p.DpD41
a.P_prev4["Dysplasiahpp", "Dysplasiahpn", 1:10]    <- p.DpDn
a.P_prev4["Dysplasiahpp", "Dysplasiahpn", 11]    <- p.p4DpDn
a.P_prev4["Dysplasiahpp", "Dysplasiahpn", 12:n.t]    <- p.DpDn
a.P_prev4["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_prev4["Dysplasiahpp", "1preclinical", ]     <- p.Dp1p
a.P_prev4["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_prev4["1preclinical", "1preclinical", 1:10]    <- 1 - p.1p1ca - p.1p2p - p.1pD3039	
a.P_prev4["1preclinical", "1preclinical", 11]    <- 1 - p.p4_1p1ca - p.1p2p - p.1pD40
a.P_prev4["1preclinical", "1preclinical", 12:n.t]    <- 1 - p.1p1ca - p.1p2p - p.1pD41
a.P_prev4["1preclinical", "1clinicala", 1:10 ]    <- p.1p1ca	
a.P_prev4["1preclinical", "1clinicala", 11 ]    <- p.p4_1p1ca
a.P_prev4["1preclinical", "1clinicala", 12:n.t ]    <- p.1p1ca	
a.P_prev4["1preclinical", "2preclinical", ]    <- p.1p2p	
a.P_prev4["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_prev4["2preclinical", "2preclinical", 1:10]    <- 1 - p.2p2ca - p.2p3p - p.2pD3039	
a.P_prev4["2preclinical", "2preclinical", 11]    <- 1 - p.p4_2p2ca - p.2p3p - p.2pD40
a.P_prev4["2preclinical", "2preclinical", 12:n.t]    <- 1 - p.2p2ca - p.2p3p - p.2pD41
a.P_prev4["2preclinical", "2clinicala", 1:10]    <- p.2p2ca	
a.P_prev4["2preclinical", "2clinicala", 11]    <- p.p4_2p2ca	
a.P_prev4["2preclinical", "2clinicala", 12:n.t]    <- p.2p2ca	
a.P_prev4["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_prev4["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_prev4["3preclinical", "3preclinical", 1:10]    <- 1 - p.3p3ca - p.3p4p - p.3pD3039	
a.P_prev4["3preclinical", "3preclinical", 11]    <- 1 - p.p4_3p3ca - p.3p4p - p.3pD40
a.P_prev4["3preclinical", "3preclinical", 12:n.t]    <- 1 - p.3p3ca - p.3p4p - p.3pD41
a.P_prev4["3preclinical", "3clinicala", 1:10]    <- p.3p3ca	
a.P_prev4["3preclinical", "3clinicala", 11]    <- p.p4_3p3ca	
a.P_prev4["3preclinical", "3clinicala", 12:n.t]    <- p.3p3ca	
a.P_prev4["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_prev4["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_prev4["4preclinical", "4preclinical", ]    <- 1 - p.4p4ca - p.4pD	
a.P_prev4["4preclinical", "4clinicala", ]    <- p.4p4ca	
a.P_prev4["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_prev4["1clinicala", "1clinicala", ]    <- 0	
a.P_prev4["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_prev4["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_prev4["2clinicala", "2clinicala", ]    <- 0	
a.P_prev4["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_prev4["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_prev4["3clinicala", "3clinicala", ]    <- 0 	
a.P_prev4["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_prev4["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_prev4["4clinicala", "4clinicala", ]    <- 0	
a.P_prev4["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_prev4["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_prev4["1clinicalb", "1clinicalb", ]    <- 0	
a.P_prev4["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_prev4["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_prev4["2clinicalb", "2clinicalb", ]    <- 0	
a.P_prev4["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_prev4["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_prev4["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_prev4["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_prev4["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_prev4["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_prev4["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_prev4["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_prev4["1clinicalc", "1clinicalc", ]    <- 0
a.P_prev4["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_prev4["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_prev4["1clinicald", "1clinicald", ]    <- 0
a.P_prev4["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_prev4["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_prev4["1clinicale", "1clinicale", ]    <- 0
a.P_prev4["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_prev4["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_prev4["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_prev4["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_prev4["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_prev4["2clinicalc", "2clinicalc", ]    <- 0
a.P_prev4["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_prev4["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_prev4["2clinicald", "2clinicald", ]    <- 0
a.P_prev4["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_prev4["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_prev4["2clinicale", "2clinicale", ]    <- 0
a.P_prev4["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_prev4["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_prev4["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_prev4["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_prev4["2clinicalf", "Dead", ]    <- p.2cfD	

### From 3 clinical c	
a.P_prev4["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_prev4["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_prev4["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_prev4["3clinicald", "3clinicald", ]    <- 0
a.P_prev4["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_prev4["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_prev4["3clinicale", "3clinicale", ]    <- 0
a.P_prev4["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_prev4["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_prev4["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_prev4["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_prev4["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_prev4["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_prev4["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_prev4["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_prev4["4clinicald", "4clinicald", ]    <- 0
a.P_prev4["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_prev4["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_prev4["4clinicale", "4clinicale", ]    <- 0
a.P_prev4["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_prev4["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_prev4["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_prev4["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_prev4["4clinicalf", "Dead", ]    <- p.4cfD

### From Dead
a.P_prev4["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_prev4)
###


### From Normal helicobacter (-)
#m.P_prev4["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p4NnNp - p.p4NnGn
#m.P_prev4["Normalhpn", "Normalhpp"]    <- p.p4NnNp
#m.P_prev4["Normalhpn", "Gastritishpn"]    <- p.p4NnGn
#m.P_prev4["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
#m.P_prev4["Normalhpp", "Normalhpp"]    <- 1 - p.p4NpGp - p.p4NpNn - p.NpD 
#m.P_prev4["Normalhpp", "Gastritishpp"]    <- p.p4NpGp
#m.P_prev4["Normalhpp", "Normalhpn"]    <- p.p4NpNn
#m.P_prev4["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
#m.P_prev4["Gastritishpn", "Gastritishpn"]    <- 1 - p.p4GnAn - p.p4GnGp - p.p4GnNn - p.GnD
#m.P_prev4["Gastritishpn", "Atrophyhpn"]    <- p.p4GnAn
#m.P_prev4["Gastritishpn", "Gastritishpp"]    <- p.p4GnGp
#m.P_prev4["Gastritishpn", "Normalhpn"]    <- p.p4GnNn
#m.P_prev4["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
#m.P_prev4["Gastritishpp", "Gastritishpp"]    <- 1 - p.p4GpAp - p.p4GpGn - p.p4GpNp - p.GpD
#m.P_prev4["Gastritishpp", "Atrophyhpp"]    <- p.p4GpAp
#m.P_prev4["Gastritishpp", "Gastritishpn"]    <- p.p4GpGn
#m.P_prev4["Gastritishpp", "Normalhpp"]    <- p.p4GpNp
#m.P_prev4["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
#m.P_prev4["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p4AnIn - p.p4AnAp - p.p4AnGn -  p.AnD
#m.P_prev4["Atrophyhpn", "Intestinalhpn"]    <- p.p4AnIn
#m.P_prev4["Atrophyhpn", "Atrophyhpp"]    <- p.p4AnAp
#m.P_prev4["Atrophyhpn", "Gastritishpn"]    <- p.p4AnGn
#m.P_prev4["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
#m.P_prev4["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p4ApIp - p.p4ApAn - p.p4ApGp - p.ApD
#m.P_prev4["Atrophyhpp", "Intestinalhpp"]    <- p.p4ApIp
#m.P_prev4["Atrophyhpp", "Atrophyhpn"]    <- p.p4ApAn
#m.P_prev4["Atrophyhpp", "Gastritishpp"]    <- p.p4ApGp
#m.P_prev4["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
#m.P_prev4["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p4InDn - p.p4InIp - p.p4InAn - p.InD 
#m.P_prev4["Intestinalhpn", "Dysplasiahpn"]    <- p.p4InDn
#m.P_prev4["Intestinalhpn", "Intestinalhpp"]    <- p.p4InIp
#m.P_prev4["Intestinalhpn", "Atrophyhpn"]    <- p.p4InAn
#m.P_prev4["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
#m.P_prev4["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p4IpDp - p.p4IpIn - p.p4IpAp - p.IpD  
#m.P_prev4["Intestinalhpp", "Dysplasiahpp"]    <- p.p4IpDp
#m.P_prev4["Intestinalhpp", "Intestinalhpn"]    <- p.p4IpIn
#m.P_prev4["Intestinalhpp", "Atrophyhpp"]    <- p.p4IpAp
#m.P_prev4["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
#m.P_prev4["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p4DnDp - p.p4DnIn -  p.DnD - p.p4_Dn1p
#m.P_prev4["Dysplasiahpn", "Dysplasiahpp"]    <- p.p4DnDp
#m.P_prev4["Dysplasiahpn", "Intestinalhpn"]    <- p.p4DnIn
#m.P_prev4["Dysplasiahpn", "1preclinical"]    <- p.p4_Dn1p
#m.P_prev4["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
#m.P_prev4["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p4DpDn - p.p4DpIp - p.DpD - p.p4_Dp1p
#m.P_prev4["Dysplasiahpp", "Dysplasiahpn"]    <- p.p4DpDn
#m.P_prev4["Dysplasiahpp", "Intestinalhpp"]    <- p.p4DpIp
#m.P_prev4["Dysplasiahpp", "1preclinical"]    <- p.p4_Dp1p
#m.P_prev4["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
#m.P_prev4["1preclinical", "1preclinical"]    <- 1 - p.p4_1p1ca - p.p4_1p2p - p.1pD	
#m.P_prev4["1preclinical", "1clinicala"]    <- p.p4_1p1ca	
#m.P_prev4["1preclinical", "2preclinical"]    <- p.p4_1p2p	
#m.P_prev4["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
#m.P_prev4["2preclinical", "2preclinical"]    <- 1 - p.p4_2p2ca - p.p4_2p3p - p.2pD	
#m.P_prev4["2preclinical", "2clinicala"]    <- p.p4_2p2ca	
#m.P_prev4["2preclinical", "3preclinical"]    <- p.p4_2p3p	
#m.P_prev4["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
#m.P_prev4["3preclinical", "3preclinical"]    <- 1 - p.p4_3p3ca - p.p4_3p4p - p.3pD	
#m.P_prev4["3preclinical", "3clinicala"]    <- p.p4_3p3ca	
#m.P_prev4["3preclinical", "4preclinical"]    <- p.p4_3p4p	
#m.P_prev4["3preclinical", "Dead"]    <- p.3pD	

### From 4 preclinical	
#m.P_prev4["4preclinical", "4preclinical"]    <- 1 - p.p4_4p4ca - p.4pD	
#m.P_prev4["4preclinical", "4clinicala"]    <- p.p4_4p4ca	
#m.P_prev4["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
#m.P_prev4["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p4_1ca1cb	
#m.P_prev4["1clinicala", "1clinicalb"] <- p.p4_1ca1cb	
#m.P_prev4["1clinicala", "Dead"]    <- p.1caD	

### From 2 clinical a	
#m.P_prev4["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p4_2ca2cb	
#m.P_prev4["2clinicala", "2clinicalb"] <- p.p4_2ca2cb	
#m.P_prev4["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
#m.P_prev4["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p4_3ca3cb 	
#m.P_prev4["3clinicala", "3clinicalb"] <- p.p4_3ca3cb	
#m.P_prev4["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
#m.P_prev4["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p4_4ca4cb	
#m.P_prev4["4clinicala", "4clinicalb"] <- p.p4_4ca4cb  	
#m.P_prev4["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
#m.P_prev4["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p4_1cb1cc	
#m.P_prev4["1clinicalb", "1clinicalc"] <- p.p4_1cb1cc	
#m.P_prev4["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
#m.P_prev4["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p4_2cb2cc	
#m.P_prev4["2clinicalb", "2clinicalc"] <- p.p4_2cb2cc	
#m.P_prev4["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
#m.P_prev4["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p4_3cb3cc 	
#m.P_prev4["3clinicalb", "3clinicalc"] <- p.p4_3cb3cc	
#m.P_prev4["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
#m.P_prev4["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p4_4cb4cc 	
#m.P_prev4["4clinicalb", "4clinicalc"] <- p.p4_4cb4cc	
#m.P_prev4["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
#m.P_prev4["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p4_1ccNn	
#m.P_prev4["1clinicalc", "Normalhpn"] <- p.p4_1ccNn	
#m.P_prev4["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
#m.P_prev4["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p4_2ccNn	
#m.P_prev4["2clinicalc", "Normalhpn"] <- p.p4_2ccNn	
#m.P_prev4["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
#m.P_prev4["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p4_3ccNn 	
#m.P_prev4["3clinicalc", "Normalhpn"] <- p.p4_3ccNn
#m.P_prev4["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
#m.P_prev4["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p4_4ccNn 	
#m.P_prev4["4clinicalc", "Normalhpn"] <- p.p4_4ccNn	
#m.P_prev4["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
#m.P_prev4["Dead", "Dead"] <- 1

# check rows add up to 1
#rowSums(m.P_prev4)

# fill in the transition probability matrix with prevention strategy 5


a.P_prev5 <- array(0, 
                   dim = c(n.s, n.s, n.t), 
                   dimnames = list( from  = v.n, 
                                    to    = v.n, 
                                    cycle = 1:n.t))

a.P_prev5["Normalhpp", "Normalhpn", 1] <- p.p5NpNn
a.P_prev5["Gastritishpp", "Gastritishpn", 1] <- p.p5GpGn
a.P_prev5["Atrophyhpp", "Atrophyhpn", 1] <- p.p5ApAn
a.P_prev5["Intestinalhpp", "Intestinalhpn", 1] <- p.p5IpIn
a.P_prev5["Dysplasiahpp", "Dysplasiahpn", 1] <- p.p5DpDn
a.P_prev5["Normalhpp", "Normalhpn", 2:n.t] <- p.NpNn
a.P_prev5["Gastritishpp", "Gastritishpn", 2:n.t] <- p.GpGn
a.P_prev5["Atrophyhpp", "Atrophyhpn", 2:n.t] <- p.ApAn
a.P_prev5["Intestinalhpp", "Intestinalhpn", 2:n.t] <- p.IpIn
a.P_prev5["Dysplasiahpp", "Dysplasiahpn", 2:n.t] <- p.DpDn

### From Normal helicobacter (-)
a.P_prev5["Normalhpn", "Normalhpn", ] <- 1 - p.NnD - p.NnNp - p.NnGn
a.P_prev5["Normalhpn", "Normalhpp", ]    <- p.NnNp
a.P_prev5["Normalhpn", "Gastritishpn", ]    <- p.NnGn
a.P_prev5["Normalhpn", "Dead", ]    <- p.NnD

### From Normal helicobacter (+) 
a.P_prev5["Normalhpp", "Normalhpp", 1 ]    <- 1 - p.NpGp - p.p5NpNn - p.NpD30
a.P_prev5["Normalhpp", "Normalhpp", 2:n.t]    <- 1 - p.NpGp - p.NpNn - p.NpD31 
a.P_prev5["Normalhpp", "Gastritishpp", ]    <- p.NpGp
##a.P_prev["Normalhpp", "Normalhpn"]    <- p.pNpNn --> Probabilidad que cambia, está arriba 
a.P_prev5["Normalhpp", "Dead", ]    <- p.NpD

### From Gastritis Helicobacter (-)
a.P_prev5["Gastritishpn", "Gastritishpn", ]    <- 1 - p.GnAn - p.GnGp - p.GnNn - p.GnD
a.P_prev5["Gastritishpn", "Atrophyhpn", ]    <- p.GnAn
a.P_prev5["Gastritishpn", "Gastritishpp", ]    <- p.GnGp
a.P_prev5["Gastritishpn", "Normalhpn", ]    <- p.GnNn
a.P_prev5["Gastritishpn", "Dead", ]    <- p.GnD

### From Gastritis Helicobacter (+)
a.P_prev5["Gastritishpp", "Gastritishpp", 1 ]    <- 1 - p.GpAp - p.p5GpGn - p.GpNp - p.GpD30
a.P_prev5["Gastritishpp", "Gastritishpp", 2:n.t ]    <- 1 - p.GpAp - p.GpGn - p.GpNp - p.GpD31
a.P_prev5["Gastritishpp", "Atrophyhpp", ]    <- p.GpAp
##a.P_prev["Gastritishpp", "Gastritishpn"]    <- p.pGpGn
a.P_prev5["Gastritishpp", "Normalhpp", ]    <- p.GpNp
a.P_prev5["Gastritishpp", "Dead", ]    <- p.GpD

### From Atrophy Helicobacter (-) 
a.P_prev5["Atrophyhpn", "Atrophyhpn", ]    <- 1 - p.AnIn - p.AnAp - p.AnGn -  p.AnD
a.P_prev5["Atrophyhpn", "Intestinalhpn", ]    <- p.AnIn
a.P_prev5["Atrophyhpn", "Atrophyhpp", ]    <- p.AnAp
a.P_prev5["Atrophyhpn", "Gastritishpn", ]    <- p.AnGn
a.P_prev5["Atrophyhpn", "Dead", ]    <- p.AnD

### From Atrophy Helicobacter (+) 
a.P_prev5["Atrophyhpp", "Atrophyhpp", 1]    <- 1 - p.ApIp - p.p5ApAn - p.ApGp - p.ApD30
a.P_prev5["Atrophyhpp", "Atrophyhpp", 2:n.t]    <- 1 - p.ApIp - p.ApAn - p.ApGp - p.ApD31
a.P_prev5["Atrophyhpp", "Intestinalhpp", ]    <- p.ApIp
##a.P_prev["Atrophyhpp", "Atrophyhpn"]    <- p.pApAn
a.P_prev5["Atrophyhpp", "Gastritishpp", ]    <- p.ApGp
a.P_prev5["Atrophyhpp", "Dead", ]    <- p.ApD

### From Intestinal Helicobacter (-)
a.P_prev5["Intestinalhpn", "Intestinalhpn", ]    <- 1 - p.InDn - p.InIp - p.InAn - p.InD 
a.P_prev5["Intestinalhpn", "Dysplasiahpn", ]    <- p.InDn
a.P_prev5["Intestinalhpn", "Intestinalhpp", ]    <- p.InIp
a.P_prev5["Intestinalhpn", "Atrophyhpn", ]    <- p.InAn
a.P_prev5["Intestinalhpn", "Dead", ]    <- p.InD

### From Intestinal Helicobacter (+)
a.P_prev5["Intestinalhpp", "Intestinalhpp", 1]    <- 1 - p.IpDp - p.p5IpIn - p.IpAp - p.IpD30
a.P_prev5["Intestinalhpp", "Intestinalhpp", 2:n.t]    <- 1 - p.IpDp - p.IpIn - p.IpAp - p.IpD31  
a.P_prev5["Intestinalhpp", "Dysplasiahpp", ]    <- p.IpDp
#a.P_prev["Intestinalhpp", "Intestinalhpn"]    <- p.pIpIn
a.P_prev5["Intestinalhpp", "Atrophyhpp", ]    <- p.IpAp
a.P_prev5["Intestinalhpp", "Dead", ]    <- p.IpD

### From Dysplasia Helicobacter (-)
a.P_prev5["Dysplasiahpn", "Dysplasiahpn", ]    <- 1 - p.DnDp - p.DnIn -  p.DnD - p.Dn1p
a.P_prev5["Dysplasiahpn", "Dysplasiahpp", ]    <- p.DnDp
a.P_prev5["Dysplasiahpn", "Intestinalhpn", ]    <- p.DnIn
a.P_prev5["Dysplasiahpn", "1preclinical", ]    <- p.Dn1p
a.P_prev5["Dysplasiahpn", "Dead", ]    <- p.DnD

### From Dysplasia Helicobacter (+)
a.P_prev5["Dysplasiahpp", "Dysplasiahpp", 1]    <- 1 - p.p5DpDn - p.DpIp - p.DpD30 - p.Dp1p
a.P_prev5["Dysplasiahpp", "Dysplasiahpp", 2:n.t]    <- 1 - p.DpDn - p.DpIp - p.DpD31 - p.Dp1p
##a.P_prev["Dysplasiahpp", "Dysplasiahpn"]    <- p.pDpDn
a.P_prev5["Dysplasiahpp", "Intestinalhpp", ]    <- p.DpIp
a.P_prev5["Dysplasiahpp", "1preclinical", ]    <- p.Dp1p
a.P_prev5["Dysplasiahpp", "Dead", ]    <- p.DpD

### From 1 preclinical 	
a.P_prev5["1preclinical", "1preclinical",]    <- 1 - p.1p1ca - p.1p2p - p.1pD	
a.P_prev5["1preclinical", "1clinicala",]    <- p.1p1ca	
a.P_prev5["1preclinical", "2preclinical",]    <- p.1p2p	
a.P_prev5["1preclinical", "Dead", ]    <- p.1pD	

### From 2 preclinical	
a.P_prev5["2preclinical", "2preclinical",]    <- 1 - p.2p2ca - p.2p3p - p.2pD	
a.P_prev5["2preclinical", "2clinicala", ]    <- p.2p2ca	
a.P_prev5["2preclinical", "3preclinical", ]    <- p.2p3p	
a.P_prev5["2preclinical", "Dead", ]    <- p.2pD	

### From 3 preclinical	
a.P_prev5["3preclinical", "3preclinical", ]    <- 1 - p.3p3ca - p.3p4p - p.3pD	
a.P_prev5["3preclinical", "3clinicala", ]    <- p.3p3ca	
a.P_prev5["3preclinical", "4preclinical", ]    <- p.3p4p	
a.P_prev5["3preclinical", "Dead", ]    <- p.3pD	

### From 4 preclinical	
a.P_prev5["4preclinical", "4preclinical", ]    <- 1 - p.4p4ca - p.4pD	
a.P_prev5["4preclinical", "4clinicala", ]    <- p.4p4ca	
a.P_prev5["4preclinical", "Dead", ]    <- p.4pD	

### From 1 clinical a	
a.P_prev5["1clinicala", "1clinicala", ]    <- 0	
a.P_prev5["1clinicala", "1clinicalb", ] <- 1 - p.1caD	
a.P_prev5["1clinicala", "Dead", ]    <- p.1caD	

### From 2 clinical a	
a.P_prev5["2clinicala", "2clinicala", ]    <- 0	
a.P_prev5["2clinicala", "2clinicalb", ] <- 1 - p.2caD	
a.P_prev5["2clinicala", "Dead", ]    <- p.2caD	

### From 3 clinical a	
a.P_prev5["3clinicala", "3clinicala", ]    <- 0 	
a.P_prev5["3clinicala", "3clinicalb", ] <- 1 - p.3caD		
a.P_prev5["3clinicala", "Dead", ]    <- p.3caD	

### From 4 clinical a 	
a.P_prev5["4clinicala", "4clinicala", ]    <- 0	
a.P_prev5["4clinicala", "4clinicalb", ] <- 1 - p.4caD  	
a.P_prev5["4clinicala", "Dead", ]    <- p.4caD	

### From 1 clinical b	
a.P_prev5["1clinicalb", "1clinicalb", ]    <- 0	
a.P_prev5["1clinicalb", "1clinicalc", ] <- 1 - p.1cbD	
a.P_prev5["1clinicalb", "Dead", ]    <- p.1cbD	

### From 2 clinical b	
a.P_prev5["2clinicalb", "2clinicalb", ]    <- 0	
a.P_prev5["2clinicalb", "2clinicalc", ] <- 1 - p.2cbD	
a.P_prev5["2clinicalb", "Dead", ]    <- p.2cbD	

### From 3 clinical b	
a.P_prev5["3clinicalb", "3clinicalb", ]    <- 0 	
a.P_prev5["3clinicalb", "3clinicalc", ] <- 1 - p.3cbD	
a.P_prev5["3clinicalb", "Dead", ]    <- p.3cbD	

### From 4 clinical b	
a.P_prev5["4clinicalb", "4clinicalb", ]    <- 0 	
a.P_prev5["4clinicalb", "4clinicalc", ] <- 1 - p.4cbD	
a.P_prev5["4clinicalb", "Dead", ]    <- p.4cbD	

### From 1 clinical c	
a.P_prev5["1clinicalc", "1clinicalc", ]    <- 0
a.P_prev5["1clinicalc", "1clinicald", ] <- 1 - p.1ccD	
a.P_prev5["1clinicalc", "Dead", ]    <- p.1ccD	

### From 1 clinical d	
a.P_prev5["1clinicald", "1clinicald", ]    <- 0
a.P_prev5["1clinicald", "1clinicale", ] <- 1 - p.1cdD	
a.P_prev5["1clinicald", "Dead", ]    <- p.1cdD	

### From 1 clinical e	
a.P_prev5["1clinicale", "1clinicale", ]    <- 0
a.P_prev5["1clinicale", "1clinicalf", ] <- 1 - p.1ceD	
a.P_prev5["1clinicale", "Dead", ]    <- p.1ceD	

### From 1 clinical f	
a.P_prev5["1clinicalf", "1clinicalf", ]    <- 1 - p.1cfNn - p.1cfD
a.P_prev5["1clinicalf", "Normalhpp", ] <- p.1cfNn	
a.P_prev5["1clinicalf", "Dead", ]    <- p.1cfD	

### From 2 clinical c	
a.P_prev5["2clinicalc", "2clinicalc", ]    <- 0
a.P_prev5["2clinicalc", "2clinicald", ] <- 1 - p.2ccD	
a.P_prev5["2clinicalc", "Dead", ]    <- p.2ccD	

### From 2 clinical d	
a.P_prev5["2clinicald", "2clinicald", ]    <- 0
a.P_prev5["2clinicald", "2clinicale", ] <- 1 - p.2cdD	
a.P_prev5["2clinicald", "Dead", ]    <- p.2cdD	

### From 2 clinical e	
a.P_prev5["2clinicale", "2clinicale", ]    <- 0
a.P_prev5["2clinicale", "2clinicalf", ] <- 1 - p.2ceD	
a.P_prev5["2clinicale", "Dead", ]    <- p.2ceD	

### From 2 clinical f	
a.P_prev5["2clinicalf", "2clinicalf", ]    <- 1 - p.2cfNn - p.2cfD
a.P_prev5["2clinicalf", "Normalhpp", ] <- p.2cfNn	
a.P_prev5["2clinicalf", "Dead", ]    <- p.2cfD

### From 3 clinical c	
a.P_prev5["3clinicalc", "3clinicalc", ]    <- 0 	
a.P_prev5["3clinicalc", "3clinicald", ] <- 1 - p.3ccD
a.P_prev5["3clinicalc", "Dead", ]    <- p.3ccD	

### From 3 clinical d	
a.P_prev5["3clinicald", "3clinicald", ]    <- 0
a.P_prev5["3clinicald", "3clinicale", ] <- 1 - p.3cdD	
a.P_prev5["3clinicald", "Dead", ]    <- p.3cdD	

### From 3 clinical e	
a.P_prev5["3clinicale", "3clinicale", ]    <- 0
a.P_prev5["3clinicale", "3clinicalf", ] <- 1 - p.3ceD	
a.P_prev5["3clinicale", "Dead", ]    <- p.3ceD	

### From 3 clinical f	
a.P_prev5["3clinicalf", "3clinicalf", ]    <- 1 - p.3cfNn - p.3cfD
a.P_prev5["3clinicalf", "Normalhpp", ] <- p.3cfNn	
a.P_prev5["3clinicalf", "Dead", ]    <- p.3cfD	

### From 4 clinical c	
a.P_prev5["4clinicalc", "4clinicalc", ]    <- 0 	
a.P_prev5["4clinicalc", "4clinicald", ] <- 1 - p.4ccD	
a.P_prev5["4clinicalc", "Dead", ]    <- p.4ccD	

### From 4 clinical d	
a.P_prev5["4clinicald", "4clinicald", ]    <- 0
a.P_prev5["4clinicald", "4clinicale", ] <- 1 - p.4cdD	
a.P_prev5["4clinicald", "Dead", ]    <- p.4cdD	

### From 4 clinical e	
a.P_prev5["4clinicale", "4clinicale", ]    <- 0
a.P_prev5["4clinicale", "4clinicalf", ] <- 1 - p.4ceD	
a.P_prev5["4clinicale", "Dead", ]    <- p.4ceD	

### From 4 clinical f	
a.P_prev5["4clinicalf", "4clinicalf", ]    <- 1 - p.4cfNn - p.4cfD
a.P_prev5["4clinicalf", "Normalhpp", ] <- p.4cfNn	
a.P_prev5["4clinicalf", "Dead", ]    <- p.4cfD	

### From Dead
a.P_prev5["Dead", "Dead", ] <- 1

# check rows add up to 1
rowSums(a.P_prev2)


# check rows add up to 1
rowSums(a.P_prev5)
###

### From Normal helicobacter (-)
#m.P_prev5["Normalhpn", "Normalhpn"] <- 1 - p.NnD - p.p5NnNp - p.p5NnGn
#m.P_prev5["Normalhpn", "Normalhpp"]    <- p.p5NnNp
#m.P_prev5["Normalhpn", "Gastritishpn"]    <- p.p5NnGn
#m.P_prev5["Normalhpn", "Dead"]    <- p.NnD

### From Normal helicobacter (+) 
#m.P_prev5["Normalhpp", "Normalhpp"]    <- 1 - p.p5NpGp - p.p5NpNn - p.NpD 
#m.P_prev5["Normalhpp", "Gastritishpp"]    <- p.p5NpGp
#m.P_prev5["Normalhpp", "Normalhpn"]    <- p.p5NpNn
#m.P_prev5["Normalhpp", "Dead"]    <- p.NpD

### From Gastritis Helicobacter (-)
#m.P_prev5["Gastritishpn", "Gastritishpn"]    <- 1 - p.p5GnAn - p.p5GnGp - p.p5GnNn - p.GnD
#m.P_prev5["Gastritishpn", "Atrophyhpn"]    <- p.p5GnAn
#m.P_prev5["Gastritishpn", "Gastritishpp"]    <- p.p5GnGp
#m.P_prev5["Gastritishpn", "Normalhpn"]    <- p.p5GnNn
#m.P_prev5["Gastritishpn", "Dead"]    <- p.GnD

### From Gastritis Helicobacter (+)
#m.P_prev5["Gastritishpp", "Gastritishpp"]    <- 1 - p.p5GpAp - p.p5GpGn - p.p5GpNp - p.GpD
#m.P_prev5["Gastritishpp", "Atrophyhpp"]    <- p.p5GpAp
#m.P_prev5["Gastritishpp", "Gastritishpn"]    <- p.p5GpGn
#m.P_prev5["Gastritishpp", "Normalhpp"]    <- p.p5GpNp
#m.P_prev5["Gastritishpp", "Dead"]    <- p.GpD

### From Atrophy Helicobacter (-) 
#m.P_prev5["Atrophyhpn", "Atrophyhpn"]    <- 1 - p.p5AnIn - p.p5AnAp - p.p5AnGn -  p.AnD
#m.P_prev5["Atrophyhpn", "Intestinalhpn"]    <- p.p5AnIn
#m.P_prev5["Atrophyhpn", "Atrophyhpp"]    <- p.p5AnAp
#m.P_prev5["Atrophyhpn", "Gastritishpn"]    <- p.p5AnGn
#m.P_prev5["Atrophyhpn", "Dead"]    <- p.AnD

### From Atrophy Helicobacter (-) 
#m.P_prev5["Atrophyhpp", "Atrophyhpp"]    <- 1 - p.p5ApIp - p.p5ApAn - p.p5ApGp - p.ApD
#m.P_prev5["Atrophyhpp", "Intestinalhpp"]    <- p.p5ApIp
#m.P_prev5["Atrophyhpp", "Atrophyhpn"]    <- p.p5ApAn
#m.P_prev5["Atrophyhpp", "Gastritishpp"]    <- p.p5ApGp
#m.P_prev5["Atrophyhpp", "Dead"]    <- p.ApD

### From Intestinal Helicobacter (-)
#m.P_prev5["Intestinalhpn", "Intestinalhpn"]    <- 1 - p.p5InDn - p.p5InIp - p.p5InAn - p.InD 
#m.P_prev5["Intestinalhpn", "Dysplasiahpn"]    <- p.p5InDn
#m.P_prev5["Intestinalhpn", "Intestinalhpp"]    <- p.p5InIp
#m.P_prev5["Intestinalhpn", "Atrophyhpn"]    <- p.p5InAn
#m.P_prev5["Intestinalhpn", "Dead"]    <- p.InD

### From Intestinal Helicobacter (+)
#m.P_prev5["Intestinalhpp", "Intestinalhpp"]    <- 1 - p.p5IpDp - p.p5IpIn - p.p5IpAp - p.IpD  
#m.P_prev5["Intestinalhpp", "Dysplasiahpp"]    <- p.p5IpDp
#m.P_prev5["Intestinalhpp", "Intestinalhpn"]    <- p.p5IpIn
#m.P_prev5["Intestinalhpp", "Atrophyhpp"]    <- p.p5IpAp
#m.P_prev5["Intestinalhpp", "Dead"]    <- p.IpD

### From Dysplasia Helicobacter (-)
#m.P_prev5["Dysplasiahpn", "Dysplasiahpn"]    <- 1 - p.p5DnDp - p.p5DnIn -  p.DnD - p.p5_Dn1p
#m.P_prev5["Dysplasiahpn", "Dysplasiahpp"]    <- p.p5DnDp
#m.P_prev5["Dysplasiahpn", "Intestinalhpn"]    <- p.p5DnIn
#m.P_prev5["Dysplasiahpn", "1preclinical"]    <- p.p5_Dn1p
#m.P_prev5["Dysplasiahpn", "Dead"]    <- p.DnD

### From Dysplasia Helicobacter (+)
#m.P_prev5["Dysplasiahpp", "Dysplasiahpp"]    <- 1 - p.p5DpDn - p.p5DpIp - p.DpD - p.p5_Dp1p
#m.P_prev5["Dysplasiahpp", "Dysplasiahpn"]    <- p.p5DpDn
#m.P_prev5["Dysplasiahpp", "Intestinalhpp"]    <- p.p5DpIp
#m.P_prev5["Dysplasiahpp", "1preclinical"]    <- p.p5_Dp1p
#m.P_prev5["Dysplasiahpp", "Dead"]    <- p.DpD

### From 1 preclinical 	
#m.P_prev5["1preclinical", "1preclinical"]    <- 1 - p.p5_1p1ca - p.p5_1p2p - p.1pD	
#m.P_prev5["1preclinical", "1clinicala"]    <- p.p5_1p1ca	
#m.P_prev5["1preclinical", "2preclinical"]    <- p.p5_1p2p	
#m.P_prev5["1preclinical", "Dead"]    <- p.1pD	

### From 2 preclinical	
#m.P_prev5["2preclinical", "2preclinical"]    <- 1 - p.p5_2p2ca - p.p5_2p3p - p.2pD	
#m.P_prev5["2preclinical", "2clinicala"]    <- p.p5_2p2ca	
#m.P_prev5["2preclinical", "3preclinical"]    <- p.p5_2p3p	
#m.P_prev5["2preclinical", "Dead"]    <- p.2pD	

### From 3 preclinical	
#m.P_prev5["3preclinical", "3preclinical"]    <- 1 - p.p5_3p3ca - p.p5_3p4p - p.3pD	
#m.P_prev5["3preclinical", "3clinicala"]    <- p.p5_3p3ca	
#m.P_prev5["3preclinical", "Dead"]    <- p.3pD	
#m.P_prev5["3preclinical", "4preclinical"]    <- p.p5_3p4p	

### From 4 preclinical	
#m.P_prev5["4preclinical", "4preclinical"]    <- 1 - p.p5_4p4ca - p.4pD	
#m.P_prev5["4preclinical", "4clinicala"]    <- p.p5_4p4ca	
#m.P_prev5["4preclinical", "Dead"]    <- p.4pD	

### From 1 clinical a	
#m.P_prev5["1clinicala", "1clinicala"]    <- 1 - p.1caD - p.p5_1ca1cb	
#m.P_prev5["1clinicala", "Dead"]    <- p.1caD	
#m.P_prev5["1clinicala", "1clinicalb"] <- p.p5_1ca1cb	

### From 2 clinical a	
#m.P_prev5["2clinicala", "2clinicalb"] <- p.p5_2ca2cb	
#m.P_prev5["2clinicala", "2clinicala"]    <- 1 - p.2caD - p.p5_2ca2cb	
#m.P_prev5["2clinicala", "Dead"]    <- p.2caD	

### From 3 clinical a	
#m.P_prev5["3clinicala", "3clinicala"]    <- 1 - p.3caD - p.p5_3ca3cb 	
#m.P_prev5["3clinicala", "3clinicalb"] <- p.p5_3ca3cb	
#m.P_prev5["3clinicala", "Dead"]    <- p.3caD	

### From 4 clinical a 	
#m.P_prev5["4clinicala", "4clinicala"]    <- 1 - p.4caD - p.p5_4ca4cb	
#m.P_prev5["4clinicala", "4clinicalb"] <- p.p5_4ca4cb  	
#m.P_prev5["4clinicala", "Dead"]    <- p.4caD	

### From 1 clinical b	
#m.P_prev5["1clinicalb", "1clinicalb"]    <- 1 - p.1cbD - p.p5_1cb1cc	
#m.P_prev5["1clinicalb", "1clinicalc"] <- p.p5_1cb1cc	
#m.P_prev5["1clinicalb", "Dead"]    <- p.1cbD	

### From 2 clinical b	
#m.P_prev5["2clinicalb", "2clinicalb"]    <- 1 - p.2cbD - p.p5_2cb2cc	
#m.P_prev5["2clinicalb", "2clinicalc"] <- p.p5_2cb2cc	
#m.P_prev5["2clinicalb", "Dead"]    <- p.2cbD	

### From 3 clinical b	
#m.P_prev5["3clinicalb", "3clinicalb"]    <- 1 - p.3cbD - p.p5_3cb3cc 	
#m.P_prev5["3clinicalb", "3clinicalc"] <- p.p5_3cb3cc	
#m.P_prev5["3clinicalb", "Dead"]    <- p.3cbD	

### From 4 clinical b	
#m.P_prev5["4clinicalb", "4clinicalb"]    <- 1 - p.4cbD - p.p5_4cb4cc 	
#m.P_prev5["4clinicalb", "4clinicalc"] <- p.p5_4cb4cc	
#m.P_prev5["4clinicalb", "Dead"]    <- p.4cbD	

### From 1 clinical c	
#m.P_prev5["1clinicalc", "1clinicalc"]    <- 1 - p.1ccD - p.p5_1ccNn	
#m.P_prev5["1clinicalc", "Normalhpn"] <- p.p5_1ccNn	
#m.P_prev5["1clinicalc", "Dead"]    <- p.1ccD	

### From 2 clinical c	
#m.P_prev5["2clinicalc", "2clinicalc"]    <- 1 - p.2ccD - p.p5_2ccNn	
#m.P_prev5["2clinicalc", "Normalhpn"] <- p.p5_2ccNn	
#m.P_prev5["2clinicalc", "Dead"]    <- p.2ccD	

### From 3 clinical c	
#m.P_prev5["3clinicalc", "3clinicalc"]    <- 1 - p.3ccD - p.p5_3ccNn 	
#m.P_prev5["3clinicalc", "Normalhpn"] <- p.p5_3ccNn
#m.P_prev5["3clinicalc", "Dead"]    <- p.3ccD	

### From 4 clinical c	
#m.P_prev5["4clinicalc", "4clinicalc"]    <- 1 - p.4ccD - p.p5_4ccNn 	
#m.P_prev5["4clinicalc", "Normalhpn"] <- p.p5_4ccNn	
#m.P_prev5["4clinicalc", "Dead"]    <- p.4ccD	

### From Dead
#m.P_prev5["Dead", "Dead"] <- 1

# check rows add up to 1
#rowSums(m.P_prev5)






#### 04 Run Markov model ####

state_membership <- array(NA_real_, dim= c(n.t, n.s), dimnames = list (cycle= 1:n.t, state =v.n))

for (i in 2:n.t) { state_membership[i, ] <- state_membership[i - 1, ] %*% a.P_noprev[, , i-1]  }


#matplot(1:n.t, state_membership, type = "l")

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_no_prev[t + 1, ] <- m.M_no_prev[t, ] %*% a.P_noprev [ , , t]       # estimate the Markov trace for cycle the next cycle (t + 1)
}

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_prev[t + 1, ] <- m.M_prev[t, ] %*% a.P_prev [, , t]        # estimate the Markov trace for cycle the next cycle (t + 1)
}        

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_prev2[t + 1, ] <- m.M_prev2[t, ] %*% a.P_prev2 [, , t]        # estimate the Markov trace for cycle the next cycle (t + 1)
}

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_prev3[t + 1, ] <- m.M_prev3[t, ] %*% a.P_prev3 [, , t]         # estimate the Markov trace for cycle the next cycle (t + 1)
}

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_prev4[t + 1, ] <- m.M_prev4[t, ] %*% a.P_prev4 [, , t]         # estimate the Markov trace for cycle the next cycle (t + 1)
}

for (t in 1:n.t){                         # loop through the number of cycles
  m.M_prev5[t + 1, ] <- m.M_prev5[t, ] %*% a.P_prev5 [, , t]        # estimate the Markov trace for cycle the next cycle (t + 1)
}

head(m.M_no_prev)  # show the first 6 lines of the matrix
head(m.M_prev)  # show the first 6 lines of the matrix
head(m.M_prev2)  # show the first 6 lines of the matrix
head(m.M_prev3)  # show the first 6 lines of the matrix
head(m.M_prev4)  # show the first 6 lines of the matrix
head(m.M_prev5)  # show the first 6 lines of the matrix



#### 05 Compute and Plot Epidemiological Outcomes ####
  ##    05.1 Cohort trace #####
col_set <- wes_palette("Darjeeling1", n.s, type = "continuous")
  viridis::turbo(n.s)

matplot(m.M_no_prev, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace - No intervention", 
        col = col_set)              # create a plot of the data
      

legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph


matplot(m.M_prev, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace - Prevention Strategy 1: Urea", 
        col = col_set)              # create a plot of the data


legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph

matplot(m.M_prev2, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace - Prevention Strategy 2: 
        Antígeno Fecal", 
        col = col_set)              # create a plot of the data


legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph

matplot(m.M_prev3, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace - Prevention Strategy 3 EDA", 
        col = col_set)              # create a plot of the data


legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph

matplot(m.M_prev4, type = 'l', 
        ylab = "Probability of state occupancy Strategy4: PS-EDA",
        xlab = "Cycle",
        main = "Cohort Trace", 
        col = col_set)              # create a plot of the data


legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph

matplot(m.M_prev5, type = 'l', 
        ylab = "Probability of state occupancy Strategy5: Serología-Antígeno Fecal",
        xlab = "Cycle",
        main = "Cohort Trace", 
        col = col_set)              # create a plot of the data


legend("topright", v.n, lty = 1:n.s, bty = "n", cex=0.35, ncol=6, col=col_set)  # add a legend to the graph


  ##    05.2 Overall Survival (OS) #####
v.os_no_prev <- 1 - m.M_no_prev[, "Dead"]       # calculate the overall survival (OS) probability for no prevention
v.os_prev <- 1 - m.M_prev[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 1
v.os_prev2 <- 1 - m.M_prev2[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 2
v.os_prev3 <- 1 - m.M_prev3[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 3 
v.os_prev4 <- 1 - m.M_prev4[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 4 
v.os_prev5 <- 1 - m.M_prev5[, "Dead"]       # calculate the overall survival (OS) probability for prevention strategy 5


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

  ##    05.2.1 Life Expectancy (LE) #####
#Without prevention
v.le <- sum(v.os_no_prev)                       # summing probablity of OS over time  (i.e. life expectancy)
v.le1 <- sum(v.os_prev)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 1
v.le2 <- sum(v.os_prev2)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 2
v.le3 <- sum(v.os_prev3)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 3
v.le4 <- sum(v.os_prev4)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 4
v.le5 <- sum(v.os_prev5)                       # summing probablity of OS over time  (i.e. life expectancy) Prevention Strategy 5

print(v.le + 30)
print(v.le1 + 30)
print(v.le2 + 30)
print(v.le3 + 30)
print(v.le4 + 30)
print(v.le5 + 30)

#### 05.3 Disease Prevalence 
  ##    05.3 H.Pylori prevalence  --------------------------------------------------

## 05.3.1 No prevention 

v.noprev <- rowSums(m.M_no_prev[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_no_prev
plot(v.noprev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence")

## 05.3.2 Prevention Strategy 1 

v.preva1 <- rowSums(m.M_prev[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev
plot(v.preva1,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence EDA-Ureasa")

## 05.3.3 Prevention Strategy 2 
v.preva2 <- rowSums(m.M_prev2[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev2
plot(v.preva2,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Antigeno")

## 05.3.4 Prevention Strategy 3 
v.preva3 <- rowSums(m.M_prev3[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev3
plot(v.preva3,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Serología")

## 05.3.5 Prevention Strategy 4 

v.preva4 <- rowSums(m.M_prev4[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev4
plot(v.preva4,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Antígeno fecal")

## 05.3.6 Prevention Strategy 5 
v.preva5 <- rowSums(m.M_prev5[, c("Normalhpp", "Gastritishpp", "Atrophyhpp", "Intestinalhpp", "Dysplasiahpp" )])/v.os_prev5
plot(v.preva5,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "H. pilory prevalence Antígeno fecal")


  ##    05.4 Cancer Stages  ----------------------------------------------------------


## 05.4.1 Cancer Stages with no No prevention
v.preva_no_prev <- rowSums(m.M_no_prev[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf")])/v.os_no_prev
plot(v.preva_no_prev*100000,
     ylim = c(0, 50),
     ylab = "Prevalence per 100000",
     xlab = "Cycle",
     main = "Disease prevalence with no prevention")

## 05.4.2 Cancer Stages with prevention

##Prevention Strategy 1 
v.preva_prev <- rowSums(m.M_prev[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf")])/v.os_prev
plot(v.preva_prev*100000,
     ylim = c(0, 20),
     ylab = "Prevalence per 100000",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 2 
v.preva_prev2 <- rowSums(m.M_prev2[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf")])/v.os_prev2
plot(v.preva_prev2*100000,
     ylim = c(0, 20),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 3
v.preva_prev3 <- rowSums(m.M_prev3[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf" )])/v.os_prev3
plot(v.preva_prev3*100000,
     ylim = c(0, 30),
     ylab = "Prevalence per 100000",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 4
v.preva_prev4 <- rowSums(m.M_prev4[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf")])/v.os_prev4
plot(v.preva_prev4 * 100000,
     ylim = c(0, 40),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

##Prevention Strategy 5
v.preva_prev5 <- rowSums(m.M_prev5[, c("1preclinical", "1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf", "2preclinical", "2clinicala", "2clinicalb", "2clinicalc", "2clinicald", "2clinicale", "2clinicalf","3preclinical", "3clinicala", "3clinicalb", "3clinicalc", "3clinicald", "3clinicale", "3clinicalf", "4preclinical", "4clinicala", "4clinicalb" , "4clinicalc", "4clinicald", "4clinicale", "4clinicalf" )])/v.os_prev5
plot(v.preva_prev5*100000,
     ylim = c(0, 20),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 05.4.2 Proportion of sick in Gastritis - No prevention
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
     main = "Proportion of sick in Gastritis- EDA", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 4 
v.prop.G4 <- rowSums(m.M_prev4[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev4
plot(0:n.t, v.prop.G4,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - PSEDA", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Gastritis - Prevention Strategy 5 
v.prop.G5 <- rowSums(m.M_prev5[, c("Gastritishpp", "Gastritishpn")]) / v.os_prev4
plot(0:n.t, v.prop.G5,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in Gastritis - Serología", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - No Prevention
v.prop.1c_noprev <- rowSums(m.M_no_prev[, c("1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf" )]) / v.preva_no_prev
plot(0:n.t, v.prop.1c_noprev,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 1 
v.prop.1c_prev <- rowSums(m.M_prev[, c("1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf" )]) / v.preva_prev
plot(0:n.t, v.prop.1c_prev,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 2 
v.prop.1c_prev2 <- rowSums(m.M_prev2[, c("1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf" )]) / v.preva_prev2
plot(0:n.t, v.prop.1c_prev2,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 3 
v.prop.1c_prev3 <- rowSums(m.M_prev3[, c("1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf" )]) / v.preva_prev3
plot(0:n.t, v.prop.1c_prev3,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 4 
v.prop.1c_prev4 <- rowSums(m.M_prev4[, c("1clinicala", "1clinicalb", "1clinicalc", "1clinicald", "1clinicale", "1clinicalf" )]) / v.preva_prev4
plot(0:n.t, v.prop.1c_prev4,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")

#### 05.4 Proportion of sick in Stage 1 - Prevention Strategy 5 
v.prop.1c_prev5 <- rowSums(m.M_prev5[, c("1clinicala", "1clinicalb", "1clinicalc" )]) / v.preva_prev5
plot(0:n.t, v.prop.1c_prev5,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in stage 1", 
     col = "black", type = "l")


#### 06 Compute Cost-Effectiveness Outcomes ####
  ##    06.1 Vector with costs and utilities  -------------------------------------

### Vectors with costs and utilities by treatment
v.u_no_prev <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)
v.u_prev    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)
v.u_prev2    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)
v.u_prev3    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)
v.u_prev4    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)
v.u_prev5    <- c(u.Nn, u.Np, u.Gn, u.Gp, u.An, u.Ap, u.In, u.Ip, u.Dn, u.Dp, u.1p, u.2p, u.3p, u.4p, u.1ca, u.1cb, u.1cc, u.1cd, u.1ce, u.1cf, u.2ca, u.2cb, u.2cc, u.2cd, u.2ce, u.2cf, u.3ca, u.3cb, u.3cc, u.3cd, u.3ce, u.3cf, u.4ca, u.4cb, u.4cc, u.4cd, u.4ce, u.4cf, u.D)


v.c_no_prev <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)

## Arrays for costs --> No funcionó 

#a.costs <- array(NA_real_, dim = c(n.str, n.s, n.t), dimnames = list( strategy  = v.names.str, states = v.n, cycle  = 1:n.t)) 
  
#a.costs[1, , ] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[2, , 1] <- c(c.Nn + c.UreaAire + c.errn.1, c.Np + c.UreaAire + c.errp.1, c.Gn + c.UreaAire + c.errn.1, c.Gp + c.UreaAire + c.errp.1, c.An + c.UreaAire + c.errn.1, c.Ap + c.UreaAire+ c.errp.1, c.In + c.UreaAire + c.errn.1, c.Ip + c.UreaAire + c.errp.1, c.Dn + c.UreaAire + c.errn.1, c.Dp + c.UreaAire + c.errp.1, c.1p , c.2p , c.3p , c.4p , c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[2, , 2:n.t ] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[3, , 1] <- c(c.Nn +  c.Antigenofecal + c.errn.2, c.Np + c.Antigenofecal + c.errp.2, c.Gn + c.Antigenofecal + c.errn.2, c.Gp + c.Antigenofecal+ c.errp.2, c.An + c.Antigenofecal+ c.errn.2, c.Ap + c.Antigenofecal + c.errp.2, c.In +  c.Antigenofecal + c.errn.2, c.Ip + c.Antigenofecal + c.errp.2, c.Dn + c.Antigenofecal + c.errn.2, c.Dp + c.Antigenofecal + c.errp.2, c.1p , c.2p , c.3p , c.4p , c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[3, , 2:n.t] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[4, , 1:15] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[4, , 16] <- c(c.Nn + c.EDABp, c.Np + c.EDABp + c.errp.3, c.Gn + c.EDABp, c.Gp+ c.EDABp + c.errp.3, c.An + c.EDABp + c.errp.3, c.Ap + c.EDABp + c.errp.3, c.In + c.EDABp , c.Ip + c.EDABp + c.errp.3, c.Dn + c.EDABp, c.Dp + c.EDABp + c.errp.3, c.1p + c.EDABp, c.2p + c.EDABp, c.3p + c.EDABp, c.4p + c.EDABp, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[4, , 17:n.t] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[5, , 1:10] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[5, , 11] <- c(c.Nn + c.PSEda, c.Np + c.PSEda, c.Gn + c.PSEda, c.Gp + c.PSEda, c.An + c.PSEda, c.Ap + c.PSEda, c.In + c.PSEda, c.Ip + c.PSEda, c.Dn + c.PSEda, c.Dp + c.PSEda, c.1p , c.2p , c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[5, , 12:n.t] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[6, , 1] <- c(c.Nn + c.Serología + c.errn.5, c.Np + c.Serología + c.errp.5, c.Gn +  c.Serología + c.errn.5, c.Gp + c.Serología + c.errp.5, c.An +  c.Serología + c.errn.5, c.Ap + c.Serología + c.errp.5, c.In + c.Serología + c.errn.5, c.Ip + c.Serología + c.errp.5, c.Dn + c.Serología + c.errn.5, c.Dp + c.Serología + c.errp.5, c.1p , c.2p , c.3p , c.4p , c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
#a.costs[6, , 2:n.t] <- c(c.Nn, c.Np, c.Gn, c.Gp, c.An, c.Ap, c.In, c.Ip, c.Dn, c.Dp, c.1p, c.2p, c.3p, c.4p, c.1ca, c.1cb, c.1cc, c.1cd, c.1ce, c.1cf, c.2ca, c.2cb, c.2cc, c.2cd, c.2ce, c.2cf, c.3ca, c.3cb, c.3cc, c.3cd, c.3ce, c.3cf, c.4ca, c.4cb, c.4cc, c.4cd, c.4ce, c.4cf, c.D)
  
#a.costs[, , 1]


  ##    06.2 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_prev <- m.M_no_prev %*% v.u_no_prev
v.tu_prev <- m.M_prev %*% v.u_prev
v.tu_prev2 <- m.M_prev2 %*% v.u_prev2
v.tu_prev3 <- m.M_prev3 %*% v.u_prev3
v.tu_prev4 <- m.M_prev4 %*% v.u_prev4
v.tu_prev5 <- m.M_prev5 %*% v.u_prev5

v.tc_no_prev <- m.M_no_prev %*% v.c_no_prev
v.tc_prev    <- m.M_prev %*% v.c_no_prev
v.tc_prev2    <- m.M_prev2 %*% v.c_no_prev
v.tc_prev3    <- m.M_prev3 %*% v.c_no_prev
v.tc_prev4    <- m.M_prev4 %*% v.c_no_prev
v.tc_prev5    <- m.M_prev5 %*% v.c_no_prev

v.tc_prev[1]    <- v.tc_prev[1] + c.UreaAire + c.errp.1*0.78 + c.errn.1*0.22 
v.tc_prev2[1]    <- v.tc_prev2[1] + c.Antigenofecal + c.errp.2*0.78* + c.errn.2*0.22
v.tc_prev3[16]    <- v.tc_prev3[16] + c.EDABp + c.errp.3*0.78  ## este número se podría afinar 
v.tc_prev4[11]    <- v.tc_prev4[11] + c.PSEda + c.errp.4 
v.tc_prev5[1]    <- v.tc_prev5[1] + c.Serología + c.errp.5*0.78 + c.errn.5*0.22 


  ##    06.3 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_prev <- t(v.tu_no_prev) %*% v.dwe  
tu.d_prev    <- t(v.tu_prev)    %*% v.dwe
tu.d_prev2    <- t(v.tu_prev2)    %*% v.dwe
tu.d_prev3    <- t(v.tu_prev3)    %*% v.dwe
tu.d_prev4    <- t(v.tu_prev4)    %*% v.dwe
tu.d_prev5    <- t(v.tu_prev5)    %*% v.dwe

tc.d_no_prev <- t(v.tc_no_prev) %*% v.dwc
tc.d_prev    <- t(v.tc_prev)    %*% v.dwc
tc.d_prev2    <- t(v.tc_prev2)    %*% v.dwc
tc.d_prev3    <- t(v.tc_prev3)    %*% v.dwc
tc.d_prev4    <- t(v.tc_prev4)    %*% v.dwc
tc.d_prev5    <- t(v.tc_prev5)    %*% v.dwc

### Vector
v.tc.d <- c(tc.d_no_prev, tc.d_prev, tc.d_prev2, tc.d_prev3, tc.d_prev4, tc.d_prev5)
v.tu.d <- c(tu.d_no_prev, tu.d_prev, tu.d_prev2, tu.d_prev3, tu.d_prev4, tu.d_prev5)

# Matrix with discounted costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.tc.d,
                   Effect   = v.tu.d)
m.ce



#### 07 Compute ICERs of FONIS CaGa model ####
m.cea <- calculate_icers(cost = m.ce$Cost,
                         effect = m.ce$Effect,
                         strategies = m.ce$Strategy)
m.cea[4,4] <- m.ce$Cost[4]-m.ce$Cost[1]

m.cea$Inc_Cost <- m.cea$Cost-m.ce$Cost[1]
m.cea$Inc_Effect <- m.cea$Effect-m.ce$Effect[1]
m.cea$ICER <- m.cea$Inc_Cost/m.cea$Inc_Effect
m.cea

plot(m.cea)
plot(m.cea, plot_frontier_only = TRUE)

n.wtp <- 387.729768

m.cea$`NMB 100` <- (m.cea$Effect*n.wtp) - m.cea$Cost
m.cea




#### 08 Plot frontier of Sick-Sicker model ####
plot(m.cea)
ggsave("figs/Markov-FONISCaGast-CEA-Frontier-gastricmucosa.png", width = 8, height = 6)


