#####################################################################################
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

#### 02 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Prevention", "Prevention1")  

## Number of strategies
n.str <- length(v.names.str)
## Markov model parameters
age     <- 20                                 # age at baseline
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles

v.n  <- c("Healthy", "1preclinical", "2preclinical", "3preclinical", "4preclinical", "1clinical", "2clinical", "3clinical", "4clinical", "Dead")    # state names
n.s  <- length(v.n)                     # number of states
n.t  <- 60                              # number of cycles

p.HD <- 0.02  # probability to die when healthy

#Probability of no prevention

p.H1p <- 0.05  # probability to become stage 1 preclinical when healthy
p.1p1c <- 0.05  # probability to become stage 1 clinical when stage 1 preclinical
p.1p2p <- 0.05  # probability to become stage 2 preclinical when stage 1 preclinical
p.2p2c <- 0.05  # probability to become stage 2 clinical when stage 2 preclinical
p.2p3p <- 0.05  # probability to become stage 3 preclinical when stage 2 preclinical
p.3p3c <- 0.05  # probability to become stage 3 clinical when stage 3 preclinical
p.3p4p <- 0.05  # probability to become stage 4 preclinical when stage 3 preclinical
p.4p4c <- 0.05  # probability to become stage 4 clinical when stage 4 preclinical
p.1cH <- 0.08 # probability to become healthy when stage 1 clinical
p.2cH <- 0.05 # probability to become healthy when stage 2 clinical 
p.3cH <- 0.04 # probability to become healthy when stage 3 clinical 
p.4cH <- 0.03 # probability to become healthy when stage 4 clinical 

# Probability of transition when prevention (llenar según efectividades)

p.pH1p <- 0.05  # probability to become stage 1 preclinical when healthy
p.p1p1c <- 0.05  # probability to become stage 1 clinical when stage 1 preclinical
p.p1p2p <- 0.05  # probability to become stage 2 preclinical when stage 1 preclinical
p.p2p2c <- 0.05  # probability to become stage 2 clinical when stage 2 preclinical
p.p2p3p <- 0.05  # probability to become stage 3 preclinical when stage 2 preclinical
p.p3p3c <- 0.05  # probability to become stage 3 clinical when stage 3 preclinical
p.p3p4p <- 0.05  # probability to become stage 4 preclinical when stage 3 preclinical
p.p4p4c <- 0.05  # probability to become stage 4 clinical when stage 4 preclinical
p.p1cH <- 0.08 # probability to become healthy when stage 1 clinical
p.p2cH <- 0.05 # probability to become healthy when stage 2 clinical 
p.p3cH <- 0.04 # probability to become healthy when stage 3 clinical 
p.p4cH <- 0.03 # probability to become healthy when stage 4 clinical 

#Hazard ratio 

hr.1p <- 1 #hazard ratio of death in stage 1 preclinical vs healthy
hr.2p <- 1 #hazard ratio of death in stage 2 preclinical vs healthy
hr.3p <- 1 #hazard ratio of death in stage 3 preclinical vs healthy
hr.4p <- 1 #hazard ratio of death in stage 4 preclinical vs healthy
hr.1c <- 3 #hazard ratio of death in stage 1 clinical vs healthy
hr.2c <- 6  #hazard ratio of death in stage 2 clinical vs healthy   
hr.3c <- 9 #hazard ratio of death in stage 3 clinical vs healthy
hr.4c <- 10 #hazard ratio of death in stage 4 clinical vs healthy

r.HD    <- - log(1 - p.HD) # rate of death in healthy
r.1pD   <- hr.1p * r.HD  	 # rate of death in 1 preclinical
r.2pD   <- hr.2p * r.HD  	 # rate of death in 2 preclinical
r.3pD   <- hr.3p * r.HD  	 # rate of death in 3 preclinical
r.4pD   <- hr.4p * r.HD  	 # rate of death in 4 preclinical
r.1cD   <- hr.1c * r.HD  	 # rate of death in 1 clinical
r.2cD   <- hr.2c * r.HD  	 # rate of death in 2 clinical
r.3cD   <- hr.1c * r.HD  	 # rate of death in 3 clinical
r.4cD   <- hr.2c * r.HD  	 # rate of death in 4 clinical 
    
p.1pD <- 1- exp(-r.1pD)  # probability to die when stage 1 preclinical 
p.2pD <- 1- exp(-r.2pD) # probability to die when stage 2 preclinical
p.3pD <- 1- exp(-r.3pD) # probability to die when stage 3 preclinical
p.4pD <- 1- exp(-r.4pD) # probability to die when stage 4 preclinical
p.1cD <- 1- exp(-r.1cD)  # probability to die when stage 1 clinical 
p.2cD <- 1- exp(-r.2cD) # probability to die when stage 2 clinical
p.3cD <- 1- exp(-r.3cD) # probability to die when stage 3 clinical
p.4cD <- 1- exp(-r.4cD) # probability to die when stage 4 clinical

# Costs and utilities  
c.H  <- 400                     # cost of remaining one cycle healthy
c.1p  <- 400                     # cost of remaining one cycle stage 1 preclinical
c.2p  <- 400                     # cost of remaining one cycle stage 2 preclinical
c.3p  <- 400                     # cost of remaining one cycle stage 3 preclinical
c.4p  <- 400                     # cost of remaining one cycle stage 4 preclinical
c.1c  <- 4000                     # cost of remaining one cycle stage 1 clinical
c.2c  <- 8000                     # cost of remaining one cycle stage 2 clinical
c.3c  <- 12000                     # cost of remaining one cycle stage 3 clinical
c.4c  <- 15000                     # cost of remaining one cycle stage 4 clinical
c.D  <- 0  # cost of remaining one cycle dead
c.prev <- 200 # cost of prevention strategy 

u.H  <- 1                     # utility when healthy 
u.1p  <- 1                     # utility when stage 1 preclinical
u.2p  <- 1                     # utility when stage 2 preclinical
u.3p  <- 0.7                     # utility when stage 3 preclinical
u.4p  <- 0.8                     # utility when stage 4 preclinical
u.1c  <- 0.5                     # utility when stage 1 clinical
u.2c  <- 0.4                     # utility when stage 2 clinical
u.3c  <- 0.3                     # utility when stage 3 clinical
u.4c  <- 0.2                     # utility when stage 4 clinical
u.D  <- 0                       # utility when dead

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

# The cohort starts as healthy

m.M_no_prev[1, ] <- m.M_prev[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)                     # initialize first cycle of Markov trace

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
### From Healthy
m.P_noprev["Healthy", "Healthy"] <- 1 - p.HD - p.H1p
m.P_noprev["Healthy", "1preclinical"]    <- p.H1p
m.P_noprev["Healthy", "Dead"]    <- p.HD

### From 1 preclinical
m.P_noprev["1preclinical", "1preclinical"]    <- 1 - p.1p1c - p.1p2p - p.1pD
m.P_noprev["1preclinical", "1clinical"]    <- p.1p1c
m.P_noprev["1preclinical", "2preclinical"]    <- p.1p2p
m.P_noprev["1preclinical", "Dead"]    <- p.1pD

### From 2 preclinical
m.P_noprev["2preclinical", "2preclinical"]    <- 1 - p.2p2c - p.2p3p - p.2pD
m.P_noprev["2preclinical", "2clinical"]    <- p.2p2c
m.P_noprev["2preclinical", "3preclinical"]    <- p.2p3p
m.P_noprev["2preclinical", "Dead"]    <- p.2pD

### From 3 preclinical
m.P_noprev["3preclinical", "3preclinical"]    <- 1 - p.3p3c - p.3p4p - p.3pD
m.P_noprev["3preclinical", "3clinical"]    <- p.3p3c
m.P_noprev["3preclinical", "4preclinical"]    <- p.3p4p
m.P_noprev["3preclinical", "Dead"]    <- p.3pD

### From 4 preclinical
m.P_noprev["4preclinical", "4preclinical"]    <- 1 - p.4p4c - p.4pD
m.P_noprev["4preclinical", "4clinical"]    <- p.4p4c
m.P_noprev["4preclinical", "Dead"]    <- p.4pD

### From 1 clinical 
m.P_noprev["1clinical", "1clinical"]    <- 1 - p.1cD - p.1cH
m.P_noprev["1clinical", "Dead"]    <- p.1cD
m.P_noprev["1clinical", "Healthy"] <- p.1cH

### From 2 clinical 
m.P_noprev["2clinical", "2clinical"]    <- 1 - p.2cD - p.2cH
m.P_noprev["2clinical", "Dead"]    <- p.2cD
m.P_noprev["2clinical", "Healthy"] <- p.2cH

### From 3 clinical 
m.P_noprev["3clinical", "3clinical"]    <- 1 - p.3cD - p.3cH
m.P_noprev["3clinical", "Dead"]    <- p.3cD
m.P_noprev["3clinical", "Healthy"] <- p.3cH

### From 4 clinical 
m.P_noprev["4clinical", "4clinical"]    <- 1 - p.4cD - p.4cH
m.P_noprev["4clinical", "Dead"]    <- p.4cD
m.P_noprev["4clinical", "Healthy"] <- p.4cH

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
### From Healthy
m.P_prev["Healthy", "Healthy"] <- 1 - p.HD - p.pH1p
m.P_prev["Healthy", "1preclinical"]    <- p.pH1p
m.P_prev["Healthy", "Dead"]    <- p.HD

### From 1 preclinical
m.P_prev["1preclinical", "1preclinical"]    <- 1 - p.p1p1c - p.p1p2p - p.1pD
m.P_prev["1preclinical", "1clinical"]    <- p.p1p1c
m.P_prev["1preclinical", "2preclinical"]    <- p.p1p2p
m.P_prev["1preclinical", "Dead"]    <- p.1pD

### From 2 preclinical
m.P_prev["2preclinical", "2preclinical"]    <- 1 - p.p2p2c - p.p2p3p - p.2pD
m.P_prev["2preclinical", "2clinical"]    <- p.p2p2c
m.P_prev["2preclinical", "3preclinical"]    <- p.p2p3p
m.P_prev["2preclinical", "Dead"]    <- p.2pD

### From 3 preclinical
m.P_prev["3preclinical", "3preclinical"]    <- 1 - p.p3p3c - p.p3p4p - p.3pD
m.P_prev["3preclinical", "3clinical"]    <- p.p3p3c
m.P_prev["3preclinical", "4preclinical"]    <- p.p3p4p
m.P_prev["3preclinical", "Dead"]    <- p.3pD

### From 4 preclinical
m.P_prev["4preclinical", "4preclinical"]    <- 1 - p.p4p4c - p.4pD
m.P_prev["4preclinical", "4clinical"]    <- p.p4p4c
m.P_prev["4preclinical", "Dead"]    <- p.4pD

### From 1 clinical 
m.P_prev["1clinical", "1clinical"]    <- 1 - p.1cD - p.p1cH
m.P_prev["1clinical", "Dead"]    <- p.1cD
m.P_prev["1clinical", "Healthy"] <- p.p1cH

### From 2 clinical 
m.P_prev["2clinical", "2clinical"]    <- 1 - p.2cD - p.p2cH
m.P_prev["2clinical", "Dead"]    <- p.2cD
m.P_prev["2clinical", "Healthy"] <- p.p2cH

### From 3 clinical 
m.P_prev["3clinical", "3clinical"]    <- 1 - p.3cD - p.p3cH
m.P_prev["3clinical", "Dead"]    <- p.3cD
m.P_prev["3clinical", "Healthy"] <- p.p3cH

### From 4 clinical 
m.P_prev["4clinical", "4clinical"]    <- 1 - p.4cD - p.p4cH
m.P_prev["4clinical", "Dead"]    <- p.4cD
m.P_prev["4clinical", "Healthy"] <- p.p4cH

### From Dead
m.P_prev["Dead", "Dead"] <- 1

# check rows add up to 1
rowSums(m.P_prev)

#### 04 Run Markov model ####
for (t in 1:n.t){                                         # loop through the number of cycles
  m.M_no_prev[t + 1, ] <- t(m.M_no_prev[t, ]) %*% m.P_noprev # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_prev[t + 1, ]    <- t(m.M_prev[t, ])    %*% m.P_prev   # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_prev)  # show the first 6 lines of the matrix

#### 05 Compute and Plot Epidemiological Outcomes ####
#### 05.1 Cohort trace #####
matplot(m.M_no_prev, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace")              # create a plot of the data
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

#### 05.2 Overall Survival (OS) #####
v.os_no_prev <- 1 - m.M_no_prev[, "Dead"]       # calculate the overall survival (OS) probability for no prevention

plot(0:n.t, v.os_no_prev, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

#### 05.2.1 Life Expectancy (LE) #####
v.le <- sum(v.os_no_prev)                       # summing probablity of OS over time  (i.e. life expectancy)

#### 05.3 Disease prevalence #####
v.preva <- rowSums(m.M_no_prev[, c("1preclinical", "1clinical", "2preclinical", "2clinical", "3preclinical", "3clinical", "4preclinical", "4clinical" )])/v.os_no_prev
plot(v.preva,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#### 05.4 Proportion of sick in 1p state #####
v.prop.1p <- m.M_no_prev[, "1preclinical"] / v.preva
plot(0:n.t, v.prop.1p,
     xlab = "Cycle", 
     ylab = "Proportion", 
     main = "Proportion of sick in 1p state", 
     col = "black", type = "l")


#### 06 Compute Cost-Effectiveness Outcomes ####
### Vectors with costs and utilities by treatment
v.u_no_prev <- c(u.H, u.1p, u.1c, u.2p, u.2c, u.3p, u.3c, u.4p, u.4c, u.D)
v.u_prev    <- c(u.H, u.1p, u.1c, u.2p, u.2c, u.3p, u.3c, u.4p, u.4c, u.D)

v.c_no_prev <- c(c.H, c.1p, c.1c, c.2p, c.2c, c.3p, c.3c, c.4p, c.4c, c.D)
v.c_prev    <- c(c.H + c.prev, c.1p, c.1c, c.2p, c.2c, c.3p, c.3c, c.4p, c.4c, c.D)

#### 06.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# estimate mean QALys and costs
v.tu_no_prev <- m.M_no_prev %*% v.u_no_prev
v.tu_prev <- m.M_prev %*% v.u_prev

v.tc_no_prev <- m.M_no_prev %*% v.c_no_prev
v.tc_prev    <- m.M_prev    %*% v.c_prev

#### 06.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_prev <- t(v.tu_no_prev) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
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

#### 07 Compute ICERs of FONIS CaGa model ####
m.cea <- calculate_icers(cost = m.ce$Cost,
                         effect = m.ce$Effect,
                         strategies = m.ce$Strategy)
m.cea

#### 08 Plot frontier of FONIS CaGa model ####
plot(m.cea, xlim = c(15.7, 16.5))
ggsave("figs/Markov-FONISCaGast-CEA-Frontier.png", width = 8, height = 6)

