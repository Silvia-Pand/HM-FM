# ---- updaload packages ----
rm(list=ls())
if(!("LMest"%in%installed.packages())) install.packages("LMest")
if(!("mvtnorm"%in%installed.packages())) install.packages("mvtnorm")
require(LMest) 
require(mvtnorm)

source("estFM2.R")
source("HM_FM_Test.R")
source("sqg.R")

# ---- load data ----
load("exampleData.Rdata")
head(Y)

#Perform the classical LR test for two different number of states 
test1 <-HM_FM_Test(data=Y,id = "id",time="time",lv=2:3,responses=NULL,type=c("single"))

# ---- display output ----
print("p-value")
print(test1$pvalueS)

#Perform the test based on consecutive triplets of observations for l=4 
test2 <-HM_FM_Test(data=Y,id = "id",time="time",lv=4,responses=NULL,type=c("triplets"))

# ---- display output ----
print("p-value")
print(test2$pvalueT)

