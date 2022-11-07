# ---- updaload packages ----
rm(list=ls())
if(!("LMest"%in%installed.packages())) install.packages("LMest")
if(!("mvtnorm"%in%installed.packages())) install.packages("mvtnorm")
require(LMest) 
require(mvtnorm)
source("estFM2.R")
source("HM_FM_Test.R")
source("sqg.R")
source("print.HM_FM_Test.R")

# ---- load data ----
load("exampleData.Rdata")
head(Y)

#Perform the classical LR test for two different number of states 
test1 <-HM_FM_Test(data=Y,id = "id",time="time",lv=2:3,responses=NULL,type="single")

#Perform the LR test based on consecutive triplets for two different number of states 
test2 <-HM_FM_Test(data=Y,id = "id",time="time",lv=2:3,responses=NULL,type="triplets")


# ---- display output ----
print(test1)

print(test2)

