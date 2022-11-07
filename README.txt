This repository contains the R functions used for the paper
"ON A CLASS OF FINITE MIXTURE MODELS THAT INCLUDES HIDDENMARKOV MODELS FOR LONGITUDINAL DATA AND RELATEDMISSPECIFICATION TESTS" by
- F. Bartolucci (University of Perugia, IT)
- S. Pandolfi (University of Perugia, IT)  
- F. Pennoni (University of Milano-Bicocca, IT)

est_FM2.R ---> estimate the constrained finite mixture model named FM2 model with suitable constraints on the conditional distribution of the response vectors

HM_FM_Test.R ---> perform the misspecification tests for hidden Markov model, allowing one to choose the classical LR test or the test based on consecutive triplets of observations (with both Bonferroni and Simes corrections)

exampleTest.R ---> example file that loads the workspace file "exampleData.RData" and shows how to perform the misspecification tests

exampleData.RData ---> workspace file containing a simulated dataset 

sqg.R ---> internal function that creates a matrix with all possible combinations of sequence of number given in input 

print.HM_FM_Test ---> internal function. Given the output, it is written in a readable form