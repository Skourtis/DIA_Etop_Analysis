#Correp
BiocManager::install("CORREP")
library("CORREP")
d0 <- NULL
for(l in 1:10)
    d0 <- rbind(d0, rnorm(100))
## The simulated data corresponds to the real-world data of 25 genes and 10 conditions, each gene expression
## profiles was replicated 4 times.
d0<- t(d0)
## This step is to make the standard deviation of each replicate equals to 1
## so that we can model the covariance matrix as correlation matrix.
d0.std <- apply(d0, 1, function(x) x/sd(x))
input <- t(d0.std)
Methods_DIA$DIA_report_file_method2
M <- cor.balance(input, m=4, G=25)
