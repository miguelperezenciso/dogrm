#!/usr/bin/env Rscript
# script to estimate genomic heritabilities using BGLR package
library(BGLR)
n=10000
niter = 50000
nburn=3000
# should contain add and dom G
filed='../dom.grm'
filea='../add.grm'
# trait analyzed
itrait = 1 
# should contain phenotypes
filey='../Phenotypes_TRN_10k.txt'
Y = read.table(filey)[,itrait+2]
Y = scale(Y[1:n])

# PC decomposition of Gs
G = matrix(scan(filea, n = n*n), n, n, byrow = TRUE)
diag(G) = diag(G)*1.05
EVD<-eigen(G)
PC1<-EVD$vectors%*%diag(sqrt(EVD$values))
save(PC1, file='PC1.RData')

G = matrix(scan(filed, n = n*n), n, n, byrow = TRUE)
diag(G) = diag(G)*1.05
EVD<-eigen(G)
PC2<-EVD$vectors%*%diag(sqrt(EVD$values))
save(PC2, file='PC2.RData')


load('../PC1.RData')
load('../PC2.RData')
ETA=list(list(X=PC1,model='BRR'), list(X=PC2,model='BRR'))
fm<-BGLR(y=Y,ETA=ETA,nIter=niter,burnIn=nburn,verbose=T)

vU<-scan('ETA_1_varB.dat')
vD<-scan('ETA_2_varB.dat')
vE<-scan('varE.dat')
# contains varA/varY and varD/varY
h2a<-vU/(vU+vD+vE)
h2d<-vD/(vU+vD+vE)

# plot
plot(density(h2a),xlim=c(0,1))
lines(density(h2d),col='blue')



