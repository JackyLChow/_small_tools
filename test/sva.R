###################################################
### load packages
###################################################
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)


###################################################
### metadata
###################################################
pheno = pData(bladderEset)


###################################################
### counts data
###################################################
edata = exprs(bladderEset)


###################################################
### model matrix
###################################################
mod = model.matrix(~as.factor(cancer), data=pheno)


###################################################
### null model
###################################################
mod0 = model.matrix(~1,data=pheno)


###################################################
### estimate number of surrogate variables
###################################################
n.sv = num.sv(edata,mod,method="leek")
n.sv


###################################################
### calculate surrogate variables
###################################################
svobj = sva(edata,mod,mod0,n.sv=n.sv)


###################################################
### model matricies with surrogate variables
###################################################
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)


###################################################
### limma pipeline
###################################################
fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")








