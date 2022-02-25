library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

fstats_obj = rowFtests(edata,as.factor(pdata$strain))
hist(fstats_obj$p.value,col=2)

edge_study = build_study(edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3)

set.seed(3333)
B = 1000
tstats_obj = rowttests(edata,pdata$strain)
tstat0 = matrix(NA,nrow=dim(edata)[1],ncol=B)
tstat = tstats_obj$statistic
strain = pdata$strain
for(i in 1:B){
    strain0 = sample(strain)
    tstat0[,i] = rowttests(edata,strain0)$statistic
}

emp_pvals = empPvals(tstat,tstat0)
hist(emp_pvals,col=2)


















