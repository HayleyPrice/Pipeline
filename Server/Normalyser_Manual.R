# #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop('At least one argument must be supplied (input file).n', call.=FALSE)
}

#setwd('E:/OneDrive/PhD/Project/Thesis/4_Pipeline/Pipeline/Server')
#file <- 'PXD004682_raw_test.txt'
#file <- 'PXD012039_input.txt'
#dataSet <- 'PXD004684'
#dataSet <- 'PXD012039'
####### ***** Log2 transformed
#library(NormalyzerDE)
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(MASS))

file <- args[1]
dataSet <- args[2]
jobdir<- paste(dataSet, "/", sep = "")

getrawdata <- as.matrix((read.table(file, header=F,sep="\t",stringsAsFactors=F,quote="")))

cols <- getrawdata[1,]
getrawdata <- getrawdata[-1,]
getrawdata[getrawdata == 0] <- 0.00001
getrawdata <- rbind(cols, getrawdata)

if(file.exists(dataSet)) {
  unlink(dataSet, recursive = TRUE)
} 

dir.create(dataSet)
dirs <- c("/NormalisedData", 
          "/Ttest", "/Ttest/Results", "/Ttest/PA", 
          "/Qmodel", "/Qmodel/subs", "/Qmodel/Results", "/Qmodel/FDR", "/Qmodel/PA", 
          "/Summary")

for(dir in dirs) {
  dir.create(paste(jobdir, dir, sep = ""))
}

setwd(paste('/mnt/hc-storage/users/hprice/Pipeline/', dataSet, "/NormalisedData", sep = ""))

#Sort the uploaded data based on replicates
b<-NULL
b<-as.factor(getrawdata[1,])
l<-levels(b)
b<-NULL
for(i in 1:length(l)){
  b<-cbind(b,getrawdata[,which(getrawdata[1,]==l[as.numeric(i)])])
}
getrawdata<-b

#Parse data for errors

checkrep<-getrawdata[1,]
repunique<-unique(checkrep)
for(i in 1:length(repunique))
{
  if(repunique[i]!=0)
  {
    if(length(grep(repunique[i],checkrep))<2)
    {
      abc<-paste("Number of replicates are less than 2 for the group ", repunique[i],sep="")
      class(abc)="try-error"
      if(inherits(abc,"try-error")){return(abc)}
      stop(paste("Number of replicates are less than 2 for the group ", repunique[i],sep=""))
    }
  }

}
#replace 0 with NA
rep0<-getrawdata[-1,]
rep0[which(rep0==0)]<-NA
getrawdata<-rbind(getrawdata[1,],rep0)


#HKflag=T
getEDdata<-((getrawdata[1,]))
filterrawdata1<-getrawdata
countna<-rowSums(!is.na(filterrawdata1[,which(checkrep>0)]))
filterrawdata1<-filterrawdata1[countna>=(1*ncol(filterrawdata1[,which(checkrep>0)])),]
filterED<-as.numeric(getEDdata[-which(getEDdata<1)])
HKflag=F

filterED<-as.numeric(getEDdata[-which(getEDdata<1)])
filterrawdata<-getrawdata[,-(1:(length(getEDdata)-length(filterED)))]
colnames(filterrawdata)<-getrawdata[2,-(1:(length(getEDdata)-length(filterED)))]
filterrawdata<-(as.matrix((filterrawdata[-(1:2),])))
class(filterrawdata)<-"numeric"


#CONVERT TO LOG

####### ***** Log
data2log<-log2((filterrawdata))
#Total intensity normalization
data2GI<-matrix(nrow=nrow(filterrawdata),ncol=ncol(filterrawdata),byrow=T)
data2ctr<-matrix(nrow=nrow(filterrawdata),ncol=ncol(filterrawdata),byrow=T)
data2med<-matrix(nrow=nrow(filterrawdata),ncol=ncol(filterrawdata),byrow=T)
data2mean<-matrix(nrow=nrow(filterrawdata),ncol=ncol(filterrawdata),byrow=T)

colsum<-colSums(filterrawdata,na.rm=T)
medofdata<-apply(filterrawdata,2,FUN="median",na.rm=T)
meanofdata<-apply(filterrawdata,2,FUN="mean",na.rm=T)

avgcolsum<-median(colsum)
for(i in 1:nrow(filterrawdata))
{
  data2GI[i,]<-unlist(sapply(1:ncol(filterrawdata),function(zd) {(filterrawdata[i,zd]/colsum[zd])*(avgcolsum)}))
  data2med[i,]<-unlist(sapply(1:ncol(filterrawdata),function(zd) {(filterrawdata[i,zd]/medofdata[zd])*mean(medofdata)}))
  data2mean[i,]<-unlist(sapply(1:ncol(filterrawdata),function(zd) {(filterrawdata[i,zd]/meanofdata[zd])*mean(meanofdata)}))
}

methodlist<-list(data2log,data2GI,data2med,data2mean)
methodnames<-c("Log2","TI-G","MedI-G","AI-G")

#Perform other norm. if the dataSet is not small
#if(nrow(filterrawdata)>100)
#{
#VSN NORMALIZATION
data2vsn<-justvsn((filterrawdata))

#QUANTILE NORMALIZATION
data2quantile<-normalize.quantiles((data2log),copy=T)

#SMAD normalization
mediandata<-apply(data2log,2,"median",na.rm=T)
maddata<-apply(data2log,2,function(x) mad(x,na.rm=T))
data2mad<-t(apply(data2log,1,function(x) ((x-mediandata)/maddata)))
data2mad<-data2mad+mean(mediandata)

#Global loess
data2loess<-normalizeCyclicLoess(data2log,method="fast")

#Global RLR
mediandata<-apply(data2log,1,"median",na.rm=T)
flag1=1
for(j in 1:ncol(data2log))
{
  LRfit<-rlm(as.matrix(data2log[,j])~mediandata,na.action=na.exclude)
  Coeffs<-LRfit$coefficients
  a<-Coeffs[2]
  b<-Coeffs[1]
  if(flag1==1)
  {
    globalfittedRLR<-(data2log[,j]-b)/a
    flag1=2
  }
  else
  {
    globalfittedRLR<-cbind(globalfittedRLR,(data2log[,j]-b)/a)
  }
}
colnames(globalfittedRLR)<-colnames(data2log)


#sink(sinkfile,type="output")
closeAllConnections()
#colnames(fittedLR)<-colnames(data2log)
colnames(data2quantile)<-colnames(data2log)

methodlist<-list(data2log,data2loess,globalfittedRLR,data2vsn,data2GI,data2med,data2mean,data2quantile)
methodnames<-c("None", "Loess-G","RLR-G","VSN-G","TI-G","MedI-G","AI-G","Quantile")

for(i in 1:length(methodlist)) {
  write.table(file=paste(dataSet, '_', methodnames[i],"-normalized_Log2",sep=""),
              cbind(getrawdata[-(1:2),(1:(length(getEDdata)-length(filterED)))],methodlist[[i]]),sep="\t",
              row.names=F,col.names=cols, quote = FALSE)
}
# 
