
library(gtools)
source('R/GENHETv3.1.R')


fj=read.table('Data/TheData.txt',header=TRUE) #read the data#
fj1=fj[,-2] #ID must be kept for GENHET to work#

#We ask R to scan the AllLoci.txt file with the names of all loci; just change the route#
locusname=scan("Data/AllLoci.txt",
               what="character",sep="\t")

#Run GENHET function in "GENHETv3.1.R" (the whole function)#
Htest=GENHET(dat=fj1,estimfreq="T",locname=locusname)
Htest
df=as.data.frame(Htest)
write.table(x=df,file="Het_All.txt")