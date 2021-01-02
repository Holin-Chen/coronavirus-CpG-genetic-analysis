library(ggplot2)
library(ggmsa)
library(seqmagick)
library(Biostrings)

#import coronovirus genetic data
filelist <- list.files(path="./raw fasta data/", pattern = "*.fasta")
fa_seq <- lapply(paste0("./raw fasta data/" ,filelist), readDNAStringSet)
fa_seq <- do.call(c, fa_seq)
fa_seq

fa_seq@ranges@NAMES[3] <- "NC_004718.3_SARS-CoV"
fa_seq@ranges@NAMES[4] <- "NC_038294.1_MERS-CoV"
fa_seq@ranges@NAMES[5] <- "NC_045512.2_SARS-CoV-2"

## get GC content
bases=alphabetFrequency(fa_seq,baseOnly=TRUE)
bases[1:5,1:4]
ntotBases=apply(bases[1:5,1:4], 1, sum)
baseFreq=bases[1:5,1:4]/ntotBases
GCcontent=baseFreq[,"C"]+baseFreq[,"G"]
ATcontent=baseFreq[,"A"]+baseFreq[,"T"]

## look at CG dinucleotide content
cg=vmatchPattern("CG", fa_seq)
ncg=lengths(cg)
## compute the observed to expected ratio
ncg/(baseFreq[,"C"]*baseFreq[,"G"]*ntotBases) ## this shows CG rarely stay together.

## compare to the observed to expected ratio of TG
tg=vmatchPattern("TG", fa_seq)
ntg=lengths(tg)
ntg/(baseFreq[,"T"]*baseFreq[,"G"]*ntotBases) ## this shows TG presented more than expected.

######### look at GC content and CG dinucleotide distribution in 1000 bp windows in whole genome.

ss_1=seq(1, lengths(fa_seq)[1], by=1000)
ss_1=ss_1[-length(ss_1)] ## remove the last one
bat_SL_CoVZC45=DNAStringSet(fa_seq$MG772933, start=ss_1, end=ss_1+999)
ff_1=alphabetFrequency(bat_SL_CoVZC45, baseOnly=TRUE)
pCG_1=(ff_1[,"C"]+ff_1[,"G"])/rowSums(ff_1)
hist(pCG_1[pCG_1>0],100, main = "GC content in 1000 bp window")

## CG occurance
nCG_1=vcountPattern("CG", bat_SL_CoVZC45)
obsExp_1=nCG_1*1000/(ff_1[,"C"]*ff_1[,"G"])
mean(obsExp_1,na.rm=TRUE)
hist(obsExp_1,100,main = "observed-to-expected CG ratio in 1000 bp window") ## see a long tail, those are CpG islands
abline(v=mean(obsExp_1,na.rm=TRUE),col="red", lwd=1, lty=2)
#####################################################################

ss_2=seq(1, lengths(fa_seq)[2], by=1000)
ss_2=ss_2[-length(ss_2)] ## remove the last one
bat_SL_CoVZXC21=DNAStringSet(fa_seq$MG772934, start=ss_2, end=ss_2+999)
ff_2=alphabetFrequency(bat_SL_CoVZXC21, baseOnly=TRUE)
pCG_2=(ff_2[,"C"]+ff_2[,"G"])/rowSums(ff_2)
hist(pCG_2[pCG_2>0],100, main = "GC content in 1000 bp window")

## CG occurance
nCG_2=vcountPattern("CG", bat_SL_CoVZXC21)
obsExp_2=nCG_2*1000/(ff_2[,"C"]*ff_2[,"G"])
mean(obsExp_2,na.rm=TRUE)
hist(obsExp_2,100,main = "observed-to-expected CG ratio in 1000 bp window") ## see a long tail, those are CpG islands
abline(v=mean(obsExp_2,na.rm=TRUE),col="red", lwd=1, lty=2)
#####################################################################

ss_3=seq(1, lengths(fa_seq)[3], by=1000)
ss_3=ss_3[-length(ss_3)] ## remove the last one
SARS_CoV=DNAStringSet(fa_seq$`NC_004718.3_SARS-CoV`, start=ss_3, end=ss_3+999)
ff_3=alphabetFrequency(SARS_CoV, baseOnly=TRUE)
pCG_3=(ff_3[,"C"]+ff_3[,"G"])/rowSums(ff_3)
hist(pCG_3[pCG_3>0],100, main = "GC content in 1000 bp window")

## CG occurance
nCG_3=vcountPattern("CG", SARS_CoV)
obsExp_3=nCG_3*1000/(ff_3[,"C"]*ff_3[,"G"])
mean(obsExp_3,na.rm=TRUE)
hist(obsExp_3,100,main = "observed-to-expected CG ratio in 1000 bp window") ## see a long tail, those are CpG islands
abline(v=mean(obsExp_3,na.rm=TRUE),col="red", lwd=1, lty=2)
#####################################################################

ss_4=seq(1, lengths(fa_seq)[4], by=1000)
ss_4=ss_4[-length(ss_4)] ## remove the last one
MERS_CoV=DNAStringSet(fa_seq$`NC_038294.1_MERS-CoV`, start=ss_4, end=ss_4+999)
ff_4=alphabetFrequency(MERS_CoV, baseOnly=TRUE)
pCG_4=(ff_4[,"C"]+ff_4[,"G"])/rowSums(ff_4)
hist(pCG_4[pCG_4>0],100, main = "GC content in 1000 bp window")

## CG occurance
nCG_4=vcountPattern("CG", MERS_CoV)
obsExp_4=nCG_4*1000/(ff_4[,"C"]*ff_4[,"G"])
mean(obsExp_4,na.rm=TRUE)
hist(obsExp_4,100,main = "observed-to-expected CG ratio in 1000 bp window") ## see a long tail, those are CpG islands
abline(v=mean(obsExp_4,na.rm=TRUE),col="red", lwd=1, lty=2)
#####################################################################

ss_5=seq(1, lengths(fa_seq)[5], by=1000)
ss_5=ss_5[-length(ss_5)] ## remove the last one
SARS_CoV2=DNAStringSet(fa_seq$`NC_045512.2_SARS-CoV-2`, start=ss_5, end=ss_5+999)
ff_5=alphabetFrequency(SARS_CoV2, baseOnly=TRUE)
pCG_5=(ff_5[,"C"]+ff_5[,"G"])/rowSums(ff_5)
hist(pCG_5[pCG_5>0],100, main = "GC content in 1000 bp window")

## CG occurance
nCG_5=vcountPattern("CG", SARS_CoV2)
obsExp_5=nCG_5*1000/(ff_5[,"C"]*ff_5[,"G"])
mean(obsExp_5,na.rm=TRUE)
hist(obsExp_5,100,main = "observed-to-expected CG ratio in 1000 bp window") ## see a long tail, those are CpG islands
abline(v=mean(obsExp_5,na.rm=TRUE),col="red", lwd=1, lty=2)
#####################################################################

## compare with genome wide distribution of GC content
d1=density(pCG_1)
d2=density(pCG_2)
d3=density(pCG_3)
d4=density(pCG_4)
d5=density(pCG_5)
plot(d5, lwd=2, main="Density Plot on GC Content in 1000 bp Window",
     xlab="Genomic GC Proportion (pCG)",xlim=c(0.3,0.55))
lines(d1, col="red",lwd=2)
lines(d2, col="blue",lwd=2)
lines(d3, col="green",lwd=2)
lines(d4, col="purple",lwd=2)
par(cex = 0.75)
legend("topright", legend=c("bat_SL_CoVZC45", "bat_SL_CoVZXC21","SARS_CoV","MERS_CoV","SARS_CoV2"), 
       lwd=2, col=c("red","blue","green","purple","black"))
## compare with genome wide distribution of GC content
d1_ICpG=density(obsExp_1)
d2_ICpG=density(obsExp_2)
d3_ICpG=density(obsExp_3)
d4_ICpG=density(obsExp_4)
d5_ICpG=density(obsExp_5)
plot(d5_ICpG,col="black", lwd=2,ylim=c(0,3.0), main="Density Plot on Observed-To-Expected CG Ratio in 1000 bp Window",
     xlab="CpG Deficiency (Observed CG/expected CG)")
lines(d1_ICpG, col="red",lwd=2)
lines(d2_ICpG, col="blue",lwd=2)
lines(d3_ICpG, col="green",lwd=2)
lines(d4_ICpG, col="purple",lwd=2)
par(cex = 0.75)
legend("topright", legend=c("bat_SL_CoVZC45", "bat_SL_CoVZXC21","SARS_CoV","MERS_CoV","SARS_CoV2"), 
       lwd=2, col=c("red","blue","green","purple","black"))

#create scatter plot
pCG <- c(mean(pCG_1), mean(pCG_2), mean(pCG_3), mean(pCG_4), mean(pCG_5))
obsExp <- c(mean(obsExp_1), mean(obsExp_2), mean(obsExp_3), mean(obsExp_4), mean(obsExp_5))

plot(x=pCG, y=obsExp)
plot(x=pCG_1, y=obsExp_1, col="red", main="Scatter Plot on Viral Genomic GC Content vs. CpG Deficiency in 1000 bp Window",
     xlab="Genomic GC Proportion (pCG)", ylab="CpG Deficiency (Observed CG/expected CG)",xlim=c(0.3,0.5),ylim=c(0,1.0))
points(x=pCG_5, y=obsExp_5, col="black")
points(x=pCG_2, y=obsExp_2,col="blue")
points(x=pCG_3, y=obsExp_3,col="green")
points(x=pCG_4, y=obsExp_4,col="purple")
points(x=mean(pCG_1),y=mean(obsExp_1),col="red",pch=17,cex=2)
points(x=mean(pCG_4),y=mean(obsExp_4),col="purple",pch=17, cex=1.5)
points(x=mean(pCG_2),y=mean(obsExp_2),col="blue",pch=17, cex=1.5)
points(x=mean(pCG_3),y=mean(obsExp_3),col="green",pch=17, cex=1.5)
points(x=mean(pCG_5),y=mean(obsExp_5),col="black",pch=17, cex=1.5)
par(cex = 0.75)
legend("topleft", legend=c("bat_SL_CoVZC45", "bat_SL_CoVZXC21","SARS_CoV","MERS_CoV","SARS_CoV2"), 
       pch=1, col=c("red","blue","green","purple","black"),bg="transparent")

