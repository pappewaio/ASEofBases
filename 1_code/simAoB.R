###########################################################################################
###### Allele specific expression
## TIIMING ~ 1 HOUR 
## Functions
## Matrix
## Create design matrix for genotype effects from genotypes
## Simulations 
##-> Simulate RNA-seq counts of single gene under ASE
##-> Simulate fpr and power for single individual, some genes with ASE
##-> FPR simulations:
## Logistic Regression 
##-> overdispersed logistic regression
##-> logistic regression
##-> logistic mixed effects regression
## Binomial
##-> snpwise binomial tests; correct for reference
##-> snpwise binomial tests correct for ref and type
##-> ROC curve simulations:
##-> POWER simulation
#-> power depth
#-> power SNPs per gene 
##-> Simulating the structure of the genes across the population
##-> Simulate RNA-seq counts of individuals with varying ASE for genes
##-> Simulate singe individual data for fpr and power purposes
##-> Simulate false positives:
##-> estimated obs level random effects
##-> population wide data
#
## NOTE: CHANGE PATHS and DEFINITIONS AS APPROPRIATE
#
#########################################################################################
################## Paths and definitions ##################
### set working directory
setwd("~/Desktop/CEU_FINAL/")
pop<-"CEU"
chr<-22
samplesList<-paste0("geuvadis.",pop,".kG.ind")
indList<-paste0("geuvadis.",pop,".kG.ind")
dataPath <- "~/Desktop/CEU_FINAL/"
plotPath <- "6_plot/"
### AoB
#' AoB.
#'
#' @name AoB
#' @docType package
NULL
############################################################################################
### Load Libraries ####
library(parallel)
library(qvalue)
library(lme4)
############################################################################################
### Functions ##################
ADAMpal <- function() {
  palette(c("#b2182b", "#2166ac", "#4DAF4A", "#FF7F00", "#F781BF","#984EA3"))
}
ADAMpal()
# 1 is (red)
# 2 is (blue)
# 3 is (green)
# 4 is (orange)
# 5 is (pink)
# 6 is (gray) 
############################################################################################
### Matrix ##################
HETS <- c("AC","CA","AG","GA","AT","TA","CG","GC","CT","TC","GT","TG")
NUMHETS <- 1:12
names(NUMHETS) <- HETS

roccurve <- function(stat,posi,clr=1,add=FALSE,...){
  pOrd <- posi[order(stat)]
  rocA <- cumsum(pOrd)/sum(pOrd)
  rocB <- cumsum(1-pOrd)/sum(1-pOrd)
  if(add){
    lines(rocB,rocA,lwd=3,col=clr,xlab="FP ratio", ylab="TP ratio")
  }else{
    plot(rocB,rocA,type="l",lwd=3,col=clr,xlab="FP ratio", ylab="TP ratio",...)
  }
}

qqp<-function(x,ci=TRUE,add=FALSE,ylab="Observed -log10(p-value)",xlab="Expected -log10(p-value)",maxLogP,...){
  x<-x[!is.na(x)]
  if(!missing(maxLogP))
    x[x<10^-maxLogP]<-10^-maxLogP
  N<-length(x)
  x<-sort(x)
  #lambda<-round(median(x)/qchisq(0.5,1),2)
  e<- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
    abline(0,1,col=2,lwd=2)
  }
  #legend("topleft",paste("lambda=",lambda))
  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1) 
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }
}

## numbers to characters that sorts correctly numerically:
n2c <- function(x,nd){
  y <- paste(x+10^nd)
  substr(y,start=2,stop=4)
}

## taken from package remix
rbind.list <- function(l){
  n <- length(l)
  results <- NULL
  for (i in 1:n) {
    results <- rbind(results, l[[i]])
  }
  results
}
#' Create design matrix for genotype effects from genotypes
#'
#' @param g genotypes in 1:12 corresponding to ("AC","CA","AG","GA","AT","TA","CG","GC","CT","TC","GT","TG")
#' @return matrix with 1 row per snp and 6 columns (AC up to TG) with e.g. value 1 if AC and -1 if CA 
gtDesign <- function(g){
  g <- HETS[g]
  z <- cbind((g=="AC")-(g=="CA"),
             (g=="AG")-(g=="GA"),
             (g=="AT")-(g=="TA"),
             (g=="CG")-(g=="GC"),
             (g=="CT")-(g=="TC"),
             (g=="GT")-(g=="TG")) # design matrix for gentype effects
  colnames(z)<-c("AC","AG","AT","CG","CT","GT")
  return(z)
}
############################################################################################
### Simulations ##################
##-> Simulate RNA-seq counts of single gene under ASE
#'
#' @param n number of heterozygous SNPs in gene
#' @param ba ASE log odds of observing allele 1
#' @param br log odds for mapping to reference (vs non-ref) [1]
#' @param bg genotype specific log odds of sequencing base 1 (vs base 2) (fixed such that AC is -CA) [6]
#' @param s overdispersion random effect standard deviation [1]
#' @param d average depth at SNPs
#' @return matrix with columns: reference r, genotype g, total t, allele 1 counts [n,4] 
#' @export
#' @examples
simCountsGene <- function(n, ba, br, bg, s,d){
  ###n=3;ba=log(0.6/0.4);br=log(0.53/0.47);bg=matrix(log(c(0.5,0.56,0.53,0.5,0.5,0.5)/(1-c(0.5,0.56,0.53,0.5,0.5,0.5))));s=0.2;d=50
  r <- sample(x=c(0,1),size=n,replace=TRUE)*2-1 # reference=1,non-reference=-1
  g <- sample(x=HETS,size=n,replace=TRUE) # phased genotype
  z <- cbind((g=="AC")-(g=="CA"),
             (g=="AG")-(g=="GA"),
             (g=="AT")-(g=="TA"),
             (g=="CG")-(g=="GC"),
             (g=="CT")-(g=="TC"),
             (g=="GT")-(g=="TG")) # design matrix for gentype effects
  eta <- ba+br*r+z%*%bg+rnorm(n=n,mean=0,sd=s) # linear predictor
  t <- rpois(n,d) # total depth at each site
  c1 <- rbinom(n,t,plogis(eta)) # counts of reads mapping to allele 1
  return(cbind(r=r,g=NUMHETS[g],t=t,c1=c1))
}

##-> Simulate fpr and power for single individual, some genes with ASE
simInd <- function(ng=300,nase=150,sd=0.3,de=50,ns=5,ind=1,pa=0.7,pr=0.6,pg=c(0.55,0.45,0.55,0.45,0.55,0.45))
{
  ge <- 1:ng
  bre <- qlogis(pr)
  bge <- qlogis(pg)
  bas <- qlogis(pa)
  bav=c(rep(0,ng-nase),rep(bas,nase))
  fun <- function(x){
    res <- simCountsGene(n=ns, ba=bav[x], br=bre, bg=bge, s=sd,d=de)
    return(cbind(ind=ind,ge=ge[x],res))
  }
  sims <- lapply(X=1:ng,FUN=fun)
  sims <- rbind.list(sims)
  rownames(sims) <- NULL
  sims <- cbind(sims, gtDesign(sims[,"g"]))
  sims <- data.frame(sims)
  sims$ge <- n2c(sims$ge,3)
  sims$ind <- as.character(sims$ind)
  sims$obs <- as.factor(1:(dim(sims)[1]))
  return(sims)
}

##-> FPR simulations:
ng=300;nase=150;pa=0.5;
sims <- simInd(ng=ng,nase=nase,pa=pa)

###### Logistic Regression  
##-> overdispersed logistic regression
system.time(fitQ <- glm(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT,data=sims,family=quasibinomial()))
sumfitQ <- summary(fitQ)
pvalsQ <- sumfitQ$coefficients[1:ng,4]
save(pvalsQ,file="5_out/fdrPvalsQ.Rdata")
pdf(file="6_plot/fdrPvalsQ.pdf")
qqp(pvalsQ,col=3,cex.lab=1.5,pch=19,main="Logistic regression with overdispersion",cex.main=1.5)
dev.off()
fun2 <- function(x){
  binom.test(x=sims$c1[x],n=sims$t[x],p=0.5+sims$r[x]*0.1+sum(gtDesign(sims$g[x])*c(0.05,-0.05,0.05,-0.05,0.05,-0.05)))$p.value
}
pvalsE <- lapply(X=1:dim(sims)[1],FUN=fun2)
pvalsE <- rbind.list(pvalsE)
##-> logistic regression
system.time(fitB <- glm(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT,data=sims,family=binomial()))
sumfitB <- summary(fitB)
pvalsB <- sumfitB$coefficients[1:ng,4]
save(pvalsB,file="5_out/fdrPvalsB.Rdata")
pdf(file="6_plot/fdrPvalsB.pdf")
qqp(pvalsB,col=4,cex.lab=1.5,pch=19,main="Logistic regression no overdispersion",cex.main=1.5)
dev.off()
##-> logistic mixed effects regression
#system.time(fit <- glmer(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT+(1|obs),data=sims,family=binomial()))
#sumfit <- summary(fit)
#qqp(sumfit$coefficients[1:ng,4])

##### Binomial #### TIMING ~
##-> snpwise binomial tests; correct for reference
fun2 <- function(x){
  binom.test(x=sims$c1[x],n=sims$t[x],p=0.5+sims$r[x]*0.1)$p.value
}
pvalsS <- lapply(X=1:dim(sims)[1],FUN=fun2)
pvalsS <- rbind.list(pvalsS)
qqp(pvalsS)
if(FALSE){### ROC curve simulations:
  ng=300;nase=150;pa=0.65;ns=5
  sims <- simInd(ng=ng,nase=nase,pa=pa)
  posiGene <- c(rep(0,(ng-nase)),rep(1,(nase))) # is gene ase
  posiSnp <- c(rep(0,(ng-nase)*ns),rep(1,(nase*ns))) # is snp ase
  ### overdispersed logistic regression
  fitQ <- glm(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT,data=sims,family=quasibinomial())
  sumfitQ <- summary(fitQ)
  pvalsQ <- sumfitQ$coefficients[1:ng,4]
  fun2 <- function(x){
    binom.test(x=sims$c1[x],n=sims$t[x],p=0.5+sims$r[x]*0.1+sum(gtDesign(sims$g[x])*c(0.05,-0.05,0.05,-0.05,0.05,-0.05)))$p.value
  }
  pvalsE <- lapply(X=1:dim(sims)[1],FUN=fun2)
  pvalsE <- rbind.list(pvalsE)
  pdf("6_plot/roc.pdf")
  roccurve(pvalsQ,posiGene,clr=3,main="ROC curves",cex.lab=1.5,cex.main=1.5)
  roccurve(pvalsE,posiSnp,add=TRUE,clr=1)
  legend("bottomright",col=c(3,1),lwd=3,legend=c("Overdispersed logistic regression", "SNP-wise binomial testing"),bty='n',cex=1.5)
  dev.off()
  ### Adding fits with missing correction for type and for ref
  fitQnt <- glm(cbind(c1,t-c1)~ge-1+r,data=sims,family=quasibinomial())
  sumfitQnt <- summary(fitQnt)
  pvalsQnt <- sumfitQnt$coefficients[1:ng,4]
  fitQntr <- glm(cbind(c1,t-c1)~ge-1,data=sims,family=quasibinomial())
  sumfitQntr <- summary(fitQntr)
  pvalsQntr <- sumfitQntr$coefficients[1:ng,4]
  fun3 <- function(x){
    binom.test(x=sims$c1[x],n=sims$t[x],p=0.5+sims$r[x]*0.1)$p.value
  }
  pvalsEnt <- lapply(X=1:dim(sims)[1],FUN=fun3)
  pvalsEnt <- rbind.list(pvalsEnt)
  fun4 <- function(x){
    binom.test(x=sims$c1[x],n=sims$t[x],p=0.5)$p.value
  }
  pvalsEntr <- lapply(X=1:dim(sims)[1],FUN=fun4)
  pvalsEntr <- rbind.list(pvalsEntr)
  pdf("6_plot/rocNoTypeRef.pdf")
  roccurve(pvalsQ,posiGene,clr=3,main="ROC curves",cex.lab=1.5,cex.main=1.5)
  roccurve(pvalsQnt,posiGene,add=TRUE,clr=3,lty=2)
  roccurve(pvalsQntr,posiGene,add=TRUE,clr=3,lty=3)
  roccurve(pvalsE,posiSnp,add=TRUE,clr=1)
  roccurve(pvalsEnt,posiSnp,add=TRUE,clr=1,lty=2)
  roccurve(pvalsEntr,posiSnp,add=TRUE,clr=1,lty=3)
  legend("bottomright",col=c(3,1,"black","black","black"),lty=c(1,1,1,2,3),lwd=3,legend=c("Overdispersed logistic regression", "SNP-wise binomial testing","Reference and type corrected","Reference corrected","No correction"),bty='n',cex=1.5)
  dev.off()
}
##-> snpwise binomial tests correct for ref and type
estRef <- (sum(sims$c1[sims$r==1])+sum((sims$t-sims$c1)[sims$r==-1]))/sum(sims$t)
fun3 <- function(x){
  binom.test(x=sims$c1[x],n=sims$t[x],p=0.5+sims$r[x]*0.1+sum(gtDesign(sims$g[x])*c(0.05,-0.05,0.05,-0.05,0.05,-0.05)))$p.value
}
pvalsT <- lapply(X=1:dim(sims)[1],FUN=fun3)
pvalsT <- rbind.list(pvalsT)
save(pvalsT,file="5_out/fdrPvalsT.Rdata")
pdf(file="6_plot/fdrPvalsT.pdf")
qqp(pvalsT,col=1,cex.lab=1.5,pch=19,main="Type- and ref-corrected binomial tests",cex.main=1.5)
dev.off()

##-> ROC curve simulations:
ng=300;nase=150;pa=0.65;ns=5
sims <- simInd(ng=ng,nase=nase,pa=pa)
posiGene <- c(rep(0,(ng-nase)),rep(1,(nase))) # is gene ase
posiSnp <- c(rep(0,(ng-nase)*ns),rep(1,(nase*ns))) # is snp ase
pdf("6_plot/roc.pdf")
roccurve(pvalsQ,posiGene,clr=3,main="ROC curves",cex.lab=1.5,cex.main=1.5)
roccurve(pvalsE,posiSnp,add=TRUE,clr=1)
dev.off()

##-> POWER simulation
pavec <- seq(0.5,0.95,by=0.025)
simPower <- function(pa,depth=50,ns=5){
  alevel <- 0.05
  ng=300;nase=150
  sims <- simInd(ng=ng,nase=nase,pa=pa,de=depth,ns=ns)
  posiGene <- c(rep(0,(ng-nase)),rep(1,(nase))) # is gene ase
  posiSnp <- c(rep(0,(ng-nase)*ns),rep(1,(nase*ns))) # is snp ase
  fitQ <- glm(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT,data=sims,family=quasibinomial())
  sumfitQ <- summary(fitQ)
  pvalsQ <- sumfitQ$coefficients[1:ng,4]
  sum(pvalsQ[as.logical(posiGene)]<alevel/nase)/nase # power
}
#pow <- sapply(X=pavec,FUN=simPower)
#plot(pavec/(1-pavec),pow) # function of odds
#plot(pavec,pow) # function of probs
pow20 <- sapply(X=pavec,FUN=function(x)simPower(x,20))
pow30 <- sapply(X=pavec,FUN=function(x)simPower(x,30))
pow50 <- sapply(X=pavec,FUN=function(x)simPower(x,50))
pow60 <- sapply(X=pavec,FUN=function(x)simPower(x,50))
pow70 <- sapply(X=pavec,FUN=function(x)simPower(x,70))
write(pow20, "6_plot/dpow2.txt")
write(pow30, "6_plot/dpow30.txt")
write(pow50, "6_plot/dpow50.txt")
write(pow60, "6_plot/dpow60.txt")
write(pow70, "6_plot/dpow70.txt")
pdf("6_plot/powerDepth.pdf")
par(mar=c(5,6,4,2)+0.1)
plot(pavec,pow20,col=1,lwd=3,type="l",xlab="Probability of allele 1", ylab="Statistical power",cex.lab=2)
lines(pavec,pow30,col=4,lwd=3)
lines(pavec,pow50,col=3,lwd=3) 
lines(pavec,pow60,col=5,lwd=3) 
lines(pavec,pow70,col=2,lwd=3)
legend("bottomright",col=c(1,4,3,5,2),lwd=3,legend=c("20 reads", "30 reads","50 reads","60 reads","70 reads"),cex=1.5)
dev.off()

pow1 <- sapply(X=pavec,FUN=function(x)simPower(x,50,1))
pow2 <- sapply(X=pavec,FUN=function(x)simPower(x,50,2))
pow3 <- sapply(X=pavec,FUN=function(x)simPower(x,50,3))
pow4 <- sapply(X=pavec,FUN=function(x)simPower(x,50,4))
pow5 <- sapply(X=pavec,FUN=function(x)simPower(x,50,5))
pow6 <- sapply(X=pavec,FUN=function(x)simPower(x,50,6))
pow7 <- sapply(X=pavec,FUN=function(x)simPower(x,50,7))
write(pow1, "6_plot/pow1.txt")
write(pow2, "6_plot/pow2.txt")
write(pow3, "6_plot/pow3.txt")
write(pow4, "6_plot/pow4.txt")
write(pow5, "6_plot/pow5.txt")
write(pow6, "6_plot/pow6.txt")
write(pow7, "6_plot/pow7.txt")
pdf("6_plot/powerSNPsPerGene.pdf")
par(mar=c(5,6,4,2)+0.1)
plot(pavec,pow3,col=1,lwd=3,type="l",xlab="Probability of allele 1", ylab="Statistical power",cex.lab=2)
lines(pavec,pow1,col=6,lwd=3) 
lines(pavec,pow2,col=5,lwd=3) 
lines(pavec,pow4,col=4,lwd=3)
lines(pavec,pow5,col=3,lwd=3) 
lines(pavec,pow6,col=2,lwd=3)
lines(pavec,pow7,col=6,lwd=3)
legend("bottomright",col=c(6,5,1,4,3,2,6),lwd=3,legend=c("1 SNP", "2 SNPs", "3 SNPs", "4 SNPs","5 SNPs","6 SNPs", "7 SNPs"),cex=1.5)
dev.off()

##-> Simulating the structure of the genes across the population
#' @param ni is number of individuals
#' @param ng is number of genes
#' @param ase is prob of seeing allele one due to ase
#' @param sd is standard deviation of overdispersion random effect
#' @param os offset for numbering of genes 
#' @return matrix with row per individual*gene and columns: n number of snps in gene, ba log odds for ase, s random eff std dev, d average gene depth, ind individual id ,ge gene id
simGenesPopV1 <- function(ni,ng,ase=0.5,sd=0.2,os=0,md=30){
  n=rpois(ng*ni,2)+1 # number of snps in genes
  ba=qlogis(rep(ase,ng*ni)) # log odds for ase
  s=rep(rep(sd,ni),each=ng) # overdispersion random effect std
  d=rep(floor(rnorm(ng,mean=md)),times=ni) # average gene depth
  ind=rep(1:ni,each=ng) # individual ids 
  ge=rep(1:ng+os,times=ni) # gene ids
  return(cbind(n=n,ba=ba,s=s,d=d,ind=ind,ge=ge))
}

#br=qlogis(0.6)
#bg=matrix(qlogis(c(0.5,0.56,0.53,0.5,0.5,0.5)))

##-> Simulate RNA-seq counts of individuals with varying ASE for genes
#' @param sg simulated info on genes in the pop [ng*ni,6]
#' @param br log odds for mapping to reference (vs non-ref) [1]
#' @param bg genotype specific log odds of sequencing base 1 (vs base 2) (fixed such that AC is -CA) [6]
#' @return matrix with columns: reference r, genotype g, total t, allele 1 counts, gene g. [n,4] 
#' @export
#' @examples
simCountsPop <- function(sg,br,bg){
  ## ni=3;ng=10;sg=simGenesPopV1(ni,ng);br=qlogis(0.6);bg=matrix(qlogis(c(0.5,0.56,0.53,0.5,0.5,0.5)))
  fun=function(x){
    gene <- simCountsGene(n=sg[x,'n'], ba=sg[x,'ba'], br=br, bg=bg, s=sg[x,'s'],d=sg[x,'d'])
    return(cbind(ind=sg[x,'ind'],ge=sg[x,'ge'],gene))
  }
  res <- sapply(X=1:dim(sg)[1],FUN=fun)
  res <- rbind.list(res)
  rownames(res) <- NULL
  return(res)
}
##-> Simulate singe individual data for fpr and power purposes
#' @param ai allelic imbalance for first 5 genes
simSingleInd <- function(ai=0.5,ng1=5,ng2=45,pr=0.6,pg=c(0.55,0.45,0.55,0.45,0.55,0.45),sd=0.2){
  ###ai=0.5
  ni <- 1
  #    ng1 <- 5
  #    ng2 <- 50
  si <- rbind(simGenesPopV1(ni,ng1,ase=ai,sd=sd),simGenesPopV1(ni,ng2,os=ng1,sd=sd))
  sc <- simCountsPop(sg=si,br=qlogis(pr),bg=matrix(qlogis(pg)))
  sc <- cbind(sc, gtDesign(sc[,"g"]))
  sc <- data.frame(sc)
  sc$ge <- as.character(sc$ge)
  sc$ind <- as.character(sc$ind)
  sc$obs <- as.factor(1:(dim(sc)[1]))
  return(sc)
}
##-> Simulate false positives:
dat <- simSingleInd(ai=0.5)
geneNames <- names(table(dat$ge))
system.time(fit <- glmer(cbind(c1,t-c1)~ge-1+r+AC+AG+AT+CG+CT+GT+(1|obs),data=dat,family=binomial()))
sumfit <- summary(fit)
alevel <- 0.05
