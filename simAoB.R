###########################################################################################
###### Simulation 
## TIMING for one chromsome ~ 
## Functions
## Matrix
## Create design matrix for genotype effects from genotypes
## Simulations 
#-> Simulate RNA-seq counts of single gene under ASE
#-> Simulate fpr and power for single individual, some genes with ASE
#-> Simulate ng genes in nInd individuals, whereof nas has alellic imbalance
#-> Simulates pvalue 
#
## NOTE: CHANGE PATHS and DEFINITIONS AS APPROPRIATE
#
#########################################################################################
################## Paths and definitions ##################
### set working directory
setwd("~/Desktop/CEU_FINAL")

############################################################################################
### Load Libraries ####
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

## taken from package remix
rbind.list <- function(l){
  n <- length(l)
  results <- NULL
  for (i in 1:n) {
    results <- rbind(results, l[[i]])
  }
  results
}
qqp<-function(x,ci=TRUE,add=FALSE,ylab="Observed -log10(p-value)",xlab="Expected
              -log10(p-value)",maxLogP,...){
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
############################################################################################
### Create design matrix for genotype effects from genotypes #################
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

## numbers to characters that sorts correctly numerically:
n2c <- function(x,nd){
  y <- paste(x+10^nd)
  substr(y,start=2,stop=4)
}
############################################################################################
### Simulation ##################
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
  sims$obs <- as.character(1:(dim(sims)[1]))
  return(sims)
}
##-> Simulate ng genes in nInd individuals, whereof nas has alellic imbalance
#'
simPop <- function(nInd=50,asePop=0.3,depth=50,nSnp=5, nas=1){
  aseVec <- runif(nInd,0.5-asePop, 0.5+asePop)
  indSim <- lapply(1:nInd,function(i)simInd(ng=50,nase=nas,ind=i,pa=aseVec[i],de=depth,ns=nSnp))
  indSim <- rbind.list(indSim)
  return(indSim)
}
##-> Simulate pvalue
#'
simPopPval <- function(effect,dep,nsnp,nase=1){
  popSim <- simPop(nInd=50,asePop=effect,depth=dep,nSnp=nsnp,nas=nase)
  popSim$obs <- paste0(popSim$ind,"o",popSim$obs)
  popGenes <- names(table(popSim$ge))
  add1model <- function(i){
    as.formula(paste(
      "cbind(c1,t-c1)~r+AC+AG+AT+CG+CT+GT+",
      paste('(as.numeric(ge=="',popGenes[i],'")-1|ind)',sep=""),
      "+(1|obs)"
    ))
  }
  form0 <- as.formula("cbind(c1,t-c1)~r+AC+AG+AT+CG+CT+GT+(1|obs)")
  system.time(logme0 <-glmer(formula=form0,
                             data=popSim,
                             family=binomial(),
                             nAGQ=0L,
                             control=glmerControl(optCtrl=list(maxfun=200000))
  ))
  testGeneAdd1 <- function(i){
    logmeAdd <-glmer(formula=add1model(i),
                     data=popSim,
                     family=binomial(),
                     nAGQ=0L,
                     control=glmerControl(optCtrl=list(maxfun=200000))
    )
    return(logmeAdd)
  }
  
  logmeAdd <- testGeneAdd1(50)
  estVarComp <- as.numeric(summary(logmeAdd)$varcor[[2]])
  pval <- anova(logmeAdd,logme0)[2,8]
  return(c(est=estVarComp,pval=pval))
}
#system.time(print(simPopPval(0.3,50,5)))
#system.time(print(simPopPval(0.3,30,5)))
#system.time(print(simPopPval(0.3,30,3)))
#system.time(print(simPopPval(0.2,50,5)))
#system.time(print(simPopPval(effect=0.3,dep=30,nsnp=3,nase=1)))
#system.time(print(simPopPval(effect=0.3,dep=30,nsnp=3,nase=10)))
#system.time(print(simPopPval(effect=0.3,dep=30,nsnp=3,nase=20)))

##-> Simulate 
simPower <- function(eff){
  pvalVec <- lapply(1:150,function(x)simPopPval(eff,30,4,nase=1))
  pvalVec <- rbind.list(pvalVec)
}
effVec <- seq(0,0.2,by=0.01)
powerSims <- parallel::mclapply(effVec, simPower, mc.cores=20)
save(powerSims,file="5_out/powerSimsPop4ns30ds.Rdata")
powVec <- lapply(powerSims,function(x)sum(x[,2]<0.05/150)/150)
pdf("powerSimsPop4ns30ds.pdf")
plot(effVec+0.5,powVec,col=3,lwd=3,type="l",xlab="Probability of allele 1", ylab="Statistical power",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()

# TIMING ~ 5 HOURS
meanEstVec <-  lapply(powerSims,function(x)mean(x[,1]))
pdf("5_out/meanEstPowerSimsPop4ns30ds.pdf")
plot(effVec+0.5,meanEstVec,col=3,lwd=3,type="l",xlab="Probability of allele 1", ylab="Mean of estimated random effect variance",ps = 16, cex = 1, cex.main = 1)

nullSims <-  parallel::mclapply(1:300,function(x)simPopPval(0,30,4,nase=1),mc.cores=30)
save(nullSims,file="nullSimsPop4ns30ds.Rdata")
nullSims <- rbind.list(nullSims)
pdf("nullSimsPop4ns30ds.pdf")
qqp(nullSims[,2],col=3,cex.lab=1.5,pch=19,cex.main=1.5)
dev.off()

if(FALSE){
  popSim <- simPop(nInd=50,asePop=0.3)
  popSim$obs <- paste0(popSim$ind,"o",popSim$obs)
  library(lme4)
  popGenes <- names(table(popSim$ge))
  add1model <- function(i){
    as.formula(paste(
      "cbind(c1,t-c1)~r+AC+AG+AT+CG+CT+GT+",
      paste('(as.numeric(ge=="',popGenes[i],'")-1|ind)',sep=""),
      "+(1|obs)"
    ))
  }
  form0 <- as.formula("cbind(c1,t-c1)~r+AC+AG+AT+CG+CT+GT+(1|obs)")
  system.time(logme0 <-glmer(formula=form0,
                             data=popSim,
                             family=binomial(),
                             nAGQ=0L,
                             control=glmerControl(optCtrl=list(maxfun=200000))
  ))
  testGeneAdd1 <- function(i){
    logmeAdd <-glmer(formula=add1model(i),
                     data=popSim,
                     family=binomial(),
                     nAGQ=0L,
                     control=glmerControl(optCtrl=list(maxfun=200000))
    )
    return(logmeAdd)
  }
  nGene <- length(popGenes)
  
  add1tests <- parallel::mclapply(1:nGene,testGeneAdd1,mc.cores=10)
  getPopPval <- function(logmeAdd){
    estVarComp <- as.numeric(summary(logmeAdd)$varcor[[2]])
    pval <- anova(logmeAdd,logme0)[2,8]
    return(c(est=estVarComp,pval=pval))
  }
  simEstPval <- rbind.list(lapply(add1tests,getPopPval))
  qqp(simEstPval[1:49,2])
}
############################################################################################
