###########################################################################################
### This program will filter the transcriptome data on user-defined filters and will preform binomal and logistic regressions tests: 
# Minimum depth per site
# Minimum depth per gene
# Minimum snps per gene
# Minimum number of individuals with gene
# Minimum number of genes per individual
#
### Contains the following sections:
# Paths and definitions
# Filtering Options
# Functions #update
# Read data
# Filtering #? update
# Filtering Check #? update
#
### Binomal and Logistic Regression Tests:
# downsample to median depth larger than 30 #? what is this doing? #that has a lot of variation in coverage, and down-sample it to 20x to reduce the variation??
# fitting log reg with dispersion for each individual # individual 
# fitting log reg with dispersion for joint for all individuals # population 
# Binomal Tests  
# Estimate reference bias
# Estimate SNP type bias
# Logistics Regression without downsampling depth #?
# Logistics Regression without downsampling #?
## for citation: Storey JD. (2002) A direct approach to false discovery rates. ##Journal of the Royal Statistical Society, Series B, 64: 479-498.
#
# Remove 1 test - see if the gene is significantly imbalances expressed
#
## NOTE: CHANGE PATHS and DEFINITIONS AS APPROPRIATE
##
#########################################################################################

################## Paths and definitions ##################
### set working directory
setwd("/home/user/ASEofBases")
pop<-"CEU" # population 
chr<-0     # chromosome
dataPath <- "/home/user/ASEofBases/4_data/" # path to data files
outPath <- paste0("5_out/chr",chr,"/") # path to output data folder, changed per population and chromosome 
plotPath <- paste0("6_plot/chr",chr,"/") # path to output plot folder, changed per population and chromosome 
samplesList<-paste0("3_raw/geuvadis.",pop,".kG.ind") # list of individuals in the population which is being called 

############################################################################################
### Filtering Options ##################
minDepthSite <- 10 ## Minimum depth per site
minDepthGene <- 10 ## Minimum depth per gene
minSNPsGene <-  2 ## Minimum snps per gene #2
minIndGene <- 5 ## Minimum number of individuals with gene
minGeneInd <- 100 ## Minimum number of genes per individual (consider changing per chr) #200
downSampleQuantile <- 0.75 #? what is this doing 
d=10 # depth
ns=2 # number of snps
############################################################################################
### Load Libraries ####
library(vcd) # Visualization techniques, data sets, summary and inference procedures aimed particularly at categorical data
library(lme4) # Fit linear and generalized linear mixed-effects models.
library(parallel) # handles running much larger chunks of computations in parallel.
library(qvalue) # takes a list of p-values resulting from the simultaneous testing of many hypotheses and estimates their q-values and local FDR values. 
                # The q-value of a test measures the proportion of false positives incurred (called the false discovery rate) when that particular test is called significant. 
                # The local FDR measures the posterior probability the null hypothesis is true given the test's p-value. Various plots are automatically generated, allowing one to make sensible significance cut-offs.
############################################################################################
### Functions ##################
ADAMpal <- function() {
  palette(c("#b2182b", "#2166ac", "#4DAF4A", "#FF7F00", "#F781BF","#984EA3")) #? colors?
}
ADAMpal()

## Function readIndividualData -> makes a table using the individual data files that will be read in
readIndividualData <- function(ind,chr,dataPath){
  filename<-paste(dataPath,"chr",chr,".ind",ind,".data",sep="")
  data<-read.table(filename,as.is=TRUE)
  n<-dim(data)[1] # Dimension, dim(X), is an integer vector giving the number of rows and columns
  colnames(data)<-c("chr","pos","ref","alt","gene","g1","g2","type","tot","A","C","G","T") 
          # type - whether allele 1 or allele 2 is the reference allele (takes the value 1 and 2)
          # tot - total counts of reads overlapping 
  base.lookup<-c(A=1,C=2,G=3,T=4) # assigning numeric value 
  data$h1<-base.lookup[data$g1] # g1, gene copy 1
  data$h2<-base.lookup[data$g2] # g2, gene copy 2
  counts<-t(as.matrix(data[,10:13])) # counting the reads for A, C, G, T 
  data$c1<-counts[4*(0:(n-1))+data$h1] # numbers corresponding to allele the (1 to allele A, 2 to C, 3 to G, 4 to T)
  data$c2<-counts[4*(0:(n-1))+data$h2]
  data$a1a2 <- paste(data$g1,data$g2,sep="") # SNP type (CT, AG, ext..)
  data$AC <- (data$a1a2=="AC")-(data$a1a2=="CA") # SNP types desgn 
  data$AG <- (data$a1a2=="AG")-(data$a1a2=="GA")
  data$AT <- (data$a1a2=="AT")-(data$a1a2=="TA")
  data$CG <- (data$a1a2=="CG")-(data$a1a2=="GC")
  data$CT <- (data$a1a2=="CT")-(data$a1a2=="TC")
  data$GT <- (data$a1a2=="GT")-(data$a1a2=="TG")
  data$type2 <- 2*(data$type-1.5) # ref factor/design 1 or -1 #? more on this 
  data$bothObs <- ((data$c1>0)&(data$c2>0)) # 
  data$tri <- (data$c1+data$c2)!=data$tot # TRUE or FALSE #? more on this 
  data$diff <- (data$c1+data$c2)-data$tot # difference in counts 
  data$ind <- ind # individual 
  return(data) 
}
# Function downSampler -> #? I'm not sure what this is doing 
downSampler <- function(ds,counts){
  if(sum(counts)<=ds){
    return(counts)
  }else{
    down <- sample(x=rep(c(1,2),times=counts),size=ds)
    return(c(sum(down==1),sum(down==2)))
  }
}

############################################################################################
### Read Data ##################
samples <- read.table(file=samplesList, as.is=TRUE) # which indiviaul IDs from the samplesList
nSamples <- dim(samples)[1] # count the number of samples (individuals)
dat <- readIndividualData(samples[1,],chr,dataPath) # makes a table called dat using the function readIndividualData, only calls the first individual 
for(k in 2:nSamples){ # for each individual(sample) run funciton readInividualData which makes a table called dat
  ind <- samples[k,]
  print(ind) # each individual(sample) that is being called
  print(k) # prints how many times is this happening/ which sample/individual is it on now
  dat <- rbind(dat,readIndividualData(ind,chr,dataPath))
}
write.csv(dat, paste0("5_out/",chr,"/chr",chr,"_",pop,"_OrginalData.csv")) # output data set for a reference, may be commented out 

NumOfSites_var<- paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumOfSites.txt") # making a variable for the defined pathway 
dimdat_var<-paste0("5_out/",chr,"/dimdat.txt") # making a variable for the defined pathway 
NumOfSites<-print("Number of sites:")
write.table(NumOfSites, NumOfSites_var, sep="\t") # write table 
dimdat<-print(dim(dat)) # dimdat, dimension of data (observed number of variables)
write.table(dimdat, dimdat_var, sep="\t")
file.append(NumOfSites_var, dimdat_var) # outputting title(number of sites) and the dimdat into one text file. 
file.remove(dimdat_var) # removing the no longer needed dimdat text file since it has been appended to the number of sites text file 

NumOfGenesWithCov_var<- paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumOfGenesWithCov.txt") # making a variable for the defined pathway 
datagene_var<-paste0("5_out/",chr,"/datagene.txt") # making a variable for the defined pathway 
NumOfGenesWithCov<-print("Number of genes with coverage:")
write.table(NumOfGenesWithCov, NumOfGenesWithCov_var, sep="\t") # write table 
datagene<-print(length(table(dat$gene))) # number of unique genes 
write.table(datagene, datagene_var, sep="\t")
file.append(NumOfGenesWithCov_var, datagene_var) # outputting title(number of genes with coverage) and the datagene into one text file. 
file.remove(datagene_var) # removing the no longer needed datagene_var text file since it has been appended to the number of genes with coverage text file 

MedianCov_var<- paste0("5_out/",chr,"/chr",chr,"_",pop,"_MedianCoverage.txt") # making a variable for the defined pathway 
Median_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_Median.txt") # making a variable for the defined pathway 
MedianCov<-print("Median coverage:")
write.table(MedianCov, MedianCov_var, sep="\t")
Median<-print(median(dat$c1+dat$c2)) #dat$c1: c1, count of reads that are allele 1     dat$c2: c2, count of reads that are allele 2
write.table(Median, Median_var, sep="\t")
file.append(MedianCov_var, Median_var) # outputting title(median coverage) and the median into one text file. 
file.remove(Median_var) # removing the no longer needed Median_var text file since it has been appended to the number of median coverage text file 

NumSites_minDepthGene<-print(paste("Number of sites minDepthGene>=",minDepthGene))
NumSites_minDepthGene_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumSites_minDepthGene>=.txt")
write.table(NumSites_minDepthGene, NumSites_minDepthGene_var, sep="\t")
minDepthGene_length<-print(with(subset(dat, c1+c2 >= minDepthGene, drop=TRUE),length(type2))) # number of sites with depth coverage equal to or greater than user defined depth d= 
minDepthGene_length_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_minDepthGene_length.txt")
write.table(minDepthGene_length, minDepthGene_length_var, sep="\t")
file.append(NumSites_minDepthGene_var, minDepthGene_length_var)
file.remove(minDepthGene_length_var)

NumSNPsPerGeneWith_minDepthGene_greater<-print(paste("Number of SNPs per gene*ind with minDepthGene>",minDepthGene))
NumSNPsPerGeneWith_minDepthGene_greater_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumSNPsPerGeneWith_minDepthGene_greater.txt")
write.table(NumSNPsPerGeneWith_minDepthGene_greater, NumSNPsPerGeneWith_minDepthGene_greater_var, sep="\t")
NumberOfSNPsPerGeneIndWith_minDepthGene_greater<-print(with(subset(dat, c1+c2 >= minDepthGene, drop=TRUE), table(table(gene,ind)))) # number of SNPs per gene with depth equal to or greater than user defined depth d=
NumberOfSNPsPerGeneIndWith_minDepthGene_greater_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumberOfSNPsPerGeneIndWith_minDepthGene_greater.txt")
write.table(NumberOfSNPsPerGeneIndWith_minDepthGene_greater, NumberOfSNPsPerGeneIndWith_minDepthGene_greater_var, sep="\t")
file.append(NumSNPsPerGeneWith_minDepthGene_greater_var, NumberOfSNPsPerGeneIndWith_minDepthGene_greater_var)
file.remove(NumberOfSNPsPerGeneIndWith_minDepthGene_greater_var)

NumGenesWith_nSNP_greaterequal_minSNPsGene<-print(paste("Number of genes*inds with nSNP>=",minSNPsGene))
NumGenesWith_nSNP_greaterequal_minSNPsGene_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumGenesWith_nSNP_greaterequal_minSNPsGene.txt")
write.table(NumGenesWith_nSNP_greaterequal_minSNPsGene, NumGenesWith_nSNP_greaterequal_minSNPsGene_var, sep="\t")
NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE<-print(with(subset(dat, c1+c2 >= minDepthGene, drop=TRUE),sum(table(table(gene,ind))[-c(1:minSNPsGene)]))) # number of genes with equal to or greater than the number of SNPs per gene which is user defined ns= 
NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE.txt")
write.table(NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE, NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE_var, sep="\t")
file.append(NumGenesWith_nSNP_greaterequal_minSNPsGene_var, NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE_var)
file.remove(NumGenesWith_nSNP_greaterequal_minSNPsGeneTABLE_var)

############################################################################################
### Matrix ##################
## taken from package remix
# remix provides remix, a quick and easy function for describing datasets. It can be view as a mix of cast (in package reshape) 
rbind.list <- function(l){ # what is this section doing? 
  n <- length(l)
  results <- NULL
  for (i in 1:n) {
    results <- rbind(results, l[[i]])
  }
  results
}

HETS <- c("AC","CA","AG","GA","AT","TA","CG","GC","CT","TC","GT","TG") # 12 differnt possible combinations of the 4 bases A,T,C,G 
NUMHETS <- 1:12 # 12 different combinations AC is 1, CA is 2 (see above)
names(NUMHETS) <- HETS

#' Create design matrix for genotype effects from genotypes
#'
#' @param g genotypes in 1:12 corresponding to ("AC","CA","AG","GA","AT","TA","CG","GC","CT","TC","GT","TG") # @ is an object
#' @return matrix with 1 row per snp and 6 columns (AC up to TG) with e.g. value 1 if AC and -1 if CA 
gtDesign <- function(g){
  g <- HETS[g] #? what is this doing? 
  z <- cbind((g=="AC")-(g=="CA"),
             (g=="AG")-(g=="GA"),
             (g=="AT")-(g=="TA"),
             (g=="CG")-(g=="GC"),
             (g=="CT")-(g=="TC"),
             (g=="GT")-(g=="TG")) # design matrix for gentype effects
  colnames(z)<-c("AC","AG","AT","CG","CT","GT")
  return(z)
}

gtDesign2 <- function(g){ # function of type in AC, AG, ...
  z <- cbind((g=="AC")-(g=="CA"),
             (g=="AG")-(g=="GA"),
             (g=="AT")-(g=="TA"),
             (g=="CG")-(g=="GC"),
             (g=="CT")-(g=="TC"),
             (g=="GT")-(g=="TG")) # design matrix for gentype effects
  colnames(z)<-c("AC","AG","AT","CG","CT","GT")
  return(z)
}

qqp<-function(x,ci=TRUE,add=FALSE,ylab="Observed -log10(p-value)",xlab="Expected -log10(p-value)",maxLogP,...){ # expected is no allelic imbalance 
  x<-x[!is.na(x)]
  if(!missing(maxLogP))
    x[x<10^-maxLogP]<-10^-maxLogP
  N<-length(x)
  x<-sort(x)
  #lambda<-round(median(x)/qchisq(0.5,1),2) # commented out from original code #? what was this for? 
  e<- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
    abline(0,1,col=2,lwd=2) # add a line
  }
 # legend("topleft",paste("lambda=",lambda)) # commented out from original code 
  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1) 
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }
}
############################################################################################
### Filtering ##################
## Generate plots to visually investigate results of filtering
dat$deep <- (dat$c1+dat$c2)>(minDepthSite-1) # depth = total allele coverage > minimum depth at that site -1 (minDepthSite was set at the beginning of the script and is user defined)
    # creates a new column in the dat of TRUE FALSE. TRUE if depth is greater than minimum depth and FALSE if depth coverage is less than minimum depth 

##  Output text files 
# Number of genes with coverage larger than minimun depth
NumGenesCovLargerThanMinDepthGene_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumGenesCovLargerThanMinDepthGene.txt") # making a variable for the defined pathway 
dataGeneDataDeep_var<-paste0("5_out/",chr,"/dataGeneDataDeep.txt") # making a variable for the defined pathway 
NumGenesCovLargerThanMinDepthGene<-print("Number of genes with coverage larger than minDepthGene:")
write.table(NumGenesCovLargerThanMinDepthGene, NumGenesCovLargerThanMinDepthGene_var, sep="\t") # write table
dataGeneDataDeep<-print(length(table(dat$gene[dat$deep]))) # number of genes with depth coverage greater than minimum
write.table(dataGeneDataDeep, dataGeneDataDeep_var, sep="\t") # write table
file.append(NumGenesCovLargerThanMinDepthGene_var, dataGeneDataDeep_var) # outputting title(Number of genes with coverage larger than minDepthGene) and the dataGeneDataDeep into one text file. 
file.remove(dataGeneDataDeep_var) # removing the no longer needed dataGeneDataDeep_var text file

## How many reads per individual per gene
readsTable_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumOfReadsPerIndividualPerGene.txt")
readsTable <- xtabs(tot~ind+gene,data=dat[dat$deep,]) # How many reads per individual per gene
write.table(readsTable, readsTable_var, sep="\t") # export readsTable

## How many SNPs per individual per gene
snpTable_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumOfSNPsPerIndividualPerGene.txt")
snpTable <- with(dat[dat$deep,],t(table(gene, ind)))
write.table(snpTable, snpTable_var, sep="\t") # export snpTable

## Post DepthGene filtering
hist(colSums(readsTable>0),main="Frequency of number of individuals 
     with gene present",xlab="Number of individuals", ylab="Number of genes", col=c("purple"))
hist(colSums(readsTable>minDepthGene),main="Frequency of number of individuals 
     with gene present",xlab="Number of individuals", ylab="Number of genes", col=c("yellow"), add =TRUE)
legend("top", c("gene present before DepthGene filter", "gene present after DepthGene filter"), fill=c("purple","yellow"), bty = "n")
DepthGenehistPath_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_GenePresentDepthGeneFilter.jpg")
#dev.copy(jpeg, DepthGenehistPath_var) # option to save histogram plot 
#dev.off()

## Post SNPsGene filtering
hist(colSums(snpTable>0),main="Frequency of number of individuals with gene present",xlab="Number of individuals", ylab="Number of genes", col=c("purple"))
hist(colSums(snpTable>minSNPsGene),main="Frequency of number of individuals with gene present", xlab="Number of individuals", ylab="Number of genes", col=c("yellow"), add=TRUE)
legend("top", c("gene present before SNPsGene filter", "gene present after SNPsGene filter"), fill=c("purple","yellow"), bty = "n")
SNPsGenehistPath_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_GenePresentSNPsGeneFilter.jpg")
#dev.copy(jpeg, SNPsGenehistPath_var) # opition to save histogram plot
#dev.off()

sum(colSums(readsTable>minDepthGene)>minIndGene) # minIndGene <- () Minimum number of individuals with gene, user defined at beginning of script
                                                 # minDepthGene <- () Minimum depth per gene, user defined at beginning of script 
                                                 # readsTable <- How many reads per individual per gene
                                                 #? colSums - function 
                                                 # = total number of genes that have at least the minimum number of individuals with gene present and at least the minimum depth per gene (taken from the reads table)
sum(colSums(snpTable>minSNPsGene)>minIndGene/2) # minIndGene/2 <- Minimum number of individuals with gene, user defined at beginning of script dividided by 2 #?why is it divided by 2? 
                                                # minSNPsGene <- minimum number of SNPs per gene 
                                                # snpTable <- table of genes and individuals and the number of SNPs per gene per individual 
sum((colSums(readsTable>minDepthGene)>minIndGene)&(colSums(snpTable>minSNPsGene)>minIndGene/2)) # overlapping genes? why & and not +?

# Genes with appropriate coverage
selectedGenes <- dimnames(readsTable)[[2]][colSums(readsTable>minDepthGene)>minIndGene&(colSums(snpTable>minSNPsGene)>minIndGene/2)] 
write.csv(selectedGenes, paste0("5_out/",chr,"/chr",chr,"_",pop,"_SelectedGenes.csv"))

# Add indicator for selected genes
dat$select <- with(dat,(gene %in% selectedGenes)) # adds to the dat table, select = True if gene is selected, False if it is not selected 
tmptab <- with(dat[dat$select & dat$deep,],table(gene,ind)) # temporary table of genes (row) and individuals (columns) and the depth coverage

# In how many individuals does genes occur
tmptabInd <- table(colSums(tmptab>0)) # genes per individual
GenesPerInd<- plot(as.numeric(dimnames(tmptabInd)[[1]]),as.numeric(tmptabInd),xlab="Genes per individual",ylab="Individuals", type='h')
GenesPerInd_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_GenesPerInd.jpg")
#dev.copy(jpeg, GenesPerInd_var)
#dev.off()

# How many genes are observed in the individuals
tmptabGene <- table(rowSums(tmptab>0)) # individuals per gene
IndPerGene<-plot(as.numeric(dimnames(tmptabGene)[[1]]),as.numeric(tmptabGene),xlab="Individuals per gene",ylab="Genes", type='h')
IndPerGene_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_IndPerGene.jpg")
#dev.copy(jpeg, IndPerGene_var)
#dev.off()

# log depth per site 
with(data=dat[dat$deep & dat$select, ],hist(log(c1+c2),main="Log depth per site", xlab="Log(depth)")) 
LogDepthPerSite_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_LogDepthPerSite.jpg")
#dev.copy(jpeg, LogDepthPerSite_var)
#dev.off()

############################################################################################
## Filtering and Binomal & Logistic Regression Tests ################## 
##-> Filtering gene to have more than:
# minSNPsGene
# min depth reads per snp
# downsample to median depth larger than 30
# fitting log reg with dispersion for each individual
# fitting log reg with dispersion for joint for all individuals

################
### add of gene names is happening here. 
# ? what is this section doing?
ind2pvals <- function(indi,downSNPs=FALSE,downDepth=TRUE){
  datSingle <- subset(dat, ind==samples[indi,] & c1+c2>=d)
  useGenes <- with(datSingle, rownames(table(gene))[table(gene)>=ns])
  if(downSNPs){# downsample to minSNPsGene SNPs per gene 
    fun1 <- function(ge){
      datGe <-subset(datSingle, gene==ge)
      nGe <- dim(datGe)[1]
      if(nGe >minSNPsGene){
        datGe <- datGe[sample(nGe,minSNPsGene),]
      }
      return(datGe)
    }
    datDown <- lapply(useGenes,fun1)
    datDown <- rbind.list(datDown)
  }else{# no dowsampling of SNPs
    datDown <- subset(datSingle, gene %in% useGenes)
  }
  if(downDepth){ # downsample to median depth larger than d per SNP
    m <- median(datDown$c1+datDown$c2)
    ## change m to d below to downsample allele to same depth
    e1e2 <- apply(X=cbind(datDown$c1,datDown$c2),1,function(x)downSampler(m,x)) 
    datDown$e1 <- e1e2[1,]
    datDown$e2 <- e1e2[2,]
    ## fit overdispersed logistic regression
    system.time(fitE <- glm(cbind(e1,e2)~gene-1+type2+AC+AG+AT+CG+CT+GT,data=datDown,family=quasibinomial()))
  }else{ # no downsamcple of SNPs per gene
    system.time(fitE <- glm(cbind(c1,c2)~gene-1+type2+AC+AG+AT+CG+CT+GT,data=datDown,family=quasibinomial()))
  }
  sumfitE <- summary(fitE)
  pvalsE <- sumfitE$coefficients[1:length(useGenes),4]
  return(data.frame(ind=samples[indi,],pvals=pvalsE,disp=sumfitE$dispersion,gene=names(pvalsE),stringsAsFactors=FALSE))
}
tryInd2pvals<-function(indi,downSNPs=FALSE,downDepth=TRUE){
  r<-try(ind2pvals(indi,downSNPs,downDepth),TRUE)
  r
}
## test while downsampling to median depth:
set.seed=0 # ? why is the set.seed=0 ?
popindPvals <-parallel::mclapply(X=1:89,tryInd2pvals,mc.cores=20) # parallel analysis of the individuals to speed up process. # changed from X=1:78 to X=1:89 beucase there are now 89 inividuals 

# Logistic Regression for each individual independently text, csv, and QQ plot 
table(l<-sapply(popindPvals,length)) # ? is this informative? should I print to a txt file?
popindPvals<-popindPvals[l==length(popindPvals[[1]])]
popindPvals <- rbind.list(popindPvals)
nPvals <- length(popindPvals$pvals)
write.table(popindPvals[order(popindPvals$pvals),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_LogRegIndividualPvals_",minSNPsGene,"minSNPsGene_",minDepthGene,"minDepthGene.txt"),quote=FALSE,row.names=FALSE)
write.csv(popindPvals[order(popindPvals$pvals),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_LogRegIndividualPvals_",minSNPsGene,"minSNPsGene_",minDepthGene,"minDepthGene.csv"),quote=FALSE,row.names=FALSE)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQLogReg_IndividualPvals_",minSNPsGene,"minSNPsGene_",minDepthGene,"minDepthGene.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(popindPvals$pvals,col=3,pch=20,cex.lab=1,cex=1)
dev.off()

#############################################################
##->fitting log reg with dispersion for joint for all individuals
## population wide 
## only keeping genes were select is TRUE 
filterInd <- function(indi,downSNPs=FALSE){
  datSingle <- subset(dat, ind==samples[indi,] & c1+c2>=d)
  useGenes <- with(datSingle, rownames(table(gene))[table(gene)>=ns])
  if(downSNPs){
    fun1 <- function(ge){
      datGe <-subset(datSingle, gene==ge)
      nGe <- dim(datGe)[1]
      if(nGe >ns){
        datGe <- datGe[sample(nGe,ns),]
      }
      return(datGe)
    }
    ## downsample to ns snps per gene
    datDown <- lapply(useGenes,fun1)
    datDown <- rbind.list(datDown)
  }else{
    datDown <- subset(datSingle, gene %in% useGenes)
  }
  ## downsample to median of depths larger than minDepthGene
  m <- median(datDown$c1+datDown$c2)
  ## change m to minDepthGene below to downsample allele to same depth
  e1e2 <- apply(X=cbind(datDown$c1,datDown$c2),1,function(x)downSampler(m,x)) 
  datDown$e1 <- e1e2[1,]
  datDown$e2 <- e1e2[2,]
  return(datDown)
}
tryFilterInd<-function(indi,downSNPs=FALSE){
  r<-try(filterInd(indi,downSNPs),TRUE)
  r
}

# only keeping genes were select is TRUE 
# Logistic Regression for joint of all individuals in the population  
set.seed=0
popfilterInd <-parallel::mclapply(X=1:89,tryFilterInd,mc.cores=20) 
table(l<-sapply(popfilterInd,length))
popfilterInd<-popfilterInd[l==length(popfilterInd[[1]])]
popfilterInd <- rbind.list(popfilterInd)
popfilterInd$geneInd <- with(popfilterInd,paste0(gene,ind))
#write.table(popfilterInd[order(popfilterInd$geneInd),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_FilterInd",minSNPsGene,"minSNPsGene",d,"d.txt"),quote=FALSE,row.names=FALSE)
write.csv(popfilterInd[order(popfilterInd$geneInd),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_FilterInd",minSNPsGene,"minSNPsGene",d,"d.csv"),quote=FALSE,row.names=FALSE)
################## 
##-> Binomial tests
#Estimate reference bias
#Estimate SNP type bias
estBias <- function(e1,e2,type2,a1a2){
  #### estimate reference bias:
  qest <-  (sum(e1[type2==1])+sum(e2[type2==-1]))/(sum(e2)+sum(e1))
  #### estimate SNP type bias:
  bAC <- (sum(e1[a1a2=="AC"])+sum(e2[a1a2=="CA"]))/((sum(e1[a1a2=="AC"])+sum(e2[a1a2=="CA"]))+(sum(e2[a1a2=="AC"])+sum(e1[a1a2=="CA"])))
  bAG <- (sum(e1[a1a2=="AG"])+sum(e2[a1a2=="GA"]))/((sum(e1[a1a2=="AG"])+sum(e2[a1a2=="GA"]))+(sum(e2[a1a2=="AG"])+sum(e1[a1a2=="GA"])))
  bAT <- (sum(e1[a1a2=="AT"])+sum(e2[a1a2=="TA"]))/((sum(e1[a1a2=="AT"])+sum(e2[a1a2=="TA"]))+(sum(e2[a1a2=="AT"])+sum(e1[a1a2=="TA"])))
  bCG <- (sum(e1[a1a2=="CG"])+sum(e2[a1a2=="GC"]))/((sum(e1[a1a2=="CG"])+sum(e2[a1a2=="GC"]))+(sum(e2[a1a2=="CG"])+sum(e1[a1a2=="GC"])))
  bCT <- (sum(e1[a1a2=="CT"])+sum(e2[a1a2=="TC"]))/((sum(e1[a1a2=="CT"])+sum(e2[a1a2=="TC"]))+(sum(e2[a1a2=="CT"])+sum(e1[a1a2=="TC"])))
  bGT <- (sum(e1[a1a2=="GT"])+sum(e2[a1a2=="TG"]))/((sum(e1[a1a2=="GT"])+sum(e2[a1a2=="TG"]))+(sum(e2[a1a2=="GT"])+sum(e1[a1a2=="TG"])))
  return(list(q=qest,b=c(bAC,bAG,bAT,bCG,bCT,bGT)))
}
biasEst <- with(popfilterInd,estBias(e1,e2,type2,a1a2))
fun2 <- function(x,dat){
  with(dat,binom.test(x=c1[x],n=c1[x]+c2[x],p=0.5+type2[x]*(biasEst$q-0.5)+sum(gtDesign2(a1a2[x])*(biasEst$b-0.5)))$p.value)
}
pvalsBinom <- lapply(X=1:dim(popfilterInd)[1],FUN=fun2,dat=popfilterInd)
pvalsBinom <- rbind.list(pvalsBinom)
pvalsBinom2 <- pvalsBinom
pvalsBinom2[pvalsBinom<1e-100] <- 1e-100
qqp(pvalsBinom2)

str(pvalsBinom2)
str(popfilterInd)
minBinomPvalueGeneNoDS <- tapply(pvalsBinom2,popfilterInd$geneInd,min)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQminBinomAndIndPvalsNoDS_",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(minBinomPvalueGeneNoDS,ylim=c(0,55),col=1,cex.lab=1.5,pch=20,cex=1)
qqp(popindPvals$pvals,add=TRUE,col=3,pch=19)
dev.off()

randBinomPvalueGeneNoDS <- tapply(pvalsBinom2,popfilterInd$geneInd,function(x)sample(x,size=1))
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQrandBinomAndIndPvalsNoDS_",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(randBinomPvalueGeneNoDS,ylim=c(0,55),col=1,cex.lab=1.5,pch=19,cex=1)
qqp(popindPvals$pvals,add=TRUE,col=3,pch=19)
dev.off()

### What is this function doing?
fun3 <- function(x,dat){
  with(dat,binom.test(x=e1[x],n=e1[x]+e2[x],p=0.5+type2[x]*(biasEst$q-0.5)+sum(gtDesign2(a1a2[x])*(biasEst$b-0.5)))$p.value)
}
pvalsBinomDS <- lapply(X=1:dim(popfilterInd)[1],FUN=fun3,dat=popfilterInd)
pvalsBinomDS <- rbind.list(pvalsBinomDS)
qqp(pvalsBinomDS)

minBinomPvalueGene <- tapply(pvalsBinomDS,popfilterInd$geneInd,min)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQminBinomAndIndPvals_",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(minBinomPvalueGene,ylim=c(0,55),col=1,cex.lab=1.5,pch=19,cex=1)
qqp(popindPvals$pvals,add=TRUE,col=3,pch=19)
dev.off()

randBinomPvalueGene <- tapply(pvalsBinomDS,popfilterInd$geneInd,function(x)sample(x,size=1))
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQrandBinomAndIndPvals_",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(randBinomPvalueGene,ylim=c(0,50),col=1,cex.lab=1.5,pch=19,cex=1)
qqp(popindPvals$pvals,add=TRUE,col=3,pch=19)
dev.off()
popindPvals_varCSV<- paste0("5_out/",chr,"/chr",chr,"_",pop,"_PopIndividualPvals.csv")
write.csv(popindPvals, popindPvals_varCSV)
popindPvals_varTxt<- paste0("5_out/",chr,"/chr",chr,"_",pop,"_PopIndividualPvals.txt")
write.table(popindPvals, popindPvals_varTxt)

randBinomPvalueGene_varCSV<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_RandBinomPvalueGene.csv")
randBinomPvalueGen_varTxt<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_RandBinomPvalueGene.txt")
write.csv(randBinomPvalueGene, randBinomPvalueGene_varCSV)
write.table(randBinomPvalueGene, randBinomPvalueGen_varTxt)

################## 
##-> logistic regression without downsampling of depth
set.seed=0
popindPvalsNoDS <-parallel::mclapply(X=1:89,tryInd2pvals,mc.cores=20,downDepth=FALSE)
table(l<-sapply(popindPvalsNoDS,length))
popindPvalsNoDS<-popindPvalsNoDS[l==length(popindPvalsNoDS[[1]])]
popindPvalsNoDS <- rbind.list(popindPvalsNoDS)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_LogRegQQindPvalsNoDS",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(popindPvalsNoDS$pvals,col=3,cex.lab=1.5,pch=19,cex=1)
dev.off()

qqp(popindPvals$pvals,col=3,cex.lab=2,pch=19)
qqp(popindPvalsNoDS$pvals,col="grey",add=TRUE)

###-> logistics Regression without downsampling
minBinomPvalueGeneNoDS <- tapply(pvalsBinom2,popfilterInd$geneInd,min)
minBinomPvalueGeneNoDS_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_MinBinomPvalueGeneNoDS.csv")
write.csv(minBinomPvalueGeneNoDS, minBinomPvalueGeneNoDS_var)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQ_LogReg_MinBinom_IndividualPvals_NoDSBoth_",minSNPsGene,"minSNPsGene",d,"d.pdf"))
par(mar=c(5.1,5.1,5.1,1.1))
qqp(minBinomPvalueGeneNoDS,ylim=c(0,55),col=1,cex.lab=1.5,pch=19,cex=1)
qqp(popindPvalsNoDS$pvals,add=TRUE,col=3,pch=19)
dev.off()
popindPvalsNoDSGene_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_PopIndividualPvalsNoDS_Gene.csv")
write.csv(popindPvalsNoDS$gene, popindPvalsNoDSGene_var)

# Filtering opitions 
#ns=4 # changed from 4 to 3
#d=30 # from 30 to 20 

tmpGenes <- with(subset(dat, c1+c2 >= minDepthGene, drop=TRUE),names(table(gene))) # changed d to minDepthGene. >= d, minDepthGene
tmpGenes_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_TmpGenes.csv")
write.csv(tmpGenes, tmpGenes_var)
tmpIndPerGene <- with(subset(dat, c1+c2 >= minDepthGene, drop=TRUE),tapply(ind,gene,table)) #changed d to minDepthGene.  >= d, minDepthGene

# how many inds with more than 3 snps
NumOfInd_MoreThan3SNP<-print("Number of Individuals with more than 3 SNPs")
NumOfInd_MoreThan3SNP_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"NumOfIndividualsWithMoreThan",minSNPsGene,"SNPsPerGene.txt")
write.table(NumOfInd_MoreThan3SNP, NumOfInd_MoreThan3SNP_var)
NumOfInd_MoreThan3SNPTABLE<-table(sapply(tmpIndPerGene,function(x)sum(x>3)))
NumOfInd_MoreThan3SNP_varTable<-paste0("5_out/",chr,"/chr",chr,"_",pop,"NumOfIndividualsWithMoreThan",minSNPsGene,"SNPsPerGeneTABLE.txt")
write.table(NumOfInd_MoreThan3SNPTABLE, NumOfInd_MoreThan3SNP_varTable, sep="\t")
file.append(NumOfInd_MoreThan3SNP_var, NumOfInd_MoreThan3SNP_varTable)
file.remove(NumOfInd_MoreThan3SNP_varTable)

# how many genes with at least 3 inds with more than 3 snp
NumGenes_withMin3inds_More3snps<-print("how many genes with at least 3 inds with more than 3 snp")
NumGenes_withMin3inds_More3snps_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumGenes_withMin_3_IndividualsAndAtLeast_3_SNPs.txt")
write.table(NumGenes_withMin3inds_More3snps, NumGenes_withMin3inds_More3snps_var, sep="\t")
NumGenes_withMin3inds_More3snpsSUM<-sum(table(sapply(tmpIndPerGene,function(x)sum(x>3)))[-c(1:3)])
NumGenes_withMin3inds_More3snpsSUM_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_NumGenes_withMin_3_IndividualsAndAtLeast_3_SNPsSUM.txt")
write.table(NumGenes_withMin3inds_More3snpsSUM, NumGenes_withMin3inds_More3snpsSUM_var, sep="\t")
file.append(NumGenes_withMin3inds_More3snps_var, NumGenes_withMin3inds_More3snpsSUM_var)
file.remove(NumGenes_withMin3inds_More3snpsSUM_var)

# gene names of those genes:
popWideGenes <- names(sapply(tmpIndPerGene,function(x)sum(x>3))[sapply(tmpIndPerGene,function(x)sum(x>3))>2])
## These are the genes that will be tested!
popWideGenes_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_PopulationWideGenes.txt")
write.table(popWideGenes, popWideGenes_var, sep="\t")


## Filter and fit the data:
## min d reads pr snp and downsample to median depth larger than 30, 
filterIndPop <- function(indi,downSNPs=FALSE){
  datSingle <- subset(dat, ind==samples[indi,] & c1+c2>=minDepthGene) 
  ##    useGenes <- with(datSingle, rownames(table(gene)))
  useGenes <- with(datSingle, rownames(table(gene))[table(gene)>=minDepthSite])
  #    useGenes <- useGenes[useGenes %in% popWideGenes]
  if(downSNPs){
    fun1 <- function(ge){
      datGe <-subset(datSingle, gene==ge)
      nGe <- dim(datGe)[1]
      if(nGe>minDepthSite){
        datGe <- datGe[sample(nGe,minDepthSite),]
      }
      return(datGe)
    }
    ## downsample to minSNPsGene snps per gene
    datDown <- lapply(useGenes,fun1)
    datDown <- rbind.list(datDown)
  }else{
    datDown <- subset(datSingle, gene %in% useGenes)
  }
  ## downsample to median of depths larger than minDepthGene
  m <- median(datDown$c1+datDown$c2)
  ## change m to minDepthGene below to downsample allele to same depth
  e1e2 <- apply(X=cbind(datDown$c1,datDown$c2),1,function(x)downSampler(m,x)) 
  datDown$e1 <- e1e2[1,]
  datDown$e2 <- e1e2[2,]
  return(datDown)
}

tryFilterIndPop<-function(indi,downSNPs=FALSE){
  r<-try(filterIndPop(indi,downSNPs),TRUE)
  r
}

# Population 
set.seed=0
popfilterPop <-parallel::mclapply(X=1:89,tryFilterIndPop,mc.cores=20)
table(l<-sapply(popfilterPop,length))
popfilterPop<-popfilterPop[l==length(popfilterPop[[1]])]
popfilterPop <- rbind.list(popfilterPop)
popfilterPop$geneInd <- with(popfilterPop,paste0(gene,ind))
popsubGenes <- subset(popfilterPop,gene %in% popWideGenes)
popsubGenes$obs <- 1:dim(popsubGenes)[1] 

nGene <- length(popWideGenes)

formGenes <- function(genes){
  as.formula(
    paste(
      "cbind(c1,c2)~type2+AC+AG+AT+CG+CT+GT+",
      paste('(as.numeric(gene=="',genes,'")-1|ind)',
            sep="",
            collapse="+"),
      "+(1|obs)"
    ))    
}

# Individuals 
set.seed=0
popfilterInd <-parallel::mclapply(X=1:89,tryFilterInd,mc.cores=20) 
table(l<-sapply(popfilterInd,length))
popfilterInd<-popfilterInd[l==length(popfilterInd[[1]])]
popfilterInd <- rbind.list(popfilterInd)
popfilterInd$geneInd <- with(popfilterInd,paste0(gene,ind))
write.table(popfilterInd[order(popfilterInd$geneInd),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_filterIndividuals",minSNPsGene,"minSNPsGene",d,"d.txt"),quote=FALSE,row.names=FALSE)
popsubGenes <- subset(popfilterInd,gene %in% popWideGenes)
popsubGenes$obs <- 1:dim(popsubGenes)[1] 
popsubGenes_var<-paste0("5_out/",chr,"/chr",chr,"_",pop,"_PopulationSubGenes.csv")
write.csv(popsubGenes, popsubGenes_var)
nGene <- length(popWideGenes)  # nGene is the number of significant imbalanced genes 

formGenes <- function(genes){
  as.formula(
    paste(
      "cbind(c1,c2)~type2+AC+AG+AT+CG+CT+GT+",
      paste('(as.numeric(gene=="',genes,'")-1|ind)',
            sep="",
            collapse="+"),
      "+(1|obs)"
    ))    
}

############# Population Based Test Add 1!!!!! ################
nGene <- length(popWideGenes) # nGene is the total number of genes that we test, so the test is repeated for each gene

form0 <- as.formula("cbind(c1,c2)~type2+AC+AG+AT+CG+CT+GT+(1|obs)")
system.time(logme0 <-glmer(formula=form0, # logme0 is the null model assuming all genes are balanced 
                           data=popsubGenes,
                           family=binomial(),
                           nAGQ=0L,
                           control=glmerControl(optCtrl=list(maxfun=200000))
))

add1model <- function(i){
  as.formula(paste(
    "cbind(c1,c2)~type2+AC+AG+AT+CG+CT+GT+",
    paste('(as.numeric(gene=="',popWideGenes[i],'")-1|ind)',sep=""),
    "+(1|obs)"
  ))
}
testGeneAdd1 <- function(i){
  logmeAdd <-glmer(formula=add1model(i), # logmeAdd is the model with one gene added
                   data=popsubGenes,
                   family=binomial(),
                   nAGQ=0L,
                   control=glmerControl(optCtrl=list(maxfun=200000))
  )
  #    return(logmeAdd)
  estVarComp <- as.numeric(summary(logmeAdd)$varcor[[2]])
  pval <- anova(logmeAdd,logme0)[2,8] # here the two models are compared (using anova) to see if the add gene has a imbalance that is actually significant
  c(est=estVarComp,pval=pval)
}

add1tests <- parallel::mclapply(1:nGene,testGeneAdd1,mc.cores=20)
add1tests <- rbind.list(add1tests)
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQpopulation_add1testsPvals",ns,"ns",d,"d.pdf"))
qqp(add1tests[,2],cex.lab=1.5,pch=19,col=3,ci=FALSE)
dev.off()
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQpopulation_add1testsPvals",ns,"ns",d,"dZOOM.pdf"))
qqp(add1tests[,2],cex.lab=1.5,pch=19,col=3,ci=FALSE,ylim=c(0,50))
dev.off()
pdf(paste0("6_plot/",chr,"/chr",chr,"_",pop,"_QQpopulation_add1testsPvals",ns,"ns",d,"dZOOM2.pdf"))
qqp(add1tests[,2],cex.lab=1.5,pch=19,col=3,ci=FALSE,ylim=c(0,10))
dev.off()
save(add1tests,file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_add1testsPvals.Rdata"))
add1testsNew <- data.frame(add1tests,genes=popWideGenes)
write.table(add1testsNew[order(add1tests[,2]),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_add1testsPvals",ns,"ns",d,"d.txt"))
write.csv(add1testsNew[order(add1tests[,2]),],file=paste0("5_out/",chr,"/chr",chr,"_",pop,"_add1testsPvals",ns,"ns",d,"d.csv"))

#########################################################################################################
