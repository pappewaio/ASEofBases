# ASEofBases: Allele Specific Expression method and pipeline for analysis

## What is ASEofBases?
Biased allele expression, refers to the imbalanced expression of the two alleles in a diploid genome. Unequal transcription of alleles may occur due to cis-regulatory element variation or allele-specific epigenetic modifications.  Allelic imbalance may be incorrectly inferred due to technical variation inherent in RNA-Seq data, including read depth, reference mapping bias, and the overdispersion of reads. To correct for technical variation we develop a logistic regression model with a mixed effects approach to combine information regarding biased allele expression from many individuals in a population, and across multiple genes. Here we describe a new method for inferring allelic imbalance that combines information from multiple SNP sites within a transcribed unit, using a logistic regression model that explicitly models the effects of reference bias and SNP type biases. Additionally, by using a mixed effects approach the method also makes it possible to combine information from many individuals in a population for each gene and to test hypothesis regarding allelic imbalance in specific genes within and across populations. 

## Citing ASEofBases
If you use ASEofBases package, please cite:

- Inference on imbalanced allelic expression from RNA sequencing data in individuals and across groups using logistic regression models. Skotte L, Olney K, Nielsen R, and Wilson Sayres M. (https://dx.doi.org/10.6084/m9.figshare.1564499.v1) Manuscript in preparation

If you use the `1000Genomes` or `Geuvidas` data sets, please also cite:

-  A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi: [10.1038/nature15393](http://www.1000genomes.org/home).
- Lappalainen et al. Nature 2013 : Transcriptome and genome sequencing uncovers functional variation in humans (http://dx.doi.org/10.1038/nature12531).

## Key features
- simulations for modeling statistcal power for various read depth and various number of SNPs per gene 
- simulations for modeling SNP-wise binomial, binomial-based logistic regression, and binomial-based logistic regression correcting for the over dispersion of read counts. 
- model of allele-specific expression in single individuals  
- model of allele-specific expression in a population 

## Installation
Download [the latest release of ASEofBases](https://github.com/WilsonSayresLab/ASEofBases/archive/master.zip), unzip and add the `bin` directory to your `PATH`. E.g.:

		wget -O ASEofBases-master.zip https://github.com/WilsonSayresLab/ASEofBases/archive/master.zip
		unzip ASEofBases-master.zip && cd ASEofBases-master
		export PATH=`pwd`/bin:$PATH

Directory overview of ASEofBases:

/ASEofBases/1_code/	# bash scripts
- 1_get.sh
- 2_run.sh
- 3_makeData.sh
- ASEofBasesAnalysis.R
- simAoB.R

/ASEofBases/2_prog/	# locally written C++ programs
- getliners
- ieatgor
- VCFmergeGTF3

/ASEofBases/3_raw/	# raw and parsed data

/ASEofBases/4_data/	# parsed RNAseq data for filtering & analysis in R

/ASEofBases/5_out/	# analysis and output from simulations made in R

/ASEofBases/6_plot/	# plots made from the output data made in R 

You will also need the following perviosuly published programs installed:

- angsd0.563				https://github.com/ANGSD/angsd
- vcftools_0.1.12b			http://sourceforge.net/projects/vcftools/files/
- zlib-1.2.8				http://www.zlib.net	
- samtools/htslib			https://github.com/samtools/htslib
- bigwig				https://genome.ucsc.edu/goldenPath/help/bigWig.html

	
# ASEofBases README
This set of bash commands downloads initial data and does some processing 
- Convert Mapability file to bed format and create file to filter on mapability
- Get protein coding gene annotations
- Download & parse variation and genotype information from 1000Genomes individuals
- Download & parse individual information from Geuvadis RNAseq data

## 1. Compile programs/code to be executable
Getliners - Merge the tmp.keys (positions) with tmp.het (heterozygous protein coding SNPs for individual) and greps the set of keys in tmp.keys in column 2 of the file tmp.het
ieatgor - Filter for alignability output chr$chr.ind$ind.data, greps any entry in tmp.data that has chromosome and position within the regions and specified in the "targetfile" (regions of the genome that are "callable")
VCFmergeGTF - This code is merging the genotype calls (from vcf) with position of protein coding genes (from gencode)

		cd /ASEofBases/2_prog/
		g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -lz
		g++ -O3 -o ieatgor ieatgorV2.cpp -lZ
		g++ -O3 -o getliners getliners.cpp -lz

## 2. Program check
Make sure each of these programs are installed

	angsd0.563			          https://github.com/ANGSD/angsd
	vcftools_0.1.12a	                  http://sourceforge.net/projects/vcftools/files/
	zlib-1.2.8			          http://www.zlib.net	
	samtools/htslib		                  https://github.com/samtools/htslib
	bigwig				          https://genome.ucsc.edu/goldenPath/help/bigWig.html

## 3. Download Mapability file, 
**** includes some point and clicking *****
The mapability file is used for ...........

Download the Mapability file:
	http://moma.ki.au.dk/genome-mirror/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability
  	wgEncodeEH000320 for wgEncodeCrgMapabilityAlign50mer.bigWig

    	mkdir ASEofBases/2_prog/BigWig/
    	mv wgEncodeCrgMapabilityAlign50mer.bigWig /ASEofBases/2_prog/BigWig/

## 4. Get raw data and conduct initial parsing
1_getRaw.sh is a bash script for ..... 

    sh /ASEofBases/bash_scripts/1_getRaw.sh


## 5. Parse and merge genotype and transcriptome data
This bash script will call another script and together will conduct filtering and parsing of the genotype and RNAseq data.
		sh /ASEofBases/1_code/2_run.sh # this program runs 3_makeData.sh

This set of bash commands downloads initial data and does some processing 
For all Individuals, chromosomes, and population:


 	ANGSD 
		Use ANGSD to get counts, mapq filter 40, remove duplicates from BAM files 


	Output of data parsing and filtering
	The output of this step will result in a file for each individual with these columns:
		chr: chromosome 
		pos: position
		ref: reference allele
		alt: alternative allele
		g: gene
		a1: observed phased allele 1
		a2: observed phased allele 2
		type: snp type (1 if allele 1 is the alternative allele and 2 if allele 2 is the alternative allele)
		n: total counts of reads overlapping
		nA: counts of As
		nC: counts of Cs
		nG: counts of Gs
		nT: counts of Ts

## 6. Run regression analysis and make plots in R
For regression analysis using data generated in pervious steps

		ASEofBases.R

For simulations

		simAoB.R

![simulations collage](https://raw.github.com/WilsonSayresLab/ASEofBases/master/doc/simulations.png)	


