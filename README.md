# ASEofBases
Allele Specific Expression method and pipeline for analysis


# OverView
---Publicly available programs-------------------------

		angsd0.563			https://github.com/ANGSD/angsd
		vcftools_0.1.12a		http://sourceforge.net/projects/vcftools/files/
		zlib-1.2.8			http://www.zlib.net	
		samtools/htslib			https://github.com/samtools/htslib
		bigwig				https://genome.ucsc.edu/goldenPath/help/bigWig.html

---In house programs-------------------------

	/ASEofBases/bash_scripts/
		1_get.sh
		2_run.sh
		3_makeData.sh

	/ASEofBases/programs/
		getliners
		ieatgor
		VCFmergeGTF3

	/ASEofBases/R_scripts/
		1_analysis.R
		2_individual.R
		3_population.R
		simAoB.R
		simpop.R
	
# ASEofBases README
# 1. Prepare directories and files 
make directories

		cd 				# move to home directory or location you would like to run the analysis
		mkdir /ASEofBases
		mkdir /ASEofBases/1_code			        # bash scripts
		mkdir /ASEofBases/2_prog 			        # locally written C++ programs
		mkdir /ASEofBases/2_prog/R_scripts 	  		# locally written R programs
		mkdir /ASEofBases/3_raw  		         	# raw and parsed data
		mkdir /ASEofBases/4_data 			        # parsed RNAseq data for filtering & analysis in R
		mkdir /ASEofBases/5_out  			        # analysis and output from simulations made in R
		mkdir /ASEofBases/6_plot 			        # plots made from the output data made in R 

move scripts and programs to appropriate locations

		mv ASEofBases_codes/bash_scripts/* ASEofBases/1_code/
		mv ASEofBases_codes/programs/* ASEofBases/2_prog/
		mv ASEofBases_codes/R_scripts/* ASEofBases/2_prog/R_scripts/

compile programs/code to be executable

		cd /ASEofBases/2_prog/
		g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -lz
		g++ -O3 -o ieatgor ieatgorV2.cpp -lZ
		g++ -O3 -o getliners getliners.cpp -lz

# 2. program check
Make sure each of these programs is loaded on the server or desktop

		angsd0.563			          https://github.com/ANGSD/angsd
		vcftools_0.1.12a	                  http://sourceforge.net/projects/vcftools/files/
		zlib-1.2.8			          http://www.zlib.net	
		samtools/htslib		                  https://github.com/samtools/htslib
		bigwig				          https://genome.ucsc.edu/goldenPath/help/bigWig.html

# 3. Download Mapability file, **** includes some point and clicking *****
Download the Mapability file:
  http://moma.ki.au.dk/genome-mirror/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability
  wgEncodeEH000320 for wgEncodeCrgMapabilityAlign50mer.bigWig

    mkdir ASEofBases/2_prog/BigWig/
    mv wgEncodeCrgMapabilityAlign50mer.bigWig /ASEofBases/2_prog/BigWig/

# 4. Get raw data and conduct initial parsing

    sh /ASEofBases/bash_scripts/1_getRaw.sh

This set of bash commands downloads initial data and does some processing 
		1. Convert Mapability file to bed format and create file to filter on mapability
		2. Get protein coding gene annotations
		3. Download & parse sample information from 1000 genomes individuals
		4. Download & parse sample information from Geuvadis RNAseq data

# 5. Parse and merge genotype and transcriptome data
This bash script will call another script and together will conduct filtering and parsing of the genotype and RNAseq data.
		sh /ASEofBases/1_code/2_run.sh # this program runs 3_makeData.sh

This set of bash commands downloads initial data and does some processing 
For all Individuals, chromosomes, and population:

	VCFmergeGTF
		This code is merging the genotype calls (from vcf) with position of 
		protein coding genes (from gencode)

 	ANGSD 
		Use ANGSD to get counts, mapq filter 40, remove duplicates from BAM files 

 	Getliners
		Merge the tmp.keys (positions) with tmp.het (heterozygous protein coding SNPs for individual)
		greps the set of keys in tmp.keys in column 2 of the file tmp.het

 	ieatgor
 		Filter for alignability output chr$chr.ind$ind.data
		greps any entry in tmp.data that has chromosome and position within the regions 
		specified in the "targetfile" (regions of the genome that are "callable")

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

# 6. Run regression analysis and make plots in R



