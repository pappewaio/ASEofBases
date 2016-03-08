// g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -std=c++11  -lz


// ./VCFmergeGTF3 TEST_REGIONS.gore ../playground/13.GBR.mac2.recode.vcf TEST_INDIVIDUALS.ind TEST_out.line.gz



#include<iostream>
#include<cstdlib>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <math.h>
#include "zlib.h"
//#include <regex>
#define LENS 2000000

int majorminor[4][4] = 
  {
    {0,1,2,3},
    {1,4,5,6},
    {2,5,7,8},
    {3,6,8,9} 
  };

int baseToNum[256];
 char bases[4]={'A','C','G','T'};


double addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}


size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

gzFile openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"%s %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  gzFile fp = gzopen(c,"w");
  delete [] c;
  return fp;
}

gzFile openFile(const char* a,const char* b,int c,int i){
  char ary[100];//dragon
  snprintf(ary,100,"%s%d.chr%d.gz",b,i,c);
  return openFile(a,ary);  
}



FILE *getFILE(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  if(0&&writeFile&&fexists(fname)){//DRAGON
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  fprintf(stderr,"\t-> opening filehandle for: %s\n",fname);
  return fp;
}

double *getData(const char *fname,int& returnSize,int nind){
  FILE *tmp = getFILE(fname,"r");
  size_t  filesize = fsize(fname);
  fprintf(stderr,"Will open list of true var sites: %lu\n",filesize);
  double *monkey = new double[filesize/sizeof(double)];
  int sizeread=fread(monkey,1,filesize,tmp);
  if(sizeread!=filesize/sizeof(double)){
    fprintf(stderr,"wrong file size\n");
  }
  returnSize=filesize/sizeof(double)/nind/10;
  return monkey;
}



int readRegion(const char* regionFilename,int *&start,int *&stop,int *&chr){
  const char *delims = " \t\n";
  FILE *infile;
  char buffer[LENS]; 

  if(NULL==(infile=fopen(regionFilename,"r"))){
    fprintf(stderr,"Error opening file: %s\n",regionFilename);
    exit(0);
  }
  int nRegions=0;
  while(fgets(buffer,LENS,infile)) 
    nRegions++;
  fclose(infile);

  infile=fopen(regionFilename,"r");
  start = new int[nRegions];
  chr = new int[nRegions];
  stop = new int[nRegions];
  //  start[0]=startR;
  //stop[0]=stopR;

  for(int i=0;i<nRegions;i++){
    if(!fgets(buffer,LENS,infile)){
      exit(0);
    }
    strtok(buffer,"r");
    //    fprintf(stderr,"throw away is: %s\n",strtok(buffer,"r"));
    //fprintf(stderr,"next is: %s\n",strtok(NULL,delims));
    //exit(0);
  
    chr[i]=atoi(strtok(NULL,delims));
    start[i]=atoi(strtok(NULL,delims));
    stop[i]=atoi(strtok(NULL,delims));

  }
  fclose(infile);
 
   return nRegions;
}



std::vector<char*> getFilenames(char * name){
  if(!fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t";
  std::vector<char*> ret;
  std::ifstream pFile(name,std::ios::in);

  char buffer[LENS];
  while(!pFile.eof()){
    pFile.getline(buffer,LENS);
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      if(tok[0]!='#')
	ret.push_back(strdup(buffer));
      tok = strtok(NULL,delims);
    }
  }
  return ret;
}

int readSNPinfo(const char* regionFilename,int *&pos,int *&chr,int *&major,int *&minor){
  const char *delims = " \t\n";
  FILE *infile;
  char buffer[LENS]; 
  if(NULL==(infile=fopen(regionFilename,"r"))){
    fprintf(stderr,"Error opening file: %s\n",regionFilename);
    exit(0);
  }
  int nRegions=0;
  while(fgets(buffer,LENS,infile)) 
    nRegions++;

  fclose(infile);
  infile=fopen(regionFilename,"r");

  pos = new int[nRegions];
  chr = new int[nRegions];
  major = new int[nRegions];
  minor = new int[nRegions];

  for(int i=0;i<nRegions;i++){
    if(!fgets(buffer,LENS,infile))
      fprintf(stdout,"no data ! \n");
    strtok(buffer,delims);
    chr[i]=atoi(strtok(NULL,delims));
    pos[i]=atoi(strtok(NULL,delims));
    major[i]=baseToNum[strtok(NULL,delims)[0]];
    minor[i]=baseToNum[strtok(NULL,delims)[0]];


  }
  fclose(infile);
 
   return nRegions;
}



gzFile openGLF(const char *fname){
  //  printf("glf[%d]=initing with name:%s\n",id,fname);
  if(!fexists(fname)) {
    fprintf(stderr,"file %s does not exits\n",fname);
    exit(0);
  }

  gzFile theFile;
  theFile = gzopen(fname,"rb");

  //  double *lk = new double [10*nInd];
  //size_t bytesRead = gzread(gz,lk,sizeof(double)*10);
  return theFile;
}



void printRegion(int nSites,int nInd, gzFile *openGLFfiles,int where, int *pos,int *chr,int *major,int *minor,gzFile outfile,gzFile markerFile){

  
  if(nSites<0){
    double **loglikes = new double*[nInd];
    for(int i=0;i<nInd;i++){
      loglikes[i] = new double [(-nSites)*10];
      //      fprintf(stderr,"reading ind %d\n",i);
      size_t bytesRead = gzread(openGLFfiles[i],loglikes[i],sizeof(double)*(-nSites)*10);
      delete[] loglikes[i];
    }
    delete[] loglikes;
    return ;
  }


  gzprintf(outfile,"marker\tallele1\tallele2");
  for(int i=0;i<nInd;i++)
    gzprintf(outfile,"\tInd%d\tInd%d\tInd%d",i,i,i);
  gzprintf(outfile,"\n");
 

  double **loglikes = new double*[nInd];
  for(int i=0;i<nInd;i++){
    loglikes[i] = new double [nSites*10];
    //fprintf(stderr,"reading ind %d\n",i);
    size_t bytesRead = gzread(openGLFfiles[i],loglikes[i],sizeof(double)*nSites*10);
  }

  for(int j=0;j<nSites;j++){

    gzprintf(outfile,"chr%d_%d\t%c\t%c",chr[j+where],pos[j+where],bases[major[j+where]],bases[minor[j+where]]);

    for(int i=0;i<nInd;i++){
      double norm=addProtect3(loglikes[i][j*10+majorminor[major[j+where]][major[j+where]]],loglikes[i][j*10+majorminor[major[j+where]][minor[j+where]]],loglikes[i][j*10+majorminor[minor[j+where]][minor[j+where]]]);
     
      gzprintf(outfile,"\t%f",exp(loglikes[i][j*10+majorminor[major[j+where]][major[j+where]]]-norm));
      gzprintf(outfile,"\t%f",exp(loglikes[i][j*10+majorminor[major[j+where]][minor[j+where]]]-norm));
      gzprintf(outfile,"\t%f",exp(loglikes[i][j*10+majorminor[minor[j+where]][minor[j+where]]]-norm));
      // fprintf(outfile,"\t%f\t%f\t%d\t%d",loglikes[i][j*10],norm,majorminor[major[j+where]][minor[j+where]],major[j+where]);

    }
  gzprintf(outfile,"\n");
  }

  //print marker file
  for(int j=0;j<nSites;j++){
    gzprintf(markerFile,"chr%d_%d\t%d\t%c\t%c\n",chr[j+where],pos[j+where],pos[j+where],bases[major[j+where]],bases[minor[j+where]]);

  }
  for(int i=0;i<nInd;i++)
    delete[] loglikes[i];
  delete[] loglikes;


}
  
int getNumberOfSNPs(int *reg,int *pos,int *chr,int where,int c){
  int i=0;
  fprintf(stderr,"region start %d\t%d where %d\n",pos[where],chr[where],where);
  while(pos[where-i]<reg[0]|chr[where-i]!=c){
    if((pos[where-i]>reg[1]&&chr[where-i]==c)||chr[where-i]>22||chr[where-i]<1)
      fprintf(stdout,"%d\t%d\t%d\n",i,chr[where-i],c);
    i--;
  }

  if(i<0){
    fprintf(stderr,"getN %d \n",i);
    return i;
  }
  while(pos[i+where]<reg[1]&&chr[i+where]==c){
    // fprintf(stderr,"region start %d\t%d\n",pos[where+i],chr[where+i]);
    i++;
  }
  fprintf(stderr,"getN %d \n",i);
  return i;

}

int main(int argc, char *argv[]){
  
  baseToNum['A'] = 0;
  baseToNum['C'] = 1;
  baseToNum['G'] = 2;
  baseToNum['T'] = 3;

  ///////////////////////////////
  //   get regions 
  ///////////////////////////////////
  int *chrRegion;
  int *start;
  int *stop;
  
  int nRegions=readRegion(argv[1],start,stop,chrRegion);
  
  fprintf(stderr,"Number of regions %d\n",nRegions);
  //  for(int i=0;i<nRegions;i++)
  //    fprintf(stderr,"%d\t%d\t%d\n",chrRegion[i],start[i],stop[i]);
  
  char** geneIDs = new char*[nRegions];    
  char** geneNames = new char*[nRegions];  
  
  FILE *gtffile;
  char gtfbuffer[LENS]; 
  
  if(NULL==(gtffile=fopen(argv[1],"r"))){
    fprintf(stderr,"Error opening file: %s\n",argv[1]);
    exit(0);
  }
  
  for(int i=0;i<nRegions;i++){// READING REGIONS
    
    if(!fgets(gtfbuffer,LENS,gtffile)){
      exit(0);
    }
    strtok(gtfbuffer,"\"");
    geneIDs[i]=strdup(strtok(NULL,"\""));
    //transcript_id: 
    strtok(NULL,"\"");
    strtok(NULL,"\"");
    //gene_type:
    strtok(NULL,"\"");
    strtok(NULL,"\"");
    //gene_status:
    strtok(NULL,"\"");   
    strtok(NULL,"\"");
    // gene_name:
    strtok(NULL,"\"");
    geneNames[i]=strdup(strtok(NULL,"\"")); 
    
  } // END OF READING REGIONS

  //  fprintf(stderr,"geneNames 1 and 2 are %s and %s \n",geneNames[1],geneNames[2]);
  fclose(gtffile);
  

  
  ////////////////////////////////////////
  // open ind file (max 2000)
  char** indNames = new char*[2000];
  
  FILE *infile;
  char bufferNames[LENS]; 
  if(NULL==(infile=fopen(argv[3],"r"))){
    fprintf(stderr,"Error opening file: %s\n",argv[3]);
    exit(0);
  }
  int nKeepInds=0;
  while(fgets(bufferNames,LENS,infile)){// READING QUERY INDS
    indNames[nKeepInds]=strdup(strtok(bufferNames,"\n"));
    
    fprintf(stdout,"%s\n",indNames[nKeepInds]);
    nKeepInds++;
  } //END OF READING QUERY INDS
  fprintf(stdout,"Number of query individuals: %d\n",nKeepInds);
  
  ///////////////////////////
  /// open outfile
  
  const char *outfiles=argv[4];
  gzFile lineFile = gzopen(outfiles,"w");

  //////////////////////////////////////////
  // open vcf files
  /////////////////////////////////
  gzFile openVCFfiles;
  char *filelist = NULL;
  filelist = strdup(argv[2]);
  openVCFfiles=openGLF(filelist);
  //   fprintf(stdout,"%s\n",filelist);
  
  int keepInd[2000];
  int nInd=0;
  const char *delims = " \t\n";
  char buffer[LENS]; 
  buffer[0]='#';
  buffer[1]='#';
  while(buffer[0]=='#'&&buffer[1]=='#'){ // REMOVING COMMENT LINES FROM VCF FILE
    gzgets(openVCFfiles,buffer,LENS);  
    //    fprintf(stdout,"%c\t%c\n",buffer[0],buffer[1]);
    //fprintf(stdout,"%s",buffer);
  }
  
  //READING INDS FROM VCF HEADER
  strtok(buffer,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  strtok(NULL,delims);
  int nMatch=0;
  
  // HERE HERE
  while(nInd<1092){
    char *next_field=strtok(NULL,delims);
    if(next_field == NULL){
      break;
    }
    // stuff below does not work if next_field is NULL
    // which happens at end of line!
    //fprintf(stdout,"current ind %s\t%c\n",next_field,next_field[0]);
    char *indN=strdup(next_field); //here is segfault
    //fprintf(stdout,"current ind %s\t%c\n",indN,indN[0]);
    //fprintf(stdout,"we reached %d\n",nInd);
    keepInd[nInd]=0;
    for(int i=0;i<nKeepInds;i++){
      
      if(strcmp(indN,indNames[i])==0){
	fprintf(stdout,"match %s\t%s\n",indN,indNames[i]);
	keepInd[nInd]=1;
	nMatch++;
      }
      //fprintf(stdout,"we reached %s\t%s\n",indN,indNames[i]);
    }
    nInd++;
    //fprintf(stdout,"we reached %d\n",nInd);
  } // END OG READING INDS FROM VCF HEADER
  fprintf(stdout,"Number of ind in vcf file %d\t number of match %d\n",nInd,nMatch);
  
  
  int posVCF=0;
  int chrVCF=0;
  
  gzprintf(lineFile,"chr\tpos");
  gzprintf(lineFile,"\tref\talt\tgene_id");
  for(int i=0;i<nKeepInds;i++)
    gzprintf(lineFile,"\t%s\t%s\ttype",indNames[i],indNames[i]); // THIS PRINTS ALL QUERY INDS, SHOULD BE CHANGED
  gzprintf(lineFile,"\n");
  

  // LOOPING trough genes (called regions)
  for(int r=0;r<nRegions;r++){
    //  nRegions=-1;
    int reg[2];
    reg[0]=start[r];
    reg[1]=stop[r];
    //    fprintf(stdout,"Region %d - %d \t chr %d\n",reg[0],reg[1],chrRegion[r]);
    
    while(gzgets(openVCFfiles,buffer,LENS)){ // LOOPING TROUGH SNPS
      
      chrVCF=atoi(strtok(buffer,delims));
      posVCF=atoi(strtok(NULL,delims));
      
      if(posVCF>stop[r]) // IS SNP AFTER CURRENT GENE
	break;
      fflush(stdout);
      if(posVCF<start[r]) // IS SNP BEFORE CURRENT GENE
	continue;
      
      // get rid of ID (rs number):
      strtok(NULL,delims); 
      char mmVCF[2];
      // takes REF:
      mmVCF[0]=strtok(NULL,delims)[0];
      // takes ALT:       
      mmVCF[1]=strtok(NULL,delims)[0];
      // get rid of QUAL:
      strtok(NULL,delims);
      // see if FILTER is PASS and throw snp if not:
      char *pass_filter = strtok(NULL,delims);
      assert(pass_filter!=NULL);
      if(strcmp(pass_filter,"PASS")==0){
	//fprintf(stdout,"pass_filter is %s\n",pass_filter);
      }
      else{
	fprintf(stdout,"pass_filter is %s and we skip it (at pos %d) \n",pass_filter,posVCF);
	continue;
      }
      
      char *info = strtok(NULL,delims); // INFO IS PROBABLY NOT USED?
      // first three columns of each snp prints here
      gzprintf(lineFile,"%d\t%d",chrVCF,posVCF);
      // also print ref and alt:
      gzprintf(lineFile,"\t%c\t%c\t%s",mmVCF[0],mmVCF[1],geneNames[r]);


      for(int i=0;i<nInd;i++){
	// first time lose FORMAT, then lose residual data for current individual
	strtok(NULL,delims);
	//fprintf(stdout,"%s\t\n",strtok(NULL,delims));
	if(keepInd[i]==0) // not individual of interest?
	  continue;
	// assumes phased data:
	char *line_er_sej = strtok(NULL,"|");
	//	 fprintf(stderr,"LINE ER SEJ: %s\n",line_er_sej); 
	//	 exit(0);
	int ind1a=atoi(line_er_sej);
	//fflush(stdout);
	int ind1b=atoi(strtok(NULL,":"));
	if((ind1a!=1&&ind1a!=0)||(ind1b!=1&&ind1b!=0)){
	  fprintf(stdout,"error prob %d\t%d \t pos %d ind %d\n",ind1a,ind1b,posVCF,i);
	  fprintf(stdout,"%s\t\n",strtok(NULL,delims));
	  exit(0);
	}
	//fflush(stdout);
	// prints the phased genotypes
	gzprintf(lineFile,"\t%c\t%c",mmVCF[ind1a],mmVCF[ind1b]);
	//fprintf(stdout,"\t%d %c %c",i,mmVCF[ind1a],mmVCF[ind1b]);
	//prints reference-phase info
	gzprintf(lineFile,"\t%d",ind1a+2*ind1b);
      } // END OF LOOP THROUGH INDIVIDUALS (i)
      //fprintf(stdout,"\n");
      gzprintf(lineFile,"\n");
    } // END OF LOOP THROUGH SNPS (while gzgets)
  } // END OF LOOP OVER REGIONS (r)

  gzclose(lineFile);
  
  //now cleanup allocated memory stuff
  for(int i=0;i<nRegions;i++){
    free(geneNames[i]);
    free(geneIDs[i]);
  }
  delete [] geneNames;
  delete [] geneIDs;
  return 0;
}
