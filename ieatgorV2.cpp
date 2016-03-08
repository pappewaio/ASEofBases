////////////////////////////////////////////////////////////////////////////
//    g++ -O3 -o ieatgor ieatgorV2.cpp -lz
//    ./getTarget targetRegion.txt temp.mafs
//   prog/getTarget target.temp ~/temp.gorout  | head -n1
//  Anders Albrechtsen
// input files needs to be sorted
// works for all files with chr and positions as the two first tab seperated columns
//////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <string>
#include "zlib.h"
#define LENS 100000 //this is the max number of bytes perline, should make bigger
int main(int argc, char **argv){
  char *saveptrMaf, *saveptrTar;
  if(argc==1){
    fprintf(stderr,"\nRUN:\tieatgor <targetFile> <file> OPTIONS\n\n");
    fprintf(stderr,"target file format:\t chrName:start-end\n\n");
    fprintf(stderr,"OPTIONS:\n");
    fprintf(stderr,"\t-offset INT\t offset the target positions by INT\n");
    exit(0);
  }
  const char* targetFileName=argv[1];
  const char* mafFileName=argv[2];
  const char *delimsTarget = "\t-:\n";
  const char *delimsMafs = "\t\n";
  FILE *targetFile = NULL;
  gzFile mafFile = NULL;
  int offset=0;
  if(NULL==(targetFile=fopen(targetFileName,"r"))){
    fprintf(stderr,"Error opening file: %s\n",targetFileName);
    exit(0);
  }
   if(NULL==(mafFile=gzopen(mafFileName,"r"))){
    fprintf(stderr,"Error opening file: %s\n",mafFileName);
    exit(0);
   }


   for(int i=3;i<argc;i++){
     if(strcmp(argv[i],"-offset")==0){
       i++;
       offset=atoi(argv[i]);
     }
     else {
       printf("\tUnknown arguments: %s\n",argv[i]);
      printf("USE -offset AFTER the target and file\n");
      return 0;
    }

   }
   char bufTarget[LENS];
   char bufMaf[LENS];

   if(fgets(bufTarget,LENS,targetFile)==0)
     fprintf(stderr,"no data in target:\n");

   char *tok = strtok_r(bufTarget,delimsTarget,&saveptrTar);//strtok_r(bufTarget,delimsTarget,&saveptrTar);

  if(tok==NULL){
    fprintf(stderr,"NULL pointer:\n");
    return 0;
  }

  char* chrTar=strdup(tok); //not really need strdub in this case
  int start=atoi(strtok_r(NULL,delimsTarget,&saveptrTar));
  int stop=atoi(strtok_r(NULL,delimsTarget,&saveptrTar));
  char *f=gzgets(mafFile,bufMaf,LENS);
  char* chrMaf=NULL;
  while((f=gzgets(mafFile,bufMaf,LENS))) {
    free(chrMaf);
    chrMaf=strdup(strtok_r(bufMaf,delimsMafs,&saveptrMaf));
    int comp=strcmp(chrMaf,chrTar);
    if(comp<0)
      continue;
    char *tok2 = strtok_r(NULL,delimsMafs,&saveptrMaf);
    int pos=atoi(tok2);

    if(comp==0&&pos<start)
      continue;
    while(comp>0||(comp==0&&pos>stop)){
      if(fgets(bufTarget,LENS,targetFile)==0){
	fprintf(stderr,"no more target file \n");
	goto gotoEndOfLoop;

      }
      free(chrTar);
      chrTar=strdup(strtok_r(bufTarget,delimsTarget,&saveptrTar)); 
      start=atoi(strtok_r(NULL,delimsTarget,&saveptrTar))+offset;
      stop=atoi(strtok_r(NULL,delimsTarget,&saveptrTar))+offset;
      comp=strcmp(chrMaf,chrTar);
      //  fprintf(stdout,"all %s %d %d\n",chrTar,start,stop);
    }
    if(comp<0||pos<start)
      continue;
    else{
      // fprintf(stdout,"realTar %s %d %d\n",chrTar,start,stop);
      fprintf(stdout,"%s\t%d\t%s",chrMaf,pos,bufMaf+strlen(chrMaf)+strlen(tok2)+2);
    }
  

  }
 gotoEndOfLoop:;
  fclose(targetFile);
  gzclose(mafFile);
  free(chrTar);
  free(chrMaf);
 
  return 0;
}


    //strcmp(chr,chr2)==0
    //chr<chr2 => negativ
