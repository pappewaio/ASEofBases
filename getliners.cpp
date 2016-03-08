//version 3
#include <iostream>
#include <fstream> //used for ifstream
#include <cstdio> //used for printf
#include <cstring> //used for strtok strcmp
#include <map>
#include <utility> //used for make_pair
#include <cstdlib> //used atoi

#define LENS 1000000 //buffer length


//comparison operator used when comparing char*
struct cmp_char {
  bool operator()(const char *first,const char* second) const {
    int tmp = std::strcmp(first, second);
    return tmp<0;
  }
};





typedef std::map <char*,int,cmp_char> aMap;


void printMap(const aMap &m){
  for(aMap::const_iterator it = m.begin(); it != m.end(); ++it){
    char *l = it->first;
    int val = it->second;
    fprintf(stderr,"%s\t%d\n",l,val) ;
  }
}





aMap build_map(const char *filename,const char *delims=" \t"){
  std::ifstream pFile(filename,std::ios::in);
  char buffer[LENS];

  
  aMap ret;
  while(pFile.getline(buffer,LENS)){
    char *key = strtok(buffer,delims);
    while(key!=NULL){
      ret.insert(std::make_pair(strdup(key),1));
      key =strtok(NULL,delims);
    }
  }
  
  return ret;
  
}



int main(int argc, char *argv[]){
  if(argc==1){
    printf("USE -c column -d delimter -k keysfile -f infile\n");
    return 0;
  }







  int extractCol=2;

  int argPos = 1;
  const char *delims = "\t";
  const char* indexfile = "keys";
  const char* datafile = "test.txt";
  while(argPos <argc){
    if(strcmp(argv[argPos],"-c")==0)
      extractCol  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-d")==0)
      delims = strdup(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-k")==0)
      indexfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-f")==0)
      datafile = argv[argPos+1];
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      printf("USE -c column -d delimter -k keysfile -f infile\n");
      return 0;
    }
    argPos+=2;
  }



  aMap asso = build_map(indexfile,delims);
  //  printMap(asso);


  std::ifstream pfile(datafile);
  char buffer[LENS];
  while(pfile.getline(buffer,LENS)){
    char *original = strdup(buffer);
    char *cmp = strtok(buffer,delims);
    for(int counter=1;counter<extractCol;counter++)
      cmp = strtok(NULL,delims);
    //cmp now contains the correct column
    aMap::iterator it = asso.find(cmp);
    
    if(it!=asso.end()){
      printf("%s\n",original);
      it->second =0;
    }
    //    else{ fprintf(stderr,"KEY:%s doesn't exist\n",cmp);        }
  }
  //  printMap(asso);
  return 0;
}
