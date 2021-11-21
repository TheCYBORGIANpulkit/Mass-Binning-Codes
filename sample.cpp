#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_linalg.h>

using namespace std;
/*
int linesinfile(char *infname)
  {
  int nlines;
  char *c,buf[1000];
  printf("Opening %s \n",infname);
  FILE *infile;
  if((infile = fopen(infname,"r")) == NULL){
    printf("File %s does not exist ... aborting \n",infname);
    perror("FILE DOES NOT EXIST ... ABORTING \n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }//if
  nlines = 0;
  while((c = fgets(buf,1000,infile)) != NULL)
    nlines++;
  rewind(infile);
  fclose(infile);
  printf("nlines = % d in %s \n",nlines,infname);
  return(nlines);
  }//linesinfile


int  main(int argc, char** argv){
    //no. of arguments expected 
    if(argc != 2){
        cout << "Expecting 2 arguements.The code itself and the txt file." << endl;
        return 0;
    }
    for(int i = 0;i< argc;i++)cout << "arguement is:" << i << " " << argv[i] << endl;

    //1. Making the file and counting the number of lines in the file. 
    ifstream infile(argv[1]);
    int n = linesinfile(argv[1]);
    double arr[n] = { 0 };
    for (int i = 0; i < n; i++){
        infile >> (arr[i]);   
    }
    cout << arr[0] << endl;
    return 0;
}
*/

     
     
     int
     main (void)
     {
       //double x = 5.0;
       //double y = gsl_sf_bessel_J0 (x);
       //printf ("J0(%g) = %.18e\n", x, y);
       return 0;
     }