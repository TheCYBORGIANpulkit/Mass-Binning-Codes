// MassBinning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;
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

int main(int argc, char** argv)
{   
    //no. of arguments expected 
    if(argc != 2){
        cout << "Expecting 2 arguements.The code itself and the txt file." << endl;
        return 0;
    }
    for(int i = 0;i< argc;i++)cout << "arguement is:" << i << " " << argv[i] << endl;

    //1. Making the file and counting the number of lines in the file. 
    //1.Function to insert the data in an array 
    ifstream file1(argv[1]);
    const int n = linesinfile(argv[1]);
    double arr[n] = { 0 };
    for (int i = 0; i < n; i++){
        file1 >> (arr[i]);   
    }

    //const int n = 7857;function to read the number of lines will be written here
    //2. Find maximum and minimum of data 
    
    double max = arr[0];
    double min = arr[0];
    for (int i = 0; i < n; i++) {
        //if(arr[i] >= 6))
        if (arr[i] >= max) max = arr[i];
        else if (arr[i] <= min) min = arr[i];
        
    } 
    cout << min << endl;
    cout << max << endl;
    
    // double max = 11.00;
    // double min = 6.00;
    //3. Make Bins
    const int nbins = 100; // it'll be later taken as an input
    ofstream file2("realizations_100.csv"); //output file
    double bins[nbins] = { 0 };
    double halfw = (max - min) /(2*nbins);
    cout << "w = " << 2*halfw << endl;
    for (int i = 0; i < nbins; i++) bins[i] = min + (2 * i + 1) * halfw;
  
   
    //4. Make counter and count
    
    int total = 0;
    double count[nbins] = { 0 };
    double normalcount[nbins];
    
    for (int i = 0; i < nbins ; i++) {
        float l = bins[i] - halfw;
        float h = bins[i] + halfw;
        
        for (int j = 0; j < n; j++) {
            if ((l < arr[j]) && (arr[j] < h)) count[i]++;
            else if(arr[j] == h)count[i]++;
            //if(i == 3 )cout << arr[j] << endl;
        }
        
        total = total + count[i];
        cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << count[i] << " total " << total <<  endl; 
        //cout << "i: " << i  << " bins[i]: " <<  bins[i] << " l: " << l << " h: " << h << endl;
        
        double Vol = 0.638 * (pow(211.7039,3)/3);
        count[i] = count[i] / (n * 2*halfw * Vol); //normalization
        file2 << bins[i] << ","<< count[i] << endl;
        cout << "i: " << i << " count:  " << count[i] << endl;
    }
    
    //cout << bins[23] << endl
    /*;
    //Method 2 for Binning
    int k = 0;
    for(int j = 0;j<n;j++){
        //if(j > 90 && j <=  95)
        //cout << bins[23] << endl;
    
       double a = (arr[j] - bins[0])/(2*halfw);
       k = (int)(a + 0.5001);
       count[k]++;
       //if(j == 95)
          //cout<<arr[j]<<" "<<a<<" "<<i<<" "<<bins[i]<<" "<<count[i]<<endl;

    }
    //cout << bins[23] << endl;
    //for(int i = 0;i<nbins;i++) cout << i << ": " << (2*i+1)*halfw <<"; "<< bins[i] << endl;
    
    for(int i = 0;i<nbins;i++){
        total = total + count[i];
        cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << count[i] << " total " << total <<  endl;
        double Vol = 0.638 * (pow(211.7039,3)/3);
        normalcount[i] = count[i] / (n * 2*halfw * Vol); //normalization
        file2 << bins[i] << ","<< count[i] << endl;
    }
    */
    file1.close(); //closing the file
    file2.close();
    
    return 0;
    
    
}


