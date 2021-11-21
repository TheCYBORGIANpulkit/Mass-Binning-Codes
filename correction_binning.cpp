#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>
#include </home/pullu1729/Desktop/cosmo/codes/MassBinning/distances.cpp>
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

//power Law function
double F1(double LW){
    double S = 0;
    if(LW <= 2.5 )S = 0.5*LW - 1.207;
    else S = LW - 2.457;
    return pow(10,S);
}

//Luminosity distance and Flux function
double F2(double S, double M){
    double f = 2.356*(pow(10,5));
    double m = pow(10, M);
    double Dsq = (m/(S*f));
    double D = pow(Dsq,0.5);
    return D;
}
int main(int argc, char** argv){   
    //no. of arguments expected 
    if(argc != 5){
        cout << "Expecting 3 arguements.The code itself and the txt file." << endl;
        return 0;
    }
    for(int i = 0;i< argc;i++)cout << "arguement is:" << i << " " << argv[i] << endl;

    //1. Making the file and counting the number of lines in the file. 
    //1.Function to insert the data in an array 
    ifstream file1(argv[1]);
    ifstream file2(argv[2]);
    const int n = linesinfile(argv[1]);
    double arrM[n] = { 0 };
    double arrW[n] = { 0 };
    for (int i = 0; i < n; i++){
        file1 >> arrM[i]; 
        file2 >> arrW[i];  
    }

    //const int n = 7857;function to read the number of lines will be written here
    //2. Find maximum and minimum of data 
    /*
    double max = arr[0]
    double min = arr[0];
    for (int i = 0; i < n; i++) {
        //if(arr[i] >= 6))
        if (arr[i] >= max) max = arr[i];
        else if (arr[i] <= min) min = arr[i];
        
    } 
    cout << min << endl;
    */
    double max = 11.00;
    double min = 6.00;
    //3. Make Bins
    const int nbins = 25; // it'll be later taken as an input
   
    double bins[nbins] = { 0 };
    double halfw = (max - min) /(2*nbins);
    cout << "w = " << 2*halfw << endl;
    for (int i = 0; i < nbins; i++) bins[i] = min + (2 * i + 1) * halfw;

    //F1 calculating S-int from log W_50
    double arrS[n] = {0};
    for(int j = 0;j< n;j++) arrS[j] = F1(arrW[j]);
    //if(j == 7058) cout << arrS[j] << endl;}
    cout << "*********************" << endl;

    //F2 calculating Dl from S
    double arrDlmax[n] = { 0 };
    for(int j = 0;j<n;j++) {
        arrDlmax[j] = F2(arrS[j], arrM[j]);
        //if(j == 7058) cout << arrDlmax[j] << endl;
    }
    cout << "******************************" << endl;

    //making a 2-D array for matching the luminosity distance to z-value
    const int nz = 1000;
    double arrNZ[nz][2];
    ifstream file3(argv[3]);
    for(int i = 0;i<nz;i++){
        for(int j = 0;j<2;j++) file3 >> arrNZ[i][j];
    }
    cout << "The last element of luminosity array is: " << arrNZ[nz - 1][1] << endl;
    //matching the arrays (z corresponding to dl_max) arrz contains z, cd_max, comov_vol and the corresponding weight.
    double Vol = 0.638 * (pow(211.7039,3)/3); //total survey volume
    double arrz[n][4] = {0};
    for(int i = 0;i < n;i++){
        //cout << i <<  "M:  " << arrM[i] << " W: " <<  arrW[i] << endl;
        double closest_dl = 0;//Checking for predecessor and successor
        double pre = 0;
        int j = 0;
        while(closest_dl < arrDlmax[i] && j < nz){
            closest_dl = arrNZ[j][1];
            pre = arrNZ[j - 1][1];
            arrz[i][0] = arrNZ[j][0];
            //cout << "closest_dl: " << closest_dl << endl;
            j++;
        }
        int a = j - 1;
        if(abs(pre - arrDlmax[i]) < abs(closest_dl - arrDlmax[i]))a = j - 2; //Checking for predecessor and successor 
        arrz[i][0] = arrNZ[a][0];
        if(i == 7058)cout << "closest_dl" << arrNZ[a][1] <<  " z " << arrz[i][0] << endl;

        //cout << "*************************" << endl;

        //calculating the comoving distance and then volume
        double dc_max = Trapezoidal(com_dis,0.3,0,0.7,0,arrz[i][0],1000);
        arrz[i][1] = dc_max;
        double comov_vol =  0.638*pow(dc_max,3)/3;
        arrz[i][2] = comov_vol;
        
        //asigning weight and counting
        if(arrz[i][2] < Vol) arrz[i][3] = Vol/arrz[i][2];
        else arrz[i][3] = 1;
        //if(i == 7058)cout  << arrNZ[a][1] << ",,," <<   arrz[i][0] <<  ",,," << arrz[i][1] << ",,,,," << arrz[i][2]<< ",,," << arrz[i][3] << endl;
    }   
    cout << "*************************" << endl;

    /*calculating the comoving distance and then volume
    for(int i =0;i < n;i++){
        double dc_max = Trapezoidal(com_dis,0.3,0,0.7,0,arrz[i][0],1000);
        arrz[i][1] = dc_max;
        double comov_vol =  0.638*pow(dc_max,3)/3;
        arrz[i][2] = comov_vol;
         if(i < 10)cout  << arrz[i][0] <<  ",,," << arrz[i][1] << ",,,,," << arrz[i][2] << endl;
    }
    */
    
    cout << "*************************" << endl;
    //4. Make counter and count
    int total = 0;
    double count[nbins][3] = { 0 };
    
    ofstream file4(argv[4]); //output file
    
    for (int i = 0; i < nbins; i++) {
        float l = bins[i] - halfw;
        float h = bins[i] + halfw;
        
        for (int j = 0; j < n; j++) {
            if ((l < arrM[j]) && (arrM[j] <= h)){
                count[i][0] = count[i][0] + arrz[j][3];
                count[i][1]++; 
                //cout << j << " weight: "  << arrz[j][3] << "  rel_count:  " << count[i][0] << "  count:  " << count[i][1] << endl;
            }
            /*
            else if(arrM[j] == h){
                count[i][0]= count[i][0] = count[i][0] + arrz[j][3];
                count[i][1]++;
                cout <<  arrz[j][3] << " ,,, " << count[i][0] << " ,,,, " << count[i][1] << endl;
            }*/
            //if(i == 17 )cout << arrz[j][3] << endl;
        }
        
        total = total + count[i][1];
        count[i][2] = 1/pow(count[i][1], 0.5);
        //cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << count[i][0] << ",,," << count[i][1] << " total " << total <<  endl; 
        //cout << "i: " << i  << " bins[i]: " <<  bins[i] << " l: " << l << " h: " << h << endl;
        
        //double Vol = 0.638 * (pow(211.7039,3)/3);
        count[i][0] = count[i][0] / (2*halfw * Vol); //normalization
       
        count[i][2] = count[i][2] * abs(count[i][0]); // poisson error in liear scale
        count[i][2] = log(count[i][0] + count[i][2]) - log(count[i][0]);
        cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << log10(count[i][0]) << " ,,, " 
        << count[i][1] << " ,,, " << count[i][2] << " total " << total <<  endl;
        
        file4 << bins[i]  << " " << log10(count[i][0]) << " " << count[i][2] << endl;
        //cout << "i: " << i << " count:  " << count[i] << endl;
    }
    
    /*
    
    //Method 2 for Binning
    int k = 0;
    for(int j = 0;j<n;j++){
        //if(j > 90 && j <=  95)
        //cout << bins[23] << endl;
    
       double a = (arrM[j] - bins[0])/(2*halfw);
       k = (int)(a + 0.5001);
       //count[k]++;
       count[k][0] = count[k][0] + arrz[j][3];
        count[k][1]++; 
        if(j == 95)cout<<arrM[j] << endl;
    }

    //cout << bins[23] << endl;
    //for(int i = 0;i<nbins;i++) cout << i << ": " << (2*i+1)*halfw <<"; "<< bins[i] << endl;
    for(int i = 0;i< nbins;i++){
        total = total + count[i][1];
        count[i][2] = 1/pow(count[i][1], 0.5);
        //cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << count[i][0] << ",,," << count[i][1] << " total " << total <<  endl; 
        //cout << "i: " << i  << " bins[i]: " <<  bins[i] << " l: " << l << " h: " << h << endl;
        
        //double Vol = 0.638 * (pow(211.7039,3)/3);
        count[i][0] = count[i][0] / (2*halfw * Vol); //normalization
       
        count[i][2] = count[i][2] * abs(count[i][0]); //final poisson error
        count[i][2] = log(count[i][0] + count[i][2]) - log(count[i][0]);
        cout << "i: " << i  << " bins[i]: " <<  bins[i] << " count: " << log10(count[i][0]) << " ,,, " 
        << count[i][1] << " ,,, " << count[i][2] << " total " << total <<  endl;
        
        file4 << bins[i] << " "<< log10(count[i][0]) << " " << count[i][2] << endl;
        //cout << "i: " << i << " count:  " << count[i] << endl;
    }
    */
    file1.close(); //closing the file
    file2.close();
    file3.close();
    file4.close();
    cout << "done" << endl;
    
    return 0;
}