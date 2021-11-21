// MassBinning.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std;
/*
int main(int argc, char** argv)
{
    int linesinfile(char* infname)
    {
        int nlines;
        char* c, buf[1000];
        printf("Opening %s \n", infname);
        FILE* infile;
        if ((infile = fopen(infname, "r")) == NULL) {
            printf("File %s does not exist ... aborting \n", infname);
            perror("FILE DOES NOT EXIST ... ABORTING \n");
            fflush(stderr);
            exit(EXIT_FAILURE);
        }//if
        nlines = 0;
        while ((c = fgets(buf, 1000, infile)) != NULL)
            nlines++;
        rewind(infile);
        fclose(infile);
        printf("nlines=%d in %s \n", nlines, infname);
        return(nlines);
    }//linesinfile

    int n, nr;
    /////////////////////////////////////////////////////////////
    // here you specify how many arguments the code should expect
    //
    if (argc != 3)
    {
        cout << "\nCommand line args\n";
        cout << "./a.out   [catalog.txt]  [#realizations]\n";
        cout << "./a.out  all_1p5sigma.d        300  \n";
        return(-1);
    }

    n = linesinfile(argv[1]);
}*/
int main()
{
    //ofstream file2("histogram.csv");

    const int n = 11942;// function to read the number of lines will be written here
    //1.Function to insert the data in an array
    double arr[n] = { 0 };
    ifstream file1;
    file1.open("logmass_1.txt");
    for (int i = 0; i < n; i++){
        file1 >> (arr[i]);
        //cout <<  arr[i] << endl;
    }
    cout << "1st element " << arr[0] << endl;
    cout << "Last element " << arr[n-1] << endl;

    //ofstream file2("histogram.csv");

    //2. Find maximum and minimum of data
    double max = arr[0];
    double min = arr[0];
    for (int i = 0; i < n; i++) {
        if (arr[i] >= max) max = arr[i];
        else if (arr[i] <= min) min = arr[i];
        }
    }
    //cout << max << endl;
    //cout << min << endl;

    //4. Make Bins
    file1.close(); //closing the file
    return 0;
}





