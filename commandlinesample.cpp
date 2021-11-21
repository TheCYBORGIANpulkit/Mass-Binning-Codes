#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;
int main(int argc, char** argv){
    int linesinfile(char* infname){
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
        while ((c = fgets(buf, 1000, infile)) != NULL) nlines++;
        rewind(infile);
        fclose(infile);
        printf("nlines=%d in %s \n", nlines, infname);
        return(nlines);
    }//linesinfile

    int n, nr;
    /////////////////////////////////////////////////////////////
    // here you specify how many arguments the code should expect
    //
    if (argc != 3){
        cout << "\nCommand line args\n";
        cout << "./a.out   [catalog.txt]  [#realizations]\n";
        cout << "./a.out  all_1p5sigma.d        300  \n";
        return(-1);
    }

    n = linesinfile(argv[1]);
    cout << "number of lines in file: " << n << endl;
    return 0;
}

