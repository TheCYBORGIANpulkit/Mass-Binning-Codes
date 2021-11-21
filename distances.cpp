#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>


using namespace std;
//co-moving distance expression for m = Omega_matter, r = Omega_radiation and L = Omega_Lambda

double com_dis(double x, double m, double r, double L){
    double s = m*(pow(1+x,3)) + r*(pow(1+x,4)) + L;
    //second term is 0 because \Omega_r = 0;
    double com_distance = 3000/(0.7 *pow(s,0.5));
    // c/H_0 = 3000 MPc
    return com_distance;
}

double Simpsons(double (*fn)(double, double, double, double), double m, double r, double L, double a, double b, int N){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i <= N;i++){
        double  x = (h/3)*((*fn)(a + i*h, m, r, L));
        if(i == 0 || i == N)s = s + x;
        else if(i%2 == 0)s = s + 2*x;
        else s = s + 4*x;
    }
    return s;
}

double Trapezoidal(double (*fn)(double, double, double, double), double m, double r, double L, double a, double b, int N){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<=N;i++){
        double x = a + i*h;
        double y = (h/2)*((*fn)(x, m, r, L));
        if(i == 0 || i == N)s = s + y;
        else s = s+2*y;
    }
    return s;
}
/*
int main(){
    ofstream outfile1;
    ofstream outfile2;
    ofstream outfile3;
    outfile1.open("co_d.csv");
    outfile2.open("ang_dia_d.csv");
    outfile3.open("lum_d.csv");
    for(int i = 0;i<200;i++){
        double s = i/10.00; //s is values for "z" fr scatter plot.
        double a = Trapezoidal(com_dis,0.3,0,0.7,0,s,1000);
        outfile1 << s << "," << a << endl; //comoving distance
        outfile2 << s << "," << a/(s + 1) << endl; //angular diameter
        outfile3 << s << "," << a*(s + 1) << endl; //luminosity distance
        //cout << i/5 <<"," << Trapezoidal(f,0,s,1000) << endl;
    }
    //cout << Simpsons(f,0,1,10000) << endl;
    outfile1.close();
    outfile2.close();
    outfile3.close();
    return 0;
}*/
