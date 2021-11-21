#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
using namespace std;

double fx(double y,double t, double omega_m, double omega_r, double omega_v, double omega_T){
    double f1 = omega_m*(pow(y,-3)) + omega_r*(pow(y,-4)) + omega_v*pow(y,2) + (1 - omega_T);
    //second term is 0 because \Omega_r = 0;
    double f = pow(f, 0.5);
    //double com_distance = 3000/(0.7 *pow(s,0.5));
    // c/H_0 = 3000 MPc
    return f;
}

double RK4_1storder(float y0,double a,double i, double h, double(*fn)(double,double, double, double, double, double), double omega_m, double omega_r, double omega_v, double omega_T){
    
    double t = a + i*h;
    double y = y0;
	double k11, k12, k13, k14; //calculatig slopes
	k11 = h * (*fn)(y0,t, omega_m, omega_r, omega_v, omega_T);

	k12 = h * (*fn)(y0 + 0.5 * k11, t + 0.5 * h,omega_m, omega_r, omega_v, omega_T);

	k13 = h * (*fn)(y0 + 0.5 * k12, t + 0.5 * h, omega_m, omega_r, omega_v, omega_T);

	k14 = h * (*fn)(y0 + k13, t + h, omega_m, omega_r, omega_v, omega_T);

	y += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
    return y;
}

int main(){
    double N = 1000;
    double a = pow(10, -25);// initial value 
    double b = 10; //final value
    double h = (b - a)/N;
    double omega_m = 0.3; 
    double omega_r = 0.00084; 
    double omega_v = 0.3;
    double omega_T = omega_m + omega_r + omega_v;
    double y0 = pow(10, -28); //initial value of a(t)
    double y = y0;
    ofstream outfile1;
    outfile1.open("a(t).txt");
    
    for(int i = 0;i<N;i++){
        cout << y0 << endl;
        y = RK4_1storder(y0,a,i,h,fx,omega_m, omega_r, omega_v, omega_T);
        outfile1 << i << " " << y << endl;
        cout << i << " " << y << endl;
        y0 = y;
    }
    return 0;
}