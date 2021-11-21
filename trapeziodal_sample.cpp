#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>


using namespace std;

double f(double x){
    double s = 0.3*(pow(1+x,3)) + 0 + 0.7;
    //second term is because \Omega_k = 0;
    double E = 3000/pow(s,0.5);
    return E;
}

int main()
{
    int n,i;        //n is for subintervals and i is for loop
    double a,b,h,sum=0,integral;
    cout<<"Enter the limits of integration,\nInitial limit,a=";    //get the limits of integration
    cin>>a;
    cout<<"Final limit, b=";
    cin>>b;
    cout<<"Enter the no. of subintervals, n=";            //get the no. of subintervals
    cin>>n;
    double x[n+1],y[n+1];
    h=(b-a)/n;                //get the width of the subintervals
    for (i=0;i<=n;i++)
    {                    //loop to evaluate x0,...xn and y0,...yn
        x[i]=a+i*h;            //and store them in arrays
        y[i]=f(x[i]);
    }
    for (i=1;i<n;i++)            //loop to evaluate h*(y1+...+yn-1)
    {
        sum=sum+h*y[i];
    }
    integral=h/2.0*(y[0]+y[n])+sum;        //h/2*[y0+yn+2(y1+y2+y3+...yn-1)]
    cout<<"The definite integral  is "<<integral<<endl;
    return 0;
}




