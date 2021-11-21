#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#define rho_c 2.8e11 // h^2 M_sun Mpc^-3
using namespace std;

int main(int argc, char **argv)
  {
  if(argc != 9){
    printf("\nCommand line args\n");
    printf("./a.out [phistar] [d_phistar] [mstar] [d_mstar] [alfa] [d_alfa] [h]  [#realizations] \n");
    printf("./a.out  0.0048     0.0002     9.96     0.02    -1.33    0.02   0.7       20000 \n");
    return(-1);
  }//if
 
  double phistar, mstar, alpha;
  double phistar_err, mstar_err, alpha_err;
  double up_bound, low_bound, h;
  double vol;
  int i, nr;
  double p, E, d1, d2, d3;;

  phistar = atof(argv[1]);
  phistar_err = atof(argv[2]);
  mstar = atof(argv[3]);
  mstar_err = atof(argv[4]);
  alpha = atof(argv[5]);
  alpha_err = atof(argv[6]);
  h = atof(argv[7]);
  nr = atoi(argv[8]);


  printf("USING GAMMA FUNCTION \n");
  printf("gamma=%g \n",gsl_sf_gamma (alpha+2));
  printf("rho_HI=%g \n",gsl_sf_gamma (alpha+2)*phistar*pow(10,mstar));
  printf("omega_HI=%g  * 1e-4 /h \n",gsl_sf_gamma (alpha+2)*phistar*pow(10,mstar)/rho_c/1e-4/h/h);
  printf("============================================== \n");


  FILE *rfile;
  rfile = fopen("realization_values.dat","w"); 

  double a_i, m_i, p_i, *Omega;
  double sigma, mean;
  mean = 0;
  sigma = 0;
  Omega = (double *)malloc(sizeof(double)*nr);
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
 
  for(i=0;i<nr;i++)
    {
    a_i = alpha + gsl_ran_gaussian (r, alpha_err);
    m_i = mstar + gsl_ran_gaussian (r, mstar_err);
    p_i = phistar + gsl_ran_gaussian (r, phistar_err);
    Omega[i] = gsl_sf_gamma (a_i+2)*p_i*pow(10,m_i)/rho_c/h/h;
    mean = mean + Omega[i];
    fprintf(rfile,"%g \n",Omega[i]);
    }//for i
  mean = 1.0*mean/nr;

  gsl_rng_free (r);

  for(i=0;i<nr;i++)
    {
    sigma = sigma + (mean - Omega[i])*(mean - Omega[i]);
    }//for i
  sigma = sqrt(1.0*sigma/nr);

  printf("nr=%d, mean=%g * 1e-4 /h, sigma=%g * 1e-4 /h \n",nr,mean/1e-4,sigma/1e-4);
  printf("============================================== \n");

  fclose(rfile);
  return(0);
  }