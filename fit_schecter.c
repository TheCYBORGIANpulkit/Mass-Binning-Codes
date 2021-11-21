#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#define WORKSPACE 1000

int linesinfile(char *);
double schecter(double,double,double,double);
void find_errors(double *,double *,int,double *,double *,int,char *);

double *chisq;
double chisqmin;

double minm,maxm,mina,maxa,minphistar,maxphistar;
double *mstar,*alpha,*phistar,dm,da,dphi;

int ndata,nbin;
FILE *infile;

double norm;
double dlogm;
double vol,*logmh1,*logphi,*counts,*logerrphi;
double *hist, *phi;


int nsum; //number of points in the tail of the counts used for normalisation  
double ntot,phi2dswmltot;

int lflag;// likelyhoodflag

int main(int argc, char **argv){
  int i,j,k,l,q;

  if(argc != 10){
    printf("\nCommand line args\n");
    printf("./a.out [minm] [maxm] [mina] [maxa] [minphi] [maxphi] [nbin] [file_index]  [lflag]\n");
    printf("./a.out   9.0   10.5  -2.0    -1.0    0.001    0.009    50       1-6          0    \n");
    printf("lflag: 0 for computing errors and 1 to skip ");
    return(-1);
  }//if


  minm = atof(argv[1]);
  maxm = atof(argv[2]);
  mina = atof(argv[3]);
  maxa = atof(argv[4]);
  minphistar = atof(argv[5]);
  maxphistar = atof(argv[6]);
  nbin = atoi(argv[7]);
  q = atoi(argv[8]);
  lflag = atoi(argv[9]);

  printf("\n\n");
  printf("Command line args ...\n");
  printf("./a.out %g %g %g %g %g %g %d %d %d\n",
         minm,maxm,mina,maxa,minphistar,maxphistar,nbin,q,lflag);
  vol = 0.63800*pow(211.7039,3)/3;  //new volume
  printf("vol: %g \n",vol);

  
 

  //read file and fix normalization
  FILE *infile;
  char *infilename;
  infilename = (char *)malloc(sizeof(char *)*1000);
  sprintf(infilename,"region%d_new.dat",q);
  ndata = linesinfile(infilename);
  printf("ndata=%d \n",ndata);
  infile = fopen(infilename,"r");
  printf("Opening %s \n",infilename);

  logmh1 = (double *)malloc(sizeof(double)*ndata);
  logphi = (double *)malloc(sizeof(double)*ndata);
  logerrphi = (double *)malloc(sizeof(double)*ndata);
  counts = (double *)malloc(sizeof(double)*ndata);
  phi = (double *)malloc(sizeof(double)*ndata);
  hist = (double *)malloc(sizeof(double)*ndata);
  
  for(i=0;i<ndata;i++){
    fscanf(infile,"%lf %lf %lf",&logmh1[i],&logphi[i],&logerrphi[i]);
    printf("%g %g %g \n",logmh1[i],logphi[i],logerrphi[i]);
    phi[i] = pow(10,logphi[i]);
  }//for...i

  fclose(infile);
  free(infilename);


  //------fix the parameter space for minimization---------------

  mstar = (double *)malloc(sizeof(double)*nbin);
  alpha = (double *)malloc(sizeof(double)*nbin);
  phistar = (double *)malloc(sizeof(double)*nbin);

  chisq = (double *)malloc(sizeof(double)*nbin*nbin*nbin);
    
  for(i=0;i<nbin;i++)
    for(j=0;j<nbin;j++)
      for(k=0;k<nbin;k++)
	chisq[k+nbin*(j+nbin*i)] = 0.0;
  
  mstar[0] = minm;
  dm = (maxm-minm)/(double)nbin;
  printf("minm=%g maxm=%g dm=%g\n",minm,maxm,dm);

  alpha[0] = mina;
  da = (maxa-mina)/(double)nbin;
  printf("mina=%g maxa=%g da=%g\n",mina,maxa,da);

  phistar[0] = minphistar;
  dphi = (maxphistar-minphistar)/(double)nbin;
  printf("minphistar=%g maxphistar=%g dphi=%g\n",minphistar,maxphistar,dphi);
  
  printf("\ndm=%g da=%g dphi=%g \n\n",dm,da,dphi);
  

  for(i=1;i<nbin;i++){
    mstar[i] = mstar[i-1]+dm;
    alpha[i] = alpha[i-1]+da;
    phistar[i] = phistar[i-1]+dphi;
    printf("m=%g a=%g phi=%g\n",mstar[i],alpha[i],phistar[i]);
  }
  //--FIXED ------parameter space 

  //---BEGIN------CHISQ MINIMIZATION------------------
  double mstarbf,alphabf,phistarbf;
  int imstar_bf,jalpha_bf,kphistar_bf;
  chisqmin = 1e80;
  
  for(i=0;i<nbin;i++){     //mstar
    for(j=0;j<nbin;j++){   //alpha
      for(k=0;k<nbin;k++){ //phistar

	for(l=0;l<ndata;l++){
	  chisq[k+nbin*(j+nbin*i)] += 
	    pow((schecter(logmh1[l],alpha[j],mstar[i],phistar[k])-logphi[l])/logerrphi[l],2.0);
	  
	}//for...l
	if(chisq[k+nbin*(j+nbin*i)] <= chisqmin){
          chisqmin = chisq[k+nbin*(j+nbin*i)];

          mstarbf = mstar[i];
          imstar_bf = i;

          alphabf = alpha[j];
          jalpha_bf = j;

          phistarbf = phistar[k];
          kphistar_bf = k;

	}//if...

      }//for...k
    }//for...j
    printf("DONE %d of %d \n",i,nbin);
    fflush(stdout);
  }//for...i

  printf("mstarbf=%g(i=%d) alphabf=%g(j=%d) phistarbf=%g(k=%d)chisqmin=%g ndata=%d chisqmin_red=%g \n",mstarbf,imstar_bf,
	 alphabf,jalpha_bf,phistarbf,kphistar_bf,chisqmin,ndata,chisqmin*1.0/(ndata-3));

  FILE *ofile;
  char *ofilename;
  ofilename = (char *)malloc(sizeof(char *)*1000);
  sprintf(ofilename,"chi_region%d.dat",q);
  ofile = fopen(ofilename,"w");
  printf("Opening %s \n",ofilename);

  double *data_chi;
  data_chi = (double *)malloc(sizeof(double)*ndata);

  
  
  for(i=0;i<nbin;i++)  //mstar
    for(j=0;j<nbin;j++){ //alpha
      for(k=0;k<nbin;k++){//phistar
	lhood[k+nbin*(j+nbin*i)] = exp(-(chisq[k+nbin*(j+nbin*i)]-chisqmin)*0.5); //could have gotten rid of 0.5 
                                                              //while calculating chisq
	lhood_tot += lhood[k+nbin*(j+nbin*i)];
      }//for ...k
    }//for...j  

  printf("\nlhood_tot=%g \n",lhood_tot);
  printf("Normalising 3d likelyhood\n");

  //normalise 3d likelyhood
  for(i=0;i<nbin;i++)  //mstar
    for(j=0;j<nbin;j++) //alpha
      for(k=0;k<nbin;k++) //alpha
	lhood[k+nbin*(j+nbin*i)] /= lhood_tot; 

  //marginalise
  double *lmstar,*lalpha,*lphistar;
  lmstar = (double *)malloc(sizeof(double)*nbin);
  lalpha = (double *)malloc(sizeof(double)*nbin);
  lphistar = (double *)malloc(sizeof(double)*nbin);

  for(i=0;i<nbin;i++){
    lmstar[i] = 0.0;
    lalpha[i] = 0.0;
    lphistar[i] = 0.0;
  }//for...i

  for(i=0;i<nbin;i++)  //mstar
    for(j=0;j<nbin;j++) //alpha
      for(k=0;k<nbin;k++){ //alpha

	lmstar[i] += lhood[k+nbin*(j+nbin*i)];
	lalpha[j] += lhood[k+nbin*(j+nbin*i)];
	lphistar[k] += lhood[k+nbin*(j+nbin*i)];
	
      }//for...k


  FILE *lfile;
  char *lfname;
  lfname = (char *)malloc(sizeof(char)*100);
  sprintf(lfname,"1dlikelyhood_%04d.dat",nbin);
  lfile = fopen(lfname,"w");
  
  for(i=0;i<nbin;i++)
    fprintf(lfile,"%g %g %g %g %g %g \n",
            mstar[i],lmstar[i],alpha[i],lalpha[i],phistar[i],lphistar[i]);
  free(lfname);
  fclose(lfile);

  
  if(lflag == 0){  //compute errors if likelyhood flag is set to zero

    //-----------------Get 1sigma errors --> crude way --------------------------
    //now find 1sigma errors this is an estimate
    double lmstar_max,lalpha_max,lphistar_max;
    int lmstar_index,lalpha_index,lphistar_index;

    lmstar_max = lalpha_max = lphistar_max = -100.;
    for(i=0;i<nbin;i++){
      //mstar
      if(lmstar_max <= lmstar[i]){
	lmstar_max = lmstar[i];
	lmstar_index = i;
      }//if

      //alpha
      if(lalpha_max <= lalpha[i]){
	lalpha_max = lalpha[i];
	lalpha_index = i;
      }//if

      //phistar
      if(lphistar_max <= lphistar[i]){
	lphistar_max = lphistar[i];
	lphistar_index = i;
      }//if
    
    }//for...i
  
    printf("\n");
    printf("lmstar_max=%g lmstar_index=%d \n",lmstar_max,lmstar_index);
    printf("lalpha_max=%g lalpha_index=%d \n",lalpha_max,lalpha_index);
    printf("lphistar_max=%g lphistar_index=%d \n",lphistar_max,lphistar_index);
    printf("\n");

    double sigma_g;
    sigma_g = 68.2/100;
    sigma_g *= 0.5;
  
    double lcum_right,lcum_left;
    int lindex;
    double mstar1,mstar2;
    double alpha1,alpha2;
    double phistar1,phistar2;

    //mstar
    lcum_right = 0.0;
    lindex = lmstar_index;
    while(lcum_right <= sigma_g){
      lcum_right += lmstar[lindex];
      lindex++;
    }
    mstar2 = mstar[lindex];
  
    lcum_left = 0.0;
    lindex = lmstar_index;
    while(lcum_left <= sigma_g){
      lcum_left += lmstar[lindex];
      lindex--;
    }
    mstar1 = mstar[lindex];

    printf("mstarbf=%g mstar1=%g mstar2=%g \n",mstarbf,mstar1,mstar2);

    //alpha
    lcum_right = 0.0;
    lindex = lalpha_index;
    while(lcum_right <= sigma_g){
      lcum_right += lalpha[lindex];
      lindex++;
    }
    alpha2 = alpha[lindex];
  
    lcum_left = 0.0;
    lindex = lalpha_index;
    while(lcum_left <= sigma_g){
      lcum_left += lalpha[lindex];
      lindex--;
    }
    alpha1 = alpha[lindex];

    printf("alphabf=%g alpha1=%g alpha2=%g \n",alphabf,alpha1,alpha2);

    //phistar
    lcum_right = 0.0;
    lindex = lphistar_index;
    while(lcum_right <= sigma_g){
      lcum_right += lphistar[lindex];
      lindex++;
    }
    phistar2 = phistar[lindex];
  
    lcum_left = 0.0;
    lindex = lphistar_index;
    while(lcum_left <= sigma_g){
      lcum_left += lphistar[lindex];
      lindex--;
    }
    phistar1 = phistar[lindex];

    printf("phistarbf=%g phistar1=%g phistar2=%g \n",phistarbf,phistar1,phistar2);

    //a better and general way to compute the sigmas 
    int nsigma;
    double *sigma,*bf;
    nsigma = 3;
    sigma = (double *)malloc(sizeof(double)*nsigma);
    sigma[0] = 68.2;
    sigma[1] = 95.4;
    sigma[2] = 99.7;
    bf = (double *)malloc(sizeof(double)*nsigma*3);

  
    find_errors(mstar,lmstar,nbin,bf,sigma,nsigma,"mstar");

    find_errors(alpha,lalpha,nbin,bf,sigma,nsigma,"alpha");

    find_errors(phistar,lphistar,nbin,bf,sigma,nsigma,"phistar");

  
    free(bf);
    free(sigma);
  }//if lflag == 0

  
  free(lmstar);
  free(lalpha);
  free(lphistar);

  free(chisq);
  free(mstar);
  free(alpha);
  free(phistar);

  free(logmh1);
  free(logphi);
  free(logerrphi);
  free(counts);
  
  return(0);
}//main()

void find_errors(double *param,double *lhood,int np,double *bfparam,double *sigma,int nsigma,char *paramfield){
  
  double dparam;
  double lmax;
  int i,imax;
  gsl_interp_accel *acc;
  gsl_interp *akima;

  double *mysigma;
  mysigma = (double *)malloc(sizeof(double)*nsigma);

 
  for(i=0;i<nsigma;i++)
    mysigma[i] = sigma[i]/200.;    //factor of 2 is for 1 side

  lmax = -100;
  for(i=0;i<np;i++)
    if(lmax <= lhood[i]){
      lmax = lhood[i];
      imax = i;
    }

  printf("\nbf[%d]=%g lhood=%g \n",imax,param[imax],lhood[imax]);
 
  dparam = param[1]-param[0];

  //left
  int npleft;
  double *lhood_left,*ltot_arr_left,*param_left;
  double ltot_left;
  double *xleft,*bfleft;
  
  npleft = imax+1;

  ltot_arr_left = (double *)malloc(sizeof(double)*npleft);
  lhood_left = (double *)malloc(sizeof(double)*npleft);
  param_left = (double *)malloc(sizeof(double)*npleft);

  for(i=0;i<npleft;i++){
    lhood_left[i] = lhood[i];
    param_left[i] = param[i];
  }

  ltot_left = 0.0;
  ltot_arr_left[0] = 0.0;
  for(i=1;i<npleft;i++){
    ltot_left +=  (lhood_left[i]+lhood_left[i-1])*dparam/2.0;
    ltot_arr_left[i] = ltot_left;
  }
  
  for(i=0;i<npleft;i++)
    ltot_arr_left[i] = 0.5*ltot_arr_left[i]/ltot_left;
 
      
  acc = gsl_interp_accel_alloc();
  akima = gsl_interp_alloc(gsl_interp_akima,npleft);
  gsl_interp_init(akima,ltot_arr_left,param_left,npleft);

  xleft = (double *)malloc(sizeof(double)*nsigma);
  bfleft = (double *)malloc(sizeof(double)*nsigma);

  for(i=0;i<nsigma;i++){
    xleft[i] = 0.5 - mysigma[i];
    bfleft[i] = gsl_interp_eval(akima,ltot_arr_left,param_left,xleft[i],acc);
  }

  free(xleft);
  free(lhood_left);
  free(param_left);
  free(ltot_arr_left);
  gsl_interp_free(akima);
  gsl_interp_accel_free(acc);
  //left


  //right
  int npright;
  double *lhood_right,*ltot_arr_right,*param_right;
  double ltot_right;
  double *xright,*bfright;
  
  npright = np-imax;

  ltot_arr_right = (double *)malloc(sizeof(double)*npright);
  lhood_right = (double *)malloc(sizeof(double)*npright);
  param_right = (double *)malloc(sizeof(double)*npright);

  for(i=imax;i<np;i++){
    lhood_right[i-imax] = lhood[i];
    param_right[i-imax] = param[i];
  }

  ltot_right = 0.0;
  ltot_arr_right[0] = 0.0;
  for(i=1;i<npright;i++){
    ltot_right +=  (lhood_right[i]+lhood_right[i-1])*dparam/2.0;
    ltot_arr_right[i] = ltot_right;
  }
  
  for(i=0;i<npright;i++)
    ltot_arr_right[i] = 0.5*ltot_arr_right[i]/ltot_right;

      
  acc = gsl_interp_accel_alloc();
  akima = gsl_interp_alloc(gsl_interp_akima,npright);
  gsl_interp_init(akima,ltot_arr_right,param_right,npright);

  xright = (double *)malloc(sizeof(double)*nsigma);
  bfright = (double *)malloc(sizeof(double)*nsigma);

  for(i=0;i<nsigma;i++){
    xright[i] = mysigma[i];
    bfright[i] = gsl_interp_eval(akima,ltot_arr_right,param_right,xright[i],acc);

  }
  free(xright);
  free(lhood_right);
  free(param_right);
  free(ltot_arr_right);
  gsl_interp_free(akima);
  gsl_interp_accel_free(acc);
  //right


  for(i=0;i<nsigma;i++){
    bfparam[0+3*i] = param[imax];
    bfparam[1+3*i] = bfleft[i];
    bfparam[2+3*i] = bfright[i];
  }
  for(i=0;i<nsigma;i++)
    printf("%s %dsigma: %g %g %g %g %g \n",
	   paramfield,i+1,bfparam[0+3*i],bfparam[1+3*i],bfparam[2+3*i],bfparam[0+3*i]-bfparam[1+3*i],bfparam[0+3*i]-bfparam[2+3*i]);


  FILE *bfparamfile;
  char *bfparamfname;
  
  for(i=0;i<nsigma;i++){
    bfparamfname = (char *)malloc(sizeof(char)*300);
    sprintf(bfparamfname,"%s_%dsigma.dat",paramfield,i+1);
    bfparamfile = fopen(bfparamfname,"w");
    
    fprintf(bfparamfile,"%g %g %g \n",
	    bfparam[0+3*i],bfparam[1+3*i],bfparam[2+3*i]);
    
    fclose(bfparamfile);
    free(bfparamfname);
  }
  //
  
  
  
  free(bfleft);
  free(bfright);
  free(mysigma);


}//find_errors




//------------------schecter()-------------------------------
//returns the schecter function curve
//mh1 is in log10
//mh1starparam is in log10
//schecter() -->  will return log10[logten*phistar*exp(-m/mstar)*(m/mstar)^(alpha+1.)]
//----------------------------------------------------------

double schecter(double mh1,double alphaparam,double mh1starparam,double phistarparam){
  double schecter_model;
  double logten = log(10.0);
  

  schecter_model = log10(logten);
  schecter_model += log10(phistarparam);
  schecter_model += (alphaparam + 1.0)*(mh1-mh1starparam);
  schecter_model += log10(exp(-1*pow(10.,mh1-mh1starparam)));

  return(schecter_model);

}//schecter


int linesinfile(char *infname){
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
  
  printf("nlines=%d in %s \n",nlines,infname);
    
  return(nlines);

}//linesinfile



