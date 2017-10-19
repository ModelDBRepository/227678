
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "atp.h"

#ifndef EFUN_
#define EFUN_
double efun(z) 
double z;
{
	if (fabs(z) < 1e-4) {
	return( 1 - z/2);
	}else{
		return( z/(exp(z) - 1));
	}
}
#endif

#ifndef BOLTZ_
#define BOLTZ_
double boltz(v,half,slope)
double v,half,slope;
{
double arg;
arg = -(v-half)/slope;
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(0.0);
else  return(1.0);}
else return(1.0/(1.0 + exp(arg)));}
#endif

double tanhsig(v,half,slope) //sigmoidal tanh function
double v, half, slope;
{
double arg;
arg = (v-half)/slope;
if(arg > 50.0) return 1.0;
else if(arg < -50.0) return 0.0;
else return 0.5*(1.0+tanh(arg));
}

int main(){

double current[C];
double STEMP;
double CATEMP;
double ATEMP,ainf,oinf; 

int check;
double tempcheck;
double step;
int dir;

int i,j;
double V0;
double iapp;
double start;
double anm,bnm;
double sinf,minf,hinf,htau,ntau,mtau,ninf,arg,dlinf,dltau,nminf;
double dminf, dhinf;
double alpm,alph,betm,beth;

double sinf2, diff, difold;

double valtemp;

double CANULL, ANULL; //calcium nullcline

for(V0=-90;V0<10;V0+=0.1){
  dlinf= boltz(V0,dloff,dlslope); //-45 7.5 //L type Calcium
  nminf= NMC+(1.0-NMC)/(1.0+(1.2/NMOFF)*exp(-V0/NMSLOPE)); //NMDA activation

  minf= boltz(V0,-28.0907,8.7264); //sodium activation -28 9
  hinf= boltz(V0,-54.0,-8.7665); //sodium inactivation -54
  ninf= boltz(V0,-20.0,7.0); //fast potassium

  current[I_NA1] = G_NA*pow(minf,3.0)*hinf*(V0 - E_NA);
  current[I_K1] = G_K*pow(ninf,3.0)*(V0 - E_K);
  current[I_L1] =  G_L*(V0 - E_L);
  current[I_LCa] =G_LCa*(V0- E_CA);
  current[I_CaL] =G_CaL*dlinf*(V0- E_CA);
  current[I_NMDA] = G_NMDA*nminf*(V0);
  

  
  //current[I_CaN] = G_CaN*pow(dminf,2.0)*dhinf*(V0- E_CA);

  //iapp = -2.0e-6*(V0+80.0);
  

  
  if(VCLAMP){
  
    CATEMP = 0;
    CATEMP = - current[I_CaL] - current[I_LCa]; //ICAP
    if(CATEMP > I_CAPMAX){
        CANULL = -1;
        ANULL = -1;
        sinf = 1;
    } else {
        ainf = ACCU*CATEMP/KM;
        if (ainf < 1){
          ANULL = ainf*A_HALF/(1.0-ainf);//is this is not [ADP]?
          sinf = 1.0/(1.0+pow(S_HALF/ANULL,spower));

        } else {
        ANULL = -1;
        sinf = 1.0;
        }
      CATEMP = CATEMP/(I_CAPMAX);
      CANULL = CHALF/pow(1.0/CATEMP-1.0,1.0/CSLOPE);
    }

    printf("%e %e %e %e\n", V0,CANULL,ANULL,sinf);
    fflush(stdout);
    continue; // this skips the voltage nullcline
  }



  
  iapp = -G_GABA*(V0 - E_GABA);
  STEMP = iapp  - current[I_NA1] - current[I_K1] - current[I_L1] - current[I_LCa] -current[I_CaL] -current[I_NMDA];
  
  diff = 1.0;
  difold = 1.0;
  CATEMP = 0;
  step = 0.1;
  while(CATEMP < 150.0){ // do a 'blind' Newton search based on sign of net current
    CATEMP += step;
    current[I_CAP] = pow(CATEMP,CSLOPE)/(pow(CATEMP,CSLOPE)+pow(CHALF,CSLOPE)); 
    current[I_CAP] = current[I_CAP]*I_CAPMAX;
    ainf = ACCU*current[I_CAP]/KM;// calculates steady state reuptake
    if(ainf >= 1.0){
         //CATEMP = -1;
         ATEMP = -1;
         current[I_KATP] = G_KATP*(V1-E_K);
         //break;
    } else {
      ATEMP = ainf*A_HALF/(1.0-ainf); //steady state ADP
      sinf = 1.0/(1.0+pow(S_HALF/ATEMP,spower)); //steady state K-ATP activation
      current[I_KATP] = G_KATP*sinf*(V0-E_K);
    }
    difold = diff;
    diff = STEMP - current[I_KATP];
    /*if(V0 > -71.0 && V0 < -69.0){
       printf("%e  %e\n", CATEMP,diff);
    }//*/
    if(diff*difold < 0&&CATEMP>step){
      printf("%e  ", V0);
      printf("%e  ", CATEMP);
      printf("%e  ", ATEMP);
      printf("%e  ", sinf);
      printf("%e\n", STEMP/(G_KATP*(V0-E_K)));
      //step = -0.1*step;
    }
    /*if(CATEMP < 0){
      CATEMP = -1;
      current[I_KATP] = 0;
      //break;
    }*/
    
  }

}

return 0;
}
