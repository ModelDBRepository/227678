//4 CALCIUM compartments with dynamic buffering

#ifndef ATP_C_
#define ATP_C_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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

float rand_gauss (void) {
  float v1,v2,s;

  do {
    v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0;
  else
    return (v1*sqrt(-2.0 * log(s) / s));
}
 
 
double  gaussian(v,a,b,c,d) 
double  v,a,b,c,d ;
{
double arg;
arg = pow(((v-c)/b),2.0);
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(d);
else  return(d+a);}
else return  (d + a*exp(-arg));
}

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



int deriv_(np,xp,Y,F,P)
double *F,*Y;
int *np;
double *xp;
double *P;
{
#ifndef G_NMDA
double G_NMDA = P[0];
#endif

#ifndef S_HALF
double S_HALF = P[1];
#endif

double SINF = 0.0;



extern double current[C]; 
int el;
double time_;
double iapp; //=I_APP;
double anm,bnm;
double minf,hinf,htau,ntau,mtau,ninf,arg,dlinf,dltau,nminf;
//atau
time_ = *xp;
double fadp,cadp,sinf;
el = *np;

if(TONIC){
  iapp = -G_GABA*(Y[V1]-E_GABA);
} else {
  iapp = 0;
}
if (!(time_<1850+WARMUP || time_>1900+WARMUP)){
  //iapp = 0;
  iapp = -G_GABA*(Y[V1]-E_GABA);
}

//*/



//HH channels
minf= boltz(Y[V1],-28.0907,8.7264); //sodium activation -28 9
hinf= boltz(Y[V1],-54.0,-8.7665); //sodium inactivation -54
ninf= boltz(Y[V1],-20.0,7.0); //fast potassium
//ninf= boltz(Y[V1],-25.0,12.0); //fast potassium


htau = 56.0*(boltz(Y[V1],-21.0,-4.5) - boltz(Y[V1],-41.0,-2.0))+1.0;
mtau= (0.01 + 1.0/((15.6504+0.4043*Y[V1])/(1.0-exp(-19.565-0.50542*Y[V1])) +3.0212*exp(-7.463e-3*Y[V1])));
if(Y[V1]>-60.0 ) ntau= 1.0 + 19.0*exp(-pow(20*log(1.0 +0.05*(Y[V1]+40.0)),2.0)/300.0); //1 19 40 the 20 does all the work
else ntau=1.0; //1.0


// L-type
dlinf= boltz(Y[V1],dloff,dlslope); //-45 7.5 //L type Calcium
dltau= 1.0/(-0.020876*(Y[V1]+39.726)/(exp(-(Y[V1]+39.726)/4.711)-1) + 0.19444*exp(-(Y[V1]+15.338)/224.21));


//NMDA

nminf = NMC+(1.0-NMC)/(1.0+(1.2/NMOFF)*exp(-Y[V1]/NMSLOPE)); //NMDA activation


F[Ca0] = -2.0e6*fsca*(current[I_LCa] + current[I_CaL] +current[I_CAP])/(D_S*0.0001*FARADAY); //calcium currents



current[I_NA1] = G_NA*pow(Y[M1],3.0)*Y[H1]*(Y[V1] - E_NA);
current[I_K1] = G_K*pow(Y[N1],3.0)*(Y[V1] - E_K);

current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_LCa] =G_LCa*(Y[V1]- E_CA);

current[I_CaL] =G_CaL*Y[DL]*(Y[V1]- E_CA);

current[I_NMDA] = G_NMDA*nminf*(Y[V1]);



current[I_CAP] = pow(Y[Ca0],CSLOPE)/(pow(Y[Ca0],CSLOPE)+pow(CHALF,CSLOPE)); 

current[I_CAP] = current[I_CAP]*I_CAPMAX;


F[A] = (ACCU*current[I_CAP] - KM*Y[A]/(A_HALF+Y[A])); //unscaled with I_CAPMAX without loss of generality

if(FIXED){ 
  sinf = SINF;
} else {
  sinf = 1.0/(1.0 + pow(S_HALF/(Y[A]),spower)); // Hill function
}

current[I_KATP] = G_KATP*sinf*(Y[V1]- E_K);

F[V1] = 1000*(iapp - current[I_NA1] - current[I_K1] - current[I_L1] - current[I_LCa] -current[I_CaL] -current[I_NMDA]- current[I_KATP])/CM;



F[M1] = (minf-Y[M1])/mtau;
F[H1] = (hinf-Y[H1])/htau;
F[N1] = (ninf-Y[N1])/ntau;
F[DL] = (dlinf-Y[DL])/dltau;


return 0;
   }


void scan_(Y,P) 
double Y[N];
double P;
{FILE *fopen(),*sp;
int i;
sp = fopen("state.data","r");
for(i=0;i<N;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);
if(FIXED){
  Y[V1] = P;
  Y[M1] = boltz(Y[V1],-28.0907,8.7264); //sodium activation -28 9
  Y[H1] = boltz(Y[V1],-54.0,-8.7665); //sodium inactivation -54
  Y[N1] = boltz(Y[V1],-20.0,7.0); //fast potassium
  Y[DL] = boltz(Y[V1],dloff,dlslope); //-45 7.5 //L type Calcium
}
}

void dump_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
if(Y[0]==Y[0]){ //sanity check on jawns
  sp = fopen("end.data","w");
  for(i=0;i<N;i++) fprintf(sp,"%.16f\n",Y[i]);
  fclose(sp);
  }
}

int mas(n,amas,l)
        int *n;
        double *amas;
        int *l;
{return 0;}

int dummy(n,t,y,ydot)
        int *n;
        double *t;
        double *y;
        double *ydot;
{return 0;}

#endif

