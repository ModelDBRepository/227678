/* euler integrator with noise that can accept same input/ouput as radau5 for minimal changes
namely - it will take as inputs n (iteration), deriv_, x (time), state, xend, state, out_
and output via out_
because everything is in c we can just pull up deriv_ and out_ form out.c and atp.c repsectively with include
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atp.h"
#include "atp.c"
#include <time.h> //for random seed
double current[C];
double iapp;

int main(int argc, char** argv){
   double dx[N];
   double y[N];
   double noise[N];
   double x;
   double isi;
   double v0=-100;
   double TIME = 0;
   double step;
   double sinf;
   double P[2];
   if(argc > 2){
      P[0] = atof(argv[1]);
      P[1] = atof(argv[2]);
   }
   int n = 0;
   int spikes = 1-PRINT_STATES;
   x = (double) n;
   iapp = 0;
   double h = 1.0e-2; //tiny step because euler, but small model so doesn't take long
   int i;
   for(i=0;i<N;i++) noise[i] = 0.0;
   noise[V1]=NOISE;
   scan_(y,P[1]);
   #ifndef S_HALF
     double S_HALF = P[1];
   #endif

   time_t tseed;
   double caavg = 0;
   while(x < ENDTIME){
      deriv_(&n, &x, y, dx,P);
      if(n%100==0&&!spikes&&x>WARMUP) printf("%e ", (x-WARMUP)*1e-3);
      step = 0.5;
      for(i=0;i<N;i++){
         if(spikes&&i==Ca0&&x>WARMUP) caavg +=y[i];
         y[i] += step*h*(dx[i]+noise[i]*rand_gauss()); //rand_gauss is defined in atp.c
         if(n%100==0&&!spikes&&x>WARMUP) printf("%e ", y[i]);       
      }
      //if(y[0]!=y[0]) break;
      if(n%100==0&&!spikes&&x>WARMUP){
         sinf = 1.0/(1.0+pow(S_HALF/y[A],spower));
         printf("%e ", sinf);
         printf("\n");
      }
      x+=h*step;
      isi+=h*step;
      n++;
      if(v0 > 0 && 0 > y[0] && isi > 10){
        if(spikes && x > WARMUP){
          printf("%e  %e  %e\n",x-WARMUP, 0.001*isi, h*step*caavg/isi);
          caavg = 0;
          fflush(stdout);
        }
        isi = 0;
      }
      v0 = y[0];
   }
   dump_(y);
   return 0;
}
