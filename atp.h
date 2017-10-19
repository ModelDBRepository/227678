#ifndef ATP_H_ 
#define ATP_H_

#define V1 0  /* mV */
#define H1 1  
#define N1 2 
#define M1 3  
#define Ca0 4  /* nM */
#define A 5  /* ADP uM */
#define DL 6  

#define N  7/* number of state variables*/

#define PI 3.14159

#define I_NA1 0  
#define I_K1 1 
#define I_L1 2
#define I_LCa 3 //calcium leak 
#define I_CaL 4 //l-type calcium
#define I_NMDA 5 
#define I_KATP 6 
#define I_CAP 7
#define I_GABA 8

#define C 9 /* number of currents*/

#define CM 1.0 /*uF/cm2*/

#define E_NA  60.0
#define E_CA  60.0 //based on voltage curve fit

#define E_K  -90.0
#define E_L  -60.0
#define G_NA 1500.0e-6 
#define G_K 440.0e-6
#define G_L   10.0e-6  

#define TONIC 1 //if 1, G_GABA = constant if 0, G_GABA = 0 except in a 100 ms window.
#define G_GABA 0.0e-6 // 10e6 for 6c, 6E
#define E_GABA -70.0

#define G_LCa  0.05e-6
#define G_CaL 6.000000e-06

#define dloff -40.0 //to match nullcline
#define dlslope 7.5

//calcium pump
#define I_CAPMAX 200*G_CaL //fixing it in terms of G_CaL preserves Ca nullcine.
#define CHALF 500.0
#define CSLOPE 1.0

//morphology
#define FARADAY  96485.0
#define D_S  1.0 //microns
#define fsca 0.02 // 0.03 for 7b

//KATP

#define G_KATP 4.50000e-05

#define S_HALF 5.000000e+03 // 5000 for 4b, 4d, 5b, 6b, 7a-c
#define spower 1.3 //hill coefficient

#define A_HALF 100.0 //small is key
#define ACCU (2.0e5/D_S) // 2.0e5 ~= 4/(2*F) ie 1 calcium out = 1 ADP in
#define KM (43.0) // 45 for 7c

#define NOISE 0.0 //1.0 for 7a and 7b


//NMDA

#define G_NMDA 4.000000e-05 //0 for 4a, 4b

#define NMSLOPE 9.0 //9
#define NMC 0.02 //baseline activation of NMDA channel
#define NMOFF 102.0 //102 larger moves left 250

#define FIXED 0

#define VCLAMP 0

#define PRINT_STATES 1 //if 0 will only return spike times, ISI since last spike
#define WARMUP 0 //time (in ms) model will run before printing to remove transients
#define ENDTIME 10000+WARMUP


#endif
