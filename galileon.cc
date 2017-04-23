// Modified by Clement Leloup 09-12-2016 :
// Added the interpolation global variables and the interpolation in arrays_
// Modified by Clement Leloup 06-12-2016 :
// Added functions grhogal and gpresgal
// Changed name of file to "galileon.cc"
// Added all functions calculating perturbations part
// Modified by Clement Leloup 04-12-2016 :
// Created header and moved everything needed there
// Modified by Clement Leloup 05-10-2016 :
// Changed name of dtauda to arrays
// Added the function handxofa
// Added three global vectors intvar, hubble and x and modified everything to take it into account
// Moved common parameters to global
// Modified by Clement Leloup 29-09-2016 :
// Made the code a little more C-like (removing strings in favor of char*)
// Modified main and dtauda
// Changed type of readParamsFile type to array of char*
// Added the dotphi parameter
// Cleaned a little
// Modified by Clement Leloup 27-09-2016 :
// Changed type of dtauda to double array (here double pointer)
// Modified by Clement Leloup 22-07-2016 :
// Added the function dtauda to be linked with CAMB
// Cleaned a little
// Added the function readParamsFile
// Modified main
// Modified by Clement Leloup 23-05-2016 :
// Added the analytic tracker solution in main function
// Added the function calcHubbleTracker
// Modified by Clement Leloup 13-04-2016 :
// Added the "coord" parameter to calcHubbleGalileon, to take into account whether we use a or z to do the calculation. In case one use a to calculate, the calculation is done backward from a=1.
// Modified the main function. Now, the behaviour is adapted whether the user prefers a or z to calculate the observables.
// Modified by Clement Leloup 07-04-2016 :
// Added the contribution of orad to the evolution of Lander system when !useacoord
// Wrote of the main function to test and compare with former results
// Commented of every function except CalcValOmC2C3C4C5CGC0 and calcHubbleGalileon

// From modified CosFitter for Galileon Cosmology file lumdist.cc by Jeremy Neveu


#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include "fstream"
#include <algorithm>
#include <math.h>
#include "string.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;


// Types of solving method
enum solvmethod {JNEVEUz, JNEVEUa, BARREIRA};

// Global vectors of the integration variable (a or z), h and x
std::vector<double> intvar;
std::vector<double> hubble;
std::vector<double> x;

// Interpolation tools
gsl_interp_accel* acc;
gsl_spline* spline_h;
gsl_spline* spline_x;

// Global galileon parameters
solvmethod coord = JNEVEUz;
double om = 0;
double c2 = 0;
double c3 = 0;
double c4 = 0;
double c5 = 0;
double cG = 0;
double c0 = 0;
double orad = 0;
double h0 = 0;

// Speed of light
double lightspeed = 2.99792458e8;


// Function to read parameters from a file
char** readParamsFile(const char* filename){

  static char* params[14]; // array of parameters
  char** param = NULL;

  for(int i = 0; i<14; i++){
    params[i] = ""; // initialization of the parameters
  }

  char* line = NULL;
  size_t len = 0;
  const char* com = "#";
  FILE* stream = fopen(filename, "r");
  if(stream != NULL){
    while(getline(&line, &len, stream) != -1){
      if(line != NULL){
	if(line[0] != com[0]){
	  char input[strlen(line)];
	  strcpy(input,line);
	  char* token = strtok(input, " =#");
	  if(string(token) == "c2"){
	    param = &params[0];
	  } else if(string(token) == "c3"){
	    param = &params[1];
	  } else if(string(token) == "c4"){
	    param = &params[2];
	  } else if(string(token) == "c5"){
	    param = &params[3];
	  } else if(string(token) == "c0"){
	    param = &params[4];
	  } else if(string(token) == "cG"){
	    param = &params[5];
	  } else if(string(token) == "ratio_rho"){
	    param = &params[6];
	  } else if(string(token) == "age"){
	    param = &params[7];
	  } else if(string(token) == "ombh2"){
	    param = &params[8];
	  } else if(string(token) == "omch2"){
	    param = &params[9];
	  } else if(string(token) == "solvingMethod"){
	    param = &params[10];
	  } else if(string(token) == "dotphi"){
	    param = &params[11];
	  } else if(string(token) == "hubble"){
	    param = &params[12];
	  } else if(string(token) == "output_root"){
	    param = &params[13];
	  }
	  token = strtok(NULL, " =");
	  if(token!=NULL){
	    
	    char* ch = token;
	    if(param != NULL){
	      *param = (char*)malloc(sizeof(char)*strlen(ch));
	      strcpy(*param,ch);
	      param = 0;
	      }
	  }
	}
      }
    }
    if(line) free(line);
  }

  // for(int i = 0; i < 13; i++){
  //   std::cout << "param " << i << " = " << params[i] << endl;;
  // }

  return params;
  
}


/*
  \param[in] z Redshift
  \param[in] y Current value of integral : y[0] is h(z), y[1] is dpi/dz and y[2] is pi
  \param[in] params Array of parameters (\f$\Omega_m, c_2, c_3, c_4, c_5, c_G, c_0, orad, useacoord\f$)
  \param[out] f The value of the integral term : f[0] is dh/dz, f[1] is d^2pi/dz^2 and f[2] is dpi/dz
*/
int calcValOmC2C3C4C5CGC0(double z, const double y[3], double f[3], void* params){

  double alpha,gamma,beta,sigma,lambda,omega,OmegaP,OmegaM,OmTest;
  const double *in_params = static_cast<const double*>(params);
  bool useacoord = (int)in_params[0];
  double *OmegaPiOverOrad = NULL;
  if(useacoord) OmegaPiOverOrad = const_cast<double*>(&in_params[8]);

  //if (useacoord) std::cout << "Using a instead of z" << std::endl; 
  double zpo;
  double a=1;
  if(!useacoord){
    zpo=(1.0+z); //df/dln(a)=zpo*f'(z or a)
  } else{
    // We consider z as a :
    a = z;
    zpo=-a;
  }
  int status = 0;

  double prod = -zpo*y[0]*y[1]; // prod=h*x
  double prod2 = prod*prod;
  double prod3 = prod*prod2;
  double prod4 = prod*prod3;
  double prod5 = prod*prod4;
  double h = y[0];
  double h2 = h*h;
  double h3 = h*h2;
  double h4 = h2*h2;

  OmegaP = ( 6.0*c0*(h*prod+y[2]*h2) + c2/2.0*prod2 - 6.0*c3*h*prod3 + 22.5*c4*h2*prod4 - 21.0*c5*h3*prod5 - 9*cG*h2*prod2  )/(3.0*h2) ;
  if(!useacoord){
    OmegaM = 1 - OmegaP - orad*pow(zpo, 4)/h2;
    OmTest = om*zpo*zpo*zpo/h2; 
  } else{
    OmegaM = 1 - OmegaP - orad/(pow(a,4)*h2);
    OmTest = om/(pow(a,3)*h2);
  }

  // The equations : 
  alpha = c0*h + c2/6.0*prod - 3*c3*h*prod2 + 15*c4*h2*prod3 - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  gamma = 2*c0*h2 + c2/3.0*h*prod - c3*h2*prod2 + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  beta = -2*c3*h3*prod + c2/6.0*h2 + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  sigma = 2.0*( 1.0 - 2*c0*y[2] )*h - 2.0*c0*prod + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  lambda = 3.0*( 1 - 2*c0*y[2] )*h2 - 2.0*c0*h*prod - 2.0*c3*h*prod3 + c2/2.0*prod2 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;
  omega = -2*c0*h2 + 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;
  if(useacoord){
    lambda += orad/pow(a,4);
  } else{
    lambda += orad*pow(zpo,4);
  }

  double denom = sigma*beta - alpha*omega;
  f[0] = (omega*gamma - lambda*beta) / (-zpo*denom);
  f[1] = (alpha*lambda - sigma*gamma) / (zpo*zpo*denom) ;
  f[2] = y[1];///(-zpo);
  if(useacoord) f[1] += -2.0*y[1]/a;

  if(a==0.001 && useacoord){
    *OmegaPiOverOrad = OmegaP / ( orad/(pow(a,4)*h2) );
  }
  if(useacoord) double dhda = f[0];
  return GSL_SUCCESS;
}



/*!
  Evaluates the h(z) in galileon cosmology

  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_{rad}\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[in] coord, whether to use a or z coordinate
  \param[out] table of a or z coordinates where h is evaluated
  \param[out] table of h
  \param[out] table of x
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghost condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : noghost condition not satisfied
		       \status = 7 : imaginary speed of sound
		       \status = 8 : failed to integrate
		       \status = 9 : one of the global vectors is empty
*/
int calcHubbleGalileon(){

  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    return 9;
  }

  double params[1];
  bool useacoord = (coord == BARREIRA || coord == JNEVEUa);
  params[0] = (int)useacoord;

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-14, 1e-14, 1, 1);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(3);
  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;
 
  if(useacoord){
    double y[3] = { 1.0, 1.0, 0.0 };  //Inital value of integral
    double z = 1;
    int nstep = min(hubble.size(),x.size());
    double h = -1e-6; //Initial step guess
    double zcurrtarg;
    int st;
    for(int i = 0; i < nstep; ++i){
      zcurrtarg = intvar[i];
      while(z > zcurrtarg){
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
	double OmegaP = (0.5*c2*pow(y[0], 2)*pow(intvar[i]*y[1], 2) - 6*c3*pow(y[0], 4)*pow(intvar[i]*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(intvar[i]*y[1], 4) - 21*c5*pow(y[0], 8)*pow(intvar[i]*y[1], 5) - 9*cG*pow(y[0], 4)*pow(intvar[i]*y[1], 2))/(3.0*pow(y[0], 2)) ;
	double OmegaM = 1 - OmegaP - orad/(pow(intvar[i],4)*pow(y[0], 2));
	double OmTest = om/(pow(intvar[i],3)*pow(y[0], 2));
	if(OmegaP<0) st = 5;
	if ( fabs(OmegaM - OmTest)>1e-4  ) {
	  printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, intvar[i], pow(y[0], 2), fabs(OmegaM - OmTest));
	  st = 4;
	}
	if(st != 0){
	  gsl_odeiv_evolve_free(e);
	  gsl_odeiv_control_free(c);
	  gsl_odeiv_step_free(s);
	  return st;
	}
      }
      if(isnan(fabs(y[0])) || isnan(fabs(y[1]))  || isnan(fabs(y[2]))){
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
	printf("\nFailure with om = %f, c2 = %f, c3 = %f, c4 = %f, c5 = %f, cG = %f and c0 = %f at z = %f, h = %f, x = %f and y = %f\n", om, c2, c3, c4, c5, cG, c0, zcurrtarg, y[0], y[1], y[2]);
	return 8;
      }
      hubble[i] = y[0];
      x[i] = y[1];
    }
  } else{
    double y[3] = {1.0, -1.0, 0.0};  //Inital value of integral
    double z = 0;
    int nstep = min(hubble.size(),x.size());
    double h = 1e-6; //Initial step guess
    double zcurrtarg;
    int st;
    for(int i = 0; i < nstep; ++i){
      zcurrtarg = intvar[i];
     while(z < zcurrtarg){
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
	double OmegaP = (0.5*c2*pow(y[0], 2)*pow(-(1+intvar[i])*y[1], 2) - 6*c3*pow(y[0], 4)*pow(-(1+intvar[i])*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(-(1+intvar[i])*y[1], 4) - 21*c5*pow(y[0], 8)*pow(-(1+intvar[i])*y[1], 5) - 9*cG*pow(y[0], 4)*pow(-(1+intvar[i])*y[1], 2))/(3.0*pow(y[0], 2)) ;
	double OmegaM = 1 - OmegaP - orad*pow(1+intvar[i], 4)/pow(y[0], 2);
	double OmTest = om*pow(1+intvar[i], 3)/pow(y[0], 2);
	if(OmegaP<0) st = 5;
	if ( fabs(OmegaM - OmTest)>1e-4  ) {
	  printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, intvar[i], pow(y[0], 2), fabs(OmegaM - OmTest));
	  st = 4;
	}
	if(st != 0){
	  gsl_odeiv_evolve_free(e);
	  gsl_odeiv_control_free(c);
	  gsl_odeiv_step_free(s);
	  return st;
	}
      }
      if(isnan(fabs(y[0])) || isnan(fabs(y[1]))  || isnan(fabs(y[2]))){
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
	printf("\nFailure with om = %f, c2 = %f, c3 = %f, c4 = %f, c5 = %f, cG = %f and c0 = %f at z = %f, h = %f, x = %f and y = %f\n", om, c2, c3, c4, c5, cG, c0, zcurrtarg, y[0], y[1], y[2]);
	return 8;
      }
      hubble[i] = y[0];
      x[i] = y[1];
    }    
  }

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return 0;

}


/*!
  Calculate analytically the h(z) and x(z) tracker solution
  
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_{rad}\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[in] coord, whether to use a or z coordinate
  \param[out] table of a or z coordinates where h is evaluated
  \param[out] table of h 
  \param[out] table of x
*/
void calcHubbleTracker(){
 
  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    return ;
  }

  bool useacoord = (coord == BARREIRA || coord == JNEVEUa);
  double op = (c2/2.0 - 6.0*c3 + 22.5*c4 - 21.0*c5 - 9*cG )/3.0; //Omega_phi

  if(useacoord){
    int nstep = min(hubble.size(),x.size());
    for(int i = 0; i < nstep; ++i){
      hubble[i] = sqrt(0.5*(om/pow(intvar[i],3)+orad/pow(intvar[i],4)-3*cG)+sqrt(op/9+3*cG+0.25*pow(3*cG-om/pow(intvar[i],3)-orad/pow(intvar[i],4),2))); // Analytical solution for H
      x[i] = intvar[i]/pow(hubble[i],2); // since for tracker, h^2*x = 1
    }
  } else{
    int nstep = min(hubble.size(),x.size());
    for(int i = 0; i < nstep; ++i){
      hubble[i] = sqrt(0.5*(om*pow(1+intvar[i],3)+orad*pow(1+intvar[i],4)-3*cG)+sqrt(op+3*cG+0.25*pow(3*cG-om*pow(1+intvar[i],3)-orad*pow(1+intvar[i],4),2)));
      x[i] = -1/((1+intvar[i])*pow(hubble[i],2));
    }    
  }

}


/*
  Evaluates the age of the universe in galileon cosmology with trapezoidal method.
  Be careful, in here, we only work with a coordinate !
  
  \param[in] cond , initial conditions array (H, xi, pi, age)
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_r\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[out] table of a coordinates where h(a) is evaluated
  \param[out] table of h(a)
  \param[out] table of x(a)
  \param[out] age of universe for given cosmology
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghosr condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : scalar noghost condition not satified
		       \status = 7 : scalar imaginary speed of sound
		       \status = 8 : failed to integrate
*/
int ageOfUniverse(double &age, double cond[4]){

  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    return 9;
  }

  const double Gyr = 3.1556926e16; //Conversion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Conversion factor from Mpc to km
  const double h0bis = h0*lightspeed/1000; // hubble constant in Mpc-1

  double params[1];
  bool useacoord = (coord == BARREIRA || coord == JNEVEUa);
  params[0] = (int)useacoord;

  // Vector y = ( dh/dz, dx/dz, dy/dz )
  double y[3] = {cond[0], cond[1], cond[2]};  //Inital value of integral
  hubble[0] = cond[0];
  x[0] = cond[1];

  // ODE System
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-16, 1e-16, 1, 1);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (3);

  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;

  double a = intvar[1];
  int nstep = intvar.size();
  double h = 1e-6; //Initial step guess
  double acurrtarg;
  double prec_age = 1/(intvar[1]*cond[0]); // initial value of 1/(a*Hbar)
  int st;
  age = cond[3];

  for (int i = 2; i < nstep; i++) {
    acurrtarg = intvar[i];
    while(a < acurrtarg){
      st = gsl_odeiv_evolve_apply(e, c, s, &sys, &a, acurrtarg, &h, y);
      double OmegaP = (0.5*c2*pow(y[0], 2)*pow(intvar[i]*y[1], 2) - 6*c3*pow(y[0], 4)*pow(intvar[i]*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(intvar[i]*y[1], 4) - 21*c5*pow(y[0], 8)*pow(intvar[i]*y[1], 5) - 9*cG*pow(y[0], 4)*pow(intvar[i]*y[1], 2))/(3.0*pow(y[0], 2));
      double OmegaM = 1 - OmegaP - orad/(pow(a,4)*pow(y[0], 2));
      double OmTest = om/(pow(a,3)*pow(y[0], 2));
      // printf("OmegaM : %.16f \t OmTest : %.16f \t |OmegaM - OmTest|/OmTest : %.16f\n", OmegaM, OmTest, fabs(OmegaM - OmTest)/OmTest);
      // printf("OmegaM : %.16f \t OmTest : %.16f \t |OmegaM - OmTest| : %.16f\n", OmegaM, OmTest, fabs(OmegaM - OmTest));
      if(OmegaP<0) st = 5;
      if ( fabs(OmegaM - OmTest)>1e-4 ) {
	printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, intvar[i], pow(y[0], 2), fabs(OmegaM - OmTest));
	st = 4;
      }
      if(st != 0){
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	return st;
      }
    }
    if(isnan(fabs(y[0])) || isnan(fabs(y[1])) || isnan(fabs(y[2]))){
      gsl_odeiv_evolve_free (e);
      gsl_odeiv_control_free (c);
      gsl_odeiv_step_free (s);
      return 8;
    }
    age += 0.5*(1.0/(intvar[i]*y[0])+prec_age)*(intvar[i]-intvar[i-1]); // trapezoidal method
    prec_age = 1/(intvar[i]*y[0]);
    hubble[i-1] = y[0];
    x[i-1] = y[1];
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  // Convert the age of the Universe in Gyr
  age *= Mpc/Gyr/h0bis;

  return 0;

}

/*
  Determine the initial condition in x from the initial conditions in rho and H at z=10^6
  cf Barreira et al. 2013 : Linear perturbations in Galileon gravity models, arXiv:1208.0600
  
  \param[in] H \f$H\f$, {value of H one step before the initial step, initial value of H}
  \param[in] ratio_rho \f$ratio_\rho\f$, initial value of rho_phi/rho_m
  \param[in] age \f$Age of Universe\f$, age of the universe in Gyr
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_r\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[out] initial value of x
  \param[out] table of a coordinates where h(a) is evaluated
  \param[out] table of h(a)
  \param[out] table of x(a)
*/
struct AgeofU{
  int i; // Position in the vector
  double a; // Age of Universe
  AgeofU(int pos, double u_age) : i(pos), a(u_age) {} // Initialization of AgeofU
};

struct less_than_key{
  // To sort by the difference to the input age of universe, while keeping the initial position
  inline bool operator() (const AgeofU& struct1, const AgeofU& struct2){
    return (fabs(struct1.a) < fabs(struct2.a));
  }
};

int initCond(double &xi, double H[2], double ratio_rho, double age){

  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    return 9;
  }

  const double Gyr = 3.1556926e16; //Convertion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Convertion factor from Mpc to km
  int status = 0;
  double coeff[6]; //array of coefficients of the polynomial
  coeff[0] = -3*ratio_rho*om/pow(intvar[1], 3); coeff[1] = 0; coeff[2] = c2/2*pow(H[1], 2)-9*cG*pow(H[1], 4);
  coeff[3] = -6*c3*pow(H[1], 4); coeff[4] = 45/2*c4*pow(H[1], 6); coeff[5] = -21*c5*pow(H[1], 8);

  size_t n = sizeof(coeff)/sizeof(double);
  double cond[4]; // initial conditions
  cond[0] = H[1]; cond[2] = 0; cond[3] = Mpc/Gyr*(1/H[1]-pow(intvar[1], 2)/2*(1/(intvar[0]*H[0])-1/(intvar[1]*H[1]))/(intvar[0]-intvar[1])); // Initial value of the age of the Universe, extrapolating linearly the integral

  AgeofU ini_vec_age = AgeofU(0, 1000);
  std::vector<AgeofU> aou(5, ini_vec_age);

  // Polynomial
  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n);
  double z[2*(n-1)]; // Array of complex numbers ordered like {Re[0],Im[0],Re[1],Im[1],...}

  status = gsl_poly_complex_solve(coeff, n, w, z); // Solve coeff[0]+coeff[1]*x+...+coeff[5]*x^5=0 and store sol in z

  int nrealsol = 0; // Number of real solutions

  for(int j = 0; j<n-1; j++){
    if(z[2*j+1]==0){
      cond[1] = z[2*j]/intvar[1]; // Initial condition for d_phi/da
      aou[j].i = j;
      int status2 = ageOfUniverse(aou[j].a, cond);
      std::cout << aou[j].a << endl;
      aou[j].a -= age;
      if(status2 != 0) return status2;
    }
  }
  
  for(int j = 0; j<n-1; j++){
    if(aou[j].a != 1000){
      nrealsol++;
    }
  }

  sort(aou.begin(), aou.end(), less_than_key()); // Sort by comparing the calculated age of universe to the one provided

  for(int j = 0; j<n-2; j++){
    aou[j].a += age;
  }

  xi = z[2*(aou[0].i)];
  if(nrealsol > 1){
    cond[1] = xi/intvar[1];
    ageOfUniverse(aou[0].a, cond);
  }
  printf("Initial condition for x is %e, it gives an age of the universe equal to %f\n", xi, aou[0].a);

  return status;

} 

void SetAcoord(double amin){

  intvar.clear();

  double da = 1e-4; //!< Step in a
  double dafactor = 0.01; //!< Factor to apply to da if a<astep
  double dafactor2 = 0.001; //!< Factor to apply to da if a<astep2
  double dafactor3 = 0.0001; //!< Factor to apply to da if a<astep3
  double dafactor4 = 0.00001; //!< Factor to apply to da if a<astep4
  double dafactor5 = 0.000001; //!< Factor to apply to da if a<astep5
  double astep;
  double astep2;
  double astep3;
  double astep4;
  double astep5;

  astep = 0.02; //!< Step in a where to change da
  astep2 = 0.0005; //!< Step in a where to change da
  astep3 = 3e-5; //!< Step in a where to change da
  astep4 = 5e-6; //!< Step in a where to change da
  astep5 = 5e-7; //!< Step in a where to change da  

  int nstep = (int)floor((1.0-astep)/da);
  int nstep1 = nstep + (int)floor((astep-astep2)/(dafactor*da));
  int nstep2 = nstep1 + (int)floor((astep2-astep3)/(dafactor2*da));
  int nstep3 = nstep2 + (int)floor((astep3-astep4)/(dafactor3*da));
  int nstep4 = nstep3 + (int)floor((astep4-astep5)/(dafactor4*da));
  int nsteptot = nstep4 + (int)floor(astep5/(dafactor5*da));

  printf("nstep : %i, nstep1 : %i, nstep2 : %i, nstep3 : %i, nstep4 : %i, nsteptot : %i\n", nstep, nstep1, nstep2, nstep3, nstep4, nsteptot);

  intvar.push_back(1.0); // Initialization of acoord

  for(int i = 1; i< nsteptot; ++i){
     if(intvar[i-1]>amin){
      if(i<nstep){
	intvar.push_back(1.0-i*da);
      } else if(nstep<=i && i<nstep1){
	intvar.push_back(astep-(i-nstep)*da*dafactor);
      } else if(nstep1<=i && i<nstep2){
	intvar.push_back(astep2-(i-nstep1)*da*dafactor2);
      } else if(nstep2<=i && i<nstep3){
	intvar.push_back(astep3-(i-nstep2)*da*dafactor3);
      } else if(nstep3<=i && i<nstep4){
	intvar.push_back(astep4-(i-nstep3)*da*dafactor4);
      } else if(nstep4<=i && i<nsteptot){
	intvar.push_back(astep5-(i-nstep4)*da*dafactor5);
      }
     } else{
      break;
    }
  }

  intvar.pop_back();
  intvar.push_back(amin);
  //for (unsigned int k=0; k<intvar.size(); ++k)
  // std::cout<<intvar[k]<<std::endl;
}


// Function that calculates dtau/da = 1/(a^2*H) (tau being confromal time)
// The function is the one linked to CAMB in order to calculate the background contribution
// extern "C" void arrays_(char* infile, char* outfile, double* omegar){
extern "C" void arrays_(char* infile, double* omegar){

  fflush(stdout);

  // Clear three global vectors
  intvar.clear();
  hubble.clear();
  x.clear();

  char** params = readParamsFile(infile);

  c2 = (*params != "") ? atof(*params) : 0;
  c3 = (*(params+1) != "") ? atof(*(params+1)) : 0;
  c4 = (*(params+2) != "") ? atof(*(params+2)) : 0;
  c5 = (*(params+3) != "") ? atof(*(params+3)) : 0;
  c0 = (*(params+4) != "") ? atof(*(params+4)) : 0;
  cG = (*(params+5) != "") ? atof(*(params+5)) : 0;
  double ratio_rho = (*(params+6) != "") ? atof(*(params+6)) : 0;
  double age = (*(params+7) != "") ? atof(*(params+7)) : 0;
  om = (*(params+8) != "" && *(params+9) != "" && *(params+12) != "") ? (atof(*(params+8))+atof(*(params+9)))/pow(atof(*(params+12)), 2)*10000 : 0;
  orad = (*omegar);
  h0 = (*(params+12) != "") ? atof(*(params+12))*1000/lightspeed : 0;
  char* solvingMethod = *(params+10);
  if(solvingMethod != "") solvingMethod[strlen(solvingMethod)-1] = '\0';
  if(strcmp(solvingMethod, "JNEVEUz") != 0 && strcmp(solvingMethod, "JNEVEUa") != 0 && strcmp(solvingMethod, "BARREIRA") != 0){
    fprintf(stderr, "WARNING : invalid integration method, will use the default one\n");
  } else if(strcmp(solvingMethod, "JNEVEUz") == 0){
    coord = JNEVEUz;
  } else if(strcmp(solvingMethod, "JNEVEUa") == 0){
    coord = JNEVEUa;
  } else if(strcmp(solvingMethod, "BARREIRA") == 0){
    coord = BARREIRA;
  }
  double dotphi = (*(params+11) != "") ? atof(*(params+11)) : 0; // careful, dotphi is the initial condition on x !!!
  char* outfile = *(params+13); // take same output file name root as camb
  if(outfile != ""){
    outfile[strlen(outfile)-1] = '\0';
    outfile = strcat(outfile, "_");
  }
  outfile = strcat(outfile, "background.dat");
  

  printf("Input file : %s\nOutput file : %s\n", infile, outfile);
  printf("OmegaM0 = %f\nOmegaR0 = %.18f\nc0 = %f\nc2 = %f\nc3 = %f\nc4 = %f\nc5 = %f\ncG = %f\nh0 = %f Mpc-1\n", om, orad, c0, c2, c3, c4, c5, cG, h0);

  // The status of the integration, 0 if everything ok
  int status = 0;

  if(coord == BARREIRA){
    printf("mode : barreira\n");

    // Where to put initial condition
    double apremin = 1e-7;
    double amin = 9.99999e-7;
    double amax = 1.;
    intvar.push_back(apremin);

    // Fill the vector of a with a geometric sequence
    double q = 1.+5./10000;
    for(int i = 0; i<log(amax/amin)/log(q); i++){
      intvar.push_back(amin*pow(q, i));
    }

    printf("Number of points : %i\n", intvar.size());

    // Calculate initial conditions in H and ratio_rho
    double xi = 0;
    double H[2];
    if(dotphi == 0){ // if dotphi is not given (i.e. initial condition is on ratio_rho)
      H[0] = sqrt(om/pow(apremin, 3)*(1+ratio_rho)+orad/pow(apremin, 4));
      H[1] = sqrt(om/pow(amin, 3)*(1+ratio_rho)+orad/pow(amin, 4));
    } else{
      H[0] = sqrt(om/pow(apremin, 3)+orad/pow(apremin, 4));
      H[1] = sqrt(om/pow(amin, 3)+orad/pow(amin, 4));
      // H[1] = sqrt(om/pow(amin, 3)*(1+ratio_rho)+orad/pow(amin, 4));
      ratio_rho = pow(intvar[1], 3)/(3*om)*(0.5*c2*pow(dotphi*H[1], 2)-6*c3*H[1]*pow(dotphi*H[1], 3)+22.5*c4*pow(H[1], 2)*pow(dotphi*H[1], 4)-21*c5*pow(H[1], 3)*pow(dotphi*H[1], 5)-9*cG*pow(H[1], 2)*pow(dotphi*H[1], 2));
    }

    hubble.resize(intvar.size()-1, 999999);
    x.resize(intvar.size()-1, 999999);

    // Calculate initial conditions in x and fill hubble and x from then to now
    status = initCond(xi, H, ratio_rho, age);

    printf("\nstatus = %i\n", status);

  } else if(coord == JNEVEUa){
    printf("mode : acoord\n");

    // Fill the array of a with points where to integrate
    double amin = 0.000908;
    SetAcoord(amin);
    sort(intvar.begin(), intvar.end(), std::greater<double>());

    hubble.resize(intvar.size(), 999999);
    x.resize(intvar.size(), 999999);

    // Integrate and fill hubble and x both when tracker and not tracker
    if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
    {
      status = calcHubbleGalileon();
    }
    else {
      calcHubbleTracker();
    }

    printf("\nstatus = %i\n", status);

  } else if(coord == JNEVEUz){
    printf("mode : zcoord\n");

    // Fill the array of z with points where to integrate
    double zmax = 1100;
    double amin = 1/(1+zmax);
    SetAcoord(amin);
    sort(intvar.begin(), intvar.end(), std::greater<double>());
    for(int i=0;i<intvar.size();i++) {intvar[i]=1/(1+intvar[i]);}

    hubble.resize(intvar.size(), 999999);
    x.resize(intvar.size(), 999999);

    double op = 1.5*c2 - 18*c3 + 67.5*c4 - 63*c5 - 27*cG;

    // Integrate and fill hubble and x both when tracker and not tracker
    if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
    {
      status = calcHubbleGalileon();
    }
    else {
      printf("Tracker\n");
      calcHubbleTracker();
    }

    printf("\nstatus = %i\n", status);
  }


  if(status != 0){
    hubble.clear();
    x.clear();
  }

  FILE* f = fopen(outfile, "w");
  for(int i = 0; i<intvar.size()-1; i++){
    double alpha = c2/6*hubble[i]*intvar[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*intvar[i+1]*x[i];
    double gamma = c2/3*pow(hubble[i], 2)*intvar[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];
    double beta = c2/6*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*intvar[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
    double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2);
    double lambda = 3*pow(hubble[i], 2) + orad/pow(intvar[i+1], 4) + c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
    double omega = 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];
	
    double x_prime = -intvar[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
    double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
	
    double rho = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
    double p = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*(h_prime*intvar[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*intvar[i+1]*x[i]+2*hubble[i]*intvar[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*intvar[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2));
	
    double weff = p/rho;
    double hubble_LCDM = sqrt(om/pow(intvar[i+1], 3)+orad/pow(intvar[i+1], 4)+(1-om-orad));
    double ratio = hubble[i]/hubble_LCDM;
    fprintf(f, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", intvar[i+1], hubble[i], x[i], hubble_LCDM, ratio, weff);
  }

  double hubble_interp[hubble.size()];
  std::copy(hubble.begin(), hubble.end(), hubble_interp);
  double x_interp[x.size()];
  std::copy(x.begin(), x.end(), x_interp);
  double intvar_interp[intvar.size()-1];
  std::copy(intvar.begin()+1, intvar.end(), intvar_interp);

  spline_h = gsl_spline_alloc(gsl_interp_cspline, hubble.size());
  spline_x = gsl_spline_alloc(gsl_interp_cspline, x.size());  
  gsl_spline_init(spline_h, intvar_interp, hubble_interp, hubble.size());
  gsl_spline_init(spline_x, intvar_interp, x_interp, x.size());

}


extern "C" double* handxofa_(double* point){
  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    exit(EXIT_FAILURE);
  }

  if((*point) < 0.0 || (*point) > 1.1){
    printf("Forbidden value of a : %.12f\n", (*point));
    exit(EXIT_FAILURE);
  }

  static double hx[2];

  // Fill the vector of a with a geometric sequence
  double q = 1.+5./10000;

  double alpha = (log(*point) - log(intvar[1]))/log(q); // Solving amin*q^alpha = a
  int i = floor(alpha)+1; // i is integer part of alpha, so that a[i] <= a < a[i+1]

  //printf("%i, %i, %f, %f, %i, %f, %f, %.12f\n", i, hubble.size(), hubble[i], hubble[i-1], intvar.size(), intvar[i+1], intvar[i], (*point));

  // // Linear interpolation
  // hx[0] = (hubble[i] - hubble[i-1])/(intvar[i+1]-intvar[i])*((*point) - intvar[i]) + hubble[i-1];
  // hx[1] = (x[i] - x[i-1])/(intvar[i+1]-intvar[i])*((*point) - intvar[i]) + x[i-1];

  // Spline interpolation
  hx[0] = gsl_spline_eval(spline_h, *point, acc);
  hx[1] = gsl_spline_eval(spline_x, *point, acc);

  // // Resolution of coupled differential equations from a[i] to *a
  // double params[1];
  // bool useacoord = (coord == BARREIRA || coord == JNEVEUa);
  // params[0] = (int)useacoord;
  // // Vector y = ( dh/dz, dx/dz, dy/dz )
  // double y[3] = {hubble[i-1], x[i-1], 0};  //Inital value of integral

  // const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  // gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 3);
  // gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-16, 1e-16, 1, 1);
  // gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(3);
  // gsl_odeiv_system sys;
  // sys.function = calcValOmC2C3C4C5CGC0;
  // sys.dimension = 3;
  // sys.params = &params;

  // double var = intvar[i]; // a or z
  // double var2;
  // double h = 1e-6; //Initial step guess
  // double currtarg = (*point);
  // double currtarg2;
  // int st;
  // double OmegaP;
  // double OmegaM;
  // double OmTest;

  // if(coord == BARREIRA || coord == JNEVEUz){
  //   var2 = var;
  //   currtarg2 = currtarg;
  // } else if(coord == JNEVEUa){
  //   var2 = -var;
  //   currtarg2 = -currtarg;
  // }

  // while(var2 < currtarg2){
  //   st = gsl_odeiv_evolve_apply(e, c, s, &sys, &var, currtarg, &h, y);
  //   if(coord == BARREIRA || coord == JNEVEUz){
  //     var2 = var;
  //   } else if(coord == JNEVEUa){
  //     var2 = -var;
  //   }
  // }
  // if(coord == BARREIRA || coord == JNEVEUa){
  //   OmegaP = (0.5*c2*pow(y[0], 2)*pow(currtarg*y[1], 2) - 6*c3*pow(y[0], 4)*pow(currtarg*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(currtarg*y[1], 4) - 21*c5*pow(y[0], 8)*pow(currtarg*y[1], 5) - 9*cG*pow(y[0], 4)*pow(currtarg*y[1], 2))/(3.0*pow(y[0], 2)) ;
  //   OmegaM = 1 - OmegaP - orad/(pow(var,4)*pow(y[0], 2));
  //   OmTest = om/(pow(var,3)*pow(y[0], 2));
  // } else if(coord == JNEVEUz){
  //   OmegaP = (0.5*c2*pow(y[0], 2)*pow(-(1+currtarg)*y[1], 2) - 6*c3*pow(y[0], 4)*pow(-(1+currtarg)*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(-(1+currtarg)*y[1], 4) - 21*c5*pow(y[0], 8)*pow(-(1+currtarg)*y[1], 5) - 9*cG*pow(y[0], 4)*pow(-(1+currtarg)*y[1], 2))/(3.0*pow(y[0], 2)) ;
  //   OmegaM = 1 - OmegaP - orad*pow(1+currtarg, 4)/pow(y[0], 2);
  //   OmTest = om*pow(1+var, 3)/pow(y[0], 2);
  // }

  // if(OmegaP<0) st = 5;
  // if ( fabs(OmegaM - OmTest)>1e-4  ) {
  //   printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, intvar[i], pow(y[0], 2), fabs(OmegaM - OmTest));
  //   st = 4;
  // }
  
  // if(st != 0){
  //   gsl_odeiv_evolve_free (e);
  //   gsl_odeiv_control_free (c);
  //   gsl_odeiv_step_free (s);
  //   printf("status = %d", st);
  //   exit(EXIT_FAILURE);
  // }
    
  // if(isnan(fabs(y[0])) || isnan(fabs(y[1])) || isnan(fabs(y[2]))){
  //   gsl_odeiv_evolve_free (e);
  //   gsl_odeiv_control_free (c);
  //   gsl_odeiv_step_free (s);
  //   //std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
  //   printf("status = %d", 8);
  //   exit(EXIT_FAILURE);
  // }
    
  // hx[0] = y[0];
  // hx[1] = y[1];

  // gsl_odeiv_evolve_free (e);
  // gsl_odeiv_control_free (c);
  // gsl_odeiv_step_free (s);

  if(i<=3){
    double hubble_LCDM = sqrt(om/pow((*point), 3)+orad/pow((*point), 4)+(1-om-orad));
    if(fabs((hx[0]-hubble_LCDM)/hx[0])>1e-3) fprintf(stderr, "Warning : no continuity between LCDM and galileon background at very early time ( i = %i, a = %f, h_LCDM = %f and h_gal = %f)", i, (*point), hubble_LCDM, hx[0]);
  }

  return hx;

}


extern "C" double grhogal_(double* point){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double grhogal = 0;

  if(*point >= 9.99999e-7) grhogal = pow(h0, 2)*pow(*point, 2)*(c2/2*pow(h, 2)*pow(xgal, 2) - 6*c3*pow(h, 4)*pow(xgal, 3) + 22.5*c4*pow(h, 6)*pow(xgal, 4) - 21*c5*pow(h, 8)*pow(xgal, 5) - 9*cG*pow(h, 4)*pow(xgal, 2));

  //printf("grhogal : %.12f \t %.12f \t %.12f \t %.12f\n", *point, h, xgal, grhogal);

  return grhogal;

}


extern "C" double gpresgal_(double* point){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c2/6*h*xgal-3*c3*pow(h, 3)*pow(xgal, 2) + 15*c4*pow(h, 5)*pow(xgal, 3) - 17.5*c5*pow(h, 7)*pow(xgal, 4) - 3*cG*pow(h, 3)*xgal;
  double gamma = c2/3*pow(h, 2)*xgal-c3*pow(h, 4)*pow(xgal, 2) + 2.5*c5*pow(h, 8)*pow(xgal, 4) - 2*cG*pow(h, 4)*xgal;
  double beta = c2/6*pow(h, 2) -2*c3*pow(h, 4)*xgal + 9*c4*pow(h, 6)*pow(xgal, 2) - 10*c5*pow(h, 8)*pow(xgal, 3) - cG*pow(h, 4);
  double sigma = 2*h + 2*c3*pow(h, 3)*pow(xgal, 3) - 15*c4*pow(h, 5)*pow(xgal, 4) + 21*c5*pow(h, 7)*pow(xgal, 5) + 6*cG*pow(h, 3)*pow(xgal, 2);
  double lambda = 3*pow(h, 2) + orad/pow(*point, 4) + c2/2*pow(h, 2)*pow(xgal, 2) - 2*c3*pow(h, 4)*pow(xgal, 3) + 7.5*c4*pow(h, 6)*pow(xgal, 4) - 9*c5*pow(h, 8)*pow(xgal, 5) - cG*pow(h, 4)*pow(xgal, 2);
  double omega = 2*c3*pow(h, 4)*pow(xgal, 2) - 12*c4*pow(h, 6)*pow(xgal, 3) + 15*c5*pow(h, 8)*pow(xgal, 4) + 4*cG*pow(h, 4)*xgal;

  double xprime = -xgal+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega); // derivative wrt ln(a)
  double hprime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (x*H)'

  double gpresgal = 0;

  if(*point >= 9.99999e-7) gpresgal = pow(h0, 2)*pow(*point, 2)*(c2/2*pow(h, 2)*pow(xgal, 2) + 2*c3*pow(h, 3)*pow(xgal, 2)*(hprime*xgal+xprime*h) - c4*(4.5*pow(h, 6)*pow(xgal, 4) + 12*pow(h, 6)*pow(xgal, 3)*xprime + 15*pow(h, 5)*pow(xgal, 4)*hprime) + 3*c5*pow(h, 7)*pow(xgal, 4)*(5*h*xprime+7*hprime*xgal+2*h*xgal) + cG*(6*pow(h, 3)*pow(xgal, 2)*hprime + 4*pow(h, 4)*xgal*xprime + 3*pow(h, 4)*pow(xgal, 2)));

  return gpresgal;

}


// Function that calculates the density perturbation of galileon
extern "C" double Chigal_(double* dgrho, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha_Z = -2*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 2) 
    + 15*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 21*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - 6*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double alpha_eta = 1.5*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 3*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta = 1./(1-0.5*alpha_Z);
  double ChitildeG = c2*h0*xgal*h*(*dphiprime) 
    - c3/pow(*point, 2)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphiprime) + 2*pow(*k, 2)*pow(xgal, 2)*pow(h, 2)*(*dphi)) 
    + c4/pow(*point, 4)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphiprime) + 12*pow(*k, 2)*pow(xgal, 3)*pow(h, 4)*(*dphi)) 
    - c5/pow(*point, 6)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphiprime) + 15*pow(*k, 2)*pow(xgal, 4)*pow(h, 6)*(*dphi)) 
    - cG/pow(*point, 2)*(18*h0*xgal*pow(h, 3)*(*dphiprime) + 4*pow(*k, 2)*xgal*pow(h, 2)*(*dphi));

  double ChiG = 0;

  if(-1e-5 < alpha_Z - 2 && alpha_Z - 2 < 1e-5) printf("WARNING : 1/beta_chi is zero");
  // if (*point >= 9.99999e-7) ChiG = beta*(ChitildeG + 0.5*alpha_Z*(*dgrho) + (alpha_Z - 2*alpha_eta)*(*k)*(*eta));
  if (*point >= 1e-4) ChiG = beta*(ChitildeG + 0.5*alpha_Z*(*dgrho) + (alpha_Z - 2*alpha_eta)*(*k)*(*eta));
  //ChiG = beta*(ChitildeG + 0.5*alpha_Z*(*dgrho) + (alpha_Z - 2*alpha_eta)*(*k)*(*eta));


  //if(*point >= 9.99999e-7) printf("%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*dgrho), (*eta), (*dphi), (*dphiprime), ChiG);

  // FILE* g = fopen("chigal/chigal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", (*point), ChiG, (*dgrho), beta*ChitildeG, 0.5*beta*alpha_Z*(*dgrho), beta*(alpha_Z - 2*alpha_eta)*(*k)*(*eta));
  // fclose(g);

  return ChiG;

}


// Perturbation of heat flux from galileon
extern "C" double qgal_(double* dgq, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 2*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - cG/(1.5*pow(*point, 2))*pow(xgal, 2)*pow(h, 2);
  double beta = 1./(1-1.5*alpha);
  double qtildeG = c2*(*k)*h0*xgal*h*(*dphi) 
    - c3/pow(*point, 2)*(*k)*(6*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) - 2*pow(xgal, 2)*pow(h, 2)*(*dphiprime)) 
    + c4/pow(*point, 4)*(*k)*(-12*pow(xgal, 3)*pow(h, 4)*(*dphiprime) + 18*h0*pow(xgal, 3)*pow(h, 5)*(*dphi)) 
    - c5/pow(*point, 6)*(*k)*(-15*pow(xgal, 4)*pow(h, 6)*(*dphiprime) + 15*h0*pow(xgal, 4)*pow(h, 7)*(*dphi)) 
    - cG/pow(*point, 2)*(*k)*(-4*xgal*pow(h, 2)*(*dphiprime) + 6*h0*xgal*pow(h, 3)*(*dphi));

  double qG = 0;

  if(-1e-5 < 1.5*alpha - 1 && 1.5*alpha - 1 < 1e-5) printf("WARNING : 1/beta_q is zero");
  // if(*point >= 9.99999e-7) qG = beta*(qtildeG + 1.5*alpha*(*dgq));
  if(*point >= 1e-4) qG = beta*(qtildeG + 1.5*alpha*(*dgq));
  //qG = beta*(qtildeG + 1.5*alpha*(*dgq));

  //if(*point >= 1) printf("%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*eta), h, xgal, qG);

  // FILE* g = fopen("qgal/qgal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", (*point), qG, (*dgq), beta*qtildeG, 1.5*beta*alpha*(*dgq));
  // fclose(g);

  return qG;

}


// Perturbation of anisotropic stress from galileon
extern "C" double Pigal_(double* dgrho, double* dgq, double* dgpi, double* eta, double* dphi, double* point, double* k){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)/pow(*point, 2)); // H = adot/a
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c2/6*hoft*xgal-3*c3*pow(hoft, 3)*pow(xgal, 2) + 15*c4*pow(hoft, 5)*pow(xgal, 3) - 17.5*c5*pow(hoft, 7)*pow(xgal, 4) - 3*cG*pow(hoft, 3)*xgal;
  double gamma = c2/3*pow(hoft, 2)*xgal-c3*pow(hoft, 4)*pow(xgal, 2) + 2.5*c5*pow(hoft, 8)*pow(xgal, 4) - 2*cG*pow(hoft, 4)*xgal;
  double beta = c2/6*pow(hoft, 2) -2*c3*pow(hoft, 4)*xgal + 9*c4*pow(hoft, 6)*pow(xgal, 2) - 10*c5*pow(hoft, 8)*pow(xgal, 3) - cG*pow(hoft, 4);
  double delta = 2*hoft + 2*c3*pow(hoft, 3)*pow(xgal, 3) - 15*c4*pow(hoft, 5)*pow(xgal, 4) + 21*c5*pow(hoft, 7)*pow(xgal, 5) + 6*cG*pow(hoft, 3)*pow(xgal, 2);
  double lambda = 3*pow(hoft, 2) + orad/pow(*point, 4) + c2/2*pow(hoft, 2)*pow(xgal, 2) - 2*c3*pow(hoft, 4)*pow(xgal, 3) + 7.5*c4*pow(hoft, 6)*pow(xgal, 4) - 9*c5*pow(hoft, 8)*pow(xgal, 5) - cG*pow(hoft, 4)*pow(xgal, 2);
  double omega = 2*c3*pow(hoft, 4)*pow(xgal, 2) - 12*c4*pow(hoft, 6)*pow(xgal, 3) + 15*c5*pow(hoft, 8)*pow(xgal, 4) + 4*cG*pow(hoft, 4)*xgal;

  double xprime = -xgal+(alpha*lambda-delta*gamma)/(delta*beta-alpha*omega); // derivative wrt ln(a)
  double hprime = (*point)*(omega*gamma-lambda*beta)/(delta*beta-alpha*omega) + h; // Careful, this is the derivative of h and not hoft
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'


  double alpha_sigprime = c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 3*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 5)*xhprime;
  double alpha_sig = c4/pow(*point, 4)*(3*pow(xgal, 4)*pow(h, 4) - 6*pow(xgal, 3)*pow(h, 3)*xhprime)
    - c5/pow(*point, 6)*(-3*pow(xgal, 5)*pow(h, 5)*hprime + 12*pow(xgal, 5)*pow(h, 6) - 15*pow(xgal, 4)*pow(h, 5)*xhprime)
    + 2*cG/pow(*point, 2)*xgal*h*xhprime;
  double alpha_phi = -c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - c5/pow(*point, 6)*(6*pow(xgal, 4)*pow(h, 5)*xhprime - 6*pow(xgal, 5)*pow(h, 6))
    + 2*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta_pi = 1./(1+0.5*alpha_phi - alpha_sigprime);
  double PitildeG = pow(*k, 2)*(c4/pow(*point, 4)*(4*pow(xgal, 3)*pow(h, 4)*(*dphi) - 6*pow(xgal, 2)*pow(h, 3)*xhprime*(*dphi))
    - c5/pow(*point, 6)*(-12*pow(xgal, 3)*pow(h, 5)*xhprime*(*dphi) + 12*pow(xgal, 4)*pow(h, 6)*(*dphi) - 3*pow(xgal, 4)*pow(h, 5)*hprime*(*dphi))
    + 2*cG/pow(*point, 2)*h*xhprime*(*dphi));

  double PiG = 0;

  if(-1e-5 < (alpha_phi - 2*alpha_sigprime+2) && (alpha_phi - 2*alpha_sigprime+2) < 1e-5) printf("WARNING : 1/beta_pi is zero");
  // if(*point >= 9.99999e-7) PiG = beta_pi*(PitildeG + (alpha_sigprime - 0.5*alpha_phi)*(*dgpi) + 0.5*(2*alpha_sigprime + alpha_sig - alpha_phi)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sig + alpha_sigprime)*(*k)*(*eta));
  if(*point >= 1e-4) PiG = beta_pi*(PitildeG + (alpha_sigprime - 0.5*alpha_phi)*(*dgpi) + 0.5*(2*alpha_sigprime + alpha_sig - alpha_phi)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sig + alpha_sigprime)*(*k)*(*eta));
  //PiG = beta_pi*(PitildeG - (0.5*alpha_phi - alpha_sigprime)*(*dgpi) + (0.5*alpha_sig - 0.5*alpha_phi + alpha_sigprime)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sig + alpha_sigprime)*(*k)*(*eta));

  //if(*point >= 9.99999e-7) printf("Pigal : %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), beta_pi, PitildeG, alpha_phi, alpha_sigprime, alpha_sig, (*dgpi), (*dgrho), (*dgq), (*eta), PiG);

  // FILE* g = fopen("pigal/pigal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", (*point), PiG, (*dgpi), beta_pi*PitildeG, beta_pi*(0.5*alpha_phi - alpha_sigprime)*(*dgpi), beta_pi*(0.5*alpha_sig - 0.5*alpha_phi + alpha_sigprime)*(*dgrho), 3*beta_pi*(0.5*alpha_sig - 0.5*alpha_phi + alpha_sigprime)*h0*h/(*k)*(*dgq), beta_pi*(alpha_sig + alpha_sigprime)*(*k)*(*eta));
  // fclose(g);

  return PiG;

}

// Perturbation of the galileon field
extern "C" double dphisecond_(double* dgrho, double* dgq, double* z, double* eta, double* dphi, double* dphiprime, double* point, double* k, double* grho, double* gpres, double* grhob, double* clxb, double* clxbdot, double* grhoc, double* clxc, double* clxcdot, double* grhor, double* clxr, double* clxrdot, double* grhog, double* clxg, double* clxgdot){
// extern "C" double dphisecond_(double* grho, double* gpres, double* grhob, double* clxb, double* clxbdot, double* grhoc, double* clxc, double* clxcdot, double* grhor, double* clxr, double* clxrdot, double* grhog, double* clxg, double* clxgdot, double* dgrho, double* dgq, double* z, double* eta, double* dphi, double* dphiprime, double* point, double* k){


  // printf("%.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f\n", (*dgrho), (*eta), (*dphi), (*dphiprime), (*point), (*k));


  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)); // H = adot/a
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c2/6*hoft*xgal-3*c3*pow(hoft, 3)*pow(xgal, 2) + 15*c4*pow(hoft, 5)*pow(xgal, 3) - 17.5*c5*pow(hoft, 7)*pow(xgal, 4) - 3*cG*pow(hoft, 3)*xgal;
  double gamma = c2/3*pow(hoft, 2)*xgal-c3*pow(hoft, 4)*pow(xgal, 2) + 2.5*c5*pow(hoft, 8)*pow(xgal, 4) - 2*cG*pow(hoft, 4)*xgal;
  double beta = c2/6*pow(hoft, 2) -2*c3*pow(hoft, 4)*xgal + 9*c4*pow(hoft, 6)*pow(xgal, 2) - 10*c5*pow(hoft, 8)*pow(xgal, 3) - cG*pow(hoft, 4);
  double delta = 2*hoft + 2*c3*pow(hoft, 3)*pow(xgal, 3) - 15*c4*pow(hoft, 5)*pow(xgal, 4) + 21*c5*pow(hoft, 7)*pow(xgal, 5) + 6*cG*pow(hoft, 3)*pow(xgal, 2);
  double lambda = 3*pow(hoft, 2) + orad/pow(*point, 4) + c2/2*pow(hoft, 2)*pow(xgal, 2) - 2*c3*pow(hoft, 4)*pow(xgal, 3) + 7.5*c4*pow(hoft, 6)*pow(xgal, 4) - 9*c5*pow(hoft, 8)*pow(xgal, 5) - cG*pow(hoft, 4)*pow(xgal, 2);
  double omega = 2*c3*pow(hoft, 4)*pow(xgal, 2) - 12*c4*pow(hoft, 6)*pow(xgal, 3) + 15*c5*pow(hoft, 8)*pow(xgal, 4) + 4*cG*pow(hoft, 4)*xgal;

  double xprime = -xgal+(alpha*lambda-delta*gamma)/(delta*beta-alpha*omega); // derivative wrt ln(a)
  double hprime = (*point)*(omega*gamma-lambda*beta)/(delta*beta-alpha*omega) + h; // Careful, this is the derivative of h and not hoft
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*h)'

  double alpha_gammasecond = c2
    - 12*c3/pow(*point, 2)*xgal*pow(h, 2)
    + 54*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)
    - 60*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)
    - 6*cG/pow(*point, 2)*pow(h, 2);
  double alpha_gammaprime = 2*c2
    - c3/pow(*point, 2)*(12*h*xhprime + 12*xgal*h*hprime)
    + c4/pow(*point, 4)*(-108*pow(xgal, 2)*pow(h, 4) + 108*xgal*pow(h, 3)*xhprime + 108*pow(xgal, 2)*pow(h, 3)*hprime)
    - c5/pow(*point, 6)*(-240*pow(xgal, 3)*pow(h, 6) + 180*pow(xgal, 2)*pow(h, 5)*xhprime + 180*pow(xgal, 3)*pow(h, 5)*hprime)
    - 12*cG/pow(*point, 2)*h*hprime;
  double alpha_gamma = c2
    - c3/pow(*point, 2)*(4*xgal*pow(h, 2) + 4*h*xhprime)
    + c4/pow(*point, 4)*(-10*pow(xgal, 2)*pow(h, 4) + 24*xgal*pow(h, 3)*xhprime + 12*pow(xgal, 2)*pow(h, 3)*hprime)
    - c5/pow(*point, 6)*(-36*pow(xgal, 3)*pow(h, 6) + 36*pow(xgal, 2)*pow(h, 5)*xhprime + 24*pow(xgal, 3)*pow(h, 5)*hprime)
    - cG/pow(*point, 2)*(2*pow(h, 2) + 4*h*hprime);
  double alpha_Z = c2*xgal
    - c3/pow(*point, 2)*(6*pow(xgal, 2)*pow(h, 2) + 4*xgal*h*xhprime)
    + c4/pow(*point, 4)*(-6*pow(xgal, 3)*pow(h, 4) + 36*pow(xgal, 2)*pow(h, 3)*xhprime + 12*pow(xgal, 3)*pow(h, 3)*hprime)
    - c5/pow(*point, 6)*(-45*pow(xgal, 4)*pow(h, 6) + 60*pow(xgal, 3)*pow(h, 5)*xhprime + 30*pow(xgal, 4)*pow(h, 5)*hprime)
    - cG/pow(*point, 2)*(6*xgal*pow(h, 2) + 4*h*xhprime + 4*xgal*h*hprime);
  double alpha_Zprime = -2*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)
    + 12*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)
    - 15*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)
    - 4*cG/pow(*point, 2)*xgal*pow(h, 2);
  double alpha_eta = c4/pow(*point, 4)*(-4*pow(xgal, 3)*pow(h, 4) + 6*pow(xgal, 2)*pow(h, 3)*xhprime)
    - c5/pow(*point, 6)*(-12*pow(xgal, 4)*pow(h, 6) + 3*pow(xgal, 4)*pow(h, 5)*hprime + 12*pow(xgal, 3)*pow(h, 5)*xhprime)
    - 2*cG/pow(*point, 2)*h*xhprime;

  double dphisecond = 0;
  
  // if (alpha_gammasecond == 0) printf("WARNING : 1/beta_gamma is zero");
  // if(*point >= 9.99999e-7) dphisecond = -(alpha_gammaprime*h0*h*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + 0.5*(alpha_Z - 2*alpha_Zprime)*(*dgrho) + (alpha_Z - alpha_Zprime - 2*alpha_eta)*(*k)*(*eta))/alpha_gammasecond;


  // Derivatives of density perturbations
  double dotdeltab = (*grhob)*((*clxbdot) - 3*h0*h*(*clxb));
  double dotdeltac = (*grhoc)*((*clxcdot) - 3*h0*h*(*clxc));
  double dotdeltar = (*grhor)*((*clxrdot) - 4*h0*h*(*clxr));
  double dotdeltag = (*grhog)*((*clxgdot) - 4*h0*h*(*clxg));
  double dotdeltaf = dotdeltab + dotdeltac + dotdeltar + dotdeltag;

  // Expression of Z'
  double chiprimehat = c2*(-2*pow(h0, 2)*xgal*pow(h, 2)*(*dphiprime) + pow(h0, 2)*xgal*h*hprime*(*dphiprime) + pow(h0, 2)*pow(h, 2)*xprime*(*dphiprime))
    - c3/pow(*point, 2)*(-72*pow(h0, 2)*pow(xgal, 2)*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphiprime) + 36*pow(h0, 2)*xgal*pow(h, 4)*xprime*(*dphiprime)
  			 + pow(*k, 2)*(2*pow(xgal, 2)*pow(h, 2)*(*dphiprime) - 8*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) + 4*h0*pow(xgal, 2)*pow(h, 2)*hprime*(*dphi) + 4*h0*xgal*pow(h, 3)*xprime*(*dphi)))
    + c4/pow(*point, 4)*(-540*pow(h0, 2)*pow(xgal, 3)*pow(h, 6)*(*dphiprime) + 450*pow(h0, 2)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphiprime) + 270*pow(h0, 2)*pow(xgal, 2)*pow(h, 6)*xprime*(*dphiprime)
  			 + pow(*k, 2)*(12*pow(xgal, 3)*pow(h, 4)*(*dphiprime) - 72*h0*pow(xgal, 3)*pow(h, 5)*(*dphi) + 48*h0*pow(xgal, 3)*pow(h, 4)*hprime*(*dphi) + 36*h0*pow(xgal, 2)*pow(h, 5)*xprime*(*dphi)))
    - c5/pow(*point, 6)*(-840*pow(h0, 2)*pow(xgal, 4)*pow(h, 8)*(*dphiprime) + 735*pow(h0, 2)*pow(xgal, 4)*pow(h, 7)*hprime*(*dphiprime) + 420*pow(h0, 2)*pow(xgal, 3)*pow(h, 8)*xprime*(*dphiprime)
  			 + pow(*k, 2)*(15*pow(xgal, 4)*pow(h, 6)*(*dphiprime) - 120*h0*pow(xgal, 4)*pow(h, 7)*(*dphi) + 90*h0*pow(xgal, 4)*pow(h, 6)*hprime*(*dphi) + 60*h0*pow(xgal, 3)*pow(h, 7)*xprime*(*dphi)))
    - cG/pow(*point, 2)*(-72*pow(h0, 2)*xgal*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*xgal*pow(h, 3)*hprime*(*dphiprime) + 18*pow(h0, 2)*pow(h, 4)*xprime*(*dphiprime)
  			 + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphiprime) - 16*h0*xgal*pow(h, 3)*(*dphi) + 8*h0*xgal*pow(h, 2)*hprime*(*dphi) + 4*h0*pow(h, 3)*xprime*(*dphi)))
    ;
  double beta_gammasecond = c2*h0*xgal*h - 18*c3*h0/pow(*point, 2)*pow(xgal, 2)*pow(h, 3) + 90*c4*h0/pow(*point, 4)*pow(xgal, 3)*pow(h, 5) - 105*c5*h0/pow(*point, 6)*pow(xgal, 4)*pow(h, 7) - 18*cG*h0/pow(*point, 2)*xgal*pow(h, 3);
  double beta_Z = -2*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 2)
    + 15*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 21*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - 6*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta_eta = 1.5*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 3*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta_Z_prime = -c3*h0/pow(*point, 2)*(-4*pow(xgal, 3)*pow(h, 3) + 4*pow(xgal, 3)*pow(h, 2)*hprime + 6*pow(xgal, 2)*pow(h, 3)*xprime)
    + c4*h0/pow(*point, 4)*(-60*pow(xgal, 4)*pow(h, 5) + 60*pow(xgal, 4)*pow(h, 4)*hprime + 60*pow(xgal, 3)*pow(h, 5)*xprime)
    - c5*h0/pow(*point, 6)*(-126*pow(xgal, 5)*pow(h, 7) + 126*pow(xgal, 5)*pow(h, 6)*hprime + 105*pow(xgal, 4)*pow(h, 7)*xprime)
    - cG*h0/pow(*point, 2)*(-12*pow(xgal, 2)*pow(h, 3) + 12*pow(xgal, 2)*pow(h, 2)*hprime + 12*xgal*pow(h, 3)*xprime);
  double beta_eta_prime = c4*h0/pow(*point, 4)*(-6*pow(xgal, 4)*pow(h, 5) + 6*pow(xgal, 4)*pow(h, 4)*hprime + 6*pow(xgal, 3)*pow(h, 5)*xprime)
    - c5*h0/pow(*point, 6)*(-18*pow(xgal, 5)*pow(h, 7) + 18*pow(xgal, 5)*pow(h, 6)*hprime + 15*pow(xgal, 4)*pow(h, 7)*xprime)
    - cG*h0/pow(*point, 2)*(-2*pow(xgal, 2)*pow(h, 3) + 2*pow(xgal, 2)*pow(h, 2)*hprime + 2*xgal*pow(h, 3)*xprime);

//   // kZ' = A*(*dphisecond) + B
//   // double A = 0.5*alpha_gammasecond/((1-0.5*alpha_Z)*(*k));
//   // double B = (*dgrho)/(*k)*(1-0.5*hprime/h) + 0.5*(*dgq)/(h0*h) - hprime/h*(*eta) + 0.5*dotdeltaf/((*k)*h0*h) + 0.5*dotdeltagal/((*k)*h0*h);
//   // double A = 0.5*alpha_gammasecond/((*k)*(1-0.5*alpha_Z));
//   // double A = 0.5*alpha_gammasecond/(1-0.5*alpha_Z);
//   // // double B = ((h0*h - 0.5*(h0*hprime + alpha_Z*h0*h) + 0.25*(alpha_Z_prime + alpha_Z*h0*hprime))*(*dgrho) + 0.5*(1 - alpha_eta)*(*k)*(*dgq) - (h0*hprime + alpha_Z*h0*h + alpha_eta_prime - 2*alpha_eta*h0*h - 0.5*(alpha_Z_prime + alpha_Z*h0*hprime))*(*k)*(*eta) + 0.5*dotdeltaf + 0.5*chiprimehat)/((*k)*h0*h*(1-0.5*alpha_Z));
//   // double B = ((h0*h - 0.5*(h0*hprime + alpha_Z*h0*h) + 0.25*(alpha_Z_prime + alpha_Z*h0*hprime))*(*dgrho) + 0.5*(1 - alpha_eta)*(*k)*(*dgq) - (h0*hprime + alpha_Z*h0*h + alpha_eta_prime - 2*alpha_eta*h0*h - 0.5*(alpha_Z_prime + alpha_Z*h0*hprime))*(*k)*(*eta) + 0.5*dotdeltaf + 0.5*chiprimehat)/(h0*h*(1-0.5*alpha_Z));

//   // double A1 = (*k)*h0*h*(2.-alpha_Z);
//   // double B1 = -alpha_gammasecond;
//   // double C1 = ((h0*h - 0.5*(h0*hprime + alpha_Z*h0*h) + 0.25*(alpha_Z_prime + alpha_Z*h0*hprime))*(*dgrho) + 0.5*(1 - alpha_eta)*(*k)*(*dgq) - (h0*hprime + alpha_Z*h0*h + alpha_eta_prime - 2*alpha_eta*h0*h - 0.5*(alpha_Z_prime + alpha_Z*h0*hprime))*(*k)*(*eta) + 0.5*dotdeltaf + 0.5*chiprimehat);
//   // double C1 = 2*h0*h*(*dgrho) + (1-alpha_eta)*(*k)*(*dgq) + dotdeltaf + (alpha_Z_prime + alpha_Z*h0*hprime - 2*alpha_Z*h0*h - 2*h0*hprime)*(*k)*h0*h*(*z) + (4*alpha_eta - 2*alpha_eta_prime)*(*k)*(*eta) + chiprimehat;
//   // double C1 = dotdeltaf + (((*grho)+(*gpres))/(2.*h0*h) + h0*h + 0.5*alpha_Z_prime + 0.5*alpha_Z*h0*hprime - alpha_Z*h0*h)*(*dgrho) + chiprimehat + (((*grho)+(*gpres))/(h0*h) - 2*h0*h + alpha_Z_prime + alpha_Z*h0*hprime - 2*alpha_Z*h0*h - 2*alpha_eta_prime + 4*alpha_eta*h0*h)*(*k)*(*eta) + (1-alpha_eta)*(*k)*(*dgq);

//   // Expression of dphisecond
//   double beta_gammasecond = c2
//     - 12*c3/pow(*point, 2)*xgal*pow(h, 2)
//     + 54*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)
//     - 60*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)
//     - 6*cG/pow(*point, 2)*pow(h, 2);
//   double beta_gammaprime = 2*c2
//     - c3/pow(*point, 2)*(12*h*xhprime + 12*xgal*h*hprime)
//     + c4/pow(*point, 4)*(-108*pow(xgal, 2)*pow(h, 4) + 108*xgal*pow(h, 3)*xhprime + 108*pow(xgal, 2)*pow(h, 3)*hprime)
//     - c5/pow(*point, 6)*(-240*pow(xgal, 3)*pow(h, 6) + 180*pow(xgal, 2)*pow(h, 5)*xhprime + 180*pow(xgal, 3)*pow(h, 5)*hprime)
//     - 12*cG/pow(*point, 2)*h*hprime;
//   double beta_gamma = c2
//     - c3/pow(*point, 2)*(4*xgal*pow(h, 2) + 4*h*xhprime)
//     + c4/pow(*point, 4)*(-10*pow(xgal, 2)*pow(h, 4) + 24*xgal*pow(h, 3)*xhprime + 12*pow(xgal, 2)*pow(h, 3)*hprime)
//     - c5/pow(*point, 6)*(-36*pow(xgal, 3)*pow(h, 6) + 36*pow(xgal, 2)*pow(h, 5)*xhprime + 24*pow(xgal, 3)*pow(h, 5)*hprime)
//     - cG/pow(*point, 2)*(2*pow(h, 2) + 4*h*hprime);
//   double beta_Z = c2*xgal
//     - c3/pow(*point, 2)*(6*pow(xgal, 2)*pow(h, 2) + 4*xgal*h*xhprime)
//     + c4/pow(*point, 4)*(-6*pow(xgal, 3)*pow(h, 4) + 36*pow(xgal, 2)*pow(h, 3)*xhprime + 12*pow(xgal, 3)*pow(h, 3)*hprime)
//     - c5/pow(*point, 6)*(-45*pow(xgal, 4)*pow(h, 6) + 60*pow(xgal, 3)*pow(h, 5)*xhprime + 30*pow(xgal, 4)*pow(h, 5)*hprime)
//     - cG/pow(*point, 2)*(6*xgal*pow(h, 2) + 4*h*xhprime + 4*xgal*h*hprime);
//   double beta_Zprime = -2*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)
//     + 12*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)
//     - 15*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)
//     - 4*cG/pow(*point, 2)*xgal*pow(h, 2);
//   double beta_eta = c4/pow(*point, 4)*(-4*pow(xgal, 3)*pow(h, 4) + 6*pow(xgal, 2)*pow(h, 3)*xhprime)
//     - c5/pow(*point, 6)*(-12*pow(xgal, 4)*pow(h, 6) + 3*pow(xgal, 4)*pow(h, 5)*hprime + 12*pow(xgal, 3)*pow(h, 5)*xhprime)
//     - 2*cG/pow(*point, 2)*h*xhprime;

//   // (*dphisecond) = C*kZ' + D
//   // double C = -beta_Zprime*(*k)/beta_gammasecond;
//   // double C = -beta_Zprime/beta_gammasecond;
//   // double D = -(beta_gammaprime*h0*h*(*dphiprime) + beta_gamma*pow(*k, 2)*(*dphi) + 0.5*beta_Z*(*dgrho) + (beta_Z - 2*beta_eta)*(*k)*(*eta))/beta_gammasecond;
//   // double A2 = (*k)*beta_Zprime;
//   // double B2 = beta_gammasecond;
//   // double C2 = -(beta_gammaprime*h0*h*(*dphiprime) + beta_gamma*pow(*k, 2)*(*dphi) + 0.5*beta_Z*(*dgrho) + (beta_Z - 2*beta_eta)*(*eta));















//   // double A1 = 2*(*k)*h0*h +
//   //   2*(*k)*c3/pow(*point, 2)*h0*pow(xgal, 3)*pow(h, 3) -
//   //   15*(*k)*c4/pow(*point, 4)*h0*pow(xgal, 4)*pow(h, 5) +
//   //   21*(*k)*c5/pow(*point, 6)*h0*pow(xgal, 5)*pow(h, 7) +
//   //   6*(*k)*cG/pow(*point, 2)*h0*pow(xgal, 2)*pow(h, 3);
//   // double B1 = -c2*h0*xgal*h +
//   //   18*c3/pow(*point, 2)*h0*pow(xgal, 2)*pow(h, 3) -
//   //   90*c4/pow(*point, 4)*h0*pow(xgal, 3)*pow(h, 5) +
//   //   105*c5/pow(*point, 6)*h0*pow(xgal, 4)*pow(h, 7) +
//   //   18*cG/pow(*point, 2)*h0*xgal*pow(h, 3);
//   // double C1 = 2*(*dgrho)*h0*h - 2*(*k)*h0*h*hprime*(*z) + (*k)*(*dgq) +
//   //   + h0*(-2*h0*h)*c2*xgal*h*(*dphiprime)
//   //   - 18*h0*(-4*h0*h)*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*(*dphiprime)
//   //   + 90*h0*(-6*h0*h)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 3)*(*dphiprime)
//   //   - 105*h0*(-8*h0*h)*c5/pow(*point, 6)*pow(xgal, 2)*pow(h, 3)*(*dphiprime)
//   //   - 18*h0*(-4*h0*h)*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*(*dphiprime)

//   //   - 2*pow(*k, 2)*(-4*h0*h)*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)*(*dphi)
//   //   + 12*pow(*k, 2)*(-6*h0*h)*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*(*dphi)
//   //   - 15*pow(*k, 2)*(-8*h0*h)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*(*dphi)
//   //   - 4*pow(*k, 2)*(-4*h0*h)*cG/pow(*point, 2)*xgal*pow(h, 2)*(*dphi)
   
//   //   - 2*h0*(*k)*(-4*h0*h)*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 3)*(*z)
//   //   + 15*h0*(*k)*(-6*h0*h)*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 5)*(*z)
//   //   - 21*h0*(*k)*(-8*h0*h)*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 7)*(*z)
//   //   - 6*h0*(*k)*(-4*h0*h)*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*(*z)

//   //   + (-3)*(*k)*(-6*h0*h)*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)*(*eta)
//   //   - (-6)*(*k)*(-8*h0*h)*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)*(*eta)
//   //   - (-2)*(*k)*(-4*h0*h)*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)*(*eta)

//   //   + c2*h0*h*xhprime*(*dphiprime)
//   //   - 36*pow(h0, 2)*c3/pow(*point, 2)*xgal*pow(h, 3)*xhprime*(*dphiprime)
//   //   - 18*pow(h0, 2)*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphiprime)
//   //   - 2*pow(*k, 2)*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)*(*dphiprime)
//   //   + 270*pow(h0, 2)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 5)*xhprime*(*dphiprime)
//   //   + 80*pow(h0, 2)*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphiprime)
//   //   + 12*pow(*k, 2)*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*(*dphiprime)
//   //   - 420*pow(h0, 2)*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 7)*xhprime*(*dphiprime)
//   //   - 315*pow(h0, 2)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 7)*hprime*(*dphiprime)
//   //   - 15*pow(*k, 2)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*(*dphiprime)
//   //   - 18*pow(h0, 2)*cG/pow(*point, 2)*pow(h, 3)*xhprime*(*dphiprime)
//   //   - 36*pow(h0, 2)*cG/pow(*point, 2)*xgal*pow(h, 3)*hprime*(*dphiprime)
//   //   - 4*pow(*k, 2)*cG/pow(*point, 2)*xgal*pow(h, 2)*(*dphiprime)
    
//   //   - 4*pow(*k, 2)*h0*c3/pow(*point, 2)*xgal*pow(h, 2)*xhprime*(*dphi)
//   //   + 36*pow(*k, 2)*h0*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)*xhprime*(*dphi)
//   //   + 12*pow(*k, 2)*h0*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*hprime*(*dphi)
//   //   - 60*pow(*k, 2)*h0*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)*xhprime*(*dphi)
//   //   - 30*pow(*k, 2)*h0*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*hprime*(*dphi)
//   //   - 4*pow(*k, 2)*h0*cG/pow(*point, 2)*pow(h, 2)*xhprime*(*dphi)
//   //   - 4*pow(*k, 2)*h0*cG/pow(*point, 2)*pow(xgal, 1)*pow(h, 2)*hprime*(*dphi)

//   //   - 6*(*k)*pow(h0, 2)*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*xhprime*(*z)
//   //   + 60*(*k)*pow(h0, 2)*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 5)*xhprime*(*z)
//   //   + 15*(*k)*pow(h0, 2)*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 5)*hprime*(*z)
//   //   - 105*(*k)*pow(h0, 2)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 7)*xhprime*(*z)
//   //   - 42*(*k)*pow(h0, 2)*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 7)*hprime*(*z)
//   //   - 12*(*k)*pow(h0, 2)*cG/pow(*point, 2)*xgal*pow(h, 3)*xhprime*(*z)
//   //   - 6*(*k)*pow(h0, 2)*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*z)

//   //   + (-12)*(*k)*h0*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*xhprime*(*eta)
//   //   - (-30)*(*k)*h0*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*xhprime*(*eta)
//   //   - (-6)*(*k)*h0*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)*hprime*(*eta)
//   //   - (-4)*(*k)*h0*cG/pow(*point, 2)*xgal*pow(h, 2)*xhprime*(*eta)

//   //   + (-1.5)*(*k)*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)*(*dgq)
//   //   - (-3)*(*k)*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)*(*dgq)
//   //   - (-1)*(*k)*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2)*(*dgq)

//   //   + (-3)*h0*h*((*grhoc)*(*clxc) + (*grhob)*(*clxb)) + (*grhoc)*(*clxcdot) + (*grhob)*(*clxbdot)
//   //   + (-4)*h0*h*((*grhog)*(*clxg) + (*grhor)*(*clxr)) + (*grhog)*(*clxgdot) + (*grhor)*(*clxrdot)
//   //   ;

//   // double A2 = -2*pow(*k, 2)*c3/(*point, 2)*pow(xgal, 2)*pow(h, 2) +
//   //   12*pow(*k, 2)*c4/(*point, 4)*pow(xgal, 3)*pow(h, 4) +
//   //   15*pow(*k, 2)*c5/(*point, 6)*pow(xgal, 4)*pow(h, 6) -
//   //   4*pow(*k, 2)*cG/(*point, 2)*xgal*pow(h, 2);
//   // double B2 = (*k)*c2 -
//   //   12*(*k)*c3/pow(*point, 2)*xgal*pow(h, 2) +
//   //   54*(*k)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4) -
//   //   60*(*k)*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6) -
//   //   6*(*k)*cG/pow(*point, 2)*pow(h, 2);
//   // double C2 = -2*(*k)*c2*h0*h*(*dphiprime)
//   //   + 12*(*k)*h0*c3/pow(*point, 2)*pow(h, 2)*xhprime*(*dphiprime)
//   //   + 12*(*k)*h0*c3/pow(*point, 2)*xgal*pow(h, 2)*hprime*(*dphiprime)
//   //   - (-108)*(*k)*h0*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 5)*(*dphiprime)
//   //   - 108*(*k)*h0*c4/pow(*point, 4)*xgal*pow(h, 4)*xhprime*(*dphiprime)
//   //   - 108*(*k)*h0*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)*hprime*(*dphiprime)
//   //   + (-240)*(*k)*h0*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 7)*(*dphiprime)
//   //   + 180*(*k)*h0*c5/pow(*point, 6)*pow(xgal, 2)*pow(h, 6)*xhprime*(*dphiprime)
//   //   + 180*(*k)*h0*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)*hprime*(*dphiprime)
//   //   + 12*(*k)*h0*cG/pow(*point, 2)*pow(h, 2)*hprime*(*dphiprime)

//   //   - pow(*k, 3)*c2*(*dphi)
//   //   + 4*pow(*k, 3)*c3/pow(*point, 2)*xgal*pow(h, 2)*(*dphi)
//   //   + 4*pow(*k, 3)*c3/pow(*point, 2)*h*xhprime*(*dphi)
//   //   - (-10)*pow(*k, 3)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)*(*dphi)
//   //   - 24*pow(*k, 3)*c4/pow(*point, 4)*xgal*pow(h, 3)*xhprime*(*dphi)
//   //   - 12*pow(*k, 3)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphi)
//   //   + (-36)*pow(*k, 3)*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)*(*dphi)
//   //   + 36*pow(*k, 3)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 5)*xhprime*(*dphi)
//   //   + 24*pow(*k, 3)*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphi)
//   //   + 2*pow(*k, 3)*cG/pow(*point, 2)*pow(h, 2)*(*dphi)
//   //   + 4*pow(*k, 3)*cG/pow(*point, 2)*h*hprime*(*dphi)

//   //   - pow(*k, 2)*h0*c2*xgal*h*(*z)
//   //   + 4*pow(*k, 2)*h0*c3/pow(*point, 2)*xgal*pow(h, 2)*xhprime*(*z)
//   //   + 6*pow(*k, 2)*h0*c3/pow(*point, 2)*pow(xgal, 2)*pow(h, 3)*(*z)    
//   //   - (-6)*pow(*k, 2)*h0*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 5)*(*z)    
//   //   - 36*pow(*k, 2)*h0*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 4)*xhprime*(*z)
//   //   - 12*pow(*k, 2)*h0*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*hprime*(*z)
//   //   + (-45)*pow(*k, 2)*h0*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 7)*(*z)    
//   //   + 60*pow(*k, 2)*h0*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 6)*xhprime*(*z)
//   //   + 30*pow(*k, 2)*h0*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*hprime*(*z)
//   //   + 6*pow(*k, 2)*h0*cG/pow(*point, 2)*xgal*pow(h, 3)*(*z)    
//   //   + 4*pow(*k, 2)*h0*cG/pow(*point, 2)*pow(h, 2)*xhprime*(*z)
//   //   + 4*pow(*k, 2)*h0*cG/pow(*point, 2)*xgal*pow(h, 2)*hprime*(*z)

//   //   - 8*pow(*k, 2)*c4/pow(*point, 4)*pow(xgal, 3)*pow(h, 4)*(*eta)
//   //   - (-12)*pow(*k, 2)*c4/pow(*point, 4)*pow(xgal, 2)*pow(h, 3)*xhprime*(*eta)
//   //   + 24*pow(*k, 2)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 6)*(*eta)
//   //   + (-24)*pow(*k, 2)*c5/pow(*point, 6)*pow(xgal, 3)*pow(h, 5)*xhprime*(*eta)
//   //   + (-6)*pow(*k, 2)*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 5)*hprime*(*eta)
//   //   + (-4)*pow(*k, 2)*cG/pow(*point, 2)*h*xhprime*(*eta)
//   //   ;


  double ksi = alpha_Zprime/(h0*h*(2-beta_Z));

  // if(*point >= 9.99999e-7) dphisecond = -(alpha_gammaprime*h0*h*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + ksi*(dotdeltaf + chiprimehat + (1-beta_eta)*(*k)*(*dgq)) + (0.5*alpha_Z + ksi*(0.5*((*grho)+(*gpres))/(h0*h) + h0*h - beta_Z*h0*h + 0.5*beta_Z_prime + 0.5*beta_Z*h0*hprime))*(*dgrho) + (alpha_Z - 2*alpha_eta + ksi*(((*grho)+(*gpres))/(h0*h) - 2*h0*h + beta_Z_prime + beta_Z*h0*hprime - 2*beta_Z*h0*h - 2*beta_eta_prime + 4*beta_eta*h0*h))*(*k)*(*eta))/(alpha_gammasecond + beta_gammasecond*ksi);
  if(*point >= 1e-4) dphisecond = -(alpha_gammaprime*h0*h*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + ksi*(dotdeltaf + chiprimehat + (1-beta_eta)*(*k)*(*dgq)) + (0.5*alpha_Z + ksi*(0.5*((*grho)+(*gpres))/(h0*h) + h0*h - beta_Z*h0*h + 0.5*beta_Z_prime + 0.5*beta_Z*h0*hprime))*(*dgrho) + (alpha_Z - 2*alpha_eta + ksi*(((*grho)+(*gpres))/(h0*h) - 2*h0*h + beta_Z_prime + beta_Z*h0*hprime - 2*beta_Z*h0*h - 2*beta_eta_prime + 4*beta_eta*h0*h))*(*k)*(*eta))/(alpha_gammasecond + beta_gammasecond*ksi);


//   // if(*point >= 9.99999e-7) dphisecond = (A2*C1-A1*C2)/(A2*B1-A1*B2);
//   // if(*point >= 1e-4) dphisecond = (A1*C2-A2*C1)/(A1*B2-A2*B1);
//   // if(*point >= 1e-4) dphisecond = (C*B+D)/(1-C*A);

//   // printf("%.16f \t %.16f \t %.16f\n", (*point), (beta_gammaprime*h0*h*(*dphiprime) + beta_gamma*pow(*k, 2)*(*dphi) + 0.5*(beta_Z - 2*beta_Zprime)*(*dgrho) + (beta_Z - beta_Zprime - 2*beta_eta)*(*k)*(*eta)) - (beta_gammaprime*h0*h*(*dphiprime) + beta_gamma*pow(*k, 2)*(*dphi) + ksi*(dotdeltaf + chiprimehat + (1-alpha_eta)*(*k)*(*dgq)) + (0.5*beta_Z + ksi*(0.5*((*grho)+(*gpres))/(h0*h) + h0*h - alpha_Z*h0*h + 0.5*alpha_Z_prime + 0.5*alpha_Z*h0*hprime))*(*dgrho) + (beta_Z - 2*beta_eta + ksi*(((*grho)+(*gpres))/(h0*h) - 2*h0*h + alpha_Z_prime + alpha_Z*h0*hprime - 2*alpha_Z*h0*h - 2*alpha_eta_prime + 4*alpha_eta*h0*h))*(*k)*(*eta)), dphisecond);
//   // printf("%.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", (*point), beta_gammasecond, alpha_gammasecond, ksi, beta_gammasecond + alpha_gammasecond*ksi);
//   // printf("%.16f \t %.16f \t %.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", (*point), beta_gammaprime*h0*h*(*dphiprime)/(beta_gammasecond + alpha_gammasecond*ksi), beta_gamma*pow(*k, 2)*(*dphi)/(beta_gammasecond + alpha_gammasecond*ksi), ksi*(dotdeltaf + chiprimehat + (1-alpha_eta)*(*k)*(*dgq))/(beta_gammasecond + alpha_gammasecond*ksi), (0.5*beta_Z + ksi*(0.5*((*grho)+(*gpres))/(h0*h) + h0*h - alpha_Z*h0*h + 0.5*alpha_Z_prime + 0.5*alpha_Z*h0*h))*(*dgrho)/(beta_gammasecond + alpha_gammasecond*ksi), (beta_Z - 2*beta_eta + ksi*(((*grho)+(*gpres))/(h0*h) - 2*h0*h + alpha_Z_prime + alpha_Z*h0*hprime - 2*alpha_Z*h0*h - 2*alpha_eta_prime + 4*alpha_eta*h0*h))*(*k)*(*eta)/(beta_gammasecond + alpha_gammasecond*ksi), dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), A, C, B, D, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), A1, A2, B1, B2, C1, C2, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), A1, A2, B1, B2, A2*B1, A1*B2, A2*B1-A1*B2, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), A1, A2, C1, C2, A2*C1, A1*C2, A2*C1-A1*C2, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), (h0*h - 0.5*(h0*hprime + alpha_Z*h0*h) + 0.25*(alpha_Z_prime + alpha_Z*h0*hprime))*(*dgrho), 0.5*(1 - alpha_eta)*(*k)*(*dgq), (h0*hprime + alpha_Z*h0*h + alpha_eta_prime - 2*alpha_eta*h0*h - 0.5*(alpha_Z_prime + alpha_Z*h0*hprime))*(*k)*(*eta), 0.5*dotdeltaf, 0.5*chiprimehat, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), A2*C1, A1*C2, A1*B2, A2*B1, A2*C1-A1*C2, A2*B1-A1*B2, dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), C*A, 0.5*alpha_gammasecond, (1-0.5*alpha_Z), dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), (*dgrho), (*eta), (*z), (*dgq), dphisecond);
//   // printf("%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f\n", (*point), chiprimehat, alpha_gammasecond, alpha_Z, alpha_eta, alpha_Z_prime, alpha_eta_prime);
//   // FILE* g = fopen("dphisecond/dphisecond_q0001.dat", "a");
//   // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", (*point), dphisecond, alpha_gammaprime/alpha_gammasecond*h0*h*(*dphiprime), alpha_gamma/alpha_gammasecond*pow(*k, 2)*(*dphi), 0.5/alpha_gammasecond*(alpha_Z - 2*alpha_Zprime)*(*dgrho), (alpha_Z - alpha_Zprime - 2*alpha_eta)/alpha_gammasecond*(*k)*(*eta), (*dphiprime), (*dphi));
//   // fclose(g);

  return dphisecond;

}

// Check if conservation equations are satisfied by perturbations
extern "C" double* conservation_(double* grho, double* gpres, double* dgrho, double* grhob, double* clxb, double* clxbdot, double* grhoc, double* clxc, double* clxcdot, double* grhor, double* clxr, double* clxrdot, double* grhog, double* clxg, double* clxgdot, double* dgq, double* qr, double* qrdot, double* qg, double* qgdot, double* dgpi, double* eta, double* dphi, double* dphiprime, double* dphisecond, double* point, double* k){

  static double cons[2];

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)); // H = adot/a
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c2/6*hoft*xgal-3*c3*pow(hoft, 3)*pow(xgal, 2) + 15*c4*pow(hoft, 5)*pow(xgal, 3) - 17.5*c5*pow(hoft, 7)*pow(xgal, 4) - 3*cG*pow(hoft, 3)*xgal;
  double gamma = c2/3*pow(hoft, 2)*xgal-c3*pow(hoft, 4)*pow(xgal, 2) + 2.5*c5*pow(hoft, 8)*pow(xgal, 4) - 2*cG*pow(hoft, 4)*xgal;
  double beta = c2/6*pow(hoft, 2) -2*c3*pow(hoft, 4)*xgal + 9*c4*pow(hoft, 6)*pow(xgal, 2) - 10*c5*pow(hoft, 8)*pow(xgal, 3) - cG*pow(hoft, 4);
  double delta = 2*hoft + 2*c3*pow(hoft, 3)*pow(xgal, 3) - 15*c4*pow(hoft, 5)*pow(xgal, 4) + 21*c5*pow(hoft, 7)*pow(xgal, 5) + 6*cG*pow(hoft, 3)*pow(xgal, 2);
  double lambda = 3*pow(hoft, 2) + orad/pow(*point, 4) + c2/2*pow(hoft, 2)*pow(xgal, 2) - 2*c3*pow(hoft, 4)*pow(xgal, 3) + 7.5*c4*pow(hoft, 6)*pow(xgal, 4) - 9*c5*pow(hoft, 8)*pow(xgal, 5) - cG*pow(hoft, 4)*pow(xgal, 2);
  double omega = 2*c3*pow(hoft, 4)*pow(xgal, 2) - 12*c4*pow(hoft, 6)*pow(xgal, 3) + 15*c5*pow(hoft, 8)*pow(xgal, 4) + 4*cG*pow(hoft, 4)*xgal;

  double xprime = -xgal+(alpha*lambda-delta*gamma)/(delta*beta-alpha*omega); // derivative wrt ln(a)
  double hprime = (*point)*(omega*gamma-lambda*beta)/(delta*beta-alpha*omega) + h; // Careful, this is the derivative of h and not hoft
  // double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*h)'

  // Derivatives of density perturbations
  double dotdeltab = (*grhob)*((*clxbdot) - 3*h0*h*(*clxb));
  double dotdeltac = (*grhoc)*((*clxcdot) - 3*h0*h*(*clxc));
  double dotdeltar = (*grhor)*((*clxrdot) - 4*h0*h*(*clxr));
  double dotdeltag = (*grhog)*((*clxgdot) - 4*h0*h*(*clxg));
  double dotdelta = dotdeltab + dotdeltac + dotdeltar + dotdeltag;

  // Derivatives of heat flux perturbations
  double dotdeltaqr = (*grhor)*((*qrdot) - 4*h0*h*(*qr));
  double dotdeltaqg = (*grhog)*((*qgdot) - 4*h0*h*(*qg));
  double dotdeltaq = dotdeltaqr + dotdeltaqg;

  // Derivative of galileon's density perturbation
  double chitildeprime = c2*(h0*xgal*h*(*dphisecond) - 2*pow(h0, 2)*xgal*pow(h, 2)*(*dphiprime) + pow(h0, 2)*xgal*h*hprime*(*dphiprime) + pow(h0, 2)*pow(h, 2)*xprime*(*dphiprime))
    - c3/pow(*point, 2)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphisecond) - 72*pow(h0, 2)*pow(xgal, 2)*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphiprime) + 36*pow(h0, 2)*xgal*pow(h, 4)*xprime*(*dphiprime)
			 + pow(*k, 2)*(2*pow(xgal, 2)*pow(h, 2)*(*dphiprime) - 8*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) + 4*h0*pow(xgal, 2)*pow(h, 2)*hprime*(*dphi) + 4*h0*xgal*pow(h, 3)*xprime*(*dphi)))
    + c4/pow(*point, 4)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphisecond) - 540*pow(h0, 2)*pow(xgal, 3)*pow(h, 6)*(*dphiprime) + 450*pow(h0, 2)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphiprime) + 270*pow(h0, 2)*pow(xgal, 2)*pow(h, 6)*xprime*(*dphiprime)
			 + pow(*k, 2)*(12*pow(xgal, 3)*pow(h, 4)*(*dphiprime) - 72*h0*pow(xgal, 3)*pow(h, 5)*(*dphi) + 48*h0*pow(xgal, 3)*pow(h, 4)*hprime*(*dphi) + 36*h0*pow(xgal, 2)*pow(h, 5)*xprime*(*dphi)))
    - c5/pow(*point, 6)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphisecond) - 840*pow(h0, 2)*pow(xgal, 4)*pow(h, 8)*(*dphiprime) + 735*pow(h0, 2)*pow(xgal, 4)*pow(h, 7)*hprime*(*dphiprime) + 420*pow(h0, 2)*pow(xgal, 3)*pow(h, 8)*xprime*(*dphiprime)
			 + pow(*k, 2)*(15*pow(xgal, 4)*pow(h, 6)*(*dphiprime) - 120*h0*pow(xgal, 4)*pow(h, 7)*(*dphi) + 90*h0*pow(xgal, 4)*pow(h, 6)*hprime*(*dphi) + 60*h0*pow(xgal, 3)*pow(h, 7)*xprime*(*dphi)))
    - cG/pow(*point, 2)*(18*h0*xgal*pow(h, 3)*(*dphisecond) - 72*pow(h0, 2)*xgal*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*xgal*pow(h, 3)*hprime*(*dphiprime) + 18*pow(h0, 2)*pow(h, 4)*xprime*(*dphiprime)
			 + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphiprime) - 16*h0*xgal*pow(h, 3)*(*dphi) + 8*h0*xgal*pow(h, 2)*hprime*(*dphi) + 4*h0*pow(h, 3)*xprime*(*dphi)))
    ;
  double alpha_Z = -c3/pow(*point, 2)*(6*pow(xgal, 2)*pow(h, 2)*xprime + 6*pow(xgal, 3)*h*hprime - 8*pow(xgal, 3)*pow(h, 2))
    + c4/pow(*point, 4)*(75*pow(xgal, 4)*(h, 3)*hprime + 60*pow(xgal, 3)*pow(h, 4)*xprime - 90*pow(xgal, 4)*pow(h, 4))
    - c5/pow(*point, 6)*(147*pow(xgal, 5)*pow(h, 5)*hprime + 105*pow(xgal, 4)*pow(h, 6)*xprime - 168*pow(xgal, 5)*pow(h, 6))
    - cG/pow(*point, 2)*(18*pow(xgal, 2)*h*hprime + 12*xgal*pow(h, 2)*xprime - 24*pow(xgal, 2)*pow(h, 2));
  double alpha_eta = c4/pow(*point, 4)*(6*pow(xgal, 4)*pow(h, 3)*hprime + 6*pow(xgal, 3)*pow(h, 4)*xprime - 9*pow(xgal, 4)*pow(h, 4))
    - c5/pow(*point, 6)*(18*pow(xgal, 5)*pow(h, 5)*hprime + 15*pow(xgal, 4)*pow(h, 6)*xprime - 24*pow(xgal, 5)*pow(h, 6))
    - cG/pow(*point, 2)*(2*pow(xgal, 2)*h*hprime + 2*xgal*pow(h, 2)*xprime - 4*pow(xgal, 2)*pow(h, 2));
  double alpha_Zprime = -2*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 2) 
    + 15*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 21*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - 6*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double alpha_etaprime = 1.5*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 3*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double dotdeltagal = chitildeprime - alpha_etaprime*(*k)*(*dgq) + (0.5*alpha_Z - alpha_Zprime)*h0*h*(*dgrho) + (alpha_Z - alpha_Zprime - 2*alpha_eta)*h0*h*(*k)*(*eta);
  
  // Derivative of galileon's heat flux perturbation
  double alphaq_zsigprime = c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 2*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - cG/(1.5*pow(*point, 2))*pow(xgal, 2)*pow(h, 2);
  double alphaprimeq_zsig = c4/pow(*point, 4)*(-6*pow(xgal, 4)*pow(h, 4) + 4*pow(xgal, 4)*pow(h, 3)*hprime + 4*pow(xgal, 3)*pow(h, 4)*xprime)
    - c5/pow(*point, 6)*(-16*pow(xgal, 5)*pow(h, 6) + 12*pow(xgal, 5)*pow(h, 5)*hprime + 10*pow(xgal, 4)*pow(h, 6)*xprime)
    - cG/(1.5*pow(*point, 2))*(-4*pow(xgal, 2)*pow(h, 2) + 2*pow(xgal, 2)*h*hprime + 2*xgal*pow(h, 2)*xprime);
  double qtildeprime = c2*(h0*xgal*h*(*dphiprime) - 2*pow(h0, 2)*xgal*pow(h, 2)*(*dphi) + pow(h0, 2)*xgal*h*hprime*(*dphi) + pow(h0, 2)*pow(h, 2)*xprime*(*dphi))
    - c3/pow(*point, 2)*(-2*pow(xgal, 2)*pow(h, 2)*(*dphisecond) + 14*h0*pow(xgal, 2)*pow(h, 3)*(*dphiprime) - 4*h0*pow(xgal, 2)*pow(h, 2)*hprime*(*dphiprime) - 4*h0*xgal*pow(h, 3)*xprime*(*dphiprime) - 24*pow(h0, 2)*pow(xgal, 2)*pow(h, 4)*(*dphi) + 18*pow(h0, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphi) + 12*pow(h0, 2)*xgal*pow(h, 4)*xprime*(*dphi))
    + c4/pow(*point, 4)*(-12*pow(xgal, 3)*pow(h, 4)*(*dphisecond) + 90*h0*pow(xgal, 3)*pow(h, 5)*(*dphiprime) - 48*h0*pow(xgal, 3)*pow(h, 4)*hprime*(*dphiprime) - 36*h0*pow(xgal, 2)*pow(h, 5)*xprime*(*dphiprime) - 108*pow(h0, 2)*pow(xgal, 3)*pow(h, 6)*(*dphi) + 90*pow(h0, 2)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphi) + 54*pow(h0, 2)*pow(xgal, 2)*pow(h, 6)*xprime*(*dphi))
    - c5/pow(*point, 6)*(-15*pow(xgal, 4)*pow(h, 6)*(*dphisecond) + 135*h0*pow(xgal, 4)*pow(h, 7)*(*dphiprime) - 90*h0*pow(xgal, 4)*pow(h, 6)*hprime*(*dphiprime) - 60*h0*pow(xgal, 3)*pow(h, 7)*xprime*(*dphiprime) - 120*pow(h0, 2)*pow(xgal, 4)*pow(h, 8)*(*dphi) + 105*pow(h0, 2)*pow(xgal, 4)*pow(h, 7)*hprime*(*dphi) + 60*pow(h0, 2)*pow(xgal, 3)*pow(h, 8)*xprime*(*dphi))
    - cG/pow(*point, 2)*(-4*xgal*pow(h, 2)*(*dphisecond) + 22*h0*xgal*pow(h, 3)*(*dphiprime) - 8*h0*xgal*pow(h, 2)*hprime*(*dphiprime) - 4*h0*pow(h, 3)*xprime*(*dphiprime) - 24*pow(h0, 2)*xgal*pow(h, 4)*(*dphi) + 18*pow(h0, 2)*xgal*pow(h, 3)*hprime*(*dphi) + 6*pow(h0, 2)*pow(h, 4)*xprime*(*dphi));
  double dotdeltaqgal = (*k)*qtildeprime - 1.5*(alphaprimeq_zsig + 2*alphaq_zsigprime)*h0*h*(*dgq) - (*k)*alphaq_zsigprime*(*dgpi);

  // printf("%.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", (*point), (*k)*qtildeprime, 1.5*(alphaprimeq_zsig + 2*alphaq_zsigprime)*h0*h*(*dgq), (*k)*alphaq_zsigprime*(*dgpi), dotdeltaqgal);

  dotdelta += dotdeltagal;
  dotdeltaq += dotdeltaqgal;
  
  // Equations 
  cons[0] = (dotdelta + ((*grho)+(*gpres))*(0.5*(*dgrho) + (*k)*(*eta))/(h0*h) + 3*h0*h*(*dgrho) + (*k)*(*dgq))/dotdelta;
  cons[1] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3*(*k)*(*dgpi))/dotdeltaq;

  // printf("conservation 1 : %.16f \t %.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", (*point), ((*grho)+(*gpres))*(0.5*(*dgrho) + (*k)*(*eta))/(h0*h)/dotdelta, 3*h0*h*(*dgrho)/dotdelta, (*k)*(*dgq)/dotdelta, dotdelta, cons[0]);
  // printf("conservation 2 : %.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", (*point), 4*h0*h*(*dgq)/dotdeltaq, 2./3*(*k)*(*dgpi)/dotdeltaq, dotdeltaq, cons[1]);
  // printf("conservation : %.16f \t %.16f\n", (*point), cons[0]);
 
  return cons;

}

// calculate the conformal time derivative of z
extern "C" void zprime_(double* grho, double* dgrho, double* dgq, double* grhob, double* clxb, double* clxbdot, double* grhoc, double* clxc, double* clxcdot, double* grhor, double* clxr, double* clxrdot, double* grhog, double* clxg, double* clxgdot, double* eta, double* dphi, double* dphiprime, double* dphisecond, double* point, double* k){

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : (*point)*sqrt(om/pow((*point), 3)+orad/pow((*point), 4));
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)); // H = adot/a
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c2/6*hoft*xgal-3*c3*pow(hoft, 3)*pow(xgal, 2) + 15*c4*pow(hoft, 5)*pow(xgal, 3) - 17.5*c5*pow(hoft, 7)*pow(xgal, 4) - 3*cG*pow(hoft, 3)*xgal;
  double gamma = c2/3*pow(hoft, 2)*xgal-c3*pow(hoft, 4)*pow(xgal, 2) + 2.5*c5*pow(hoft, 8)*pow(xgal, 4) - 2*cG*pow(hoft, 4)*xgal;
  double beta = c2/6*pow(hoft, 2) -2*c3*pow(hoft, 4)*xgal + 9*c4*pow(hoft, 6)*pow(xgal, 2) - 10*c5*pow(hoft, 8)*pow(xgal, 3) - cG*pow(hoft, 4);
  double delta = 2*hoft + 2*c3*pow(hoft, 3)*pow(xgal, 3) - 15*c4*pow(hoft, 5)*pow(xgal, 4) + 21*c5*pow(hoft, 7)*pow(xgal, 5) + 6*cG*pow(hoft, 3)*pow(xgal, 2);
  double lambda = 3*pow(hoft, 2) + orad/pow(*point, 4) + c2/2*pow(hoft, 2)*pow(xgal, 2) - 2*c3*pow(hoft, 4)*pow(xgal, 3) + 7.5*c4*pow(hoft, 6)*pow(xgal, 4) - 9*c5*pow(hoft, 8)*pow(xgal, 5) - cG*pow(hoft, 4)*pow(xgal, 2);
  double omega = 2*c3*pow(hoft, 4)*pow(xgal, 2) - 12*c4*pow(hoft, 6)*pow(xgal, 3) + 15*c5*pow(hoft, 8)*pow(xgal, 4) + 4*cG*pow(hoft, 4)*xgal;

  double xprime = -xgal+(alpha*lambda-delta*gamma)/(delta*beta-alpha*omega); // derivative wrt ln(a)
  double hprime = (*point)*(omega*gamma-lambda*beta)/(delta*beta-alpha*omega) + h; // Careful, this is the derivative of h and not hoft
  // double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*h)'

  // Derivatives of density perturbations
  double dotdeltab = (*grhob)*((*clxbdot) - 3*h0*h*(*clxb));
  double dotdeltac = (*grhoc)*((*clxcdot) - 3*h0*h*(*clxc));
  double dotdeltar = (*grhor)*((*clxrdot) - 4*h0*h*(*clxr));
  double dotdeltag = (*grhog)*((*clxgdot) - 4*h0*h*(*clxg));
  double dotdelta = dotdeltab + dotdeltac + dotdeltar + dotdeltag;

  // Derivative of galileon's density perturbation
  double chitildeprime = c2*(h0*xgal*h*(*dphisecond) - 2*pow(h0, 2)*xgal*pow(h, 2)*(*dphiprime) + pow(h0, 2)*xgal*h*hprime*(*dphiprime) + pow(h0, 2)*pow(h, 2)*xprime*(*dphiprime))
    - c3/pow(*point, 2)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphisecond) - 72*pow(h0, 2)*pow(xgal, 2)*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*pow(xgal, 2)*pow(h, 3)*hprime*(*dphiprime) + 36*pow(h0, 2)*xgal*pow(h, 4)*xprime*(*dphiprime)
			 + pow(*k, 2)*(2*pow(xgal, 2)*pow(h, 2)*(*dphiprime) - 8*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) + 4*h0*pow(xgal, 2)*pow(h, 2)*hprime*(*dphi) + 4*h0*xgal*pow(h, 3)*xprime*(*dphi)))
    + c4/pow(*point, 4)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphisecond) - 540*pow(h0, 2)*pow(xgal, 3)*pow(h, 6)*(*dphiprime) + 450*pow(h0, 2)*pow(xgal, 3)*pow(h, 5)*hprime*(*dphiprime) + 270*pow(h0, 2)*pow(xgal, 2)*pow(h, 6)*xprime*(*dphiprime)
			 + pow(*k, 2)*(12*pow(xgal, 3)*pow(h, 4)*(*dphiprime) - 72*h0*pow(xgal, 3)*pow(h, 5)*(*dphi) + 48*h0*pow(xgal, 3)*pow(h, 4)*hprime*(*dphi) + 36*h0*pow(xgal, 2)*pow(h, 5)*xprime*(*dphi)))
    - c5/pow(*point, 6)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphisecond) - 840*pow(h0, 2)*pow(xgal, 4)*pow(h, 8)*(*dphiprime) + 735*pow(h0, 2)*pow(xgal, 4)*pow(h, 7)*hprime*(*dphiprime) + 420*pow(h0, 2)*pow(xgal, 3)*pow(h, 8)*xprime*(*dphiprime)
			 + pow(*k, 2)*(15*pow(xgal, 4)*pow(h, 6)*(*dphiprime) - 120*h0*pow(xgal, 4)*pow(h, 7)*(*dphi) + 90*h0*pow(xgal, 4)*pow(h, 6)*hprime*(*dphi) + 60*h0*pow(xgal, 3)*pow(h, 7)*xprime*(*dphi)))
    - cG/pow(*point, 2)*(18*h0*xgal*pow(h, 3)*(*dphisecond) - 72*pow(h0, 2)*xgal*pow(h, 4)*(*dphiprime) + 54*pow(h0, 2)*xgal*pow(h, 3)*hprime*(*dphiprime) + 18*pow(h0, 2)*pow(h, 4)*xprime*(*dphiprime)
			 + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphiprime) - 16*h0*xgal*pow(h, 3)*(*dphi) + 8*h0*xgal*pow(h, 2)*hprime*(*dphi) + 4*h0*pow(h, 3)*xprime*(*dphi)))
    ;
  double alpha_Z = -c3/pow(*point, 2)*(6*pow(xgal, 2)*pow(h, 2)*xprime + 6*pow(xgal, 3)*h*hprime - 8*pow(xgal, 3)*pow(h, 2))
    + c4/pow(*point, 4)*(75*pow(xgal, 4)*(h, 3)*hprime + 60*pow(xgal, 3)*pow(h, 4)*xprime - 90*pow(xgal, 4)*pow(h, 4))
    - c5/pow(*point, 6)*(147*pow(xgal, 5)*pow(h, 5)*hprime + 105*pow(xgal, 4)*pow(h, 6)*xprime - 168*pow(xgal, 5)*pow(h, 6))
    - cG/pow(*point, 2)*(18*pow(xgal, 2)*h*hprime + 12*xgal*pow(h, 2)*xprime - 24*pow(xgal, 2)*pow(h, 2));
  double alpha_eta = c4/pow(*point, 4)*(6*pow(xgal, 4)*pow(h, 3)*hprime + 6*pow(xgal, 3)*pow(h, 4)*xprime - 9*pow(xgal, 4)*pow(h, 4))
    - c5/pow(*point, 6)*(18*pow(xgal, 5)*pow(h, 5)*hprime + 15*pow(xgal, 4)*pow(h, 6)*xprime - 24*pow(xgal, 5)*pow(h, 6))
    - cG/pow(*point, 2)*(2*pow(xgal, 2)*h*hprime + 2*xgal*pow(h, 2)*xprime - 4*pow(xgal, 2)*pow(h, 2));
  double alpha_Zprime = -2*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 2) 
    + 15*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 21*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - 6*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double alpha_etaprime = 1.5*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4)
    - 3*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6)
    - cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double dotdeltagal = chitildeprime - alpha_etaprime*(*k)*(*dgq) + (0.5*alpha_Z - alpha_Zprime)*h0*h*(*dgrho) + (alpha_Z - alpha_Zprime - 2*alpha_eta)*h0*h*(*k)*(*eta);

  dotdelta += dotdeltagal;

  double dotz = - hprime*((*eta)/h + (*dgrho)/(2.*(*k)*h)) + dotdelta/(2.*(*k)*h0*h) + (*dgrho)/(*k) + (*dgq)/(2.*h0*h);
  double zprime = - (*eta) - (*dgrho)/(*k);

  // FILE* g = fopen("zprime_comparison_bis.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f\n", (*point), zprime, dotz);
  // fclose(g);

}



int test(){

  fflush(stdout);

  double orad = 8.2987687251764e-5;

  arrays_("params.ini", &orad);
  // FILE* f = fopen("full_integration.txt", "w");

  // for(int i = 2; i< intvar.size()-1; i++){
  //   double point = (intvar[i]+intvar[i+1])/2;
  //   double* hx = handxofa_(&point);
  //   fprintf(f, "%.16f ; %.16f ; %.16f\n", point, (*hx), (*(hx+1)));
  // }

  gsl_spline_free(spline_h);
  gsl_spline_free(spline_x);
  gsl_interp_accel_free(acc);  

  return 0;

}
