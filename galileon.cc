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

  static char* params[14];
  char** param = NULL;

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

  double prod = -zpo*y[0]*y[1]; // prod=h(z)x(z)
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
  alpha = -3*c3*h*prod2 + 15*c4*h2*prod3 + c0*h + c2/6.0*prod - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  gamma = 2*c0*h2 - c3*h2*prod2 + c2/3.0*h*prod + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
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
	//std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
	// Careful, can't throw exceptions when linking to fortran with ifort
	//throw CosFitterExcept("lumdist","calcValOmC2C3C4C5CGC0", "NaN value!",1);
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
	//std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
	// Careful, can't throw exceptions when linking to fortran with ifort
	//throw CosFitterExcept("lumdist","calcValOmC2C3C4C5CGC0", "NaN value!",1);
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
      x[i] = intvar[i]/pow(hubble[i],2);
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
  //const double h0 = 71; // Present-day Hubble constant (WMAP7)

  double params[1];
  bool useacoord = (coord == BARREIRA || coord == JNEVEUa);
  params[0] = (int)useacoord;

  // Vector y = ( dh/dz, dx/dz, dy/dz )
  double y[3] = {cond[0], cond[1], cond[2]};  //Inital value of integral
  hubble[0] = cond[0];
  x[0] = cond[1];

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
      double OmegaP = (0.5*c2*pow(y[0], 2)*pow(intvar[i]*y[1], 2) - 6*c3*pow(y[0], 4)*pow(intvar[i]*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(intvar[i]*y[1], 4) - 21*c5*pow(y[0], 8)*pow(intvar[i]*y[1], 5) - 9*cG*pow(y[0], 4)*pow(intvar[i]*y[1], 2))/(3.0*pow(y[0], 2)) ;
      double OmegaM = 1 - OmegaP - orad/(pow(a,4)*pow(y[0], 2));
      double OmTest = om/(pow(a,3)*pow(y[0], 2));
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
  age *= Mpc/Gyr/h0;

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
  //std::cout << "Initial value : " << cond[3] << endl;

  AgeofU ini_vec_age = AgeofU(0, 1000);
  std::vector<AgeofU> aou(5, ini_vec_age);

  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n);
  double z[2*(n-1)]; // Array of complex numbers ordered like {Re[0],Im[0],Re[1],Im[1],...}

  status = gsl_poly_complex_solve(coeff, n, w, z); // Solve coeff[0]+coeff[1]*x+...+coeff[5]*x^5=0 and store sol in z

  int nrealsol = 0; // Number of real solutions

  for(int j = 0; j<n-1; j++){
    //printf("Root number %i is %e + %e*i\n", j, z[2*j], z[2*j+1]);
    //std::cout << "Root number " << j << " is " << z[2*j] << " + " << z[2*j+1] << "i" << endl;
    if(z[2*j+1]==0){
      cond[1] = z[2*j]/intvar[1]; // Initial condition for d_phi/da
      aou[j].i = j;
      int status2 = ageOfUniverse(aou[j].a, cond);
      aou[j].a -= age;
      if(status2 != 0) return status2;
    }
  }
  
  for(int j = 0; j<n-1; j++){
    if(aou[j].a != 1000){
      nrealsol++;
      //printf("Root number %i gives age of universe equal to %f\n", aou[j].i, aou[j].a + age);
      //std::cout << "Root number " << aou[j].i << " gives age of universe equal to " << aou[j].a + age << endl; 
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
  //printf("Initial condition for x is %e, it is the root number %i and it gives an age of the universe equal to %f\n", xi, aou[0].i, aou[0].a);
  //std::cout << "Initial condition for x is " << xi << ", it is the root number " << aou[0].i << " and it gives an age of the universe equal to " << aou[0].a << endl;

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
  // char* default_name = "output.txt";

  c2 = (*params != "") ? atof(*params) : 0;
  c3 = (*(params+1) != "") ? atof(*(params+1)) : 0;
  c4 = (*(params+2) != "") ? atof(*(params+2)) : 0;
  c5 = (*(params+3) != "") ? atof(*(params+3)) : 0;
  c0 = (*(params+4) != "") ? atof(*(params+4)) : 0;
  cG = (*(params+5) != "") ? atof(*(params+5)) : 0;
  double ratio_rho = (*(params+6) != "") ? atof(*(params+6)) : 0;
  double age = (*(params+7) != "") ? atof(*(params+7)) : 0;
  // om = (*(params+8) != "") ? atof(*(params+8)) : 0;
  // printf("(%f + %f)/%f^2\n", atof(*(params+8)), atof(*(params+9)), atof(*(params+12))/100);
  om = (*(params+8) != "" && *(params+9) != "" && *(params+12) != "") ? (atof(*(params+8))+atof(*(params+9)))/pow(atof(*(params+12)), 2)*10000 : 0;
  // orad = (*(params+9) != "") ? atof(*(params+9)) : 0;
  orad = (*omegar);
  h0 = (*(params+12) != "") ? atof(*(params+12))*1000/lightspeed : 0;
  //char* name = (*(params+10) != "") ? *(params+10) : default_name;
  // if(name != "") name[strlen(name)-1] = '\0';
  // if(strlen(name) != 0) name[strlen(name)-1] = '\0';
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
  double dotphi = (*(params+11) != "") ? atof(*(params+11)) : 0;
  char* outfile = *(params+13);
  if(outfile != ""){
    outfile[strlen(outfile)-1] = '\0';
    outfile = strcat(outfile, "_");
  }
  outfile = strcat(outfile, "background.dat");
  

  printf("Input file : %s\nOutput file : %s\n", infile, outfile);
  printf("OmegaM0 = %f\nOmegaR0 = %.18f\nc0 = %f\nc2 = %f\nc3 = %f\nc4 = %f\nc5 = %f\ncG = %f\nh0 = %f Mpc-1\n", om, orad, c0, c2, c3, c4, c5, cG, h0);

  // The status of the integration, 0 = it went well
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
    //intvar.push_back(1.);

    printf("Number of points : %i\n", intvar.size());

    // Calculate initial conditions in H and ratio_rho
    double xi = 0;
    double H[2];
    if(dotphi == 0){
      H[0] = sqrt(om/pow(apremin, 3)*(1+ratio_rho)+orad/pow(apremin, 4));
      H[1] = sqrt(om/pow(amin, 3)*(1+ratio_rho)+orad/pow(amin, 4));
    } else{
      H[0] = sqrt(om/pow(apremin, 3)+orad/pow(apremin, 4));
      H[1] = sqrt(om/pow(amin, 3)+orad/pow(amin, 4));
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

  if(i<=3){
    double hubble_LCDM = sqrt(om/pow((*point), 3)+orad/pow((*point), 4)+(1-om-orad));
    if(fabs((hx[0]-hubble_LCDM)/hx[0])>1e-3) fprintf(stderr, "Warning : no continuity between LCDM and galileon background at very early time ( i = %i, a = %f, h_LCDM = %f and h_gal = %f)", i, (*point), hubble_LCDM, hx[0]);
  }

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
  
  return hx;

}


extern "C" double grhogal_(double* point){

  double* hx = handxofa_(point);
  double h = (*hx);
  double xgal = (*point)*(*(hx+1)); // here take x as a function of ln(a)

  double grhogal = 0;

  //grhogal = pow(h0, 2)*1e6*(pow(*point, 2)/3*(c2/2*pow(h, 2)*pow(xgal, 2) - 6*c3*pow(h, 4)*pow(xgal, 3) + 22.5*c4*pow(h, 6)*pow(xgal, 4) - 21*c5*pow(h, 8)*pow(xgal, 5) - 9*cG*pow(h, 4)*pow(xgal, 2)))/pow(lightspeed, 2);
  grhogal = pow(h0, 2)*pow(*point, 2)/3*(c2/2*pow(h, 2)*pow(xgal, 2) - 6*c3*pow(h, 4)*pow(xgal, 3) + 22.5*c4*pow(h, 6)*pow(xgal, 4) - 21*c5*pow(h, 8)*pow(xgal, 5) - 9*cG*pow(h, 4)*pow(xgal, 2));

  //printf("grhogal : %.12f \t %.12f \t %.12f \t %.12f\n", *point, h, xgal, grhogal);

  return grhogal;

}


extern "C" double gpresgal_(double* point){

  double* hx = handxofa_(point);
  double h = (*hx);
  double xgal = (*point)*(*(hx+1)); // here take x as a function of ln(a)
  
    
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

  //gpresgal = pow(h0, 2)*1e6*(pow(*point, 2)*(c2/2*pow(h, 2)*pow(xgal, 2) + 2*c3*pow(h, 3)*pow(xgal, 2)*(hprime*xgal+xprime*h) - c4*(4.5*pow(h, 6)*pow(xgal, 4) + 12*pow(h, 6)*pow(xgal, 3)*xprime + 15*pow(h, 5)*pow(xgal, 4)*hprime) + 3*c5*pow(h, 7)*pow(xgal, 4)*(5*h*xprime+7*hprime*xgal+2*h*xgal) + cG*(6*pow(h, 3)*pow(xgal, 2)*hprime + 4*pow(h, 4)*xgal*xprime + 3*pow(h, 4)*pow(xgal, 2))))/pow(lightspeed, 2);
  gpresgal = pow(h0, 2)*pow(*point, 2)*(c2/2*pow(h, 2)*pow(xgal, 2) + 2*c3*pow(h, 3)*pow(xgal, 2)*(hprime*xgal+xprime*h) - c4*(4.5*pow(h, 6)*pow(xgal, 4) + 12*pow(h, 6)*pow(xgal, 3)*xprime + 15*pow(h, 5)*pow(xgal, 4)*hprime) + 3*c5*pow(h, 7)*pow(xgal, 4)*(5*h*xprime+7*hprime*xgal+2*h*xgal) + cG*(6*pow(h, 3)*pow(xgal, 2)*hprime + 4*pow(h, 4)*xgal*xprime + 3*pow(h, 4)*pow(xgal, 2)));

  return gpresgal;

}


// int main(){

//   fflush(stdout);

//   // Clear three global vectors
//   intvar.clear();
//   hubble.clear();
//   x.clear();

//   char** params = readParamsFile("galileon.ini");

//   char* default_name = "output.txt";

//    c2 = (*params != "") ? atof(*params) : 0;
//    c3 = (*(params+1) != "") ? atof(*(params+1)) : 0;
//    c4 = (*(params+2) != "") ? atof(*(params+2)) : 0;
//    c5 = (*(params+3) != "") ? atof(*(params+3)) : 0;
//    c0 = (*(params+4) != "") ? atof(*(params+4)) : 0;
//    cG = (*(params+5) != "") ? atof(*(params+5)) : 0;
//    double ratio_rho = (*(params+6) != "") ? atof(*(params+6)) : 0;
//    double age = (*(params+7) != "") ? atof(*(params+7)) : 0;
//    om = (*(params+8) != "") ? atof(*(params+8)) : 0;
//    orad = (*(params+9) != "") ? atof(*(params+9)) : 0;
//    char* name = (*(params+10) != "") ? *(params+10) : default_name;
//    if(name != "") name[strlen(name)-1] = '\0';
//    // double barreira = (*(params+11) != "") ? atof(*(params+11)) : 0;
//   char* solvingMethod = *(params+11);
//   if(solvingMethod != "") solvingMethod[strlen(solvingMethod)-1] = '\0';
//   if(strcmp(solvingMethod, "JNEVEUz") != 0 && strcmp(solvingMethod, "JNEVEUa") != 0 && strcmp(solvingMethod, "BARREIRA") != 0){
//     fprintf(stderr, "WARNING : invalid integration method, will use the default one\n");
//   } else if(strcmp(solvingMethod, "JNEVEUz") == 0){
//     coord = JNEVEUz;
//   } else if(strcmp(solvingMethod, "JNEVEUa") == 0){
//     coord = JNEVEUa;
//   } else if(strcmp(solvingMethod, "BARREIRA") == 0){
//     coord = BARREIRA;
//   }
//    double dotphi = (*(params+12) != "") ? atof(*(params+12)) : 0;

//    //printf("OmegaM0 = %.2f\nOmegaR0 = %.2e\nc0 = %.2f\nc2 = %.2f\nc3 = %.2f\nc4 = %.2f\nc5 = %.2f\ncG = %.2f\n", om, orad, c0, c2, c3, c4, c5, cG);

//   // double coord = 0;
//   int status = 0;

//   // File to store h and x as a function of z
//   FILE* f;
  
//   /*
//   // Set of parameters (Best fit 2015 BAO+Planck+Lya from Jeremy)
//   double h = 0.762;
//   double OmegaM0 = 0.261175;
//   double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
//   double om = OmegaM0;
//   double orad = OmegaR0;
//   double c2 = -5.25803;
//   double c3 = -1.13137;
//   double cG = 0;
//   double c4 = 1.0/27*(30*(om+orad-1)-9*c2+24*c3-6*cG);
//   double c5 = 1.0/3*(4*(om+orad-1)-c2+2*c3-2*cG);
//   double c0 = 0;
//   char* name = "tracker_test.dat";
//   */
  

//   /*
//   // Best fit Unc 2015 All+JLA+noH0prior from Jeremy
//   double h = 0.736;
//   double OmegaM0 = 0.275;
//   double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
//   double om = OmegaM0;
//   double orad = OmegaR0;
//   double c2 = -4.145;
//   double c3 = -1.545;
//   double cG = 0;
//   double c4 = -0.776;
//   double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
//   double c0 = 0;
//   char* name = "5params_all_2015_combined_noH0prior_bis_bis.dat";
//   */
  

//   /*
//   // Best fit + cG +1 sigma 2015 All+JLA+noH0prior from Jeremy
//   double h = 0.727;
//   double OmegaM0 = 0.280;
//   double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
//   double om = OmegaM0;
//   double orad = OmegaR0;
//   double c2 = -3.434;
//   double c3 = -1.062;
//   double cG = 0.235;
//   double c4 = -0.610;
//   double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
//   double c0 = 0;
//   char* name = "6params_all_2015_combined_noH0prior.dat";
//   */

  
//   // Galileon parameters from Barreira
//   // double c3 = 12.8;
//   // double c4 = -1.7;
//   // double c5 = 1.0;
//   // double cG = 0;
//   // double ratio_rho = 1e-4;
//   // double c2 = -27.0;
//   // double age = 13.978;
//   // double c0 = 0;
//   // double orad = 8e-5;
//   // double om = 0.265;
//   // char* name = "barreira_galileon_1_1e-4.dat";
//   // barreira = 1;
  


//   // if(barreira)
//   if(coord == BARREIRA)
//   {
//     printf("barreira\n");

//     f = fopen(name, "w");

//     double apremin = 1e-7;
//     double amin = 9.99999e-7;
//     // std::vector<double> acoord;
//     // acoord.push_back(apremin);
//     intvar.push_back(apremin);
//     //SetAcoord(acoord, amin, barreira);
//     //sort(acoord.begin(), acoord.end());

//     double q = 1.+5./10000;
//     for(int i = 0; i<log(1/amin)/log(q) ;i++){
//       // acoord.push_back(amin*pow(q, i));
//       intvar.push_back(amin*pow(q, i));
//     }

//     // printf("%i\n", acoord.size());
//     printf("%i\n", intvar.size());

//     double xi = 0;
//     double H[2];
//     if(dotphi == 0){
//       H[0] = sqrt(om/pow(apremin, 3)*(1+ratio_rho)+orad/pow(apremin, 4));
//       H[1] = sqrt(om/pow(amin, 3)*(1+ratio_rho)+orad/pow(amin, 4));
//     } else{
//       H[0] = sqrt(om/pow(apremin, 3)+orad/pow(apremin, 4));
//       H[1] = sqrt(om/pow(amin, 3)+orad/pow(amin, 4));
//       // ratio_rho = pow(acoord[1], 3)/(3*om)*(0.5*c2*pow(dotphi*H[1], 2)-6*c3*H[1]*pow(dotphi*H[1], 3)+22.5*c4*pow(H[1], 2)*pow(dotphi*H[1], 4)-21*c5*pow(H[1], 3)*pow(dotphi*H[1], 5)-9*cG*pow(H[1], 2)*pow(dotphi*H[1], 2));
//       ratio_rho = pow(intvar[1], 3)/(3*om)*(0.5*c2*pow(dotphi*H[1], 2)-6*c3*H[1]*pow(dotphi*H[1], 3)+22.5*c4*pow(H[1], 2)*pow(dotphi*H[1], 4)-21*c5*pow(H[1], 3)*pow(dotphi*H[1], 5)-9*cG*pow(H[1], 2)*pow(dotphi*H[1], 2));
//       std::cout << ratio_rho << endl;
//     }

//     // std::vector<double> hubble(acoord.size()-1, 999999);
//     // std::vector<double> x(acoord.size()-1, 999999);
//     hubble.resize(intvar.size()-1, 999999);
//     x.resize(intvar.size()-1, 999999);

//     // int status = initCond(hubble, x, xi, acoord, H, ratio_rho, age, om, orad, c2, c3, c4, c5, cG, c0);
//     int status = initCond(xi, H, ratio_rho, age);

//     printf("OmegaM0 = %f\nOmegaR0 = %f\nc0 = %f\nc2 = %f\nc3 = %f\nc4 = %f\nc5 = %f\ncG = %f\n\n\nstatus = %i\n", om, orad, c0, c2, c3, c4, c5, cG, status);

//     if(status!=0) return 0;


//     for(int i=0; i<hubble.size(); i++){
//         // double alpha = c2/6*hubble[i]*acoord[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*acoord[i+1]*x[i];
// 	// double gamma = c2/3*pow(hubble[i], 2)*acoord[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];
// 	// double beta = c2/6*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*acoord[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
// 	// double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2);
// 	// double lambda = 3*pow(hubble[i], 2) + orad/pow(acoord[i+1], 4) + c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
// 	// double omega = 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];
	
// 	// double x_prime = -acoord[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
// 	// double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
	
// 	// double rho = c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
// 	// double p = c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*(h_prime*acoord[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*acoord[i+1]*x[i]+2*hubble[i]*acoord[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*acoord[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2));
	
// 	// double OmegaM = om/(pow(hubble[i], 2)*pow(acoord[i], 3));
// 	// double OmegaR = orad/(pow(hubble[i], 2)*pow(acoord[i], 4));
// 	// double OmegaP = rho*3*pow(hubble[i], 2);
	
// 	// double w = p/rho;
// 	// double hubble_LCDM = sqrt(om/pow(acoord[i+1], 3)+orad/pow(acoord[i+1], 4)+(1-om-orad));

// 	// fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", acoord[i+1], hubble[i], hubble_LCDM, w, x[i], OmegaP, OmegaM, OmegaR);
//         double alpha = c2/6*hubble[i]*intvar[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*intvar[i+1]*x[i];
// 	double gamma = c2/3*pow(hubble[i], 2)*intvar[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];
// 	double beta = c2/6*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*intvar[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
// 	double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2);
// 	double lambda = 3*pow(hubble[i], 2) + orad/pow(intvar[i+1], 4) + c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
// 	double omega = 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];

// 	double x_prime = -intvar[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
// 	double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
	
// 	double rho = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
// 	double p = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*(h_prime*intvar[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*intvar[i+1]*x[i]+2*hubble[i]*intvar[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*intvar[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2));
	
// 	double OmegaM = om/(pow(hubble[i], 2)*pow(intvar[i], 3));
// 	double OmegaR = orad/(pow(hubble[i], 2)*pow(intvar[i], 4));
// 	double OmegaP = rho*3*pow(hubble[i], 2);
	
// 	double w = p/rho;
// 	double hubble_LCDM = sqrt(om/pow(intvar[i+1], 3)+orad/pow(intvar[i+1], 4)+(1-om-orad));

// 	fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", intvar[i+1], hubble[i], hubble_LCDM, w, x[i], OmegaP, OmegaM, OmegaR);
//     }

//   } else if(coord == JNEVEUa)
//   {
//     printf("acoord\n");

//     f = fopen(name, "w");

//     double amin = 0.000908;
//     // std::vector<double> acoord;
//     // SetAcoord(acoord, amin, barreira);
//     SetAcoord(amin);

//     // sort(acoord.begin(), acoord.end(), std::greater<double>());
//     sort(intvar.begin(), intvar.end(), std::greater<double>());

//     // std::vector<double> hubble(acoord.size(), 999999);
//     // std::vector<double> x(acoord.size(), 999999);
//     hubble.resize(intvar.size(), 999999);
//     x.resize(intvar.size(), 999999);

//     if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
//     {
//       // status = calcHubbleGalileon(hubble, x, acoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
//       status = calcHubbleGalileon();
//     }
//     else {
//       // calcHubbleTracker(hubble, x, acoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
//       calcHubbleTracker();
//     }


//     printf("OmegaM0 = %f\nOmegaR0 = %f\nc0 = %f\nc2 = %f\nc3 = %f\nc4 = %f\nc5 = %f\ncG = %f\n\n\nstatus = %i\n", om, orad, c0, c2, c3, c4, c5, cG, status);

//     if(status!=0) return 0;

//     for(int i=0; i<hubble.size(); i++)
//     {
//       // double OmegaM = om/(pow(hubble[i], 2)*pow(acoord[i], 3));
//       // double OmegaR = orad/(pow(hubble[i], 2)*pow(acoord[i], 4));
//       // double OmegaP = (c2/2.0*pow(acoord[i]*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(acoord[i]*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(acoord[i]*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

//       // double hubble_LCDM = sqrt(om/pow(acoord[i+1], 3)+orad/pow(acoord[i+1], 4)+(1-om-orad));

//       // fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", 1/acoord[i]-1, hubble[i], hubble_LCDM, x[i], OmegaP, OmegaM, OmegaR);
//       double OmegaM = om/(pow(hubble[i], 2)*pow(intvar[i], 3));
//       double OmegaR = orad/(pow(hubble[i], 2)*pow(intvar[i], 4));
//       double OmegaP = (c2/2.0*pow(intvar[i]*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(intvar[i]*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(intvar[i]*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(intvar[i]*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(intvar[i]*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

//       double hubble_LCDM = sqrt(om/pow(intvar[i+1], 3)+orad/pow(intvar[i+1], 4)+(1-om-orad));

//       fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", 1/intvar[i]-1, hubble[i], hubble_LCDM, x[i], OmegaP, OmegaM, OmegaR);
//     }

// } else if(coord == JNEVEUz)
//   {
//     printf("zcoord\n");

//     f = fopen(name, "w");

//     double zmax = 1100;
//     double amin = 1/(1+zmax);
//     // std::vector<double> zcoord;
//     // SetAcoord(zcoord, amin, barreira);
//     SetAcoord(amin);

//     // sort(zcoord.begin(), zcoord.end(), std::greater<double>());
//     sort(intvar.begin(), intvar.end(), std::greater<double>());

//     // for(int i=0;i<zcoord.size();i++) {zcoord[i]=1/(1+zcoord[i]);}
//     // std::vector<double> hubble(zcoord.size(), 999999);
//     // std::vector<double> x(zcoord.size(), 999999);
//     for(int i=0;i<intvar.size();i++) {intvar[i]=1/(1+intvar[i]);}
//     hubble.resize(intvar.size(), 999999);
//     x.resize(intvar.size(), 999999);

//     double op = 1.5*c2 - 18*c3 + 67.5*c4 - 63*c5 - 27*cG;

//     if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
//     {
//       // status = calcHubbleGalileon(hubble, x, zcoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
//       status = calcHubbleGalileon();
//     }
//     else {
//       printf("Tracker\n");
//       // calcHubbleTracker(hubble, x, zcoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
//       calcHubbleTracker();
//     }

//     printf("OmegaM0 = %f\nOmegaR0 = %f\nOmegaPo = %f\nc0 = %f\nc2 = %f\nc3 = %f\nc4 = %f\nc5 = %f\ncG = %f\n\n\nstatus = %i\n", om, orad, op, c0, c2, c3, c4, c5, cG, status);

//     if(status!=0) return 0;

//     // for(int i=0; i<zcoord.size(); i++)
//     // {
//     //   double OmegaM = om*pow(1+zcoord[i], 3)/pow(hubble[i], 2);
//     //   double OmegaR = orad*pow(1+zcoord[i], 4)/pow(hubble[i], 2);
//     //   double OmegaP = (c2/2.0*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(-(1+zcoord[i])*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(-(1+zcoord[i])*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

//     //   double hubble_LCDM = sqrt(om*pow(1+zcoord[i+1], 3)+orad*pow(1+zcoord[i+1], 4)+(1-om-orad));

//     //   fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", zcoord[i], hubble[i], hubble_LCDM, x[i], OmegaP, OmegaM, OmegaR);
//     // }
//     for(int i=0; i<intvar.size(); i++)
//     {
//       double OmegaM = om*pow(1+intvar[i], 3)/pow(hubble[i], 2);
//       double OmegaR = orad*pow(1+intvar[i], 4)/pow(hubble[i], 2);
//       double OmegaP = (c2/2.0*pow(-(1+intvar[i])*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(-(1+intvar[i])*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(-(1+intvar[i])*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(-(1+intvar[i])*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(-(1+intvar[i])*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

//       double hubble_LCDM = sqrt(om*pow(1+intvar[i+1], 3)+orad*pow(1+intvar[i+1], 4)+(1-om-orad));

//       fprintf(f, "%.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f ; %.12f\n", intvar[i], hubble[i], hubble_LCDM, x[i], OmegaP, OmegaM, OmegaR);
//     }

//   }

//   return 0;

// }


extern "C" double Chigal_(double* dgrho, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  //printf("OmegaM0 = %f\n OmegaR0 = %.18f\n c0 = %f\n c2 = %f\n c3 = %f\n c4 = %f\n c5 = %f\n cG = %f\n", om, orad, c0, c2, c3, c4, c5, cG);

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : sqrt(om/pow((*point), 2)+orad/pow((*point), 3)+(1-om-orad));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = -2*c3/pow(*point, 2)*pow(xgal, 3)*pow(h, 2) 
    + 15*c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 21*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - 6*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta = 1./(1-0.5*alpha);
  double ChitildeG = c2*h0*xgal*h*(*dphiprime) 
    - c3/pow(*point, 2)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphiprime) + 2*pow(*k, 2)*pow(xgal, 2)*pow(h, 2)*(*dphi)) 
    + c4/pow(*point, 4)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphiprime) - 3*(*k)*pow(xgal, 4)*pow(h, 4)*(*eta) + 12*pow(*k, 2)*pow(xgal, 3)*pow(h, 4)*(*dphi)) 
    - c5/pow(*point, 6)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphiprime) - 6*(*k)*pow(xgal, 5)*pow(h, 6)*(*eta) + 15*pow(*k, 2)*pow(xgal, 4)*pow(h, 6)*(*dphi)) 
    - cG/pow(*point, 2)*(18*h0*xgal*pow(h, 3)*(*dphiprime) - 2*(*k)*pow(xgal, 2)*pow(h, 2) + 4*pow(*k, 2)*xgal*pow(h, 2)*(*dphi));

  double ChiG = 0;

  if(alpha == 2) printf("WARNING : 1/beta_chi is zero");
  if (*point >= 9.99999e-7) ChiG = beta*(ChitildeG + 0.5*alpha*((*dgrho) - pow(*k, 2)*(*eta)));

  // double coeff = -c3/pow(*point, 4)*h0*pow(xgal, 3)*pow(h, 2) + 7.5/pow(*point, 6)*c4*h0*pow(xgal, 4)*pow(h, 4) - 10.5/pow(*point, 8)*c5*h0*pow(xgal, 5)*pow(h, 6) - 3/pow(*point, 4)*cG*h0*pow(xgal, 2)*pow(h, 2);
  // double denom = 1-coeff;
  // double Znogal = (*dgrho) - pow(*k, 2)*(*eta);

  // double ChiG = 0; // 8*pi*G*drho_galileon

  // // Don't forget that parameters are different from barreira
  // // ChiG = c2/pow(*a, 2)*h0*xgal*h*(*dphiprime) - c3/pow(*a, 4)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphiprime) + (*k)*2*h0*pow(xgal, 3)*pow(h, 3)*(*Z) + pow(*k, 2)*2*pow(xgal, 2)*pow(h, 2)*(*dphi)) + c4/pow(*a, 6)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphiprime) + (*k)*15*h0*pow(xgal, 4)*pow(h, 5)*(*Z) + pow(*k, 2)*(12*pow(xgal, 3)*pow(h, 4)*(*dphi) + 1.5*pow(xgal, 4)*pow(h, 4)*(*eta))) - c5/pow(*a, 8)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphiprime) + (*k)*21*h0*pow(xgal, 5)*pow(h, 7)*(*Z) + pow(*k, 2)*(15*pow(xgal, 4)*pow(h, 6)*(dphi) + 3*pow(xgal, 5)*pow(h, 6)*(*eta))) - cG/pow(*a, 4)*(18*h0*xgal*pow(h, 3)*(*dphiprime) + (*k)*6*h0*pow(xgal, 2)*pow(h, 3)*(*Z) + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphi) + pow(xgal, 2)*pow(h, 2)*(*eta)));


  // // No explicit dependance in Z
  // ChiG = (Znogal*coeff + c2*h0*xgal*h*(*dphiprime) - c3/pow(*point, 2)*(18*h0*pow(xgal, 2)*pow(h, 3)*(*dphiprime) + pow(*k, 2)*2*pow(xgal, 2)*pow(h, 2)*(*dphi)) + c4/pow(*point, 4)*(90*h0*pow(xgal, 3)*pow(h, 5)*(*dphiprime) + pow(*k, 2)*(12*pow(xgal, 3)*pow(h, 4)*(*dphi) + 1.5*pow(xgal, 4)*pow(h, 4)*(*eta))) - c5/pow(*point, 6)*(105*h0*pow(xgal, 4)*pow(h, 7)*(*dphiprime) + pow(*k, 2)*(15*pow(xgal, 4)*pow(h, 6)*(*dphi) + 3*pow(xgal, 5)*pow(h, 6)*(*eta))) - cG/pow(*point, 2)*(18*h0*xgal*pow(h, 3)*(*dphiprime) + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphi) + pow(xgal, 2)*pow(h, 2)*(*eta))))/denom;

  //printf("chigal : %.12f \t %.12f \t %.12f \t %.12f\n", *point, h, xgal, ChiG);
  //if(*point >= 9.99999e-7) printf("%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*dgrho), (*eta), (*dphi), (*dphiprime), ChiG);
  if(*point >= 1) printf("%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*eta), h, xgal, ChiG);

  return ChiG;

}


extern "C" double qgal_(double* dgq, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  //printf("OmegaM0 = %f\n OmegaR0 = %.18f\n c0 = %f\n c2 = %f\n c3 = %f\n c4 = %f\n c5 = %f\n cG = %f\n", om, orad, c0, c2, c3, c4, c5, cG);

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : sqrt(om/pow((*point), 2)+orad/pow((*point), 3)+(1-om-orad));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  double alpha = c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    - 2*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    - 2*cG/(3*pow(*point, 2))*pow(xgal, 2)*pow(h, 2);
  double beta = 1./(1-1.5*alpha);
  double qtildeG = c2*(*k)*h0*xgal*h*(*dphi) 
    - c3/pow(*point, 2)*(*k)*(6*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) - 2*pow(xgal, 2)*pow(h, 2)*(*dphiprime)) 
    + c4/pow(*point, 4)*(*k)*(-12*pow(xgal, 3)*pow(h, 4)*(*dphiprime) + 18*h0*pow(xgal, 3)*pow(h, 5)*(*dphi)) 
    - c5/pow(*point, 6)*(*k)*(-15*pow(xgal, 4)*pow(h, 6)*(*dphiprime) + 15*h0*pow(xgal, 4)*pow(h, 7)*(*dphi)) 
    - cG/pow(*point, 2)*(*k)*(-4*xgal*pow(h, 2)*(*dphiprime) + 6*h0*xgal*pow(h, 3)*(*dphi));

  double qG = 0;

  if(1.5*alpha == 1) printf("WARNING : 1/beta_q is zero");
  if(*point >= 9.99999e-7) qG = beta*(qtildeG + 1.5*alpha*(*dgq));

  // double coeff = 1.5/pow(*point, 6)*c4*h0*pow(xgal, 4)*pow(h, 4) - 3./pow(*point, 8)*c5*h0*pow(xgal, 5)*pow(h, 6) - cG/pow(*point, 4)*h0*pow(xgal, 2)*pow(h, 2);
  // double denom = 1-coeff;
  // double qnogal = (*dgq);

  // double qG = 0; // 8*pi*G*dgq_galileon

  // // qG = c2/pow(*a, 2)*h0*xgal*h*(*dphi) - c3/pow(*a, 4)*(*k)*(6*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) - 2*pow(xgal, 2)*pow(h, 2)*(*dphiprime)) + c4/pow(*a, 6)*((*k)*(18*h0*pow(xgal, 3)*pow(h, 5)*(*dphi) - 12*pow(xgal, 3)*pow(h, 4)*(*dphiprime)) + pow(*k, 2)*pow(xgal, 4)*pow(h, 4)*((*sigma) - (*Z))) - c5/pow(*a, 8)*((*k)*(15*h0*pow(xgal, 4)*pow(h, 7)*(*dphi) - 15*pow(xgal, 4)*pow(h, 6)*(*dphiprime)) + 2*pow(*k, 2)*pow(xgal, 5)*pow(h, 6)*((*sigma) - (*Z))) - cG/pow(*a, 4)*((*k)*(6*h0*xgal*pow(h, 3)*(*dphi) - 4*xgal*pow(h, 2)*(*dphiprime)) + 2./3.*pow(*k, 2)*pow(xgal, 2)*pow(h, 2)*((*sigma) - (*Z)));


  // // No exgalplicit dependance in sigma or Z
  // qG = (qnogal*coeff + c2*(*k)*h0*xgal*h*(*dphi) - c3/pow(*point, 2)*(*k)*(6*h0*pow(xgal, 2)*pow(h, 3)*(*dphi) - 2*pow(xgal, 2)*pow(h, 2)*(*dphiprime)) + c4/pow(*point, 4)*(*k)*(18*h0*pow(xgal, 3)*pow(h, 5)*(*dphi) - 12*pow(xgal, 3)*pow(h, 4)*(*dphiprime)) - c5/pow(*point, 6)*(*k)*(15*h0*pow(xgal, 4)*pow(h, 7)*(*dphi) - 15*pow(xgal, 4)*pow(h, 6)*(*dphiprime)) - cG/pow(*point, 2)*(*k)*(6*h0*xgal*pow(h, 3)*(*dphi) - 4*xgal*pow(h, 2)*(*dphiprime)))/denom;

  //if(*point >= 1) printf("%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*eta), h, xgal, qG);

  return qG;

}


extern "C" double Pigal_(double* dgrho, double* dgq, double* dgpi, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  //printf("OmegaM0 = %f\n OmegaR0 = %.18f\n c0 = %f\n c2 = %f\n c3 = %f\n c4 = %f\n c5 = %f\n cG = %f\n", om, orad, c0, c2, c3, c4, c5, cG);

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : sqrt(om/pow((*point), 2)+orad/pow((*point), 3)+(1-om-orad)); // \mathcal{H} = a'/a = a*H
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)+(1-om-orad)); // H = adot/a
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


  double alpha_sigprime = -c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    + 3*c5/pow(*point, 6)*pow(xgal, 4)*pow(h, 5)*xhprime;
  double alpha_sig = c4/pow(*point, 4)*(3*pow(xgal, 4)*pow(h, 4) - 6*pow(xgal, 3)*pow(h, 3)*xhprime) 
    - c5/pow(*point, 6)*(-3*pow(xgal, 5)*pow(h, 5)*hprime + 12*pow(xgal, 5)*pow(h, 6) - 15*pow(xgal, 4)*pow(h, 5)*xhprime) 
    + 2*cG/pow(*point, 2)*xgal*h*xhprime;
  double alpha_phi = -c4/pow(*point, 4)*pow(xgal, 4)*pow(h, 4) 
    + 6*c5/pow(*point, 6)*pow(xgal, 5)*pow(h, 6) 
    + 2*cG/pow(*point, 2)*pow(xgal, 2)*pow(h, 2);
  double beta_pi = 1./(1+0.5*(alpha_phi + 2*alpha_sigprime));
  double PitildeG = c4/pow(*point, 4)*pow(*k, 2)*(4*pow(xgal, 3)*pow(h, 4)*(*dphi) - 6*pow(xgal, 2)*pow(h, 3)*xhprime*(*dphi)) 
    - c5/pow(*point, 6)*pow(*k, 2)*(-12*pow(xgal, 3)*pow(h, 5)*xhprime*(*dphi) + 12*pow(xgal, 4)*pow(h, 6)*(*dphi) - 3*pow(xgal, 4)*pow(h, 5)*hprime*(*dphi)) 
    + 2*cG/pow(*point, 2)*(*k, 2)*h*xhprime*(*dphi);

  double PiG = 0;

  if((alpha_phi + 2*alpha_sigprime) == -1) printf("WARNING : 1/beta_pi is zero");
  if(*point >= 9.99999e-7) PiG = beta*(PitildeG - 0.5*(alpha_phi + 2*alpha_sigprime)*(*dgpi) + 0.5*(alpha_sig-alpha_phi-2*alpha_sigprime)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sigprime+alpha_sig)*(*k)*(*eta));

  // double ksigmaprime = (*k)*h0*h*(*sigma) + 0.5*(*Pi); // equation 41 from arXiv:1208.0600
  // double k2weyl = ((*k)*((*Pi) + (*dgrho)) + 3*h0*h*(*dgq))/(2*(*k)); // Here it's rho_tot and q_tot
  // double coeff_sigmaprime = -c4/pow(*point, 6)*pow(xgal, 4)*pow(h, 4) + 3*c5/pow(*point, 8)*pow(xgal, 4)*pow(h, 5)*xhprime;
  // double coeff_weyl = -c4/pow(*point, 6)*pow(xgal, 4)*pow(h, 4) - c5/pow(*point, 8)*(6*pow(xgal, 4)*pow(h, 5)*xhprime - 6*pow(xgal, 5)*pow(h, 6)) + cG/pow(*point, 4)*pow(xgal, 4)*pow(h, 4);
  // double denom = 1 + 0.5*(2*coeff_sigmaprime + coeff_weyl);

  // double PiG = 0;  

  // PiG = (c4/pow(*point, 4)*((*k)*(3*h0*pow(xgal, 4)*pow(h, 5)*(*sigma) - 6*h0*pow(xgal, 3)*pow(h, 4)*xhprime*(*sigma)) + pow(*k, 2)*(4*pow(xgal, 3)*pow(h, 4)*(*dphi) - 6*pow(xgal, 2)*pow(h, 3)*xhprime*(*dphi))) - c5/pow(*point, 6)*((*k)*(12*h0*pow(xgal, 5)*pow(h, 7)*(*sigma) - 15*h0*pow(xgal, 4)*pow(h, 6)*xhprime*(*sigma) -3*pow(xgal, 5)*pow(h, 6)*hprime*(*sigma)) + pow(*k, 2)*(-12*pow(xgal, 3)*pow(h, 5)*xhprime*(*dphi) + 12*pow(xgal, 4)*pow(h, 6)*(*dphi) - 3*pow(xgal, 4)*pow(h, 5)*hprime*(*dphi))) - cG/pow(*point, 2)*(-2*(*k)*h0*xgal*pow(h, 2)*xhprime*(*sigma) + pow(*k, 2)*(4*pow(xgal, 3)*pow(h, 4)*(*dphi) - 6*pow(xgal, 2)*pow(h, 3)*xhprime*(*dphi))) - coeff_sigmaprime*ksigmaprime - coeff_weyl*k2weyl)/denom;

  return PiG;

}


extern "C" double dphisecond_(double* dgrho, double* eta, double* dphi, double* dphiprime, double* point, double* k){

  //printf("OmegaM0 = %f\n OmegaR0 = %.18f\n c0 = %f\n c2 = %f\n c3 = %f\n c4 = %f\n c5 = %f\n cG = %f\n", om, orad, c0, c2, c3, c4, c5, cG);

  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : sqrt(om/pow((*point), 2)+orad/pow((*point), 3)+(1-om-orad)); // \mathcal{H} = a'/a = a*H
  double hoft = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/pow((*point), 3)+orad/pow((*point), 4)+(1-om-orad)); // H = adot/a
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
  double alpha_gammaprime = 2*c2*h0*h 
    - c3/pow(*point, 2)*(12*h0*pow(h, 2)*xhprime + 12*h0*xgal*pow(h, 2)*hprime) 
    + c4/pow(*point, 4)*(-108*h0*pow(xgal, 2)*pow(h, 5) + 108*h0*xgal*pow(h, 4)*xhprime + 108*h0*pow(xgal, 2)*pow(h, 4)*hprime) 
    - c5/pow(*point, 6)*(-240*h0*pow(xgal, 3)*pow(h, 7) + 180*h0*pow(xgal, 2)*pow(h, 6)*xhprime + 180*h0*pow(xgal, 3)*pow(h, 6)*hprime) 
    - 12*cG/pow(*point, 2)*h0*pow(h, 2)*hprime;
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
  
  if (alpha_gammasecond == 0) printf("WARNING : 1/beta_gamma is zero");
  if(*point >= 9.99999e-7) dphisecond = (alpha_gammaprime*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + 0.5*(alpha_Z - 2*alpha_Zprime)*(*dgrho) - (alpha_Zprime - alpha_Z + 2*alpha_eta)*(*k)*(*eta))/alpha_gammasecond;

  // double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  // double h = (*point >= 9.99999e-7) ? (*point)*(*hx) : sqrt(om/pow((*point), 2)+orad/pow((*point), 3)+(1-om-orad));
  // double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)
  
  // double alpha = c2/6*h*xgal-3*c3*pow(h, 3)*pow(xgal, 2) + 15*c4*pow(h, 5)*pow(xgal, 3) - 17.5*c5*pow(h, 7)*pow(xgal, 4) - 3*cG*pow(h, 3)*xgal;
  // double gamma = c2/3*pow(h, 2)*xgal-c3*pow(h, 4)*pow(xgal, 2) + 2.5*c5*pow(h, 8)*pow(xgal, 4) - 2*cG*pow(h, 4)*xgal;
  // double beta = c2/6*pow(h, 2) -2*c3*pow(h, 4)*xgal + 9*c4*pow(h, 6)*pow(xgal, 2) - 10*c5*pow(h, 8)*pow(xgal, 3) - cG*pow(h, 4);
  // double sigma = 2*h + 2*c3*pow(h, 3)*pow(xgal, 3) - 15*c4*pow(h, 5)*pow(xgal, 4) + 21*c5*pow(h, 7)*pow(xgal, 5) + 6*cG*pow(h, 3)*pow(xgal, 2);
  // double lambda = 3*pow(h, 2) + orad/pow(*point, 4) + c2/2*pow(h, 2)*pow(xgal, 2) - 2*c3*pow(h, 4)*pow(xgal, 3) + 7.5*c4*pow(h, 6)*pow(xgal, 4) - 9*c5*pow(h, 8)*pow(xgal, 5) - cG*pow(h, 4)*pow(xgal, 2);
  // double omega = 2*c3*pow(h, 4)*pow(xgal, 2) - 12*c4*pow(h, 6)*pow(xgal, 3) + 15*c5*pow(h, 8)*pow(xgal, 4) + 4*cG*pow(h, 4)*xgal;

  // double xprime = -xgal+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega); // derivative wrt ln(a)
  // double hprime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
  // double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'

  // double dZ = -h0*h*(*Z) - 0.5/pow(*point, 2)*(*dgrho)/(*k); // This is the derivative of Z wrt tau, not sure about the a^2 coeff though
  // double denom = c2 - 12/pow(*point, 2)*c3*xgal*pow(h, 2) + 54/pow(*point, 4)*c4*pow(xgal, 2)*pow(h, 4) - 60/pow(*point, 6)*c5*pow(xgal, 3)*pow(h, 6) - 6/pow(*point, 2)*cG*pow(h, 2);

  // double dphisecond = 0;

  // if(denom != 0){
  //   dphisecond = ( c2*(2*h0*h*(*dphiprime) + (*k)*xgal*h*(*Z) + pow(*k, 2)*(*dphi)) - c3/pow(*point, 2)*(12*h0*h*xhprime*(*dphiprime) + 12*h0*xgal*pow(h, 2)*hprime*(*dphiprime) + (*k)*(6*h0*pow(xgal, 2)*pow(h, 3)*(*Z) + 4*h0*xgal*pow(h, 2)*xhprime*(*Z) + 2*pow(xgal, 2)*pow(h, 2)*dZ) + pow(*k, 2)*(4*xgal*pow(h, 2)*(*dphi) + 4*h*xhprime*(*dphi))) + c4/pow(*point, 4)*(-108*h0*pow(xgal, 2)*pow(h, 5)*(*dphiprime) + 108*h0*pow(xgal, 2)*pow(h, 4)*hprime*(*dphiprime) + 108*xgal*pow(h, 3)*xhprime*(*dphiprime) + (*k)*(-6*h0*pow(xgal, 3)*pow(h, 5)*(*Z) + 36*h0*pow(xgal, 2)*pow(h, 4)*xhprime*(*Z) + 12*h0*pow(xgal, 3)*pow(h, 4)*hprime*(*Z) + 12*pow(xgal, 3)*pow(h, 4)*dZ) + pow(*k, 2)*(-10*pow(xgal, 2)*pow(h, 4)*(*dphi) - 4*pow(xgal, 3)*pow(h, 4)*(*eta) + 24*xgal*pow(h, 3)*xhprime*(*dphi) + 12*pow(xgal, 2)*pow(h, 3)*hprime*(*dphi) + 6*pow(xgal, 2)*pow(h, 3)*xhprime*(*eta))) - c5/pow(*point, 6)*(h0*(-240*pow(xgal, 3)*pow(h, 7)*(*dphiprime) + 180*pow(xgal, 2)*pow(h, 6)*xhprime*(*dphiprime) + 180*pow(xgal, 3)*pow(h, 6)*hprime*(*dphiprime)) + (*k)*(-45*h0*pow(xgal, 4)*pow(h, 7)*(*Z) + 60*h0*pow(xgal, 3)*pow(h, 6)*xhprime*(*Z) + 30*h0*pow(xgal, 4)*pow(h, 6)*hprime*(*Z) + 15*pow(xgal, 4)*pow(h, 6)*dZ) + pow(*k, 2)*(-36*pow(xgal, 3)*pow(h, 6)*(*dphi) - 12*pow(xgal, 4)*pow(h, 6)*(*eta) + 3*pow(xgal, 4)*pow(h, 5)*hprime*(*eta) + 36*pow(xgal, 2)*pow(h, 5)*xhprime*(*dphi) + 24*pow(xgal, 3)*pow(h, 5)*hprime*(*dphi) + 12*pow(xgal, 3)*pow(h, 5)*xhprime*(*eta))) - cG/pow(*point, 2)*(12*h0*pow(h, 2)*hprime*(*dphiprime) + (*k)*(6*h0*xgal*pow(h, 3)*(*Z) + 4*h0*pow(h, 2)*xhprime*(*Z) + 4*h0*xgal*pow(h, 2)*hprime*(*Z) + 4*xgal*pow(h, 2)*dZ) + pow(*k, 2)*(2*h*(*dphi) + 4*h*hprime*(*dphi) + 2*h*xhprime*(*eta))) ) / denom;
  // } else{
  //   fprintf(stderr, "ERROR : invalid set of parameters\n");
  // }

  //if(*point >= 9.99999e-7) printf("dphisecond : %.10f \t %.10f \t %.10f \t %.10f %.10f \t %.10f \t %.10f \t %.10f %.10f \t %.10f \t %.10f\n", (*point), alpha_gammaprime, alpha_gamma, alpha_Z, alpha_Zprime, alpha_eta, (*dgrho), (*eta), (*dphi), (*dphiprime), dphisecond);
  
  //if(*point >= 9.99999e-7) printf("dphisecond : %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*k), (*dgrho), (*dphi), (*dphiprime), dphisecond);
  //if(*point >= 9.99999e-7) printf("%.10f \t %.10f\n", (*point), dphisecond);
  //printf("dphisecond : %.10f \t %.10f \t %.10f \t %.10f\n", (*point), (*dgrho), (*dphiprime), dphisecond);

  return dphisecond;

}




int test(){

  fflush(stdout);

  double orad = 8.2987687251764e-5;

  // arrays_("galileon.ini", "barreira_background_sample/barreira_galileon_4_5e-6_test.dat", &orad);
  //arrays_("galileon.ini", &orad);
  arrays_("params.ini", &orad);
  FILE* f = fopen("full_integration.txt", "w");

  for(int i = 2; i< intvar.size()-1; i++){
    double point = (intvar[i]+intvar[i+1])/2;
    double* hx = handxofa_(&point);
    fprintf(f, "%.16f ; %.16f ; %.16f\n", point, (*hx), (*(hx+1)));
  }

  gsl_spline_free(spline_h);
  gsl_spline_free(spline_x);
  gsl_interp_accel_free(acc);  

  return 0;

}
