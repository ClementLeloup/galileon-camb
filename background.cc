// Modified by Clement Leloup 22-07-2016 :
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
#include <iostream>
#include <iomanip>
#include <vector>
#include "fstream"
#include <algorithm>
#include <math.h>
#include "string.h"
#include <stdio.h>

using namespace std;

std::vector<std::string> readParamsFile(std::string filename){

  std::vector<std::string> params(12, "");
  std::string* param = 0;
  
  std::ifstream file;
  std::string line;
  const char* com = "#";
  file.open(filename.c_str());
  if(file.is_open()){
    while(!file.eof()){
      getline(file,line);      
      if(line != ""){	
	if(strncmp(line.c_str(), com, strlen(com)) != 0){
	  char input[strlen(line.c_str())];
	  strcpy(input,line.c_str());
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
	  } else if(string(token) == "om"){
	    param = &params[8];
	  } else if(string(token) == "orad"){
	    param = &params[9];
	  }  else if(string(token) == "name"){
	    param = &params[10];
	  } else if(string(token) == "barreira"){
	    param = &params[11];
	  }
	  token = strtok(NULL, " =");
	  if(token!=NULL){
	    
	    std::string str = token;
	    if(param != NULL){
	      *param = str;
	      param = 0;
	    }
	  }
	}
      }
    }
  }
  
  //std::cout << "c2 = " << params[0] << "\n" << "c3 = " << params[1] << "\n" << "c4 = " << params[2] << "\n" << "c5 = " << params[3] << "\n" << "c0 = " << params[4] << "\n" << "cG = " << params[5] << "\n" << "ratio_rho = " << params[6] << "\n" << "age = " << params[7] << "\n" << "om = " << params[8] << "\n" << "orad = " << params[9] << "\n" << "name = " << params[10] << "\n" << "barreira = " << params[11] << "\n";

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
  double om = in_params[0];
  double c2 = in_params[1];
  double c3 = in_params[2];
  double c4 = in_params[3];
  double c5 = in_params[4];
  double cG = in_params[5];
  double c0 = in_params[6];
  double orad = in_params[7];
  bool useacoord = (int)in_params[8];
  double *OmegaPiOverOrad = NULL;
  if(useacoord) OmegaPiOverOrad = const_cast<double*>(&in_params[9]);
  if (useacoord) std::cout << "Using a instead of z" << std::endl; 
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
  if(fabs(OmegaM - OmTest)>1e-5){
    std::cout<<z<<" "<<om<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<cG<<" "<<c0<<" "<<OmTest<<" "<<OmegaM<<" "<<OmegaP<<" "<<fabs(OmegaM - OmTest)<<std::endl;
    return 4;
  }
  if(OmegaP<0) return 5;
  // if ( useacoord && OmegaP>0.1 && abs(1/z-1000.0)<1e-8  ){
  //    std::cout<<"Warning : OmegaP is more than 10% of Omega_tot at z=1000 : OmegaP="<<OmegaP<<"\n\t with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<std::endl;
  //    std::cout<<"\tBeware of theoretical validity of some computations!"<<std::endl;
  //  }

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
*/
int calcHubbleGalileon(std::vector<double>& hubble, std::vector<double>& x, const std::vector<double>& zcoord, double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0, double coord){

  double params[10];
  params[0] = om; params[1] = c2; params[2] = c3;
  params[3] = c4; params[4] = c5; params[5] = cG;
  params[6] = c0; params[7] = orad; params[8] = coord; params[9] = 0; 
  bool useacoord = (int)params[8];

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
      zcurrtarg = zcoord[i];
      while(z > zcurrtarg){
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
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
	std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
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
      zcurrtarg = zcoord[i];
      while(z < zcurrtarg){
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
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
	std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
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
void calcHubbleTracker(std::vector<double>& hubble, std::vector<double>& x, const std::vector<double>& zcoord, double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0, double coord){
 
  bool useacoord = (int)coord;
  double op = (c2/2.0 - 6.0*c3 + 22.5*c4 - 21.0*c5 - 9*cG )/3.0; //Omega_phi

  if(useacoord){
    int nstep = min(hubble.size(),x.size());
    for(int i = 0; i < nstep; ++i){
      hubble[i] = sqrt(0.5*(om/pow(zcoord[i],3)+orad/pow(zcoord[i],4)-3*cG)+sqrt(op/9+3*cG+0.25*pow(3*cG-om/pow(zcoord[i],3)-orad/pow(zcoord[i],4),2))); // Analytical solution for H
      x[i] = zcoord[i]/pow(hubble[i],2);
    }
  } else{
    int nstep = min(hubble.size(),x.size());
    for(int i = 0; i < nstep; ++i){
      hubble[i] = sqrt(0.5*(om*pow(1+zcoord[i],3)+orad*pow(1+zcoord[i],4)-3*cG)+sqrt(op+3*cG+0.25*pow(3*cG-om*pow(1+zcoord[i],3)-orad*pow(1+zcoord[i],4),2)));
      x[i] = -1/((1+zcoord[i])*pow(hubble[i],2));
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
int ageOfUniverse(std::vector<double> &hubble, std::vector<double> &x, double &age, std::vector<double> &acoord, double cond[4], double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0){

  const double Gyr = 3.1556926e16; //Conversion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Conversion factor from Mpc to km
  const double h0 = 71; // Present-day Hubble constant (WMAP7)

  double params[10];
  params[0] = om; params[1] = c2; params[2] = c3;
  params[3] = c4; params[4] = c5; params[5] = cG;
  params[6] = c0; params[7] = orad; params[8]= 1; params[9] = 0.0;	

  // Vector y = ( dh/dz, dx/dz, dy/dz )
  double y[3] = {cond[0], cond[1], cond[2]};  //Inital value of integral
  hubble[0] = cond[0];
  x[0] = cond[1];

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-14, 1e-14, 1, 1);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (3);

  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;

  double a = acoord[1];
  int nstep = acoord.size();
  double h = 1e-6; //Initial step guess
  double acurrtarg;
  double prec_age = 1/(acoord[1]*cond[0]); // initial value of 1/(a*Hbar)
  int st;
  age = cond[3];

  for (int i = 2; i < nstep; i++) {
    acurrtarg = acoord[i];
    while(a < acurrtarg){
      st = gsl_odeiv_evolve_apply(e, c, s, &sys, &a, acurrtarg, &h, y);
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
      //std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
      return 8;
    }
    age += 0.5*(1.0/(acoord[i]*y[0])+prec_age)*(acoord[i]-acoord[i-1]); // trapezoidal method
    prec_age = 1/(acoord[i]*y[0]);
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

int initCond(std::vector<double> &hubble, std::vector<double> &x, double &xi, std::vector<double> &acoord, double H[2], double ratio_rho, double age, double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0){

  const double Gyr = 3.1556926e16; //Convertion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Convertion factor from Mpc to km
  int status = 0;
  double coeff[6]; //array of coefficients of the polynomial
  coeff[0] = -3*ratio_rho*om/pow(acoord[1], 3); coeff[1] = 0; coeff[2] = c2/2*pow(H[1], 2)-9*cG*pow(H[1], 4);
  coeff[3] = -6*c3*pow(H[1], 4); coeff[4] = 45/2*c4*pow(H[1], 6); coeff[5] = -21*c5*pow(H[1], 8);

  size_t n = sizeof(coeff)/sizeof(double);
  double cond[4]; // initial conditions
  cond[0] = H[1]; cond[2] = 0; cond[3] = Mpc/Gyr*(1/H[1]-pow(acoord[1], 2)/2*(1/(acoord[0]*H[0])-1/(acoord[1]*H[1]))/(acoord[0]-acoord[1])); // Initial value of the age of the Universe, extrapolating linearly the integral
  //std::cout << "Initial value : " << cond[3] << endl;

  AgeofU ini_vec_age = AgeofU(0, 1000);
  std::vector<AgeofU> aou(5, ini_vec_age);

  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n);
  double z[2*(n-1)]; // Array of complex numbers ordered like {Re[0],Im[0],Re[1],Im[1],...}

  status = gsl_poly_complex_solve(coeff, n, w, z); // Solve coeff[0]+coeff[1]*x+...+coeff[5]*x^5=0 and store sol in z

  int nrealsol = 0; // Number of real solutions

  for(int j = 0; j<n-1; j++){
    std::cout << "Root number " << j << " is " << z[2*j] << " + " << z[2*j+1] << "i" << endl;
    if(z[2*j+1]==0){
      cond[1] = z[2*j]/acoord[1]; // Initial condition for d_phi/da
      aou[j].i = j;
      int status2 = ageOfUniverse(hubble, x, aou[j].a, acoord, cond, om, orad, c2, c3, c4, c5, cG, c0);
      aou[j].a -= age;
      if(status2 != 0) return status2;
    }
  }
  
  for(int j = 0; j<n-1; j++){
    if(aou[j].a != 1000){
      nrealsol++;
      std::cout << "Root number " << aou[j].i << " gives age of universe equal to " << aou[j].a + age << endl; 
    }
  }

  sort(aou.begin(), aou.end(), less_than_key()); // Sort by comparing the calculated age of universe to the one provided

  for(int j = 0; j<n-2; j++){
    aou[j].a += age;
  }

  xi = z[2*(aou[0].i)];
  if(nrealsol > 1){
    cond[1] = xi/acoord[1];
    ageOfUniverse(hubble, x, aou[0].a, acoord, cond, om, orad, c2, c3, c4, c5, cG, c0);
  }

  std::cout << "Initial condition for x is " << xi << ", it is the root number " << aou[0].i << " and it gives an age of the universe equal to " << aou[0].a << endl;

  return status;

} 


void SetAcoord(std::vector<double>& acoord, double amin ){

  double da = 1e-5; //!< Step in a
  double dafactor = 0.01; //!< Factor to apply to da if a<astep
  double dafactor2 = 0.001; //!< Factor to apply to da if a<astep2
  double dafactor3 = 0.0001; //!< Factor to apply to da if a<astep3
  double dafactor4 = 0.00001; //!< Factor to apply to da if a<astep4
  double dafactor5 = 0.000001; //!< Factor to apply to da if a<astep5
  double astep = 0.02; //!< Step in a where to change da
  double astep2 = 0.0005; //!< Step in a where to change da
  double astep3 = 3e-5; //!< Step in a where to change da
  double astep4 = 5e-6; //!< Step in a where to change da
  double astep5 = 5e-7; //!< Step in a where to change da  

  if (acoord.size() != 0) acoord.clear();
  int nstep = (int)floor((1.0-astep)/da);
  int nstep1 = nstep + (int)floor((astep-astep2)/(dafactor*da));
  int nstep2 = nstep1 + (int)floor((astep2-astep3)/(dafactor2*da));
  int nstep3 = nstep2 + (int)floor((astep3-astep4)/(dafactor3*da));
  int nstep4 = nstep3 + (int)floor((astep4-astep5)/(dafactor4*da));
  int nsteptot = nstep4 + (int)floor(astep5/(dafactor5*da));

  //std::cout << nstep << " " << nstep1 << " " << nstep2 << " " << nstep3 << " " << nstep4 << " " << nsteptot << endl;

  acoord.push_back(1.0); // Initialization of acoord

  for(int i = 1; i< nsteptot; ++i){
    if(acoord[i-1]>amin){
      if(i<nstep){
	acoord.push_back(1.0-i*da);
      } else if(nstep<=i && i<nstep1){
	acoord.push_back(astep-(i-nstep)*da*dafactor);
      } else if(nstep1<=i && i<nstep2){
	acoord.push_back(astep2-(i-nstep1)*da*dafactor2);
      } else if(nstep2<=i && i<nstep3){
	acoord.push_back(astep3-(i-nstep2)*da*dafactor3);
      } else if(nstep3<=i && i<nstep4){
	acoord.push_back(astep4-(i-nstep3)*da*dafactor4);
      } else if(nstep4<=i && i<nsteptot){
	acoord.push_back(astep5-(i-nstep4)*da*dafactor5);
      }
    } else{
      break;
    }
  }

  acoord.pop_back();
  //for (unsigned int k=0; k<acoord.size(); ++k)
  // std::cout<<acoord[k]<<std::endl;
}





int main(){

  std::vector<std::string> params = readParamsFile("galileon.ini");

  double c2 = atof(params[0].c_str());
  double c3 = atof(params[1].c_str());
  double c4 = atof(params[2].c_str());
  double c5 = atof(params[3].c_str());
  double c0 = atof(params[4].c_str());
  double cG = atof(params[5].c_str());
  double ratio_rho = atof(params[6].c_str());
  double age = atof(params[7].c_str());
  double om = atof(params[8].c_str());
  double orad = atof(params[9].c_str());
  char* name = const_cast<char*>(params[10].c_str());
  double barreira = atof(params[11].c_str());

  //double barreira = 0;
  double coord = 0;
  int status = 0;

  // File to store h and x as a function of z
  ofstream f;
  //f.open("5params_all_2015_combined_noH0prior.dat", ios::out );
  
  /*
  // Set of parameters (Best fit 2015 BAO+Planck+Lya from Jeremy)
  double h = 0.762;
  double OmegaM0 = 0.261175;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -5.25803;
  double c3 = -1.13137;
  double cG = 0;
  double c4 = 1.0/27*(30*(om+orad-1)-9*c2+24*c3-6*cG);
  double c5 = 1.0/3*(4*(om+orad-1)-c2+2*c3-2*cG);
  double c0 = 0;
  char* name = "tracker_test.dat";
  */
  

  /*
  // Best fit Unc 2015 All+JLA+noH0prior from Jeremy
  double h = 0.736;
  double OmegaM0 = 0.275;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -4.145;
  double c3 = -1.545;
  double cG = 0;
  double c4 = -0.776;
  double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
  double c0 = 0;
  char* name = "5params_all_2015_combined_noH0prior_bis_bis.dat";
  */
  

  /*
  // Best fit + cG +1 sigma 2015 All+JLA+noH0prior from Jeremy
  double h = 0.727;
  double OmegaM0 = 0.280;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -3.434;
  double c3 = -1.062;
  double cG = 0.235;
  double c4 = -0.610;
  double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
  double c0 = 0;
  char* name = "6params_all_2015_combined_noH0prior.dat";
  */

  
  // Galileon parameters from Barreira
  // double c3 = 12.8;
  // double c4 = -1.7;
  // double c5 = 1.0;
  // double cG = 0;
  // double ratio_rho = 1e-4;
  // double c2 = -27.0;
  // double age = 13.978;
  // double c0 = 0;
  // double orad = 8e-5;
  // double om = 0.265;
  // char* name = "barreira_galileon_1_1e-4.dat";
  // barreira = 1;
  


  if(barreira)
  {
    cout << "barreira" << endl;
    
    f.open(name, ios::out );

    double apremin = 1e-7;
    double amin = 9.99999e-7;
    std::vector<double> acoord;
    acoord.push_back(apremin);
    SetAcoord(acoord, amin);

    sort(acoord.begin(), acoord.end());

    double xi = 0;
    double H[] = {sqrt(om/pow(apremin, 3)+orad/pow(apremin, 4)), sqrt(om/pow(amin, 3)+orad/pow(amin, 4))};

    std::vector<double> hubble(acoord.size()-1, 999999);
    std::vector<double> x(acoord.size()-1, 999999);

    int status = initCond(hubble, x, xi, acoord, H, ratio_rho, age, om, orad, c2, c3, c4, c5, cG, c0);

    std::cout << "OmegaM0 = " << om << std::endl;
    std::cout << "OmegaR0 = " << orad << std::endl;
    std::cout << "c0 = " << c0 << std::endl;
    std::cout << "c2 = " << c2 << std::endl;
    std::cout << "c3 = " << c3 << std::endl;
    std::cout << "c4 = " << c4 << std::endl;
    std::cout << "c5 = " << c5 << std::endl;
    std::cout << "cG = " << cG << std::endl;
    std::cout << "\n\nstatus = " << status << "\n" << std::endl;

    if(status!=0) return 0;


    for(int i=0; i<hubble.size(); i++)
      {
        double alpha = c2/6*hubble[i]*acoord[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*acoord[i+1]*x[i];
	double gamma = c2/3*pow(hubble[i], 2)*acoord[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];
	double beta = c2/6*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*acoord[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
	double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2);
	double lambda = 3*pow(hubble[i], 2) + orad/pow(acoord[i+1], 4) + c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
	double omega = 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];
	
	double x_prime = -acoord[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
	double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
	
	double rho = c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
	double p = c2/2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*(h_prime*acoord[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*acoord[i+1]*x[i]+2*hubble[i]*acoord[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*acoord[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2));
	
	double OmegaM = om/(pow(hubble[i], 2)*pow(acoord[i], 3));
	double OmegaR = orad/(pow(hubble[i], 2)*pow(acoord[i], 4));
	double OmegaP = rho*3*pow(hubble[i], 2);
	
	double w = p/rho;
	double hubble_LCDM = sqrt(om/pow(acoord[i+1], 3)+orad/pow(acoord[i+1], 4)+(1-om-orad));

	f << std::setprecision(12) << acoord[i+1] << " ; " << hubble[i] << " ; " << hubble_LCDM << " ; " << w << " ; " << x[i] << " ; " << OmegaP << " ; " << OmegaM << " ; " << OmegaR << std::endl;
      }

  } else if(coord)
  {
    cout << "acoord" << endl;

    f.open(name, ios::out );

    double amin = 0.000908;
    std::vector<double> acoord;
    SetAcoord(acoord, amin);

    sort(acoord.begin(), acoord.end(), std::greater<double>());

    std::vector<double> hubble(acoord.size(), 999999);
    std::vector<double> x(acoord.size(), 999999);

    if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
    {
      status = calcHubbleGalileon(hubble, x, acoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
    }
    else {
      calcHubbleTracker(hubble, x, acoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
    }

    std::cout << "OmegaM0 = " << om << std::endl;
    std::cout << "OmegaR0 = " << orad << std::endl;
    std::cout << "c0 = " << c0 << std::endl;
    std::cout << "c2 = " << c2 << std::endl;
    std::cout << "c3 = " << c3 << std::endl;
    std::cout << "c4 = " << c4 << std::endl;
    std::cout << "c5 = " << c5 << std::endl;
    std::cout << "cG = " << cG << std::endl;
    std::cout << "\n\nstatus = " << status << "\n" << std::endl;

    if(status!=0) return 0;

    for(int i=0; i<hubble.size(); i++)
    {
      double OmegaM = om/(pow(hubble[i], 2)*pow(acoord[i], 3));
      double OmegaR = orad/(pow(hubble[i], 2)*pow(acoord[i], 4));
      double OmegaP = (c2/2.0*pow(acoord[i]*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(acoord[i]*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(acoord[i]*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

      double hubble_LCDM = sqrt(om/pow(acoord[i+1], 3)+orad/pow(acoord[i+1], 4)+(1-om-orad));

      f << 1/acoord[i]-1 << " ; " << hubble[i] << " ; " << hubble_LCDM << " ; " << x[i] << " ; " << OmegaP << " ; " << OmegaM << " ; " << OmegaR << std::endl;
    }

  } else
  {
    cout << "zcoord" << endl;

    f.open(name, ios::out );

    double zmax = 1100;
    double amin = 1/(1+zmax);
    std::vector<double> zcoord;
    SetAcoord(zcoord, amin);

    sort(zcoord.begin(), zcoord.end(), std::greater<double>());

    for(int i=0;i<zcoord.size();i++) {zcoord[i]=1/(1+zcoord[i]);}
    std::vector<double> hubble(zcoord.size(), 999999);
    std::vector<double> x(zcoord.size(), 999999);

    double op = 1.5*c2 - 18*c3 + 67.5*c4 - 63*c5 - 27*cG;

    if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
    {
      status = calcHubbleGalileon(hubble, x, zcoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
    }
    else {
      std::cout << "Tracker " << std::endl;
      calcHubbleTracker(hubble, x, zcoord, om, orad, c2, c3, c4, c5, cG, c0, coord);
    }

    std::cout << "OmegaM0 = " << om << std::endl;
    std::cout << "OmegaR0 = " << orad << std::endl;
    std::cout << "OmegaP0 = " << op << std::endl;
    std::cout << "c0 = " << c0 << std::endl;
    std::cout << "c2 = " << c2 << std::endl;
    std::cout << "c3 = " << c3 << std::endl;
    std::cout << "c4 = " << c4 << std::endl;
    std::cout << "c5 = " << c5 << std::endl;
    std::cout << "cG = " << cG << std::endl;
    std::cout << "\n\nstatus = " << status << "\n" << std::endl;

    if(status!=0) return 0;

    for(int i=0; i<zcoord.size(); i++)
    {
      double OmegaM = om*pow(1+zcoord[i], 3)/pow(hubble[i], 2);
      double OmegaR = orad*pow(1+zcoord[i], 4)/pow(hubble[i], 2);
      double OmegaP = (c2/2.0*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(-(1+zcoord[i])*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(-(1+zcoord[i])*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

      double hubble_LCDM = sqrt(om*pow(1+zcoord[i+1], 3)+orad*pow(1+zcoord[i+1], 4)+(1-om-orad));

      f << std::setprecision(12) << zcoord[i] << " ; " << hubble[i] << " ; " << hubble_LCDM << " ; " << x[i] << " ; " << OmegaP << " ; " << OmegaM << " ; " << OmegaR << std::endl;
    }

  }

  return 0;

}


