#include "Index.h"
#include "Eff.h"
#include "Smear.h"
#include "CRK.h"
#include <cmath>

using namespace std;
double lambda;

CRK::CRK(Index *index, double RadiatorLength)
{
  n = index;
  L = RadiatorLength;
}

void CRK::AddEff(Eff *e)
{
  effs.push_back(e);
}

void CRK::AddSmear(Smear *s)
{
  smears.push_back(s);
}

void CRK::Simulate_Something() //Simulate_Something provides outputs of efficiencies, theta_c, Npe, etc.
{
  // initialize integration parameters
  double lambdamax = 180.0; // maximum wavelength of cherenkov light considered
  double lambdamin = 100.0; // minimum wavelength of cherenkov light considered (once all efficiencies are fully implemented, this should be able to be set to 0)
  double stepsize = .1; // step size for numerical integration
  //
  
  // initializing other variables
  double Weighted_Efficiency = 0.0; 
  double normalization = 0.0;
  double i = 0.0;
  double Mult_E;
  double dNdxIntegral = 0;
  //

  //Find normalization constant
  for (i = lambdamin; i <= lambdamax; i=i+stepsize) 
    {      
      normalization += pow(i,-2.0);
      //cout <<"normalization constant is: "<< normalization << endl;
    }
  //

  
  //Loop over wavelengths between lambdamin nm and lambdamax nm to find quantities at each wavelength in steps of stepsize nm
  for (lambda = lambdamin; lambda <= lambdamax; lambda=lambda+stepsize)   
    {


      //Efficiency calculations
      Mult_E = (effs[0]->e(lambda))*(effs[1]->e(lambda)); // Multiplying the available efficiencies by each other, this needs to change when new efficiences are added by multiplying by (eff[2]->e(lambda)
      Weighted_Efficiency += Mult_E*pow(lambda,-2)/normalization; // metric of the full survival/measurement rate of all photons
      dNdxIntegral += Mult_E*stepsize*pow(lambda,-2); // integral for use in dN/dx calculation


      // Index of refraction/Theta_C calculations
      double mPi = 0.13957; // GeV                                                                                                                                                                         
      double mKa = 0.49368;                                                                                                                                                                                
      double mPr = 0.93827;                                                                                                                                                                                
      double me = .000511;
      double PID [4] = {mPi, mKa, mPr, me};
      double p = 30; // Momentum of particle in GeV
      double p_m = p/PID[0]; //insert desired particle here, 0 = Pion, 1 = Kaon, 2 = Proton, 3 = Electron
      double beta = sqrt(p_m*p_m/(1 + p_m*p_m));
      
      cout <<"Index of Refraction at "<< lambda << " nm = " << n->n(lambda) << endl;
      double ncount = n->n(lambda);
      double ThetaCth = acos(1/(ncount*beta));
      cout << "Theta_C at " << lambda << " nm = " << ThetaCth << endl;

 
      }
  //
  
     



  
  
  //double nmin = *n(62);
  //double ThetaCthmin = acos(1/(nmin*beta));
  //cout << "Theta_C at lambdamin = " << ThetaCthmin << endl;

  //double nmax = *n(120);
  //double ThetaCthmax = acos(1/(nmax*beta));
  //cout << "Theta_C at lambdamax = " << ThetaCthmin << endl;



  //Calculate Npe
  double ThetaCth = .03; // 30 mrad
  double AlphaEM = 1/137.035999;
  double dNdx = 2*M_PI*AlphaEM*pow(sin(ThetaCth),2)*dNdxIntegral; // number of photons produced AND measured per unit length
  double Npe;
  Npe = dNdx*L*10000000; // large number is there to convert from nm to cm
  //cout << "Sin^2(thetaC) = " << pow(sin(ThetaCth),2);
  //cout << "dN/dx = " << dNdx;
  //


  // Insert smears
  double ThetaCms = ThetaCth;
  for (int k=0; k<smears.size(); k++)
    {
      ThetaCms = smears[k]->smr(ThetaCms,500);
    }
  //



  //Print Values
  cout << "Number of photons detected between "<< lambdamin << " nm and " << lambdamax <<" nm = " << Npe << endl;
  cout << "Total Photon Efficiency = " << Weighted_Efficiency << endl;
  //cout << "ThetaC(theory) = " << ThetaCth << endl;
  //cout << "ThetaC(measured) = " << ThetaCms << endl;
  //
}
