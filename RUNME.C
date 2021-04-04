
#include "CONSTANT_Smear.h"
#include "Momentum_Dep_Smear.h"
#include "CONSTANT_Efficiency.h"
#include "Eff.h"
#include "Efftextfile.h"
#include "CRK.h"
#include "Eff.C"
#include "Efftextfile.C"
#include "CRK.C"



void RUNME()
{  
  Index *n   = new Sellmeier_Index(1);
  Eff   *e1  = new Efftextfile("SiPM_Eff.txt",1);
  Eff   *e2  = new Efftextfile("Mirror_Eff.txt",1);
  Eff   *e3  = new Efftextfile("SiPM_Geom_Eff.txt",1);
  
  // Eff   *e5 = new CONSTANT_Efficiency(.88);//Glass/Mesh optical transparency 


  //   Eff   *e4  = new Efftextfile("CF4Transmission_Eff.txt",1);
  // Eff *e6 = new CONSTANT_Efficiency(.66);//Collection Efficiency
  // Eff *e7 = new CONSTANT_Efficiency(.98);//single photon detection Efficiency
  Smear *s1 = new CONSTANT_Smear(.003); // Quadrature sum of all smears, 4.3 mrad is DELPHI error per photon
  Smear *s2 = new Momentum_Dep_Smear();
  CRK *detector = new CRK(n,100); // index of refraction and length of radiator in cm
  detector->AddEff(e1);
  detector->AddEff(e2);
  detector->AddEff(e3);
  // detector->AddEff(e4);
  // detector->AddEff(e5);
  // detector->AddEff(e6);
  // detector->AddEff(e7);
  detector->AddSmear(s1);
  detector->AddSmear(s2);
 
  detector->Simulate_Momentum(); // Simulate N-Sigma separation and Theta_C distribution as a function of momenta
  // detector->Simulate_Alpha(); // Simulate N-Sigma and Theta_C distribution as a function of tracking angle error
  //  detector->MultEffPlot(); // Plot Efficiencies and spectrum weighted efficiencies as a function of wavelength
}
