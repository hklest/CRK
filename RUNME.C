#include "CONSTANT_Index.h"
#include "Sellmeier_Index.h"
#include "CONSTANT_Efficiency.h"
#include "CONSTANT_Smear.h"
#include "Photocathode_Eff.h"
#include "Mirror_Eff.h"
#include "Geom_Surface_Eff.h"
//#include "WL_Dep.h"
#include "CRK.C"


void RUNME()
{
  
  Index *n   = new Sellmeier_Index(1);
  Eff   *e1  = new Geom_Surface_Eff(1);
  Eff   *e2  = new Photocathode_Eff(1);
  Eff   *e3  = new Mirror_Eff(1);
  Smear *s1 = new CONSTANT_Smear(0.001);
  Smear *s2 = new CONSTANT_Smear(0.002);
  Smear *s3 = new CONSTANT_Smear(0.003);
  
  
  CRK *detector = new CRK(n,100); // index of refraction and length of radiator in cm
  detector->AddEff(e1);
  detector->AddEff(e2);
  detector->AddEff(e3);
  detector->AddSmear(s1);
  detector->AddSmear(s2);
  detector->AddSmear(s3);
  //  detector->Normalize();
  //detector->dNdx();
  //detector->AverageThetaC();
  detector->Simulate_Something();
}
