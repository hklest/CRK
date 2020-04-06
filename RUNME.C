#include "CONSTANT_Index.h"
#include "CONSTANT_Efficiency.h"
#include "CONSTANT_Smear.h"
#include "Mirror_Refl.h"
#include "CRK.C"

void RUNME()
{
  Index *n   = new CONSTANT_Index(1.00056);
  Eff   *e1  = new CONSTANT_Efficiency(0.95);
  Eff   *e2  = new CONSTANT_Efficiency(0.95);
  //  Eff   *e3  = new Mirror_Refl(105);
  Smear *s1 = new CONSTANT_Smear(0.001);
  Smear *s2 = new CONSTANT_Smear(0.002);
  Smear *s3 = new CONSTANT_Smear(0.003);
  
  
  CRK *detector = new CRK(n,100);
  detector->AddEff(e1);
  detector->AddEff(e2);
  // detector->AddEff(e3);
  detector->AddSmear(s1);
  detector->AddSmear(s2);
  detector->AddSmear(s3);
  
  detector->Simulate_Something();
}
