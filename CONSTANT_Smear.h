#ifndef __CONSTANT_SMEAR_H__
#define __CONSTANT_SMEAR_H__

#include "Smear.h"
#include "TRandom.h"

class CONSTANT_Smear : public Smear
{
public:
  CONSTANT_Smear(double SMEAR) {smear = SMEAR;}
  virtual ~CONSTANT_Smear() {}

  // smears the input angle in a gaussian manner, independent of wavelength;
  double smr(double ThetaCherenkov, double lambda, double p, double B, double r_in,double r_out,double ringphi,double rapidity) { return ThetaCherenkov + Randy.Gaus(0.0,smear); }

protected:
  double smear;
  TRandom Randy;
};
	
#endif /* __CONSTANT_SMEAR_H__ */
