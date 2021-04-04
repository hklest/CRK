#ifndef __EMISSION_SMEAR_H__
#define __EMISSION_SMEAR_H__

#include "Smear.h"
#include "TRandom.h"

class Emission_Smear : public Smear
{
public:
  Emission_Smear(double SMEAR);
  virtual ~Emission_Smear() {}

  // smears the input angle in a gaussian manner, independent of wavelength;
  double smr(double ThetaCherenkov, double lambda);

protected:
  double smear;
  TRandom Randy;
};
	
#endif /* __EMISSION_SMEAR_H__ */
