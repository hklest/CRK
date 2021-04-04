#ifndef __PULSE_SMEAR_H__
#define __PULSE_SMEAR_H__

#include "Smear.h"
#include "TRandom.h"

class Pulse_Smear : public Smear
{
public:
  Pulse_Smear(double SMEAR);
  virtual ~Pulse_Smear() {}

  // smears the input angle in a gaussian manner, independent of wavelength;
  double smr(double ThetaCherenkov, double lambda);

protected:
  double smear;
  TRandom Randy;
};
	
#endif /* __PULSE_SMEAR_H__ */
