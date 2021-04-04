#ifndef __PAD_SMEAR_H__
#define __PAD_SMEAR_H__

#include "Smear.h"
#include "TRandom.h"

class Pad_Smear : public Smear
{
public:
  Pad_Smear(double SMEAR);
  virtual ~Pad_Smear() {}

  // smears the input angle in a gaussian manner, independent of wavelength;
  double smr(double ThetaCherenkov, double lambda);

protected:
  double smear;
  TRandom Randy;
};
	
#endif /* __PAD_SMEAR_H__ */
