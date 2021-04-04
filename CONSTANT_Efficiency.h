#ifndef __CONSTANT_EFFICIENCY_H__
#define __CONSTANT_EFFICIENCY_H__

#include "Eff.h"

class CONSTANT_Efficiency : public Eff
{
public:
  CONSTANT_Efficiency(double EFF) {eff = EFF;}
  virtual ~CONSTANT_Efficiency() {}

  // This efficiency has no wavelength dependence
  double e(double lambda) {return eff;}

protected:
  double eff;
};
	
#endif /* __CONSTANT_EFFICIENCY_H__ */
