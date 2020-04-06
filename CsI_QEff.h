#ifndef __CSI_QEFF_H__
#define __CSI_QEFF_H__

#include "Eff.h"

class CsI_QEff : public Eff
{
public:
  CsI_QEff(double EFF) {eff = EFF;}
  virtual ~CsI_QEff() {}

  // This efficiency has a wavelength dependence
  double e(double lambda)
  {
    return eff;
  }

protected:
  double eff;
};
	
#endif /* __CSI_QEFF_H__ */
