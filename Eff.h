#ifndef __EFF_H__
#define __EFF_H__

class Eff
{
public:
  Eff() {}
  virtual ~Eff() {}

  virtual double e(double lambda)=0; // Efficiency at any wavelength
};
	
#endif /* __EFF_H__ */
