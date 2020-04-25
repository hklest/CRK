#ifndef __INDEX_H__
#define __INDEX_H__

class Index
{
public:
  Index() {}
  virtual ~Index() {}
  string GasName;
  // Index of refraction at any wavelength
  virtual double n(double lambda)=0;
};
	
#endif /* __INDEX_H__ */

