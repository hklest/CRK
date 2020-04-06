#ifndef __CONSTANT_INDEX_H__
#define __CONSTANT_INDEX_H__

#include "Index.h"

class CONSTANT_Index : public Index
{
public:
  CONSTANT_Index(double N) {index = N;}
  virtual ~CONSTANT_Index() {}

  // Index of refraction at all wavelengths
  double n(double lambda) {return index;}

protected:
  double index;
};
	
#endif /* __CONSTANT_INDEX_H__ */
