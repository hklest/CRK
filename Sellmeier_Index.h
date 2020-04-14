#ifndef __SELLMEIER_INDEX_H__
#define __SELLMEIER_INDEX_H__

#include "Index.h"
#include <cmath>



class Sellmeier_Index : public Index
{
public:
  Sellmeier_Index(double lambda) {}
  virtual ~Sellmeier_Index() {}

  // Index of refraction at all wavelengths as defined by Sellmeier equation, default is CF4 with A = .12 and Lambda0 = 2.61154E-4 with a cutoff wavelength of 62 nm
  double n(double lambda)
  {
    double A = .120; //Average of https://cds.cern.ch/record/600182/files/ep-2002-099 and Nucl. Instrum. Methods Phys. Res., A : 292(1990)593
    double Lambda0 = 61.81;
    if (lambda < 62)
      {
	cout << "Wavelength entered is below UV cutoff" << endl;
      }
    else
      {
	index = 1+(A*pow(10.0,-6.0))/((pow(Lambda0,-2.0) - pow(lambda,-2)));
	//cout << index << endl;
      }
    return index;
  }

protected:
  double index;
};
	
#endif /* __SELLMEIER_INDEX_H__ */
