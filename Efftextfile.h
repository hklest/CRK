#ifndef __EFFTEXTFILE_H__
#define __EFFTEXTFILE_H__

#include "Eff.h"
#include <vector>       // std::vector
#include <utility>
#include <string>


// Efftextfile takes a CSV "PhotocathodeEff.txt" (as produced by DataThief) and produces the quantum efficiency of a photocathode (Default CsI) at a specific wavelength in nm. This is accomplished by splitting the columns of the CSV into two vectors, one containing efficiences (CsIEffs) and one containing wavelengths (CsILambdas). The indices of the vectors have a one-to-one correspondence, so finding the index of the nearest wavelength to the one entered produces also the index of the closest efficiency.

class Efftextfile : public Eff
{
public:
  Efftextfile(std::string filename, int sparse = 1); 
  virtual ~Efftextfile() {}

  // This efficiency has a wavelength dependence
  double e(double lambda);
  
  double LMin(); // Values chosen to avoid singularities in Sellmeier coefficients. 
 
  double LMax();
 
  
protected:
  std::vector<double> Lambdas , Effs; 
};
	
#endif /* __EFFTEXTFILE_H__ */
