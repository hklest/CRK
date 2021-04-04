#ifndef __EFF_H__
#define __EFF_H__

#include <string>

class TH1;


class Eff
{
 public:
  Eff();
  virtual ~Eff() {}

  virtual double e(double lambda)=0; // Efficiency at any wavelength
  virtual double LMin(){return 100;}
  virtual double LMax(){return 400;}
  TH1* EffPlot();
  std::string EffName()
  {
    return Name;
  }

  
 protected:
  TH1* WavelengthDependence;
  std::string Name;
  void FillWavelengthDependence();

};
	
#endif /* __EFF_H__ */
