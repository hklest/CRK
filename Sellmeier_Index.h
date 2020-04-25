#ifndef __SELLMEIER_INDEX_H__
#define __SELLMEIER_INDEX_H__

#include "Index.h"
#include <cmath>



class Sellmeier_Index : public Index
{
public:
  Sellmeier_Index(double lambda) {}
  virtual ~Sellmeier_Index() {}
  bool do_CF4 = true; // Please only make one of these true....
  bool do_C2F6 = false;
  bool do_C4F10 = false;

  // Index of refraction at all wavelengths as defined by Sellmeier equation, default is CF4 with A = .12 and Lambda0 = 61.81 
  // Other options: C2F6 with A = .1746 and Lambda0 = 66.75, C4F10 with A = .2375, Lambda0 = 73.63, info scraped from https://cds.cern.ch/record/600182/files/ep-2002-099
  double n(double lambda)
  {

    double A; 
    double Lambda0;
   
    if(do_CF4)
      {
	A = .12;//Average of https://cds.cern.ch/record/600182/files/ep-2002-099 and Nucl. Instrum. Methods Phys. Res., A : 292(1990)593
	Lambda0 = 61.81;
	Index::GasName = "CF4";
	//cout << "Radiator is CF4" << endl;
      }
    
    if(do_C2F6)
      {
	A = .1746;
	Lambda0 = 66.75;
	Index::GasName = "C2F6";
	//cout << "Radiator is C2F6" << endl;
      }
    
    if(do_C4F10)
      {
	A = .2375;
	Lambda0 = 73.63;
	Index::GasName = "C4F10";
	//cout << "Radiator is C4F10" << endl;
      }
       
    if (lambda < Lambda0)
      {
	//cout << "Wavelength entered is below UV cutoff, gas and photocathode operate at difference wavelengths." << endl;
      }
    else
      {
	index = 1+(A*pow(10.0,-6.0))/((pow(Lambda0,-2.0) - pow(lambda,-2)));
	//cout << index << endl;
      }
    return index;
  }
  string Gas()
  {
    if(do_CF4)
      {
	Index::GasName = "CF4";
	//cout << "Radiator is CF4" << endl;
      }
    
    if(do_C2F6)
      {
	Index::GasName = "C2F6";
	//	cout << "Radiator is C2F6" << endl;
      }
    
    if(do_C4F10)
      {
	Index::GasName = "C4F10";
	//cout << "Radiator is C4F10" << endl;
      }
    
    return Index::GasName;
  }





 protected:
  double index;
};
	
#endif /* __SELLMEIER_INDEX_H__ */
