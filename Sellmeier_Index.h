#ifndef __SELLMEIER_INDEX_H__
#define __SELLMEIER_INDEX_H__

#include "Index.h"
#include <cmath>



class Sellmeier_Index : public Index
{
public:
  Sellmeier_Index(double lambda) {}
  virtual ~Sellmeier_Index() {}

  //
  bool fixed_index = false; // sets a fixed index of refraction (unphysical)
  double fixed_index_value = 1.004; // Whatever index you'd like
  //
  double pressure = 3; // pressure in bar, ie 2 bar = 1 atm overpressure
  //
  bool do_CF4 = false; // Please only make one of these true....
  bool do_C2F6 = false;
  bool do_C4F10 = false;
  bool do_Ar = false;
  bool do_H2O = false;
  bool do_C3F8 = true;
  // Index of refraction at all wavelengths as defined by Sellmeier equation, default is CF4 with A = .12 and Lambda0 = 61.81 
  // Other options: C2F6 with A = .1746 and Lambda0 = 66.75, C4F10 with A = .2375, Lambda0 = 73.63, info scraped from https://cds.cern.ch/record/600182/files/ep-2002-099
  double n(double lambda)
  {
    
    double A; 
    double Lambda0;
    
    if (do_Ar)
      {
	index = 1+pressure*(.00250141/(91.012-pow(.001*lambda,-2))+.000500283/(87.892-pow(.001*lambda,-2))+.0522343/(214.02-pow(.001*lambda,-2)));
	return index;
      }
    else
      {
	if(do_CF4)
	  {
	    A = .12;//Average of https://cds.cern.ch/record/600182/files/ep-2002-099 and Nucl. Instrum. Methods Phys. Res., A : 292(1990)593
	    Lambda0 = 61.81;
	    Index::GasName = "CF4";
	    
	    //cout << "Radiator is CF4" << endl;
	  }
	if(do_C2F6)
	  {
	    A = 0.1746;
	    Lambda0 = 66.75;
	    Index::GasName = "C2F6";
	    //cout << "Radiator is C2F6" << endl;
	  }
	
	if(do_C3F8)
	  {
	    A = .2305;
	    Lambda0 = 67.90;
	    Index::GasName = "C3F8";
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
	    index = 1000000000000; 
	    //cout << "Wavelength entered is below UV cutoff, gas and photocathode operate at difference wavelengths." << endl;
	  }
	else
	  {
	    index =1+pressure*(A*pow(10.0,-6.0))/((pow(Lambda0,-2.0) - pow(lambda,-2)));
	    //cout << index << endl;
	  }
	if(do_H2O)
	  {
	    index = sqrt(1+(.75821*(pow((.001*lambda),2))/(pow((.001*lambda),2)-.01007)+.08495*(pow((.001*lambda),2))/(pow((.001*lambda),2)-8.91377)));
	  }
	if (fixed_index)
	  {
	    index = fixed_index_value;
	  }
	return index;
      }
  }
  string Gas()
  {
    if(do_CF4)
      {
	Index::GasName = "CF4";
	//cout << "Radiator is CF4" << endl;
      }
    if(do_Ar)
      {
	Index::GasName = "Pressurized Argon";
	//cout << "Radiator is CF4" << endl;
      }
    
    if(do_C2F6)
      {
	Index::GasName = "C2F6";
	//	cout << "Radiator is C2F6" << endl;
      }
    if(do_C3F8)
      {
	Index::GasName = "C3F8";
	//	cout << "Radiator is C2F6" << endl;
      }
    
    if(do_C4F10)
      {
	Index::GasName = "C4F10";
	//cout << "Radiator is C4F10" << endl;
      }
    
    return Index::GasName;
  }
  
   double SMpressure()
   {
   return pressure;
   }

  



 protected:
  double index;
};
	
#endif /* __SELLMEIER_INDEX_H__ */
