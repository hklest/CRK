#ifndef __PHOTOCATHODE_EFF_H__
#define __PHOTOCATHODE_EFF_H__

#include "Eff.h"
#include <fstream>
#include <iostream>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort                                                                                                                                   
#include <vector>       // std::vector
#include <string>


// PhotoCathode_Eff takes a CSV "PhotocathodeEff.txt" (as produced by DataThief) and produces the quantum efficiency of a photocathode (Default CsI) at a specific wavelength in nm. This is accomplished by splitting the columns of the CSV into two vectors, one containing efficiences (CsIEffs) and one containing wavelengths (CsILambdas). The indices of the vectors have a one-to-one correspondence, so finding the index of the nearest wavelength to the one entered produces also the index of the closest efficiency.                    



class Photocathode_Eff : public Eff
{
public:
  Photocathode_Eff(double lambda) {}
  virtual ~Photocathode_Eff() {}

  double lambda_min, lambda_max;
  bool do_CsI = true;
  bool do_LAPPD = false;
  bool do_H12700A = false;
  string PC;
  // Options are LAPPD_Eff.txt, CsI_Eff.txt, H12700A.txt others TBD ...
  // Due to difference in wavelength coverage for the various options, the time the simulation takes can vary.

  // This efficiency has a wavelength dependence
  double e(double lambda)
  {
    if (do_H12700A)
      {
	PC = "H12700A_Eff.txt";
      }
    if(do_LAPPD)
      {
	PC = "LAPPD_Eff.txt";
      }
    if(do_CsI)
      {
	PC = "CsI_Eff.txt";
      }
   
    ifstream inFile;
    double eff;
    vector<double> Holder; //Vector for holding values from DataThief CSV                                                                                                                                   
    vector<double> PhotocathodeLambdas; // Vector that holds first column of CSV (wavelengths)                                                                                                                    
    vector<double> PhotocathodeEffs;// Vector that holds second column of CSV (efficiencies)                                                                                                                      
    string str; //Placeholder string                                                                                                                                                                        
    char * value; //pointer to placeholder value                                                                                                                                                            
    double hold; // placeholder                                                                                                                                                                             
    int i = 0; // ints for looping over                                                                                                                                                                     
    int j = 0;

    inFile.open(PC.c_str()); // Open DataThief CSV, Options are LAPPD_Eff.txt, CsI_Eff.txt, H12700A.txt others TBD ...                                                                                                                                                     
    if (!inFile)
      {
        cerr << "End of File or Unable to open file";
        exit(1);   // call system to stop if unable to open                                                                                                                                                 
      }
    while (std::getline(inFile, str)&&!str.empty()) //while sorting through lines of CSV                                                                                                                    
      {
        char cstr[str.size() + 1];
        strcpy(cstr, str.c_str());// copy string into a character array for strtok to oprate on UNIT TEST ->//cout << cstr << endl;                                                                         
        char * value = strtok(cstr,", ");// Split on "," and line break, store double in "value", strtok returns pointer to token so value needs be a pointer, if no delimiter found returns NULL UNIT TEST\
 ->//cout << value2 << endl;                                                                                                                                                                                
        while (value != NULL)
          {
            double val = atof(value); //Convert value from character pointer to double                                                                                                                      
            Holder.push_back(val);// push doubles from CSV into holder vector                                                                                                                               
            value = strtok(NULL,", /n");
          }
      }

    while (i < Holder.size()/2) // looking to pull wavelengths out of CSV (even indicies of Holder)                                                                                                         
      {
        PhotocathodeEffs.push_back(Holder.at(2*i+1)); // push odd values of holder onto PhotocathodeEffs vector                                                                                                         
        PhotocathodeLambdas.push_back(Holder.at(2*i)); // push even values of holder onto PhotocathodeLambdas vector                                                                                                   
        i++;
      }
    //   cout << "Your wavelength input is: " << lambda << "nm" <<endl;

    int index;
    double Max_Lambda = *max_element(PhotocathodeLambdas.begin(),PhotocathodeLambdas.end());
    double Min_Lambda = *min_element(PhotocathodeLambdas.begin(),PhotocathodeLambdas.end());

    
    // Determine element of wavelength vector closest to the input wavelength, then pull out the corresponding efficiency                                                                                  
    if (lambda < Min_Lambda) // Set wavelengths below min in data = 0                                                                                                                                       
      {
        eff = PhotocathodeEffs.at(0);
	// cout << "Wavelength entered is below minimum in dataset, efficiency set to 0" << endl;
      }
    else if (lambda > Max_Lambda) // set wavelengths above max in data = max                                                                                                                                
      {
        eff = PhotocathodeEffs.at(std::distance(PhotocathodeLambdas.begin(),PhotocathodeLambdas.end())-1);
        //cout << "Wavelength entered is above maximum in dataset, efficiency set to maximum value = " << eff << endl;
      }
    else
      {
        //std::vector<double>::iterator low; // find the point where                                                                                                                                        
        //low = std::lower_bound(MirrorLambdas.begin(),MirrorLambdas.end(), lambda);                                                                                                                        
        //      index = std::distance(MirrorLambdas.begin(), low);                                                                                                                                          
        //}                                                                                                                                                                                                 

        double difference = 0.0;
        double smallestdifference = 10000.00; // value to be overwritten                                                                                                                                    
        for(int k = 0; k < PhotocathodeLambdas.size(); k++) // loop for finding index of value in PhotocathodeLambdas that is closest to input wavelength.

	  for(int k = 0; k < PhotocathodeLambdas.size(); k++) // loop for finding index of value in PhotocathodeLambdas that is closest to input wavelength.                                                              
          {
            difference = lambda - PhotocathodeLambdas[k];
            difference = sqrt(pow(difference,2));//make sure difference is positive                                                                                                                         
            //cout << difference << endl;                                                                                                                                                                   
             // cout << smallestdifference << endl;                                                                                                                                                         
            if(difference < smallestdifference)
              {
                smallestdifference = difference;
                //cout << smallestdifference<< endl;                                                                                                                                                        
                //cout << k << endl;                                                                                                                                                                        
                index = k;
              }
            //cout << index << endl;                                                                                                                                                                        
            eff =  PhotocathodeEffs[index];
                                                                                                                                                                                        
          }
	
      }
    //    cout <<"photocathode efficiency is: " << eff*.01 <<endl;


    
    inFile.close(); //close io stream                                                                                                                                                                       
    return eff*.01; // return efficiency in decimal form, e.g .876             
  }


  string PCName()
  {
    string PCn;
    if (do_H12700A)
      {
	PCn = "H12700A";
      }
    if(do_LAPPD)
      {
	PCn = "LAPPD";
      }
    if(do_CsI)
      {
	PCn = "CsI GEM";
      }
    return PCn;
  }
 double LMin()
  {
    if (do_H12700A)
      {
	lambda_min = 200;
      }
    if(do_LAPPD)
      {
	lambda_min = 90;
      }
    if(do_CsI)
      {
	lambda_min = 60;
      }
    return lambda_min;
  }
 double LMax()
  {
    if (do_H12700A)
      {
	lambda_max = 650;
      }
    if(do_LAPPD)
      {
	lambda_max = 630;
      }
    if(do_CsI)
      {
	lambda_max = 200;
      }
    return lambda_max;
  }


  

protected:
  double eff;
};
	
#endif /* __PHOTOCATHODE_EFF_H__ */
