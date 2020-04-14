#ifndef __MIRROR_EFF_H__
#define __MIRROR_EFF_H__

#include "Eff.h"
#include <fstream>
#include <iostream>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector


// Mirror_Eff takes a CSV "MirrorEff.txt" (as produced by DataThief) and produces the efficiency of reflection for a mirror such as the one used at the 2015 BNL-SBU RICH test beam (https://arxiv.org/pdf/1501.03530.pdf) at a specific wavelength in nm. This is accomplished by splitting the columns of the CSV into two vectors, one containing efficiences (MirrorEffs) and one containing wavelengths(MirrorLambdas). The indices of the vectors have a one-to-one correspondence, so finding the index of the nearest wavelength to the one entered produces also the index of the closest efficiency.



class Mirror_Eff : public Eff
{
 public:
  Mirror_Eff(double lambda) {}
  virtual ~Mirror_Eff() {}
  
  // This efficiency has a wavelength dependence, goal is to take an input wavelength and return an efficiency at that specific wavelength
  double e(double lambda)
  {
    ifstream inFile;
    double eff;
    vector<double> Holder; //Vector for holding values from DataThief CSV
    vector<double> MirrorLambdas; // Vector that holds first column of CSV (wavelengths)
    vector<double> MirrorEffs;// Vector that holds second column of CSV (efficiencies)
    string str; //Placeholder string
    char * value; //pointer to placeholder value
    double hold; // placeholder
    int i = 0; // ints for looping over
    int j = 0;
    
    inFile.open("MirrorEff.txt"); // Open DataThief CSV
    if (!inFile) 
      {
	cerr << "End of File or Unable to open file";
	exit(1);   // call system to stop if unable to open
      }
    while (std::getline(inFile, str)&&!str.empty()) //while sorting through lines of CSV
      {
	char cstr[str.size() + 1];
	strcpy(cstr, str.c_str());// copy string into a character array for strtok to oprate on UNIT TEST ->//cout << cstr << endl;
	char * value = strtok(cstr,", ");// Split on "," and line break, store double in "value", strtok returns pointer to token so value needs be a pointer, if no delimiter found returns NULL UNIT TEST ->//cout << value2 << endl; 
	while (value != NULL) 
	  { 
	    double val = atof(value); //Convert value from character pointer to double
	    Holder.push_back(val);// push doubles from CSV into holder vector
	    value = strtok(NULL,", /n");
	  }
      }
    
    while (i < Holder.size()/2) // looking to pull wavelengths out of CSV (even indicies of Holder)
      {
	MirrorEffs.push_back(Holder.at(2*i+1)); // push odd values of holder onto MirrorEffs vector 
	MirrorLambdas.push_back(Holder.at(2*i)); // push even values of holder onto Mirror Lambdas vector
	i++;
      }
    
    
    //cout << "Your wavelength input is: " << lambda << "nm" <<endl;
    
    int index;    
    double Max_Lambda = *max_element(MirrorLambdas.begin(),MirrorLambdas.end());
    double Min_Lambda = *min_element(MirrorLambdas.begin(),MirrorLambdas.end());
    //  cout << Max_Lambda << "   " << Min_Lambda << endl;
    //  cout << lambda << endl;
    if (lambda < Min_Lambda) // Set wavelengths below min in data = 0
      {
	eff = 0;
	//cout << "Wavelength entered is below minimum in dataset, efficiency set to 0" << endl;
      }
    else if (lambda > Max_Lambda) // set wavelengths above max in data = max
      {
	eff = MirrorEffs.at(std::distance(MirrorLambdas.begin(),MirrorLambdas.end())-1);
	//cout << "Wavelength entered is above maximum in dataset, efficiency set to maximum value = " << eff << endl;
      }		  
    else
      {
	//std::vector<double>::iterator low; // find the point where
	//low = std::lower_bound(MirrorLambdas.begin(),MirrorLambdas.end(), lambda);
	//	index = std::distance(MirrorLambdas.begin(), low);
	//}
      	
	double difference = 0.0;
	double smallestdifference = 10000.00; // value to be overwritten
	for(int k = 0; k < MirrorLambdas.size(); k++) // loop for finding index of value in MirrorLambdas that is closest to input wavelength.
	  {
	    difference = lambda - MirrorLambdas[k];
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
	    eff =  MirrorEffs[index];
	    //cout <<eff<<endl;
	  }
      }
    // cout << "mirror efficiency is: " << eff*.01 << endl;	
    inFile.close(); //close io stream
    return eff*.01; // return efficiency in decimal form e.g .876
  }
 protected:
  double eff; // make sure eff stays here
};

#endif /*__MIRROR_EFF_H__ */

  
