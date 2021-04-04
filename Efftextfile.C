#include "Efftextfile.h"

using namespace std;

Efftextfile::Efftextfile(string filename, int sparse = 1)
{
  Name = filename;
  ifstream inFile;
  vector<double> Holder; //Vector for holding values from DataThief CSV                                                                                                                                 
  string str; //Placeholder string                                                                                                                                                                     
  char * value; //pointer to placeholder value                                                                                                                                                            

  inFile.open(filename); // Open DataThief CSV, Options are Mirror_Eff.txt, LAPPD_Eff.txt, CsI_Eff.txt, H12700A.txt others TBD ...                                                                                                                                                     
  if (!inFile)
    {
      cerr << "End of File or Unable to open file " << Name << endl;
      exit(1);   // call system to stop if unable to open                                                                                                                                               
    }
  while (std::getline(inFile, str)&&!str.empty()) //while sorting through lines of CSV                                                                                                                  
    {
      Holder.clear();
      char cstr[str.size() + 1];
      strcpy(cstr, str.c_str());// copy string into a character array for strtok to oprate on UNIT TEST ->//cout << cstr << endl;                                                                       
      char * value = strtok(cstr,", ");// Split on "," and line break, store double in "value", strtok returns pointer to token so value needs be a pointer, if no delimiter found returns NULL UNIT TEST \
 ->//cout << value2 << endl;                                                                                                                                                                              
      while (value != NULL)
	{
	  double val = atof(value); //Convert value from character pointer to double                                                                                                                    
	  Holder.push_back(val);// push doubles from CSV into holder vector                                                                                                                             
	  value = strtok(NULL,", /n");
	}
      if(Holder.size()!=2)
	{
	  cerr << "Wrong number of entries" << Holder.size() << endl;
	  exit(1);
	}
      Effs.push_back(.01*Holder.at(1)); 
      Lambdas.push_back(Holder.at(0));                                                                                                   
    }

  //Sort in wavelength order
  for(int i = 0; i < Lambdas.size() - 1;i++)
    {
      for (int j = i; j<Lambdas.size(); j++)
	{
	  if(Lambdas[j]<Lambdas[i])
	    {
	      double tempL = Lambdas[i];
	      double tempE = Effs[i];
	      Lambdas[i] = Lambdas[j];
	      Effs[i] = Effs[j];
	      Lambdas[j] = tempL;
	      Effs[j] = tempE;
	    }
	}
    }
  inFile.close(); //close io stream
}

double Efftextfile::e(double lambda) //User should make sure the CSV input file has reasonable endpoints
{
  double Max_Lambda = Lambdas[Lambdas.size()-1];
  double Min_Lambda = Lambdas[0];
  
  // Determine element of wavelength vector closest to the input wavelength, then pull out the corresponding efficiency                                                                                  
  if (lambda <= Min_Lambda) // Set wavelengths below min in data = efficiency at leftmost point                                                                                                        
    {
      return Effs[0];
    }
  
  if (lambda >= Max_Lambda) // set wavelengths above max in data = efficiency at rightmost point                                                                                                                                
    {
      return Effs[Effs.size()-1];
    }
  
  for(int k = 1; k < Lambdas.size(); k++) // loop for finding index of value in PhotocathodeLambdas that is closest to input wavelength.                                                              
    {
      if(Lambdas[k]>=lambda)
	{
	  return Effs[k-1]+(lambda-Lambdas[k-1])*(Effs[k]-Effs[k-1])/(Lambdas[k]-Lambdas[k-1]);
	}
    }
  
  return Effs[Effs.size()-1];
  
}


double Efftextfile::LMin() // Values chosen to avoid singularities in Sellmeier coefficients. 
{
  return Lambdas[0];
}
double Efftextfile::LMax()
{
  return Lambdas[Lambdas.size()-1];
}

	

