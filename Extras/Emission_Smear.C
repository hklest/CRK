#include "Emission_Smear.h"

using namespace std;

Emission_Smear::Emission_Smear(double SMEAR)
{
}
double Emission_Smear::smr(double ThetaCherenkov, double lambda) //User should make sure the CSV input file has reasonable endpoints
{
  double smearedThC = ThetaCherenkov;
  return smearedThC;
}
	

