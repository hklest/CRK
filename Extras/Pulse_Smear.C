#include "Pulse_Smear.h"

using namespace std;

Pulse_Smear::Pulse_Smear(double SMEAR)
{
}
double Pulse_Smear::smr(double ThetaCherenkov, double lambda) //User should make sure the CSV input file has reasonable endpoints
{
  double  smearedThC = ThetaCherenkov;
  return smearedThC;
}
	

