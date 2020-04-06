#include "Index.h"
#include "Eff.h"
#include "Smear.h"

#include "CRK.h"

using namespace std;
	
CRK::CRK(Index *index, double RadiatorLength)
{
  n = index;
  L = RadiatorLength;
}

void CRK::AddEff(Eff *e)
{
  effs.push_back(e);
}

void CRK::AddSmear(Smear *s)
{
  smears.push_back(s);
}

void CRK::Simulate_Something()
{
  cout << "You should do something here." << endl;
  double Total_Efficiency=1.0;
  for (int i=0; i<effs.size(); i++)
    {
      Total_Efficiency *= effs[i]->e(500);
    }

  double ThetaCth = 0.050; // about 3 degrees.
  double ThetaCms = ThetaCth;
  for (int i=0; i<smears.size(); i++)
    {
      ThetaCms = smears[i]->smr(ThetaCms,500);
    }

  cout << " Efficiency= " << Total_Efficiency;
  cout << " ThetaC(theory)= " << ThetaCth;
  cout << " ThetaC(measured)= " << ThetaCms;
  cout << endl;

}

