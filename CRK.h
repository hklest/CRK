#ifndef __CRK_H__
#define __CRK_H__

#include <vector>

class Index;
class Eff;
class Smear;
	
class CRK
{
public:
CRK(Index *n, double RadiatorLength);
virtual ~CRK() {}

void AddEff  (Eff   *);
void AddSmear(Smear *);

void Simulate_Something();

protected:
Index *n;
double L;

std::vector<Eff *> effs;
std::vector<Smear *> smears;
};
	
#endif /* __CRK_H__ */
