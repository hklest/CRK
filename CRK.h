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

  double p_m, beta;
  double lambdamin,lambdamax,p,ThetaCMax,Calpha,DeltaP,avgthetac,summedWL,summedthetac,alpha,NpMC,NAlphaMC,meanThC,sigmaThC,NP,PMax,PMin,NAlpha,AlphaMax,CentralP;
  double pi_mean = 0 , ka_mean = 0, pr_mean = 0, e_mean = 0, pi_RMS = 0, ka_RMS = 0, pr_RMS = 0, e_RMS = 0;
  double normalization = 0.0;
  double stepsize;
  double k = 0;
  double AlphaEM = 1/137.035999;
  double dNdxValue, Npe;
  
  double i = 0.0;
  double Mult_E;
  double lambda;
  double Weighted_Efficiency;
  std::vector<double> Eff_Dist;
  std::vector<double> Lambdas;
  double dNdxIntegral = 0;
  int nrolls;
  string GasName,PCname;
  unsigned seed; //= std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time                                                          
  std::default_random_engine generator;  // Initialize RNG                                                                                                                 
  std::discrete_distribution<int> dist;
  



  
  void AddEff  (Eff   *);
  void AddSmear(Smear *);
  void Simulate_Something();
  double Normalize();
  void dNdx();
  void AverageThetaC(double beta);
  void AssignPIDAlpha(int k, double alpha);
   void AssignPIDMomentum(int k, double p);
  void AlphaMCPlots();
  void MomentumMCPlots();


  // Root Tools/Histograms
    TH2 *ThetaCvsP;
   TH2 *DeltaThCvsP;
   TH2 *pi_ka_vsP;
   TH2 *ka_pr_vsP;
  TH2 *e_pi_vsP;
  TH2 *ThetaCvsAlpha;
  TH2 *DeltaThCvsAlpha;
  TH2 *pi_ka_vsAlpha;
  TH2 *ka_pr_vsAlpha;
  TH2 *e_pi_vsAlpha;

 protected:
  Index *n;
  double L;
  
   std::vector<Eff *> effs;
  std::vector<Smear *> smears;
};

#endif /* __CRK_H__ */
