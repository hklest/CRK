#ifndef __CRK_H__
#define __CRK_H__

#include <vector>

class Index;
class Eff;
class Smear;
class Efftextfile;

class CRK
{
 public:
  CRK(Index *n, double RadiatorLength);
  virtual ~CRK() {}

  
  
  double p_m, beta;
  double rapidity,pixel_size,lambdamin,lambdamax,sumstraythetas,sumstraysigmas,p,avgphotons,ThetaCMax,Calpha,DeltaP,avgthetac,summedWL,summedthetac,alpha,NpMC,NAlphaMC,meanThC,sigmaThC,NP,PMax,PMin,NAlpha,AlphaMax,CentralP,ThCreco;
  double pi_mean = 0 , ka_mean = 0, pr_mean = 0, e_mean = 0, pi_RMS = 0, ka_RMS = 0, pr_RMS = 0, e_RMS = 0, pi_photons = 0, ka_photons=0,pr_photons=0,e_photons=0;
  double normalization = 0.0;
  double stepsize;
  double k = 0;
  double AlphaEM = 1/137.035999;
  double dNdxValue, Npe;
  double Indiv_Eff,SUMDEV;
  
  double i = 0.0;
  double Mult_E;
  double lambda;
  double Weighted_Efficiency;
  std::vector<double> Eff_Dist, ThetaCs;
  std::vector<double> Lambdas;
  std::vector<double> pi_means,ka_means,pr_means,e_means,pi_RMSs,ka_RMSs,pr_RMSs,e_RMSs,alphas,Ps,pi_photonsvec,ka_photonsvec,e_photonsvec,pr_photonsvec;
  double dNdxIntegral = 0;
  int nrolls,NStrayPhotons,Aincrement;
  string GasName,PCname;
  double r_in,r_out,B,pressure;
  unsigned seed; //= std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time                                                          
  unsigned seed2,seed3,seed4,seed5;
  std::default_random_engine LambdaGen;  // Initialize RNG
  std::default_random_engine strayxgen;  // Initialize RNG
  std::default_random_engine strayygen;  // Initialize RNG                                                                                                                 
  std::discrete_distribution<int> LambdaDist;
  



  
  void AddEff  (Eff   *);
  void AddSmear(Smear *);
  void Simulate_Momentum();
  void Simulate_Alpha();
  double Normalize();
  double dNdx();
  double AverageThetaC(double beta, std::vector<double> Eff_Dist, std::vector<double> Lambdas);
  void AssignPIDAlpha(int k, double alpha);
  void AssignPIDMomentum(int k, double p);
  void AlphaMCPlots();
  void MomentumMCPlots();
  void FillAlphas();
  void ProduceEffDist();
  void AssignWLMinMax();
  void MultEffPlot();
  void PixelError(double ThCreco, double pixelsize); 
  void EmissionError(double Npe, double L);
  
  // Root Tools/Histograms
  TFile *SimResults;
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
  TH1 *MultEffs;
  TH1 *NvThC;
  TH1 *NvsWL;
 protected:
  Index *n;
  double L;
  
   std::vector<Eff *> effs;
  std::vector<Smear *> smears;
};

#endif /* __CRK_H__ */
