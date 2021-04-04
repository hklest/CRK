#include "Index.h"
#include "Eff.h"
#include "Smear.h"
#include "CRK.h"
#include "Sellmeier_Index.h"
#include "Photocathode_Eff.h"
#include "Efftextfile.h"
#include <cmath>
#include <ctime>
#include "Functions.h"
using namespace std;

void CRK::Simulate_Momentum() //Simulate_Something provides outputs of efficiencies, theta_c, Npe, etc.
{
  //Input Parameters
  r_in = 50;
  r_out = 70;
  L = r_out-r_in;
  B = 1.5; // magnetic field in Tesla
  pixel_size = .1; // in cm
  PMax = 40; // Maximum momentum in MC in GeV 
  PMin = 0;  // Minimum momentum in MC                                                                                                                                                             
  NP = 80; // number of P points along x-axis
  Calpha = 0; // Error in Alpha (in radians)
  ThetaCMax = .14; // Max of y-axis in plot
  NpMC = 5000; // Number of times to run Monte Carlo
  double Pincrement = PMax/NP;

  //
  bool do_stray_photons = false; // Currently does not work
  bool do_NvsWL = false; 
  double time_resolution = 100; // ps
  double dark_rate = 100000; // hz
  double MeanStrayPhotons = .01;
  //
  
  // Histograms
  ThetaCvsP = new TH2D("ThetaCvsP","ThetaC vs. P",NP,0.0,PMax,500,0,ThetaCMax);
  pi_ka_vsP = new TH2D("pi_ka_vsP","N Sigma Pi-Ka vs. P",NP,0.0,PMax,600,.5,100);
  ka_pr_vsP = new TH2D("ka_pr_vsP","N Sigma Ka-Pr vs. P",NP,0.0,PMax,600,.5,100);
  e_pi_vsP = new TH2D("e_pi_vsP","N Sigma e-Pi vs. P",NP,0.0,PMax,600,.5,100);
 
  //
  
  
  // Print Max and Min Wavelengths for each Efficiency
  //   for(int i = 0; i < effs.size();i++)
  // {
  //   cout << "Eff " << i + 1<< "LMin is " << effs[i]->LMin() << endl;
  //  cout << "Eff " << i + 1<< "LMax is " << effs[i]->LMax() << endl;
  // }
   
   
   // initialize integration parameters
   AssignWLMinMax();
   lambdamin = 300; // Hardcoded cutoffs
   lambdamax = 1000;
   double lambdaref = 400;
   stepsize = 1; // step size for numerical integration 
   cout << "lambdamin is " << lambdamin << endl;
   cout << "lambdamax is " << lambdamax << endl;
   if(do_NvsWL)
     {
       NvsWL = new TH1D("NvsWL","Wavelengths",1000,lambdamin,lambdamax);  
     }
   // Particle Masses
   double mPi = 0.13957; // GeV                                                                                                                                                                         
   double mKa = 0.49368;                                                                                                                                                                                
   double mPr = 0.93827;                                                                                                                                                                                
   double me = .0005110;

   
   // ThetaC Parameters
   double PID [4] = {mPi, mKa, mPr, me};
   string PIDname [4] = {"Pion","Kaon","Proton","Electron"};
   double pi_mean = 0 , ka_mean = 0, pr_mean = 0, e_mean = 0, pi_RMS = 0, ka_RMS = 0, pr_RMS = 0, e_RMS = 0;
   for (double p = PMin; p <= PMax; p = p+Pincrement)
     {
       for (int k = 0; k <= 3; k++) // edit loop to select for subsets of particles
	 {
	   p_m = p/PID[k];
	   beta = sqrt(p_m*p_m/(1 + p_m*p_m));
	   //cout << PIDname[k] << " threshold is " << beta/(n->n(lambdaref)) << endl; 
	   //  cout << "lambdaref is " << lambdaref << " n at lref is " << n->n(lambdaref) << endl;  
	   //Procedure for skipping if below cherenkov threshold, saves time
	   double ThCPrelim = acos(1/((n->n(lambdaref))*beta)); // check if below cherenkov threshold
	   cout << PIDname [k] << " Preliminary ThetaC is " << ThCPrelim << " at index of " << n->n(lambdaref) << endl;
	   if(isnan(ThCPrelim))
	     {
	       continue;
	     }
	   //
	   
	   Normalize(); // These functions can be found in Functions.h
	   dNdxIntegral = dNdx();
	   ProduceEffDist();	 
	   avgthetac = AverageThetaC(beta,Eff_Dist,Lambdas);
	   cout << "avgthetac is " << avgthetac << endl;
	   AlphaEM = 1/137.035999;
	   dNdxValue = 2*M_PI*AlphaEM*pow(sin(avgthetac),2)*dNdxIntegral; // number of photons produced AND measured per unit length
	   Npe = dNdxValue*L*10000000; // large number is there to convert from nm to cm;
	   Npe = Npe*.7;
	   cout << "Npe = " << Npe << endl;
	   
	   //RNG for Monte Carlo
	   seed = std::chrono::system_clock::now().time_since_epoch().count();
	   std::default_random_engine PoisGen (seed);
	   std::poisson_distribution<int> Poisson(Npe);
	   seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	   std::default_random_engine PoisGen2 (seed2);
	   std::poisson_distribution<int> PoissonStray(MeanStrayPhotons);
	   TRandom Randy;	 
	   seed3 = std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time 
	   std::default_random_engine LambdaGen3 (seed3);  // Initialize RNG
	   std::discrete_distribution<int> LambdaDist3(Eff_Dist.begin(), Eff_Dist.end()); // returns an integer number with probabilities defined by the Eff_Dist vector, used in finding avgth
	   //
	   
	   // Momentum Monte Carlo


	   
	   cout << "Running Momentum Monte Carlo for " << PIDname[k] << " at " << p << " GeV" << endl;		 
	   
	   meanThC = 0;
	   sigmaThC = 0;
	   NStrayPhotons = 0;
	   avgphotons = 0;
	   int i = 0;
	   while(i<=NpMC)
	    {
	      int photons = Poisson(PoisGen);
	      // cout << "Photons = " << photons << endl;
	      avgphotons += photons;
	      //  cout << "N photons is " << photons << endl;
	      double sum = 0;
	      double sumsigma = 0;
	      for (int f=0; f<photons; f++)
		{
		  int WLindex = LambdaDist3(LambdaGen3);
		  double WL = Lambdas[WLindex];
		  if(do_NvsWL){NvsWL->Fill(WL);}
		  double ThetaC = acos(1/((n->n(WL))*beta));
		  double ringphi = 2*TMath::Pi()*Randy.Rndm();
		  double ThCreco = acos(sin(Calpha)*sin(ThetaC)*cos(ringphi) + cos(Calpha)*cos(ThetaC));
		  double ThCsmeared,xSmeared,ySmeared,TotalSmear;
		  ThCsmeared = ThCreco;
		  //cout << "Initial ThCSmeared is " << ThCsmeared << endl;
		  ThCsmeared = ThCsmeared + Randy.Gaus(0,pixel_size/(sqrt(12)*L)); // smear due to pixels
		  // cout << "gaussian smear from pixel size" << Randy.Gaus(0,pixel_size/L) << endl;
		  // cout << "First ThCSmeared is " << ThCsmeared << endl;
		  ThCsmeared = smears[0]->smr(ThCsmeared,WL,p,B,r_in,r_out,ringphi); //Additional smear
		  //cout << "Second ThCSmeared is " << ThCsmeared << endl;
		  ThCsmeared = smears[1]->smr(ThCsmeared,WL,p,B,r_in,r_out,ringphi); // uniform bending smear
		  // cout << "Third ThCSmeared is " << ThCsmeared << endl;
		  // for(i=0;i<smears.size();i++)
		  //   {
		  //    ThCsmeared = (smears[i]->smr(ThCsmeared,WL));
		  //    return;
		  //  }
		  // cout << "thetac smeared is " << ThCsmeared << endl;
		  sum += ThCsmeared;
		  sumsigma += pow(abs(avgthetac-ThCsmeared),2); // is avgthetac the right thing to use here?
		  //cout << "sigma is " << pow(abs(avgthetac-ThCsmeared),2) << endl;
		}
	      if (photons > 0)
		{
		  sumstraythetas=0;
		  sumstraysigmas=0;
		  NStrayPhotons = 0;
		  // if(do_stray_photons)
		  // {
		  //   NStrayPhotons = PoissonStray(PoisGen2);
		  //   int j=0;
		  //   while(j<NStrayPhotons)
		  //{
		  //	  double strayX = 1.5*Randy.Rndm()*avgthetac;
		  //	  double strayY = 1.5*Randy.Rndm()*avgthetac;
		  //	  double strayTheta = sqrt(pow(strayX,2)+pow(strayY,2))/L;
		  //	  sumstraythetas += strayTheta;
		  //	  sumstraysigmas += pow(abs(avgthetac-strayTheta),2);
			  // cout << straythetas<< endl;
		  //	  j++;
		  //	  avgphotons++;
		  //	}
		  //  sigmaThC += pow(sumstraysigmas/NStrayPhotons,.5);
		  //}
		  double thetaApparent = (sum+sumstraythetas)/((double)photons+NStrayPhotons);// adding in scintillation photons distributed randomly, if do_stray_photons is false, sumstray,nscint=0
		  sigmaThC += pow(sumsigma/photons,.5); //sum and divide by N outside photons loop to find mean sigma at a certain P
		  //cout << "SigmaThC = " << sigmaThC << endl;
		  meanThC += thetaApparent; // sum and divide by N at the end to find mean apparent theta at a P
		  //cout << "MeanThC = " << meanThC << endl;
		  ThetaCvsP->Fill(p,thetaApparent);
		  i++;
		}
	      
	    }
	  AssignPIDMomentum(k,p);   
	 }
    }      
  // Extract information for graphic
  string PCName = effs[0]->EffName();
  cout << effs[0]->EffName() << endl;
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  pressure = SMn.SMpressure();
  cout << "inner radius is " << r_in << endl;
  cout << "outer radius is " << r_out << endl;
  
  cout << "Radiator length is " << L << endl;

  cout << "Gas is " << GasName << endl;
  cout << "Pressure is " << pressure << " bar (absolute)" << endl;
  cout << "B-field is " << B << endl;
  //
  
  //Plot
  MomentumMCPlots();

  if(do_NvsWL){NvsWL->Draw();}
  //
}

void CRK::Simulate_Alpha() //Simulate_Something provides outputs of efficiencies, theta_c, Npe, etc.
{
  //Alpha MC Parameters
  double r_in = 50;
  double r_out = 75;
  double B = 1.5; 
  AlphaMax = 0.010; // radians
  NAlpha = 10; // number of alpha points along x-axis
  p = 33; // Momentum value where the alpha MC is run at
  NAlphaMC = 1000; // number of points in MC
  //

  // Whether or not to add in background
  bool do_stray_photons = false;
  double MeanStrayPhotons = .1;
  //
  
  //Histograms
  ThetaCMax = .06;
  pi_ka_vsAlpha = new TH2D("N Sigma Pi-Ka vs. Alpha","N Sigma Pi-Ka vs. Alpha",NAlpha,0.0,AlphaMax,1000,.5,50);
  ka_pr_vsAlpha = new TH2D("N Sigma Ka-Pr vs. Alpha","N Sigma Ka-Pr vs. Alpha",NAlpha,0.0,AlphaMax,1000,.5,50);
  e_pi_vsAlpha = new TH2D("N Sigma e-Pi vs. Alpha","N Sigma e-Pi vs. Alpha",NAlpha,0.0,AlphaMax,1000,.5,50);
  ThetaCvsAlpha = new TH2D("ThetaCvAlpha","ThetaC vs. Alpha",NAlpha,0.0,AlphaMax,1000,0,ThetaCMax);
  //

  //Determine Max and Min wavelengths to consider
  for(int i = 0; i < effs.size();i++)
    {
      cout << "Eff " << i + 1<< "LMin is " << effs[i]->LMin() << endl;
      cout << "Eff " << i + 1<< "LMax is " << effs[i]->LMax() << endl;
    }
  AssignWLMinMax();
  //
  // initialize integration parameters
  stepsize = 1; // step size for numerical integration 
  cout << "lambdamin is " << lambdamin << endl;
  cout << "lambdamax is " << lambdamax << endl;
  //
  
  // Particle Masses
  double mPi = 0.13957; // GeV                                                                                                                                                                         
  double mKa = 0.49368;                                                                                                                                                                                
  double mPr = 0.93827;                                                                                                                                                                                
  double me = .0005110;
  
  // ThetaC Parameters
  double PID [4] = {mPi, mKa, mPr, me};
  string PIDname [4] = {"Pion","Kaon","Proton","Electron"};
  double pi_mean = 0 , ka_mean = 0, pr_mean = 0, e_mean = 0, pi_RMS = 0, ka_RMS = 0, pr_RMS = 0, e_RMS = 0;
  for (int k = 0; k <= 3; k++) // edit loop to select for subsets of particles
    {
      p_m = p/PID[k];
      beta = sqrt(p_m*p_m/(1 + p_m*p_m));
      
      //Procedure for skipping if below cherenkov threshold, saves time
      double ThCPrelim = acos(1/(n->n((lambdamin+lambdamax)/2)*beta)); // check if below cherenkov threshold
      cout << PIDname [k] << "Preliminary ThetaC is " << ThCPrelim << endl;
      if(isnan(ThCPrelim))
	{
	  continue;
	}
      //
      Normalize();
      dNdxIntegral = dNdx();
      ProduceEffDist();
      avgthetac = AverageThetaC(beta,Eff_Dist,Lambdas);
      cout << "avgthetac is " << avgthetac << endl;
      AlphaEM = 1/137.035999;
      dNdxValue = 2*M_PI*AlphaEM*pow(sin(avgthetac),2)*dNdxIntegral; // number of photons produced AND measured per unit length
      // cout << "dNdxValue is " << dNdxValue<<endl;
      Npe = dNdxValue*L*10000000; // large number is there to convert from nm to cm;
      Npe = Npe*.7;
      cout << "Npe = " << Npe << endl;
      
      // RNG for Monte Carlo
      seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine PoisGen (seed);
      std::poisson_distribution<int> Poisson(Npe);
      seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine PoisGen2 (seed);
      std::poisson_distribution<int> PoissonStray(MeanStrayPhotons);
      TRandom Randy;
      seed = std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time 
      std::default_random_engine LambdaGen2 (seed);  // Initialize RNG
      std::discrete_distribution<int> LambdaDist2(Eff_Dist.begin(), Eff_Dist.end()); // returns an integer number with probabilities defined by the Eff_Dist vector, used in finding avgthetac and in MCs
      //
      
     
      if(do_stray_photons)
	{
	  // cout << NScintPhotons << endl;
	  // avgphotons = MeanStrayPhotons*NpMC;
	}
      else
	{
	  avgphotons=0;
	  // MeanStrayPhotons=0;
	}
      double alpha = 0;
      for (int i=1; i<=ThetaCvsAlpha->GetNbinsX(); i++)
	{
	  alpha = ThetaCvsAlpha->GetXaxis()->GetBinCenter(i);
	  meanThC = 0;
	  sigmaThC = 0;
	  avgphotons = 0;
	  cout << "Running Alpha Monte Carlo for " << PIDname[k] << " at p = " << p << " and " << "alpha = " << alpha << " rad" << endl;
	  //cout << "alpha = "<< alpha << endl;
	  for (int j = 0; j<NAlphaMC; j++)
	    {
	      
	      int photons = Poisson(PoisGen);
	      avgphotons += photons;
	      // cout << "N photons" << photons << endl;
	      double sum = 0;
	      double sumsigma = 0;
	      for (int f=0; f<photons; f++)
		{
		  int WLindex = LambdaDist2(LambdaGen2);
		  double WL = Lambdas[WLindex];
		  double ThetaC = acos(1/((n->n(WL))*beta));
		  double ringphi = 2*TMath::Pi()*Randy.Rndm();
		  double ThCreco = acos(sin(alpha)*sin(ThetaC)*cos(ringphi) + cos(alpha)*cos(ThetaC));
		  double ThCsmeared,TotalSmear;
		  ThCsmeared = ThCreco;
		  ThCsmeared = smears[3]->smr(ThCreco,WL,p,B,r_in,r_out,ringphi);
		  // for(i=0;i<smears.size();i++)
		  //  {
		  //  ThCsmeared = (smears[i]->smr(ThCsmeared,WL));
		  //    return;
		  // }
		  // cout << "thetac smeared is " << ThCsmeared << endl;
		  sum += ThCsmeared;
		  sumsigma += pow(abs(avgthetac - ThCsmeared),2);
		  //cout << sumsigma << endl;
		  
		}
	      if (photons > 0)
		{
		  sumstraythetas=0;
		  sumstraysigmas=0;
		  NStrayPhotons=0;
		  if(do_stray_photons)
		    {
		      
		      int j=0;
		      NStrayPhotons = PoissonStray(PoisGen2);
		      while(j<NStrayPhotons)
			{
			  double strayX = 1.5*Randy.Rndm()*avgthetac;
			  double strayY = 1.5*Randy.Rndm()*avgthetac;
			  double strayTheta = sqrt(pow(strayX,2)+pow(strayY,2))/L;
			  sumstraythetas += strayTheta;
			  sumstraysigmas += pow(abs(avgthetac-strayTheta),2);
			  // cout << straythetas<< endl;
			  j++;
			  avgphotons++;
			}
		      sigmaThC += pow(sumstraysigmas/NStrayPhotons,.5);
		    }
		  double thetaApparent = (sum+sumstraythetas)/((double)photons+NStrayPhotons);// adding in stray photons distributed randomly
		  sigmaThC += pow(sumsigma/photons,.5); //sum and divide by N outside photons loop to find mean sigma at a certain P
		  meanThC += thetaApparent;
		  // cout << meanThC << endl;		  
		  ThetaCvsAlpha->Fill(alpha,thetaApparent);
		}
	    }
	  AssignPIDAlpha(k,alpha);
	}  
    }
  
  // Extract info for plots
  string PCName = effs[0]->EffName();
  cout << effs[0]->EffName() << endl;
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  cout << "Radiator length is " << L << endl;
  cout << "Gas is " << GasName << endl;
  //
  
  //Fill and Plot
  FillAlphas();
  AlphaMCPlots();      
  //
}
