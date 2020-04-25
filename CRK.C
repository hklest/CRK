#include "Index.h"
#include "Eff.h"
#include "Smear.h"
#include "CRK.h"
#include "Sellmeier_Index.h"
#include "Photocathode_Eff.h"
//#include "WL_Dep.h"
#include <cmath>
#include <ctime>
using namespace std;
double lambda = 0.0;



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
double CRK::Normalize()
{
 //Find normalization constant
  normalization = 0;
  Photocathode_Eff PCE(lambda); 
  // initialize integration parameters
  lambdamin = PCE.LMin(); // Pulling mimumum wavelength observed by photocathode
  lambdamax = PCE.LMax(); // Pulling maximum wavelength observed by photocathode
  stepsize = 1; // step size for numerical integration 
 
    for (int i = lambdamin; i <= lambdamax; i=i+stepsize)
      {
	// cout << lambdamin << lambdamax<<endl;
        normalization += pow(i,-2.0);                                                                                                                                     
      }
    //    cout << normalization << endl;
    return normalization;
}
void CRK::dNdx() // Function that produces an efficiency distribution vector and produces the integral term "dNdxIntegral"
{
  dNdxIntegral = 0;
  Photocathode_Eff PCE(lambda);
  lambdamin = PCE.LMin(); // Pulling mimumum wavelength observed by photocathode
  lambdamax = PCE.LMax(); // Pulling maximum wavelength observed by photocathode
  // cout << lambdamin << lambdamax << endl;
//Loop over wavelengths between lambdamin nm and lambdamax nm to find quantities at each wavelength in steps of stepsize nm
  for (lambda = lambdamin; lambda <= lambdamax; lambda=lambda+stepsize)   
    {
      //cout << "normalization is "<<normalization<< endl;
      //Efficiency calculations
      Mult_E = (effs[0]->e(lambda))*(effs[1]->e(lambda))*(effs[2]->e(lambda));
      // cout << "Mult_E INSIDE is " << Mult_E << endl;
      // Multiplying the available efficiencies by each other, the above needs to change when new efficiences are added by multiplying by (eff[3]->e(lambda)      
      double Indiv_Eff = Mult_E*pow(lambda,-2); // store efficiency and likelihood of emission at specific lambda
      Eff_Dist.push_back(floor(1000000*Indiv_Eff)); // create a vector of efficiencies
      Lambdas.push_back(lambda); // create a vector of lambdas
      Weighted_Efficiency += Mult_E*pow(lambda,-2)/normalization; // metric of the full survival/measurement rate of all photons
      dNdxIntegral += Mult_E*stepsize*pow(lambda,-2); // integral for use in dN/dx calculation 
    }
  cout << "dndx integral is " << dNdxIntegral<< endl;  
  return;
}



void CRK::AverageThetaC(double beta) // Function that gives an average value of theta_C for seeding the monte carlo
{
  cout << beta<< endl;
  //Use a realistic distribution of photons to determine average theta_C detected for the detector configuration
  nrolls = 100; // number of lambdas to generate 
  seed = std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time 
  std::default_random_engine generator (seed);  // Initialize RNG
  std::discrete_distribution<int> dist(Eff_Dist.begin(), Eff_Dist.end()); // returns an integer number with probabilities defined by the Eff_Dist vector.
  summedWL = 0; // summed wavelengths
  summedthetac = 0; // summed theta_c
  for (int A=0; A<nrolls; ++A) // Loop for determining average wavelength and thetaC
    {
      int WLindex = dist(generator); // Picks a random index 
      double WL = Lambdas[WLindex];
      //      cout << WL << endl;
      // cout << n->n(WL)<<endl;// picks out the value of the wavelength at that index
      double ThetaCWL = acos(1/((n->n(WL))*beta)); // Finds the thetaC at that wavelength
      if (!isnan(ThetaCWL))
	{
	  //summedWL += WL;
	  summedthetac += ThetaCWL;
	}
    }
  avgthetac = summedthetac/nrolls;
  cout << "avgtheta = " << avgthetac << endl;
  return;
}


void CRK::AssignPIDAlpha(int k, double alpha)
{
  
  if(k==0)// if pion
    {
      pi_mean = meanThC/NAlphaMC;
      cout << "pi_mean = " << pi_mean << endl;
      pi_RMS = sigmaThC/NAlphaMC;
      cout << "pi_RMS = " << pi_RMS << endl;
    }
  if(k==1) // if kaon
    {
      ka_mean = meanThC/NAlphaMC;
      cout << "ka_mean = " << ka_mean << endl;
      ka_RMS = sigmaThC/NAlphaMC;
      cout << "ka_RMS = " << ka_RMS << endl;
      pi_ka_vsAlpha->Fill(alpha,(abs(pi_mean-ka_mean)/((pi_RMS+ka_RMS)/2)));
      cout << "N sigma pi-k = " << (abs(pi_mean-ka_mean)/((pi_RMS+ka_RMS)/2)) << endl;
    }
  if(k==2)// if proton
    {
      pr_mean = meanThC/NAlphaMC;
      cout << "pr_mean = " << pr_mean << endl;
      pr_RMS = sigmaThC/NAlphaMC;
      cout << "pr_RMS = " << pr_RMS << endl;
      ka_pr_vsAlpha->Fill(alpha,(abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)));
      cout << "N sigma k-pr = " << (abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)) << endl;
    }
  if(k==3)// if electron
    {
      e_mean = meanThC/NAlphaMC;
      cout << "e_mean = " << e_mean << endl;
      e_RMS = sigmaThC/NAlphaMC;
      cout << "e_RMS = " << e_RMS << endl;
      e_pi_vsAlpha->Fill(alpha,(abs(e_mean-pi_mean)/((e_RMS+pi_RMS)/2)));
      cout << "N sigma e-pi = " << (abs(pi_mean-e_mean)/((pi_RMS+e_RMS)/2)) << endl;
    }
  cout << "Moving on....." << endl << endl;

 return;
 }


void CRK::AssignPIDMomentum(int k, double p)
{
if(k==0)// if pion
  {
    pi_mean = meanThC/NpMC;
    cout << "pi_mean = " << pi_mean << endl;
    pi_RMS = sigmaThC/NpMC;
    cout << "pi_RMS = " << pi_RMS << endl;
  }
 if(k==1) // if kaon
   {
     ka_mean = meanThC/NpMC;
     cout << "ka_mean = " << ka_mean << endl;
     ka_RMS = sigmaThC/NpMC;
     cout << "ka_RMS = " << ka_RMS << endl;
     pi_ka_vsP->Fill(p,(abs(pi_mean-ka_mean)/((pi_RMS+ka_RMS)/2)));
     cout << "N sigma k-pi = " << (abs(ka_mean-pi_mean)/((pi_RMS+ka_RMS)/2)) << endl;
   }
 if(k==2)// if proton
   {
     pr_mean = meanThC/NpMC;
     cout << "pr_mean = " << pr_mean << endl;
     pr_RMS = sigmaThC/NpMC;
     cout << "pr_RMS = " << pr_RMS << endl;
     ka_pr_vsP->Fill(p,(abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)));
     cout << "N sigma k-pr = " << (abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)) << endl;
   }
 if(k==3)// if electron
   {
     e_mean = meanThC/NpMC;
     cout << "e_mean = " << e_mean << endl;
     e_RMS = sigmaThC/NpMC;
     cout << "e_RMS = " << e_RMS << endl;
     e_pi_vsP->Fill(p,(abs(e_mean-pi_mean)/((e_RMS+pi_RMS)/2)));
     cout << "N sigma e-pi = " << (abs(pi_mean-e_mean)/((pi_RMS+e_RMS)/2)) << endl;
   }
 cout << "Moving on....." << endl << endl;
 
 
 return;
}

void CRK::MomentumMCPlots()
{
   Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  cout << "Gas is " << GasName << endl;
   Photocathode_Eff PCE(lambda);
  string PCname = PCE.PCName();
 gStyle->SetOptStat(0);
      TCanvas *c1 = new TCanvas("c1","ThetaC vs. P",100,300,800,1000);
      //c1->Divide(1,2);
      //c1->cd(1);
      ThetaCvsP->Draw("colz");
      // c1->cd(2);
      // DeltaThCvsP->Draw("colz");
      TCanvas *sigmas = new TCanvas("sigmas", "Sigmas",100,800,1200,800);
      pi_ka_vsP->SetMarkerStyle(kFullSquare);
      ka_pr_vsP->SetMarkerStyle(kFullSquare);
      e_pi_vsP->SetMarkerStyle(kFullSquare);
      sigmas->Divide(3,1);
      sigmas->cd(1);
      pi_ka_vsP->Draw("PLC PMC");
       TLine *l1=new TLine(0,3.0,PMax,3.0);
      l1->SetLineColor(kBlue);
      l1->Draw();  
      TPaveText *pt = new TPaveText(0.15,0.7,0.8,0.85,"NDC");
      pt->SetTextSize(0.03);
      pt->SetFillColor(0);
      pt->SetTextAlign(12);
      pt->AddText(Form("Gas is %s, Alpha is %g rad",GasName.c_str(),Calpha));
      pt->AddText(Form("Photocathode is %s",PCname.c_str()));
      pt->AddText(Form("Lambdamin is %g nm, lambdamax is %g nm",lambdamin,lambdamax));
      pt->Draw();
    
     
      sigmas->cd(2);
      ka_pr_vsP->Draw("PLC PMC");
       TLine *l2=new TLine(0,3.0,PMax,3.0);
      l2->SetLineColor(kBlue);
      l2->Draw();
     
      sigmas->cd(3);
      e_pi_vsP->Draw("PLC PMC");
      TLine *l3=new TLine(0,3.0,PMax,3.0);
      l3->SetLineColor(kBlue);
      l3->Draw();
  
  return;
}

void CRK::AlphaMCPlots()
{
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  cout << "Gas is " << GasName << endl;
  Photocathode_Eff PCE(lambda);
  string PCname = PCE.PCName();
gStyle->SetOptStat(0);
      TCanvas *c1 = new TCanvas("c1","",100,300,800,1000);
      c1->Divide(1,2);
      c1->cd(1);
      ThetaCvsAlpha->Draw("colz");
      c1->cd(2);
      DeltaThCvsAlpha->Draw("colz");
      TCanvas *sigmas = new TCanvas("sigmas", "Sigmas",500,800,1200,800);
      pi_ka_vsAlpha->SetMarkerStyle(kFullSquare);
      ka_pr_vsAlpha->SetMarkerStyle(kFullSquare);
      e_pi_vsAlpha->SetMarkerStyle(kFullSquare);
     
      sigmas ->Divide(3,1);
      sigmas->cd(1);
      pi_ka_vsAlpha->Draw("PLC PMC");
      TLine *l1=new TLine(0,3.0,AlphaMax,3.0);
      l1->SetLineColor(kBlue);
      l1->Draw();
      TPaveText *pt = new TPaveText(0.15,0.7,0.8,0.85,"NDC");
      pt->SetTextSize(0.03);
      pt->SetFillColor(0);
      pt->SetTextAlign(12);
      pt->AddText(Form("Gas is %s, P is %g GeV",GasName.c_str(),CentralP));
      pt->AddText(Form("Delta P is %g GeV",DeltaP));
      pt->AddText(Form("Photocathode is %s",PCname.c_str()));
      pt->AddText(Form("Lambdamin is %g nm, lambdamax is %g nm",lambdamin,lambdamax));
      pt->Draw();

     
      sigmas->cd(2);
      ka_pr_vsAlpha->Draw("PLC PMC");
      TLine *l2=new TLine(0,3.0,PMax,3.0);
      l2->SetLineColor(kBlue);
      l2->Draw();
     
      sigmas->cd(3);
      e_pi_vsAlpha->Draw("PLC PMC");
      TLine *l3=new TLine(0,3.0,PMax,3.0);
      l3->SetLineColor(kBlue);
      l3->Draw();
  return;
}


void CRK::Simulate_Something() //Simulate_Something provides outputs of efficiencies, theta_c, Npe, etc.
{
  bool do_MomentumMC = true; // Whether or not to perform a montecarlo for cherenkov angle as a function of momentum.
 
  //Momentum MC Parameters
  // TH2 *ThetaCvsP;
  // TH2 *DeltaThCvsP;
  // TH2 *pi_ka_vsP;
  // TH2 *ka_pr_vsP;
  // TH2 *e_pi_vsP;
    
  PMax = 70; // Maximum momentum in MC in GeV 
  PMin = 0;  // Minimum momentum in MC                                                                                                                                                             
  NP = 70; // number of P points along x-axis
  Calpha = 0.005; // Error in Alpha (in radians)
  ThetaCMax = .09; // Max of y-axis in plot
  ThetaCvsP = new TH2D("ThetaC vs. P","ThetaC vs. P",NP,0.0,PMax,500,0,ThetaCMax);
  // DeltaThCvsP = new TH2D("DeltaThetaC vs. P","DeltaThetaC vs. P",NP,0.0,PMax,50,0,.15*ThetaCMax);
  pi_ka_vsP = new TH2D("N Sigma Pi-Ka vs. P","N Sigma Pi-Ka vs. P",NP,0.0,PMax,60,0,20);
  ka_pr_vsP = new TH2D("N Sigma Ka-Pr vs. P","N Sigma Ka-Pr vs. P",NP,0.0,PMax,60,0,20);
  e_pi_vsP = new TH2D("N Sigma e-Pi vs. P","N Sigma e-Pi vs. P",NP,0.0,PMax,60,0,20);
  double Pincrement = PMax/NP;
  //

  
  bool do_AlphaMC = false; // Whether or not to perform a montecarlo for cherenkov angle as a function of alpha. If do_MomentumMC is also true, plot of ThetaC vs. Alpha will include all momenta.
  
  
  //Alpha MC Parameters
  
  // TH2 *ThetaCvsAlpha;
  // TH2 *DeltaThCvsAlpha;
  // TH2 *pi_ka_vsAlpha;
  // TH2 *ka_pr_vsAlpha;
  // TH2 *e_pi_vsAlpha;
  
  //  cout << GasName << endl;
  // cout << n->getGasName();
  AlphaMax = 0.010; // radians
  NAlpha = 100; // number of alpha points along x-axis
  CentralP = 30; // Momentum value where the alpha MC is run at
  DeltaP = .01*CentralP; // momentum resolution
  pi_ka_vsAlpha = new TH2D("N Sigma Pi-Ka vs. Alpha","N Sigma Pi-Ka vs. Alpha",NAlpha,0.0,AlphaMax,500,0,20);
  ka_pr_vsAlpha = new TH2D("N Sigma Ka-Pr vs. Alpha","N Sigma Ka-Pr vs. Alpha",NAlpha,0.0,AlphaMax,500,0,20);
  e_pi_vsAlpha = new TH2D("N Sigma e-Pi vs. Alpha","N Sigma e-Pi vs. Alpha",NAlpha,0.0,AlphaMax,500,0,20);
  ThetaCvsAlpha = new TH2D("ThetaC vs. Alpha","ThetaC vs. Alpha",NAlpha,0.0,AlphaMax,500,0,ThetaCMax);
  DeltaThCvsAlpha = new TH2D("DeltaThetaC vs. Alpha","ThetaC vs. Alpha",NAlpha,0.0,AlphaMax,500,0,ThetaCMax);
  if(do_AlphaMC && !do_MomentumMC)
    {
      PMax = CentralP + DeltaP;
      PMin = CentralP - DeltaP;
      Pincrement = DeltaP;
    }
  
  
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  cout << "Gas is " << GasName << endl;
  Photocathode_Eff PCE(lambda);
  string PCname = PCE.PCName();
  // initialize integration parameters
  double lambdamin = PCE.LMin(); // Pulling mimumum wavelength observed by photocathode
  double lambdamax = PCE.LMax(); // Pulling maximum wavelength observed by photocathode
  double stepsize = 1; // step size for numerical integration 
   
  // cout << "PCname is " << PCname << lambdamin << lambdamax << endl;
   
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

	  //Procedure for skipping if below cherenkov threshold, saves a TON of time
	  double ThCPrelim = acos(1/(n->n((lambdamax+lambdamin)/2)*beta)); // check if below cherenkov threshold
	  cout << PIDname [k] << " ThCPrelim is " << ThCPrelim << endl;
	  if(isnan(ThCPrelim))
	    {
	      continue;
	    }
	  //
	  Normalize();
	  dNdx();
	  AverageThetaC(beta);
	 
	  //Outputs for non-MC option
	  if (!do_AlphaMC && !do_MomentumMC)
	    {
	      cout << "Average Theta_C detected = "<< avgthetac << endl; // Returns NaN if below cherenkov limit
	      cout << "Average photon wavelength detected = " << summedWL/nrolls << endl;
	    }
	  //
	  cout << "Mult_E OUTSIDE is " << Mult_E << endl;
	  AlphaEM = 1/137.035999;
	  dNdxValue = 2*M_PI*AlphaEM*pow(sin(avgthetac),2)*dNdxIntegral; // number of photons produced AND measured per unit length
	  cout << "dNdxValue is " << dNdxValue<<endl;
	  Npe = dNdxValue*L*10000000; // large number is there to convert from nm to cm;
	  cout << "Npe = " << floor(Npe) << endl;
	  
	  //Print Values
	  if (!do_AlphaMC && !do_MomentumMC)
	    {
	      cout << "Number of photons detected between "<< lambdamin << " nm and " << lambdamax <<" nm = " << Npe << endl;
	      cout << "Total Photon Efficiency = " << Weighted_Efficiency << endl;
	    }
 
	  
	  TRandom Randy;
	  
	  // Monte Carlo for determining spread of ThetaC with respect to an error in tracking angle (alpha) and momentum resolution
	  if (do_AlphaMC)
	    {
	       double alpha = 0;
	       int NAlphaMC = 10000; // number of points in MC
	       for (int i=0; i<=NAlpha; i++)
		 {
		     meanThC = 0;
		     sigmaThC = 0;
		   alpha = ThetaCvsAlpha->GetXaxis()->GetBinCenter(i)+AlphaMax/(2*NAlpha);
		   cout << "Running Alpha Monte Carlo for " << PIDname[k] << " at alpha = " << alpha << " rad" << endl;
		   //cout << "alpha = "<< alpha << endl;
		   for (int j = 0; j<NAlphaMC; j++)
		     {
		       //Npe = 18;
		       int WLindex = dist(generator);
		       double WL = Lambdas[WLindex];
		       int photons = Randy.Poisson(floor(Npe));
		       double sum = 0;
		       double sumsigma = 0;
		       double ThetaC = acos(1/((n->n(WL))*beta));
		       for (int f=0; f<photons; f++)
			 {
			   double phi = 2*TMath::Pi()*Randy.Rndm();
			   double ThCreco = acos(sin(alpha)*sin(ThetaC)*cos(phi) + cos(alpha)*cos(ThetaC)); 
			   sum += ThCreco;
			   sumsigma += pow(abs(avgthetac - ThCreco),2);
			   //cout << sumsigma << endl;	  
			}
		      if (photons > 0)
			{
			  double thetaApparent = sum/(double)photons;
			  
			  //cout << "thetaApparent = " << thetaApparent << endl;
    
			  sigmaThC += pow(abs(sumsigma/photons),.5); //sum and divide by N outside photons loop to find mean sigma at a certain P
			  meanThC += thetaApparent;
			  // cout << meanThC << endl;		  
			  ThetaCvsAlpha->Fill(alpha,thetaApparent);
			  //	  DeltaThCvsAlpha->Fill(alpha, sigmaThC);
			}
		     }

		   AssignPIDAlpha(k,alpha);


		 }
	       
	    }
	  if(do_MomentumMC)   // Monte Carlo for determining spread of ThetaC with respect to an error in tracking angle (alpha)
	    {
	      cout << "Running Momentum Monte Carlo for " << PIDname[k] << " at " << p << " GeV" << endl;		 
	      NpMC = 1000;
	      meanThC = 0;
	      sigmaThC = 0;	 
	      // cout << Npe << endl;
	       // for (int i=1; i<=ThetaCvsP->GetNbinsX(); i++)
	       //{
		  for (int j = 0; j<NpMC; j++)
		    {
		    
		      int WLindex = dist(generator);
		      double WL = Lambdas[WLindex];
		      //Npe = 200;
		      int photons = Randy.Poisson(floor(Npe));
		      double sum = 0;
		      double sumsigma = 0;
		      double ThetaC = acos(1/((n->n(WL))*beta));
		      for (int f=0; f<photons; f++)
			{
			  double phi = 2*TMath::Pi()*Randy.Rndm();
			  double ThCreco = acos(sin(Calpha)*sin(ThetaC)*cos(phi) + cos(Calpha)*cos(ThetaC));
			  sum += ThCreco;
			  sumsigma += pow(abs(avgthetac-ThCreco),2);
			  // cout << sumsigma << endl;
			
			}
		      if (photons > 0)
			{
			  double thetaApparent = sum/(double)photons;
			  // cout << "Theta Apparent = " << thetaApparent << endl;
			  sigmaThC += pow(abs(sumsigma/photons),.5); //sum and divide by N outside photons loop to find mean sigma at a certain P
			  meanThC += thetaApparent; // sum and divide by N at the end to find mean apparent theta at a P
			  //cout << PIDname[k] << "Theta_C std dev is: " << sigmaThC << endl;
			   // cout << "Theta Apparent = " << thetaApparent << endl;
			  ThetaCvsP->Fill(p,thetaApparent);
			  // DeltaThCvsP->Fill(p,pow(abs(sumsigma/photons),.5));
			}
		      
		    } 
		  AssignPIDMomentum(k,p);
	    }
	  
	}      
    }
  // Produce Momentum MC Plots
  if (do_MomentumMC)
    {
      MomentumMCPlots();
     
      }
  //

  // Produce Alpha MC Plots
    if (do_AlphaMC)
    {
      AlphaMCPlots();
      
    }
    // 
}
