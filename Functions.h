#include "Index.h"
#include "Eff.h"
//#include "Eff.C"
#include "Efftextfile.h"
//#include "Efftextfile.C"
#include "Smear.h"
#include "CRK.h"
//#include "CRK.C"
#include "Sellmeier_Index.h"
//#include "Photocathode_Eff.h"
//#include "WL_Dep.h"                                                                                                                                                                                       
#include <cmath>
#include <ctime>



CRK::CRK(Index *index, double RadiatorLength)
{
  n = index;
  L = RadiatorLength;
}

double CRK::dNdx() // Function that produces an efficiency distribution vector and produces the integral term "dNdxIntegral"                                                                               
{
  dNdxIntegral = 0;
  // cout << lambdamin << lambdamax << endl;                                                                                                                                                                
  //Loop over wavelengths between lambdamin nm and lambdamax nm to find quantities at each wavelength in steps of stepsize nm                                                                               
  for (lambda = lambdamin; lambda <= lambdamax; lambda=lambda+stepsize)
    {
      int i = 0;
      Mult_E = 1;
      while (i < effs.size())
      {
        Mult_E *= (effs[i]->e(lambda));
        //cout << Mult_E << endl;                                                                                                                                                                           
        i++;
     }
      //  cout << "Mult_E INSIDE is " << Mult_E << endl;                                                                                                                                                    
      Weighted_Efficiency += Mult_E*pow(lambda,-2)/normalization; // metric of the full survival/measurement rate of all photons                                                                            
      dNdxIntegral += Mult_E*stepsize*pow(lambda,-2); // integral for use in dN/dx calculation                                                                                                              
    }
  //cout << "dndx integral is " << dNdxIntegral<< endl;                                                                                                                                                     
  return dNdxIntegral;
}



void CRK::ProduceEffDist()
{
for (lambda = lambdamin; lambda <= lambdamax; lambda=lambda+stepsize)
  {
    //    cout << "Effs size " << effs.size() << endl;                                                                                                                                                      
    // cout << "Eff 1 " << effs[1]->e(lambda) << endl;                                                                                                                                                      
    // cout << "Effs size " << effs.size() << endl;                                                                                                                                                         
    Mult_E = 1;
    int i = 0;
    while (i < effs.size())
      {
       //      cout << "pre mult_e effs = " << effs[i]->e(lambda) << endl;                                                                                                                                 
        Mult_E *= (effs[i]->e(lambda));
        //      cout << Mult_E << endl;                                                                                                                                                                     
        i++;
      }
    //cout << "Mult_E post loop is "<< Mult_E << endl;                                                                                                                                                      

    Eff_Dist.push_back(Mult_E*pow(lambda,-2)); // create a vector of efficiencies                                                                                                                           
    Lambdas.push_back(lambda); // create a vector of lambdas                                                                                                                                                    //   for(int i=0; i< Eff_Dist.size();i++)                                                                                                                                                                  // {                                                                                                                                                                                                       //  cout << "Eff_Dist[i] = "<<Eff_Dist[i] << endl;                                                                                                                                                         //  }                                                                                                                                                                                                      // for(int i=0; i< Lambdas.size();i++)                                                                                                                                                                     // {                                                                                                                                                                                                       //  cout << "Lambdas[i] = "<<Lambdas[i] << endl;
    Eff_Dist.push_back(Mult_E*pow(lambda,-2)); // create a vector of efficiencies
    Lambdas.push_back(lambda); // create a vector of lambdas                                                                                                                                                   //   for(int i=0; i< Eff_Dist.size();i++)                                                                                                                                                                   // {                                                                                                                                                                                                       //  cout << "Eff_Dist[i] = "<<Eff_Dist[i] << endl;                                                                                                                                                         //  }                                                                                                                                                                                                      // for(int i=0; i< Lambdas.size();i++)                                                                                                                                                                     // {                                                                                                                                                                                                       //  cout << "Lambdas[i] = "<<Lambdas[i] << endl;                                                                                                                                                           //  }                                                                                                                                                                                                 
  }
}

double CRK::AverageThetaC(double beta,vector<double> Eff_Dist,vector<double> Lambdas) // Function that gives an average value of theta_C for seeding the monte carlo                                        
{
  seed = std::chrono::system_clock::now().time_since_epoch().count(); //Set a seed for the RNG so new numbers generated each time                                                                           
  std::default_random_engine LambdaGen1 (seed);  // Initialize RNG                                                                                                                                          
  std::discrete_distribution<int> LambdaDist1(Eff_Dist.begin(), Eff_Dist.end()); // returns an integer number with probabilities defined by the Eff_Dist vector, used in finding avgthetac and in MCs       

  // cout << beta<< endl;                                                                                                                                                                                   
  //Use a realistic distribution of photons to determine average theta_C detected for the detector configuration                                                                                            
  nrolls = 1000; // number of lambdas to generate                                                                                                                                                           
  summedWL = 0; // summed wavelengths                                                                                                                                                                       
  summedthetac = 0; // summed theta_c                                                                                                                                                                       

  // double SUMDEV=0;                                                                                                                                                                                       

  int A = 0;
  while (A < nrolls) // Loop for determining average wavelength and thetaC                                                                                                                                  
    {
      int WaveLindex = LambdaDist1(LambdaGen1); // Picks a random index according to LambdaDist
      double WaveL = Lambdas[WaveLindex];
      // cout << "WaveL = " << WaveL << endl;
      // cout << " n @ wavel = "<< n->n(WaveL)<<endl;// picks out the value of the wavelength at that index
      double ThetaCWL = acos(1/((n->n(WaveL))*beta)); // Finds the thetaC at that wavelength                                                                                                               
      // cout << "nbins is " << nbins << endl;
      if (!isnan(ThetaCWL) && ThetaCWL < 1)
        {
          summedWL += WaveL;
          summedthetac += ThetaCWL;
          // cout << "summedthetac is " << summedthetac << endl;                                                                                                                                            
          A++;
        }

      // NvThC->Fill(ThetaCWL);                                                                                                                                                                             
    }
  avgthetac = summedthetac/nrolls;
  // for(int i = 0; i<ThetaCs.size();i++)                                                                                                                                                                    // {                                                                                                                                                                                                       //  SUMDEV += pow(abs(ThetaCs[i]-avgthetac),2);                                                                                                                                                                // cout << "ThetaC[i] is " << ThetaCs[i] << endl;                                                                                                                                                          // cout << SUMDEV << endl;                                                                                                                                                                             //  }                                                                                                                                                                                                      // double sigmachrom = sqrt(SUMDEV/ThetaCs.size());                                                                                                                                                        // cout << "Chromaticity contribution to resolution is " << sigmachrom << endl;                                                                                                                            //cout << "avgtheta = " << avgthetac << endl;                                                                                                                                                              //cout << "avgWL = " << summedWL/nrolls << endl;
  return avgthetac;
}




void CRK::AssignWLMinMax()
{
 lambdamin = effs[0]->LMin(); // Want lambdamin to be the smallest of the LMin()'s                                                                                                                          
 for (int i = 1; i<effs.size();i++)
   {
     if(lambdamin > effs[i]->LMin())
       {
         lambdamin = effs[i]->LMin();
	 // cout << "new lambdamin is " << lambdamin << endl;
       }
   }
 lambdamax = effs[0]->LMax(); // Want Lambdamax to be largest of the LMax()'s                                                                                                                               
for (int i = 1; i<effs.size()-3;i++)
   {
     if(lambdamax < effs[i]->LMax())
       {
         lambdamax = effs[i]->LMax();
	 // cout << "new lambdamax is " << lambdamax << endl;
       }
   }
}




void CRK::MultEffPlot() // Function that produces an efficiency distribution vector and produces the integral term "dNdxIntegral"                                                                               
{
  // cout << "Lambdamin inside multeffplot is " << lambdamin << endl;
  // cout << "Lambdamax inside multeffplot is " << lambdamax << endl;
  //  lambdamin = 100;
  // lambdamax = 200;
  // cout << "stepsize inside multeffplot is " << stepsize << endl; 
  int nbins = ceil(lambdamax-lambdamin)/stepsize;
  // cout << "nbins is " << nbins << endl;
  TCanvas *MultEff = new TCanvas("MultEff","Multiplied Efficiencies",200,500,900,1200);
  MultEff->Divide(1,2);
  TH1* MultEffs = new TH1D("Multiplied Efficiencies","Multiplied Efficiencies", nbins, lambdamin, lambdamax);
  TH1* WeightEffs = new TH1D("Multiplied Efficiencies","Spectrum-Weighted Efficiencies", nbins, lambdamin, lambdamax);
  for (lambda = lambdamin; lambda <= lambdamax; lambda=lambda+stepsize)
    {
      int i = 0;
      Mult_E = 1;
      while (i < effs.size())
      {
        Mult_E *= (effs[i]->e(lambda));
	//  cout << Mult_E << endl;
	
	i++;
      }
      WeightEffs->SetBinContent(ceil(lambda-lambdamin),Mult_E/pow(lambda,2));
      MultEffs->SetBinContent(ceil(lambda-lambdamin),Mult_E);
      //   cout << "lambda is " << lambda << "Mult_E is " << Mult_E << endl;
    }
  MultEff->cd(1);
  MultEffs->Draw("LF2");
  MultEff->cd(2);
  WeightEffs->Draw("LF2");
  return;
}


double CRK::Normalize()
{
 //Find normalization constant                                                                                                                                                                             
  normalization = 0;
  // initialize integration parameters

  // lambdamin ; // Pulling mimumum wavelength observed by photocathode
 
  //  lambdamax; // Pulling maximum wavelength observed by photocathode                                                                                                                            
  stepsize = 1; // step size for numerical integration                                                                                                                                                    

    for (int i = lambdamin; i <= lambdamax; i=i+stepsize)
      {
        // cout << lambdamin << lambdamax<<endl;                                                                                                                                                            
        normalization += pow(i,-2.0);
      }
    //    cout << normalization << endl;                                                                                                                                                                    
    return normalization;
}

void CRK::AssignPIDAlpha(int k, double alpha)
{

  if(k==0)// if pion                                                                                                                                                                                        
    {
      alphas.push_back(alpha);
      // cout << " In assign pid alpha for pion, meanThC is " << meanThC << endl;
      // cout << " In assign pid alpha for pion, sigmathc is " << sigmaThC << endl;
      pi_mean = meanThC/NAlphaMC;
      pi_means.push_back(pi_mean);
      pi_photons = avgphotons/NAlphaMC;
       cout << "Pi Photons is "<< pi_photons << endl;
      cout << "pi_mean = " << pi_mean << endl;
      pi_photonsvec.push_back(pi_photons);
      pi_RMS = pow(sigmaThC,.5)/NAlphaMC;
      pi_RMSs.push_back(pi_RMS);
      cout << "pi_RMS = " << pi_RMS << endl;
    }
  if(k==1) // if kaon                                                                                                                                                                                       
    {
      //   cout << " In assign pid alpha for kaon, meanThC is " << meanThC << endl;
      // cout << " In assign pid alpha for kaon, NAlphaMC is " << NAlphaMC << endl;
      ka_mean = meanThC/NAlphaMC;
      ka_means.push_back(ka_mean);
      cout << "ka_mean = " << ka_mean << endl;
      ka_RMS = pow(sigmaThC,.5)/NAlphaMC;
      ka_RMSs.push_back(ka_RMS);
      ka_photons = avgphotons/NAlphaMC;
      cout << "Ka Photons is "<< ka_photons << endl;
       cout << "Pi Photons inside Ka loop is "<< pi_photons << endl;
      ka_photonsvec.push_back(ka_photons);
      cout << "ka_RMS = " << ka_RMS << endl;
      cout << "N sigma ka-pi = " << sqrt((pi_photons+ka_photons)/2)*abs(ka_mean-pi_mean)/((ka_RMS+pi_RMS)/2) << endl;
    
    }
  if(k==2)// if proton                                                                                                                                                                                      
    {
      pr_mean = meanThC/NAlphaMC;
      pr_means.push_back(pr_mean);
      cout << "pr_mean = " << pr_mean << endl;
      pr_RMS = pow(sigmaThC,.5)/NAlphaMC;
      pr_RMSs.push_back(pr_RMS);
      pr_photons = avgphotons/NAlphaMC;
      pr_photonsvec.push_back(pr_photons);
      cout << "pr_RMS = " << pr_RMS << endl;
      cout << "N sigma ka-pr = " << sqrt((pr_photons+ka_photons)/2)*abs(ka_mean-pr_mean)/((ka_RMS+pr_RMS)/2) << endl;
    }
if(k==3)// if electron                                                                                                                                                                                    
    {
      e_mean = meanThC/NAlphaMC;
      e_means.push_back(e_mean);
      cout << "e_mean = " << e_mean << endl;
      e_RMS = pow(sigmaThC,.5)/NAlphaMC;
      e_photons=avgphotons/NAlphaMC;
      e_RMSs.push_back(e_RMS);
      e_photonsvec.push_back(e_photons);
      cout << "e_RMS = " << e_RMS << endl;
      cout << "N sigma e-pi = " <<sqrt((pi_photons+e_photons)/2)*abs(e_mean-pi_mean)/((pi_RMS+e_RMS)/2) << endl;
	
     
    }
  cout << "Moving on....." << endl << endl;

 return;
 }
void CRK::FillAlphas()
{
  cout <<"mean sizes are " << pi_means.size() << ka_means.size() << pr_means.size() << e_means.size() << endl;
  cout <<"RMS sizes are " << pi_RMSs.size() << ka_RMSs.size() << pr_RMSs.size() << e_RMSs.size() << endl;
  cout <<"photonsvec sizes are " << pi_photonsvec.size() << ka_photonsvec.size() << pr_photonsvec.size() << e_photonsvec.size() << endl;
  // cout <<"mean sizes are " << pi_means.size() << ka_means.size() << pr_means.size() << e_means.size() << endl;
  for(int m = 0; m <=pi_means.size(); m++)
    {
      // if(!isnan(pi_RMSs.at(m)+ka_RMSs.at(m)) && !isnan(pi_RMSs.at(m)+ka_RMSs.at(m)) && !isnan(pr_RMSs.at(m)+ka_RMSs.at(m)) && !isnan(pi_RMSs.at(m)+e_RMSs.at(m)))
      //	{
      // cout << "N sigma pi-k is " << (sqrt((pi_photonsvec[m]+ka_photonsvec[m])/2)*abs(pi_means[m]-ka_means[m])/((pi_RMSs[m]+ka_RMSs[m])/2)) << "at alpha " << alphas[m] << endl;
      if(pi_means.size()==ka_means.size() && ka_means.size() != 0)
	{
	  pi_ka_vsAlpha->Fill(alphas[m],(sqrt((pi_photonsvec[m]+ka_photonsvec[m])/2)*abs(pi_means[m]-ka_means[m])/((pi_RMSs[m]+ka_RMSs[m])/2)));
	}
      if(pr_means.size()==ka_means.size() && ka_means.size() != 0)
	{
	  ka_pr_vsAlpha->Fill(alphas[m],(sqrt((pr_photonsvec[m]+ka_photonsvec[m])/2)*abs(pr_means[m]-ka_means[m])/((pr_RMSs[m]+ka_RMSs[m])/2)));
	}
      if(pi_means.size()==e_means.size() &&  pi_means.size() != 0)
	{
	  e_pi_vsAlpha->Fill(alphas[m],(sqrt((e_photonsvec[m]+pi_photonsvec[m])/2)*abs(pi_means[m]-e_means[m])/((pi_RMSs[m]+e_RMSs[m])/2)));
	}
      //	}
    }
  
}

void CRK::AddEff(Eff *e)
{
  effs.push_back(e);
}

void CRK::AddSmear(Smear *s)
{
  smears.push_back(s);
}

void CRK::AssignPIDMomentum(int k, double p)
{
 if(k==0)// if pion                                                                                                                                                                           
   {
     pi_mean = meanThC/NpMC;
     cout << "pi_mean = " << pi_mean << endl;
     pi_photons = avgphotons/NpMC;
     cout << "pi_photons " << pi_photons<< endl;     
     pi_RMS = sigmaThC/NpMC;
     cout << "pi_RMS = " << pi_RMS << endl;
   }
 if(k==1) // if kaon                                                                                                                                                                          
   {
     ka_mean = meanThC/NpMC;
     cout << "ka_mean = " << ka_mean << endl;
     ka_RMS =  sigmaThC/NpMC;
     ka_photons = avgphotons/NpMC;
     cout << "ka_photons " << ka_photons<< endl;
     cout << "ka_RMS = " << ka_RMS << endl;
     //  cout << "pi mean in k loop = " << pi_mean << endl;
     pi_ka_vsP->Fill(p,sqrt((pi_photons+ka_photons)/2)*((abs(pi_mean-ka_mean)/((pi_RMS+ka_RMS)/2))));
     cout << "N sigma k-pi = " << sqrt((pi_photons+ka_photons)/2)*(abs(ka_mean-pi_mean)/((pi_RMS+ka_RMS)/2)) << endl;
   }
 if(k==2)// if proton                                                                                                                                                                         
   {
     pr_mean = meanThC/NpMC;
     cout << "pr_mean = " << pr_mean << endl;
     pr_RMS =  sigmaThC/NpMC;
     cout << "pr_RMS = " << pr_RMS << endl;
     pr_photons = avgphotons/NpMC;
     cout << "pr_photons " << pr_photons<< endl;
     
     ka_pr_vsP->Fill(p,sqrt((pr_photons+ka_photons)/2)*(abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)));
     cout << "N sigma k-pr = " << sqrt((pr_photons+ka_photons)/2)*(abs(ka_mean-pr_mean)/((pr_RMS+ka_RMS)/2)) << endl;
   }
 if(k==3)// if electron                                                                                                                                                                       
   {
     e_mean = meanThC/NpMC;
     cout << "e_mean = " << e_mean << endl;
     e_RMS =  sigmaThC/NpMC;
     cout << "e_RMS = " << e_RMS << endl;
     e_photons = avgphotons/NpMC;
     cout << "e_photons " << e_photons<< endl;
     if(pi_mean >0)
       {
	 e_pi_vsP->Fill(p,sqrt((pi_photons+e_photons)/2)*(abs(e_mean-pi_mean)/((e_RMS+pi_RMS)/2)));
       }
     cout << "N sigma e-pi = " << sqrt((pi_photons+e_photons)/2)*(abs(pi_mean-e_mean)/((pi_RMS+e_RMS)/2)) << endl;
   }
 cout << "Moving on....." << endl << endl;


 return;
}


void CRK::MomentumMCPlots()
{
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  string PCName = effs[0]->EffName();//Name.c_str();
  gStyle->SetOptStat(0);
  TCanvas *momenta = new TCanvas("momenta","ThetaC vs. P",100,300,800,800);
  //c1->Divide(1,2);                                                                                                                                                                              
      momenta->cd(1);
     
      ThetaCvsP->Draw("colz");
      momenta->cd(2);
      TCanvas *sigmasP = new TCanvas("sigmasP", "Sigmas",400,100,1300,800);
    
      pi_ka_vsP->SetMarkerStyle(kPlus);
      ka_pr_vsP->SetMarkerStyle(kPlus);
      e_pi_vsP->SetMarkerStyle(kPlus);
      sigmasP->Divide(3,1);
      sigmasP->cd(1);
      gPad->SetLogy(1);
      pi_ka_vsP->Draw("PLC PMC");
      TLine *l1p=new TLine(0,3.0,PMax,3.0);
      l1p->SetLineColor(kBlue);
      l1p->Draw();
      TPaveText *pt = new TPaveText(0.15,0.7,0.87,0.85,"NDC");
      pt->SetTextSize(0.03);
      pt->SetFillColor(0);
      pt->SetTextAlign(12);
      pt->AddText(Form("Gas is %s, Alpha is %g rad",GasName.c_str(),Calpha));
      pt->AddText(Form("Photocathode efficiency from %s",PCName.c_str()));
      pt->AddText(Form("Radiator length is %g cm",L));
      pt->AddText(Form("Detector inner radius is %g cm, outer radius is %g cm",r_in,r_out));
      pt->AddText(Form("B-field is %g T",B));
      pt->AddText(Form("Pixel size is %g x %g",pixel_size,pixel_size));
      pt->AddText(Form("Pressure %g Bar absolute",pressure));
      // pt->AddText(Form("Lambdamin is %g nm ",lambdamin));
      // pt->AddText(Form("Lambdamax is %g nm ",lambdamax));
      // pt->AddText(Form("Theta_C smear is %g",smears[3]->smr(1)));
      // pt->AddText(Form("Number of stray photons is %d",NScintPhotons));
      if(floor(Npe)==Npe)
	{
	  pt->AddText(Form("Npe is forced to %g",Npe));
	}
      pt->Draw();
      //pi_ka_vsP->GetXaxis()->SetLabelSize(5.2323);                                                                                                                                                        
      // pi_ka_vsP->GetXaxis()->SetTitle("GeV");                                                                                                                                                            
      // pi_ka_vsP->GetXaxis()->SetTitleOffset(1);                                                                                                                                                          
      // pi_ka_vsP->GetYaxis()->SetTitle("Number of Sigma");                                                                                                                                                
      // pi_ka_vsP->GetYaxis()->SetTitleSize(80);                                                                                                                                                           


      sigmasP->cd(2);
      gPad->SetLogy(1);
      ka_pr_vsP->Draw("PLC PMC");
      TLine *l2p=new TLine(0,3.0,PMax,3.0);
      l2p->SetLineColor(kBlue);
      l2p->Draw();

      sigmasP->cd(3);
      gPad->SetLogy(1);
      e_pi_vsP->Draw("PLC PMC");
      TLine *l3p=new TLine(0,3.0,PMax,3.0);
      l3p->SetLineColor(kBlue);
      l3p->Draw();

  return;
}


void CRK::AlphaMCPlots()
{
  Sellmeier_Index SMn(lambda);
  string GasName = SMn.Gas();
  Photocathode_Eff PCE(lambda);
  string PCname = PCE.PCName();
  gStyle->SetOptStat(0);
  TCanvas *Alphas = new TCanvas("Alphas","Theta_C vs. Alpha",100,300,800,1000);
  // c1->Divide(1,2);
  // c1->cd(1);
  ThetaCvsAlpha->Draw("colz");
  // c1->cd(2);
  //DeltaThCvsAlpha->Draw("colz");
  TCanvas *sigmasAlpha = new TCanvas("sigmasAlpha", "Sigmas",500,800,1200,800);
  pi_ka_vsAlpha->SetMarkerStyle(kPlus);
  ka_pr_vsAlpha->SetMarkerStyle(kPlus);
  e_pi_vsAlpha->SetMarkerStyle(kPlus);
  sigmasAlpha ->Divide(3,1);
  sigmasAlpha->cd(1);
  // gPad->SetLogy(1);
  pi_ka_vsAlpha->Draw("PLC PMC");
  TLine *l1=new TLine(0,3.0,AlphaMax,3.0);
  l1->SetLineColor(kBlue);
  l1->Draw();
  TPaveText *pt = new TPaveText(0.15,0.68,0.8,0.85,"NDC");
  pt->SetTextSize(0.03);
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->AddText(Form("Gas is %s, Momentum is %g GeV",GasName.c_str(),p));
  // pt->AddText(Form("Delta p is %g GeV",DeltaP));
  pt->AddText(Form("Photocathode is %s",PCname.c_str()));
  pt->AddText(Form("Radiator length is %g cm",L));
  //  pt->AddText(Form("Lambdamin is %g nm, lambdamax is %g nm",lambdamin,lambdamax));
  // pt->AddText(Form("Number of stray photons is %d",NScintPhotons));
  if(floor(Npe)==Npe)
    {
      pt->AddText(Form("Npe is forced to %g",Npe));
    }
  pt->Draw();
  
  
  sigmasAlpha->cd(2);
  // gPad->SetLogy(1);
  ka_pr_vsAlpha->Draw("PLC PMC");
  TLine *l2=new TLine(0,3.0,AlphaMax,3.0);
  l2->SetLineColor(kBlue);
  l2->Draw();
  
  sigmasAlpha->cd(3);
  //  gPad->SetLogy(1);
  e_pi_vsAlpha->Draw("PLC PMC");
  TLine *l3=new TLine(0,3.0,AlphaMax,3.0);
  l3->SetLineColor(kBlue);
  l3->Draw();
  return;
}
