#include "Eff.h"
#include "TH2D.h"

Eff::Eff()
{
  WavelengthDependence = 0;
}

TH1* Eff::EffPlot()
{
  if(!WavelengthDependence)
    {
      FillWavelengthDependence();
    }
  return WavelengthDependence;
}

void Eff::FillWavelengthDependence()
{
  WavelengthDependence = new TH1D(Name.c_str(),Name.c_str(), 301, 99.5 , 300.5);

  for (int bin = 1; bin <= WavelengthDependence->GetNbinsX(); bin++)
    {
      double lambda = WavelengthDependence->GetBinCenter(bin);
      WavelengthDependence->SetBinContent(bin,e(lambda));
    }
}
