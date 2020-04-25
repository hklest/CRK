===================================                                  
                             ,--. 
  ,----..  ,-.----.      ,--/  /| 
 /   /   \ \    /  \  ,---,': / ' 
|   :     :;   :    \ :   : '/ /  
.   |  ;. /|   | .\ : |   '   ,   
.   ; /--` .   : |: | '   |  /    
;   | ;    |   |  \ : |   ;  ;    
|   : |    |   : .  / :   '   \   
.   | '___ ;   | |  \ |   |    '  
'   ; : .'||   | ;\  \'   : |.  \ 
'   | '/  ::   ' | \.'|   | '_\.' 
|   :    / :   : :-'  '   : |     
 \   \ .'  |   |.'    ;   |,'     
  `---`    `---'      '---'       
                                  
===================================

By Tom Hemmick and Henry Klest

CRK is a fast simulation framework for simulating efficiencies and various other quantities and properties of ring-imaging cherenkov detectors.

How to: CRK takes as inputs CSV files of efficiencies, radiator properties, or smears with respect to wavelength.
Lower efficiencies decrease the number of detected photons, which result in a less accurate measurement.
Smears distort quantities such as Theta_C in a manner predetermined by the user.

The default detector configuration is the BNL/SBU 2015 testbeam setup, with 1-bar CF4 as the radiator.

To run, run RUNME.C

To input other efficiencies, create a new header file "name_Eff.h" in the same format as the other efficiencies, edit RUNME.C by adding a line in the format of  "Eff   *e4  = new name_Eff(1);"
and "detector->AddEff(e4);"

To input the indicies of refraction of your materials, open Sellmeier_Index.h, and edit the sellmeier coefficients

To play around simulating various other quantities, edit the Simulate_Something function in CRK.C

Simulate_Something currently contains two monte carlo routines (do_ThetaMC and do_MomentumMC) that produces a plot of ThetaC as a function of the error in angle
from a tracking system (known here as alpha).

 
