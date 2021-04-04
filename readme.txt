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

TL;DR:

The primary file to edit is RUNME.C, What to simulate, and with what efficiencies/smears is decided there. In CRK.C, one can decide things like how many momenta/alpha points to simulate, etc.
Functions.h includes all the functions necessary to run a CRK.

To run a simulation, do "root -l RUNME.C"


More detailed how to: CRK takes as inputs CSV files of efficiencies, radiator properties, or smears with respect to wavelength.
Lower efficiencies decrease the number of detected photons, which result in a less accurate measurement.
Smears distort quantities such as Theta_C in a manner predetermined by the user in RUNME.C.

The default detector configuration is the BNL/SBU 2015 testbeam setup, with 1-bar CF4 as the radiator.

To input the indicies of refraction of your materials, open Sellmeier_Index.h, and edit the sellmeier coefficients


 
