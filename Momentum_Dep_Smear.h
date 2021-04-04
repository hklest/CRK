#ifndef __MOMENTUM_DEP_SMEAR_H__
#define __MOMENTUM_DEP_SMEAR_H__

#include "Smear.h"
#include "TRandom.h"

class Momentum_Dep_Smear : public Smear
{
public:
  Momentum_Dep_Smear() {}
  virtual ~Momentum_Dep_Smear() {}
  //
  // Smears the cherenkov angle based on th bending of the track in a transverse magnetic field. P_t = .3*B*R,
  double smr(double ThetaCherenkov, double lambda, double p, double B, double r_in, double r_out, double ringphi, double rapidity)
  {
    double L = r_out-r_in;
    //    double polar_angle = 2*atan2(exp(-rapidity),1);
    //double p_t = p*cos(polar_angle - M_PI/2);
    //    cout << "p_t is " << p_t << endl;
    // cout << "polar angle is " << polar_angle*180/M_PI << endl;
    //cout << "r_in is " << r_in << "r_out is " << r_out << endl;
    double radius = p/(.3*B); //radius of bending
    // cout << "radius of bending is " << radius << endl;
    // cout << "2*pow(radius,2) - pow(r_out,2)/(2*pow(radius,2)) is " << (2*pow(radius,2) - 2*pow(r_out,2))/(2*pow(radius,2)) << endl;
    //cout << "2*pow(radius,2) - pow(r_in,2)/(2*pow(radius,2)) is " << (2*pow(radius,2) - 2*pow(r_in,2))/(2*pow(radius,2)) << endl;
    double phi = acos((2*pow(radius,2) - pow(.01*r_out,2))/(2*pow(radius,2)))-acos((2*pow(radius,2)-pow(.01*r_in,2))/(2*pow(radius,2)));
    //cout << "bending phi subtended is " << phi << endl;
    // double psmear = cos(.01*L/(2*radius*sin(.3*B*.01*L/p)));
    double psmear = phi/2;
    //    cout << "2*asin(.01*L/(2*(p/(.3*B)))) is " << 2*asin(.01*L/(2*(p/(.3*B)))) << endl;
    //cout << "L,B,p,ThC is " << L <<" "<< B <<" "<< p << " " << ThetaCherenkov <<" "<< "psmear is " << psmear << endl;
    //    cout << "psmear is " << psmear << endl;
    double ThCSmearedUniform = ThetaCherenkov + cos(ringphi)*Randy.Uniform(0,psmear);
    double ThCSmearedGaussian = ThetaCherenkov + cos(ringphi)*Randy.Gaus(0,psmear);
    //    cout << "Mom. dep. gauss Smear is " << ThCSmearedGaussian - ThetaCherenkov << endl;
    // cout << "Mom. dep. uniform Smear is " << ThCSmearedUniform - ThetaCherenkov << endl; 
    return ThCSmearedUniform;
  }
  
protected:
  double smear;
  TRandom Randy;
};
	
#endif /* __MOMENTUM_DEP_SMEAR_H__ */
