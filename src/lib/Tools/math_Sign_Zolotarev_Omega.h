/*!
        @file    math_Sign_Zolotarev_Omega.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#ifndef MATH_SIGN_ZOLOTAREV_OMEGA_INCLUDED
#define MATH_SIGN_ZOLOTAREV_OMEGA_INCLUDED

#include <vector>

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

//! Determination of Zolotarev coefficients.
/*!
   Temporary: for optimal domain-wall.
                                    [22 Jan 2012 H.Matsufuru]

   This class determines the Zolotarev's optimal coefficients
   of partial fractional approximation to 1/sqrt(x).
   Present implementation makes use of the code in Numerical
   Recipes, and thus cannot be put public.
   To be replaced with public implementation.
 */
class Math_Sign_Zolotarev_Omega{

 private:
  int    d_Ns;
  double d_bmax;
  double d_D0;
  double d_delta;
  std::vector<double> d_omega;
  std::vector<double> d_cl;

 protected:
  Bridge::VerboseLevel m_vl;

 public:
  Math_Sign_Zolotarev_Omega(int Ns, double bmax)
    : m_vl(CommonParameters::Vlevel())
  {
    d_Ns = Ns;
    d_bmax = bmax;
    set_sign_prms();
  }

  void get_sign_prms(std::vector<double>&, double&);

  double sign(double x);

  double sign_cl(double x);

 private:

  void set_sign_prms();

  void rescale_cl();

  void set_omega();

  double binary_search(double, double);

  void poly_Zolotarev(double bmax, double& UK);

  void Jacobi_elliptic(double uu, double emmc,
	       double &sn, double &cn, double &dn);

};

#endif
