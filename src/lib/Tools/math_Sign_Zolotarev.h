/*!
        @file    math_Sign_Zolotarev.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-02-14 17:05:33 #$

        @version $LastChangedRevision: 1944 $
*/

#ifndef MATH_SIGN_ZOLOTAREV_INCLUDED
#define MATH_SIGN_ZOLOTAREV_INCLUDED

#include <cassert>

#include "bridge_defs.h"
#include "Parameters/commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Determination of Zolotarev coefficients.

/*!
    This class determines the Zolotarev's optimal coefficients
    of partial fractional approximation to 1/sqrt(x).
    The present implementation makes use of the code in Numerical
    Recipes, and thus cannot be public. To be replaced.
                                     [28 Dec 2011 H.Matsufuru]
    Replaced with a GSL based code.  [19 Jun 2013 S.Ueda]
 */

class Math_Sign_Zolotarev
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_Np;
  double m_bmax;
  std::vector<double> m_cl;
  std::vector<double> m_bl;

 public:
  Math_Sign_Zolotarev(const int Np, const double bmax)
    : m_vl(CommonParameters::Vlevel()),
    m_Np(Np), m_bmax(bmax)
  {
    set_sign_parameters();
  }

  void get_sign_parameters(std::vector<double>& cl, std::vector<double>& bl);

  double sign(const double x);

 private:
  void set_sign_parameters();
  void poly_Zolotarev(const double bmax, double& UK);
  void Jacobi_elliptic(const double uu, const double emmc,
                       double& sn, double& cn, double& dn);
};
#endif
