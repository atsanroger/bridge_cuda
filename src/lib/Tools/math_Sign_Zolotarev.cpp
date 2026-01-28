/*!
        @file    math_Sign_Zolotarev.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2020-11-20 20:26:10 #$

        @version $LastChangedRevision: 2123 $
*/

#include "math_Sign_Zolotarev.h"

#include <cstdlib>
#include <cfloat>

using Bridge::vout;
using std::abs;

const std::string Math_Sign_Zolotarev::class_name = "Math_Sign_Zolotarev";

//====================================================================
void Math_Sign_Zolotarev::get_sign_parameters(std::vector<double>& cl,
                                              std::vector<double>& bl)
{
  //- Zolotarev coefficient defined

  if (cl.size() != 2 * m_Np) {
    vout.crucial(m_vl, "Error at %s: size of cl is not correct\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if (bl.size() != m_Np) {
    vout.crucial(m_vl, "Error at %s: size of bl is not correct\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < 2 * m_Np; ++i) {
    cl[i] = m_cl[i];
  }

  for (int i = 0; i < m_Np; ++i) {
    bl[i] = m_bl[i];
  }
}


//====================================================================
void Math_Sign_Zolotarev::set_sign_parameters()
{
  // Zolotarev coefficient defined
  m_cl.resize(2 * m_Np);
  m_bl.resize(m_Np);

  double UK = 0.0;
  vout.general(m_vl, " bmax = %12.4e\n", m_bmax);
  vout.general(m_vl, " UK   = %12.4e\n", UK);

  poly_Zolotarev(m_bmax, UK);

  //  for(int i=0; i<m_Np; i++){
  //    vout.general(m_vl, " %3d %12.4e %12.4e %12.4e\n", i, m_cl[i],
  //                                      m_cl[i+m_Np], m_bl[i]);
  //  }
}


//====================================================================
double Math_Sign_Zolotarev::sign(const double x)
{
  //-  cl[2*Np], bl[Np]: coefficients of rational approx.

  double x2R = 0.0;

  for (int l = 0; l < m_Np; l++) {
    x2R += m_bl[l] / (x * x + m_cl[2 * l]);
  }
  x2R = x2R * (x * x + m_cl[2 * m_Np - 1]);

  return x * x2R;
}


//====================================================================
void Math_Sign_Zolotarev::poly_Zolotarev(const double bmax, double& UK)
{
//    Return the coefficients of Zolotarev's approximation
//    of sign function to rational function.
//
//    Np: number of poles (2*Np is number of cl)
//    bmax: range of argument [1,bmax]
//    UK: complete elliptic integral of the 1st kind
//    cl[2*Np],bl[Np]: coefficients of Zolotarev rational approx.

  const int Nprec = 14;

  const double rk   = sqrt(1.0 - 1.0 / (bmax * bmax));
  const double emmc = 1.0 - rk * rk;

  vout.general(m_vl, "emmc = %22.14e\n", emmc);

  double sn, cn, dn;

  //- Determination of K
  double Dsr = 10.0;
  double u   = 0.0;
  for (int i_prec = 0; i_prec < Nprec + 1; i_prec++) {
    Dsr = Dsr * 0.1;

    for (int i = 0; i < 20; i++) {
      u = u + Dsr;
      Jacobi_elliptic(u, emmc, sn, cn, dn);
      vout.detailed(m_vl, " %22.14e %22.14e %22.14e\n", u, sn, cn);
      if (cn < 0.0) goto succeeded;
    }
    vout.general(m_vl, "Something wrong in setting Zolotarev\n");
    return;

succeeded:
    u = u - Dsr;
  }

  UK = u;

  Jacobi_elliptic(UK, emmc, sn, cn, dn);
  vout.general(m_vl, " %22.14e %22.14e %22.14e\n", UK, sn, cn);


  //- Determination of c_l
  const double FK = UK / (2.0 * m_Np + 1.0);

  for (int l = 0; l < 2 * m_Np; l++) {
    u = FK * (l + 1.0);
    Jacobi_elliptic(u, emmc, sn, cn, dn);
    m_cl[l] = sn * sn / (1.0 - sn * sn);
  }

  //- Determination of b_l
  double d0 = 1.0;
  for (int l = 0; l < m_Np; l++) {
    d0 = d0 * (1.0 + m_cl[2 * l]) / (1.0 + m_cl[2 * l + 1]);
  }

  for (int l = 0; l < m_Np; l++) {
    m_bl[l] = d0;
    for (int i = 0; i < m_Np; i++) {
      if (i < m_Np - 1) m_bl[l] = m_bl[l] * (m_cl[2 * i + 1] - m_cl[2 * l]);
      if (i != l) m_bl[l] = m_bl[l] / (m_cl[2 * i] - m_cl[2 * l]);
    }
  }

  //- Correction
  const int    Nj = 10000;
  const double Dx = (bmax - 1.0) / Nj;

  double d1 = 0.0;
  double d2 = 2.0;

  for (int jx = 0; jx <= Nj; jx++) {
    double x    = 1.0 + Dx * jx;
    double sgnx = sign(x);
    if (fabs(sgnx) > d1) d1 = sgnx;
    if (fabs(sgnx) < d2) d2 = sgnx;
  }
  vout.general(m_vl, " |sgnx|_up  = %8.4e \n", d1 - 1.0);
  vout.general(m_vl, " |sgnx|_dn  = %8.4e \n", 1.0 - d2);

  /*
  double d0c = 0.5*(d1+d2);
  for(int ip=0; ip<m_Np; ip++){
    m_bl[ip] = m_bl[ip]/d0c;
  };

  d1 = 0.0;
  d2 = 2.0;
  for(int jx=0; jx <= Nj; jx++){
    double x = 1.0 + Dx*jx;
    double sgnx = sign(x);
    if(fabs(sgnx) > d1) d1 = sgnx;
    if(fabs(sgnx) < d2) d2 = sgnx;
  }
  vout.general(m_vl, " |sgnx|_up  = %8.4e \n", d1-1.0);
  vout.general(m_vl, " |sgnx|_dn  = %8.4e \n", 1.0-d2);
  */

  double dmax = d1 - 1.0;
  if (dmax < 1.0 - d2) dmax = 1.0 - d2;
  vout.general(m_vl, " Amp(1-|sgn|) = %16.8e \n", dmax);
}


//====================================================================
void Math_Sign_Zolotarev::Jacobi_elliptic(const double uu, const double emmc,
                                          double& sn, double& cn, double& dn)
{
//  Return the Jacobi elliptic functions sn(u,kc),
//  cn(u,kc), and dn(u,kc), where kc = 1 - k^2, and
//  arguments:
//    uu: u
//    emmc: kc.
//
//  This routine is based on GSL. [19 Jun 2013 S.Ueda]

  // The accuracy
  const double epsilon = DBL_EPSILON;
  // Max number of iteration
  const int N = 16;

  const double emmc1 = 1 - emmc;

  if (abs(emmc1) > 1.0) {
    sn = cn = dn = 0;
    vout.crucial("Error at %s: parameter m = %f must be smaller than 1.\n", class_name.c_str(), emmc1);
    exit(EXIT_FAILURE);
  } else if (abs(emmc1) < 2.0 * epsilon) {
    sn = sin(uu);
    cn = cos(uu);
    dn = 1.0;
  } else if (abs(emmc1 - 1) < 2.0 * epsilon) {
    sn = tanh(uu);
    cn = 1.0 / cosh(uu);
    dn = cn;
  } else {
    double mu[N];
    double nu[N];

    int n = 0;

    mu[0] = 1.0;
    nu[0] = sqrt(1.0 - emmc1);

    while (abs(mu[n] - nu[n]) > epsilon * abs(mu[n] + nu[n]))
    {
      mu[n + 1] = 0.5 * (mu[n] + nu[n]);
      nu[n + 1] = sqrt(mu[n] * nu[n]);
      ++n;
      if (n >= N - 1) {
        vout.general(m_vl, "max iteration\n");
        break;
      }
    }

    double sin_umu = sin(uu * mu[n]);
    double cos_umu = cos(uu * mu[n]);

    /* Since sin(u*mu(n)) can be zero we switch to computing sn(K-u),
       cn(K-u), dn(K-u) when |sin| < |cos| */
    if (fabs(sin_umu) < fabs(cos_umu)) {
      double cot = sin_umu / cos_umu;
      double cs  = mu[n] * cot;
      dn = 1.0;

      while (n > 0)
      {
        --n;
        double r = cs * cs / mu[n + 1];
        cs *= dn;
        dn  = (r + nu[n]) / (r + mu[n]);
      }

      double sqrt_m1 = sqrt(1.0 - emmc1);
      dn = sqrt_m1 / dn;
      // Is hypot(1, cs) better sqrt(1 + cs * cs)?
      cn = dn * (cos_umu > 0 ? 1 : -1) / sqrt(1 + cs * cs);
      sn = cn * cs / sqrt_m1;
    } else {
      double cot = cos_umu / sin_umu;
      double cs  = mu[n] * cot;
      dn = 1.0;

      while (n > 0)
      {
        --n;
        double r = cs * cs / mu[n + 1];
        cs *= dn;
        dn  = (r + nu[n]) / (r + mu[n]);
      }
      // Is hypot(1, cs) better than sqrt(1 + cs * cs)
      sn = (sin_umu > 0 ? 1 : -1) / sqrt(1 + cs * cs);
      cn = cs * sn;
    }
  }
}


//====================================================================
//============================================================END=====
