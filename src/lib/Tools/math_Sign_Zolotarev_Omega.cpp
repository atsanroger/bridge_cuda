/*!
        @file    math_Sign_Zolotarev_Omega.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#include "Tools/math_Sign_Zolotarev_Omega.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
using namespace std;

//====================================================================
void Math_Sign_Zolotarev_Omega::get_sign_prms(vector<double>& omega,
                                                    double& delta){

  // Zolotarev coefficient defined

  if(omega.size() != d_Ns){
    vout.crucial(m_vl,
            "Math_Sign_Zolotarev_Omega: size of omega incorrect\n");
    abort();
  }

  for(int i=0; i<d_Ns; ++i){
    omega[i] = d_omega[i];
  };

  delta = d_delta;

}

//====================================================================
void Math_Sign_Zolotarev_Omega::set_sign_prms(){

  // Zolotarev coefficient defined

  d_omega.resize(d_Ns);
  d_cl.resize(d_Ns);

  double UK=0.0;
  vout.general(m_vl, " bmax = %12.4e\n", d_bmax);
  vout.general(m_vl, " UK   = %12.4e\n", UK);

  poly_Zolotarev(d_bmax, UK);
  d_D0 = 1.0;

  rescale_cl();
  vout.general(m_vl, "  rescaled: D0 = %16.8f\n",d_D0);

  rescale_cl();
  vout.general(m_vl, "  rescaled: D0 = %16.8f\n",d_D0);

  set_omega();

  // determine maximum deviation
  int Nj = 10000;
  double Dx = (d_bmax-1.0)/Nj;

  double d1 = 0.0;
  double d2 = 1.0e20;
  for(int jx=0; jx <= Nj; jx++){
    double x = 1.0 + Dx*jx;
    double sgnx = sign(x);
    if(fabs(sgnx) > d1) d1 = sgnx;
    if(fabs(sgnx) < d2) d2 = sgnx;
  }

  vout.general(m_vl, "  |sgnx-1|_up  = %14.8e \n", d1-1.0);
  vout.general(m_vl, "  |sgnx-1|_dn  = %14.8e \n", 1.0-d2);

  d_delta = d1-1.0;
  if(1.0-d2 > d_delta) d_delta = 1.0-d2;

}

//====================================================================
void Math_Sign_Zolotarev_Omega::set_omega(){

  int Ndiv = 10000;
  double dx = (d_bmax-1.0)/Ndiv;
  int nslv = 0;

  double x0, x1, v0, v1, omega2R;

  vout.general(m_vl, "  setting omega:\n");

  for(int div = 0; div < Ndiv+1; ++div){
    x0 = 1.0 + div*dx;
    x1 = x0 + dx;
    v0 = sign_cl(x0)-1.0;
    v1 = sign_cl(x1)-1.0;
    if(v1*v0 < 0.0){
      omega2R = binary_search(x0,x1);
      double v_zero = sign_cl(omega2R)-1.0;
      vout.general(m_vl, "    %2d %20.14f %20.14f\n",nslv,omega2R,v_zero);
      d_omega[nslv] = 1.0/omega2R;
      nslv++;
    }
  }

  if(nslv > d_Ns){
    vout.crucial(m_vl, "  too many zeros were found.\n");
    abort();
  }
}

//====================================================================
double Math_Sign_Zolotarev_Omega::sign(double x){

  double t = 1.0;
  for(int l=0; l<d_Ns; l++){
    t *= (1.0 - d_omega[l]*x)/(1.0 + d_omega[l]*x);
  };
  double s = (1-t)/(1+t);

  return s;

}

//====================================================================
double Math_Sign_Zolotarev_Omega::sign_cl(double x){

  int Np = d_Ns/2;

  double x2 = x*x;
  double sgn = d_D0;
  for(int l = 0; l < Np-1; ++l){
    sgn *= (x2 + d_cl[2*l+1])/(x2 + d_cl[2*l]);
  };
  sgn *= x/(x2 + d_cl[2*(Np-1)]);

  return sgn;

}

//====================================================================
double Math_Sign_Zolotarev_Omega::binary_search(double x1, double x2){

  double Eprec = 1.e-14;
  double Esep = x2-x1;

  double v_zero = x1;
  double v0 = sign_cl(v_zero)-1.0;
  double v1;

  while(Esep > Eprec){
    Esep *= 0.5;
    v1 = sign_cl(v_zero + Esep)-1.0;
    if(v1*v0 > 0.0){
      v_zero += Esep;
      v0 = v1;
    }
  }

  return v_zero;

}

//====================================================================
void Math_Sign_Zolotarev_Omega::poly_Zolotarev(double bmax, double& UK){

//    Return the coefficients of Zolotarev's approximation
//    of sign function to rational function.
//
//    Np: number of poles (2*Np is number of cl)
//    bmax: range of argument [1,bmax]
//    UK: complete elliptic integral of the 1st kind
//    cl[2*Np],bl[Np]: coefficients of Zolotarev rational approx.

  int Nprec = 14;

  double rk = sqrt( 1.0 - 1.0/(bmax*bmax) );
  double emmc = 1.0 - rk*rk;

  vout.general(m_vl, "emmc = %22.14e\n", emmc);

  double sn, cn, dn;

  // Determination of K

  double Dsr = 10.0;
  double u = 0.0;
  for(int iprec=0; iprec<Nprec+1; iprec++){
    Dsr = Dsr * 0.1;

    for(int i=0; i<20; i++){
      u = u + Dsr;
      Jacobi_elliptic(u,emmc,sn,cn,dn);
      vout.general(m_vl, " %22.14e %22.14e %22.14e\n", u, sn, cn);
      if(cn < 0.0) goto succeeded;
    }
    vout.general(m_vl, "Something wrong in setting Zolotarev\n");
    return;

    succeeded: 
      u = u - Dsr;

  };

  UK = u;

  Jacobi_elliptic(UK,emmc,sn,cn,dn);
  vout.general(m_vl, " %22.14e %22.14e %22.14e\n", UK, sn, cn);


// Determination of c_l

  double FK = UK/d_Ns;

  for(int l = 0; l < d_Ns-1; ++l){
    u = FK * (l + 1.0);
    Jacobi_elliptic(u,emmc,sn,cn,dn);
    d_cl[l] = sn*sn/(1.0-sn*sn);
  };

}

//====================================================================
void Math_Sign_Zolotarev_Omega::rescale_cl(){

  int Nj = 10000;
  double Dx = (d_bmax-1.0)/Nj;

  double d1 = 0.0;
  double d2 = 1.0e20;
  for(int jx=0; jx <= Nj; jx++){
    double x = 1.0 + Dx*jx;
    double sgnx = sign_cl(x);
    if(fabs(sgnx) > d1) d1 = sgnx;
    if(fabs(sgnx) < d2) d2 = sgnx;
  }

    vout.general(m_vl, "  |sgnx-1|_up  = %14.8e \n", d1-1.0);
    vout.general(m_vl, "  |sgnx-1|_dn  = %14.8e \n", 1.0-d2);

  d_D0 = d_D0/(0.5*(d1+d2));

}

//====================================================================
void Math_Sign_Zolotarev_Omega::
Jacobi_elliptic(double uu, double emmc,
                     double &sn, double &cn, double &dn){

//  Return the Jacobi elliptic functions sn(u,kc),
//  cn(u,kc), and dn(u,kc), where kc = 1 - k^2, and
//  arguments:
//    uu: u, emmc: kc.
//
//  This routine is from
//  W.H.Press et al., Numerical Recipes in Fortran 2nd ed.
//  (Cambridge Univ. Press, 1986,1992).
//                         Typed by H. Matsufuru 1 May 2005

  double const CA = 0.00000001;
  //          (The accuracy is the square of CA.)

  double a,b,c,d,em[13],en[13];

  int bo;
  //  LOGICAL bo

// main

  double emc = emmc;
  double u = uu;

  if(emc != 0.0){

    //    bo = (emc.LT.0.D0)   // how to do ?
    int bo = 0;
    if(emc < 0.0) bo = 1;

    if(bo != 0){
      d = 1.0 - emc;
      emc = - emc/d;
      d = sqrt(d);
      u = d*u;
    }

    a  = 1.0;
    dn = 1.0;

    int l;
    for(int i=0; i<13; i++){
      l = i;
      em[i] = a;
      emc = sqrt(emc);
      en[i] = emc;
      c = 0.5*(a + emc);
      if(fabs(a-emc) <= CA*a) goto continued;
      emc = a * emc;
      a = c;
    }

    continued:

    u = c*u;
    sn = sin(u);
    cn = cos(u);

    if(sn == 0.0) goto continued2;

    a = cn/sn;
    c = a * c;

    for(int ii=l; ii >= 0; --ii){
      b = em[ii];
      a = c * a;
      c = dn * c;
      dn = (en[ii] + a)/(b + a);
      a = c/b;
    }

    a = 1.0/sqrt(c*c + 1.0);

    if(sn < 0.0){
      sn = -a;
    }else{
      sn = a;
    }

    cn = c * sn;


  continued2:

   if(bo != 0){
     a = dn;
     dn = cn;
     cn = a;
     sn = sn/d;
   }

  }else{

    cn = 1.0/cosh(u);
    dn = cn;
    sn = tanh(u);

  }

}
//====================================================================
//============================================================END=====

