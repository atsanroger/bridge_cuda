/*!
        @file    fermionFlow.h

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FERMIONFLOW_INCLUDED
#define FERMIONFLOW_INCLUDED

#include "Action/action.h"
#include "Field/shiftField_lex.h"
#include "Field/field.h"
#include "Field/field_F.h"
#include "Field/field_F_1spinor.h"

#include "Measurements/Gauge/gradientFlow.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#include <iostream>
#include <fstream>
#include <vector>
using std::vector;

/*!
  This class implements fermion flow introduced in https://arxiv.org/abs/1302.5246 .
  This is written by Y. Taniguchi and K. Kanaya for ver.1.2.2.
  [April 2015 Y.Taniguchi]
  Ported to trunk by Y. Taniguchi.
  [29 May 2017 Y.Taniguchi]
*/

class FermionFlow
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_order_RK;
  double m_Estep;
  //! The Number of flow steps shoud be set by calling set_parameters(Nstep).
  int m_Nstep;
  int m_Nprec;

  int m_Ndim;
  int m_Nvol;

  Action *m_action;
  Field_G m_iP, m_iPP;
  Field_G m_U_0, m_W_1;

  std::vector<int> m_boundary;
  ShiftField_lex m_shift;
  GradientFlow m_gflow;

 public:
  FermionFlow(Action* action)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()),
    m_Nvol(CommonParameters::Nvol()),
    m_action(action),
    m_gflow(action)
  {
    m_iP.reset(m_Nvol, m_Ndim);
    m_iPP.reset(m_Nvol, m_Ndim);

    m_U_0.reset(m_Nvol, m_Ndim);
    m_W_1.reset(m_Nvol, m_Ndim);

    m_Nprec = 0;
    m_Nstep = 0;
    m_Estep = 0.0;
  }

  FermionFlow(Action* action, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()),
    m_Nvol(CommonParameters::Nvol()),
    m_action(action),
    m_gflow(action)
  {
    m_iP.reset(m_Nvol, m_Ndim);
    m_iPP.reset(m_Nvol, m_Ndim);

    m_U_0.reset(m_Nvol, m_Ndim);
    m_W_1.reset(m_Nvol, m_Ndim);

    m_Nprec = 0;
    m_Nstep = 0;
    m_Estep = 0.0;

    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const int order_RK,
                      const double Estep, const int Nprec, const std::vector<int> bc);
  void set_parameters(const int Nstep);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void get_parameters(Parameters& params) const;


  //! v=Laplacian(U)f=D_\mu(U)D_\mu(U)f
  void laplacian(Field_F& v, const Field_F& f, Field_G& U);
  void laplacian(Field_F_1spinor& v, const Field_F_1spinor& f, Field_G& U);

  //! v=U_\mu(x)f(x+\mu)-U^\dagger_\mu(x-\mu)f(x-\mu)
  void del_symmetric(int mu, Field_F& v, const Field_F& f, Field_G& U);
  void del_symmetric(int mu, Field_F_1spinor& v, const Field_F_1spinor& f, Field_G& U);

  //! Not implemented
  void fermionFlow_1st(Field_F& f, Field_G& U);

  //! Not implemented
  void fermionFlow_2nd(Field_F& f, Field_G& U);

  //! Adjoint fermion flow using 3rd oder of Runge-Kutta.
  void fermionFlow_3rd(Field_F& f, Field_G& U);
  void fermionFlow_3rd(Field_F_1spinor& f, Field_G& U);

  //! Not implemented
  void fermionFlow_4th(Field_F& f, Field_G& U);

  //! Normal order fermion flow using 3rd oder of Runge-Kutta. 1-spinor fermion field f(x,t) and 4-spinor fermion field sf(x,t)=S(x,y)f(y,t) is flowed to t+epsilon.
  void fermionFlow_normal_order_3rd(Field_F_1spinor& f, Field_F& sf, Field_G& U);

  //! Flow a fermion field f(t) to t=0. Gauge flow is evaluated at every step.
  void evolve(Field_F& f, Field_G& U);
  void evolve(Field_F_1spinor& f, Field_G& U);

  //! Flow a fermon field D_mu(U(t))f(t) to t=0. Gauge flow is evaluated at every step.
  void evolve(int mu, Field_F& f, Field_G& U);

  /*!
    Flow fermion fields f(t) and D_\mu(U(t))f(t) simultaneously.
    In the field "f" both f(t) and D_\mu(U(t))f(t) are stored.
    In U0 guage links are stored at each flow time separation.
  */
  void evolve2(Field_F_1spinor& f, const Field_G& U0, int measure, int gauge_store_interval);

  //! Flow a 1-spinor fermion field f(x,0) and 4-spinor fermion field sf(x,t)=S(x,y)f(y,0) to t in normal order. Gauge flow is performed simultaneously.
  void evolve_normal_order(Field_F_1spinor& f, Field_F& sf, Field_G& U);
};
#endif
