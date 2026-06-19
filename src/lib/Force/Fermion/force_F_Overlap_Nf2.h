/*!
        @file    force_F_Overlap_Nf2.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/


#ifndef FORCE_F_OVERLAP_NF2_INCLUDED
#define FORCE_F_OVERLAP_NF2_INCLUDED

#include "force_F_Wilson_Nf2.h"

#include "Fopr/fopr_Overlap.h"
#include "Field/index_lex.h"
#include "Solver/solver_CG.h"
#include "Solver/shiftsolver_CG.h"
#include "Tools/math_Sign_Zolotarev.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for overlap ferimon.

/*!
    At present, recursive calculation of smeqared fermion
    force is not implemented.
                              [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.      [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Overlap_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_mq, m_M0;
  int m_Np;
  double m_x_min, m_x_max;

  int m_Niter_ms;
  double m_Stop_cond_ms;

  std::vector<int> m_boundary;

  double m_kappa;

  Fopr_Wilson *m_fopr_w;
  Force_F_Wilson_Nf2 *m_force_w;

  // Zolotarev coefficients
  std::vector<double> m_cl;
  std::vector<double> m_bl;
  std::vector<double> m_sigma;

 public:
  Force_F_Overlap_Nf2()
    : m_vl(CommonParameters::Vlevel()) {}

  Force_F_Overlap_Nf2(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  ~Force_F_Overlap_Nf2()
  {
    delete m_fopr_w;
    delete m_force_w;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const double M0, const int Np,
                      const double x_min, const double x_max,
                      const int Niter_ms, const double Stop_cond_ms,
                      const std::vector<int> bc);
  void set_parameters();

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
    m_force_w->set_config(U);
  }

  //  void set_eigenmodes();

  void force_core(Field&, const Field&);
  void force_core1(Field&, const Field&, const Field&);

  void force_udiv(Field&, const Field&);
  void force_udiv1(Field&, const Field&, const Field&);

  void force_core1_impl(Field_G&, const Field_F&, const Field_F&);
};
#endif
