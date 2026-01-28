/*!
        @file    fprop_Overlap_5d.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FPROP_OVERLAP_5D_INCLUDED
#define FPROP_OVERLAP_5D_INCLUDED

#include "Fopr/fopr_Overlap.h"
#include "Fopr/fopr_Overlap_5d.h"

#include "Field/index_lex.h"
#include "Field/index_eo.h"

#include "Solver/solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Calculation of overlap propagater with 5d solver.

/*!
   This class calculates overlap quark propagator using
   5-dimensional solver.
                                 [25 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                 [21 Mar 2015 Y.Namekawa]
    Add flop_count.              [ 8 Aug 2016 Y.Namekawa]
 */


class Fprop_Overlap_5d
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  //- parameters common to overlap fermion
  double m_mq;                         // quark mass
  double m_M0;                         // domain-wall height
  int m_Np;                            // number of poles in rational approx.
  double m_x_min, m_x_max;             // valid range of approximate sign function
  int m_Niter_ms;                      // max iteration of shiftsolver (dummy)
  double m_Stop_cond_ms;               // stopping condition of shift solver (dummy)
  std::vector<int> m_boundary;

  std::string m_kernel_type;
  std::string m_repr;

  //- solver parameter
  Parameters m_params_solver;

  Fopr_Overlap *m_fopr_ov;
  Field_G *m_U;

  Fopr_Overlap_5d *m_fopr_ov5d;
  Fopr_Wilson_eo *m_fopr_w_eo;
  Field_G *m_Ueo;
  Solver *m_solver;

  Index_lex m_index;
  Index_eo m_index_eo;

  //- low-mode subtraction
  int m_Nsbt;
  std::vector<double> *m_ev;
  std::vector<Field> *m_vk;

  bool m_is_initialized;

 public:
  Fprop_Overlap_5d(Fopr_Overlap *fopr_ov, Field_G *U)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr_ov(fopr_ov),
    m_U(U),
    m_is_initialized(false)
  {}

  ~Fprop_Overlap_5d()
  { finalize(); }

 private:
  // non-copyable
  Fprop_Overlap_5d(const Fprop_Overlap_5d&);
  Fprop_Overlap_5d& operator=(const Fprop_Overlap_5d&);

 public:
  void set_parameters(const Parameters& params_fprop_overlap_5d);
  void set_parameters(const Parameters& params_overlap,
                      const Parameters& params_solver);
  void set_parameters(const double mq, const double M0, const int Np,
                      const double x_min, const double x_max,
                      const int Niter_ms, const double Stop_cond_ms,
                      const std::vector<int> bc,
                      const std::string& repr
                      );

  void get_parameters(Parameters& params) const;

  void invert_D(Field&, const Field&);

  void calc_H2inv(Field&, const Field&);

  void solve_overlap_5D_1(Field&, const Field&, int&, double&);

  void set_lowmodes(const int Nsbt, std::vector<double> *ev,
                    std::vector<Field> *vk)
  {
    m_Nsbt = Nsbt;
    m_ev   = ev;
    m_vk   = vk;
  }

  double flop_count();

 private:
  void initialize();
  void finalize();
};
#endif
