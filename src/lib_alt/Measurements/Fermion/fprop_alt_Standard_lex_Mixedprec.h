/*!
      @file    fprop_alt_Standard_lex_Mixedprec.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef FPROP_ALT_STANDARD_LEX_MIXEDPREC_INCLUDED
#define FPROP_ALT_STANDARD_LEX_MIXEDPREC_INCLUDED

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib/Fopr/afopr.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Tools/timer.h"
#include "lib/Smear/director_Smear.h"

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Solver/aprecond.h"

//! Get quark propagator for Fopr with lexical site index: alternative version.

/*!
    This is temporary implementation.
                                        [31 May 2017 H.Matsufuru]
 */

template<typename AFIELD_d, typename AFIELD_f>
class Fprop_alt_Standard_lex_Mixedprec : public Fprop_alt<AFIELD_d>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD_d::real_t real_t;
  using Fprop_alt<AFIELD_d>::m_vl;
  using Fprop_alt<AFIELD_d>::m_mode;

 private:
  AFopr<AFIELD_d> *m_kernel_d;
  AFopr<AFIELD_f> *m_kernel_f;
  AFopr<AFIELD_d> *m_fopr_d;
  AFopr<AFIELD_f> *m_fopr_f;
  ASolver<AFIELD_f> *m_solver_prec;
  APrecond<AFIELD_d> *m_precond;
  ASolver<AFIELD_d> *m_solver;

  Director_Smear *m_dr_smear;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;

 public:
  Fprop_alt_Standard_lex_Mixedprec(const Parameters& params_fopr,
                                   const Parameters& params_solver)
    : Fprop_alt<AFIELD_d>()
  { init(params_fopr, params_solver); }

  Fprop_alt_Standard_lex_Mixedprec(const Parameters& params_fopr,
                                   const Parameters& params_solver,
                                   Director_Smear *dr_smear)
    : Fprop_alt<AFIELD_d>()
  { init(params_fopr, params_solver, dr_smear); }

  ~Fprop_alt_Standard_lex_Mixedprec()
  { tidyup(); }

  void set_config(Field *);

  void invert(Field&, const Field&, int&, double&);

  void invert_D(Field&, const Field&, int&, double&);

  void invert_DdagD(Field&, const Field&, int&, double&);

  void invert(AFIELD_d&, const AFIELD_d&, int&, double&);

  void invert_D(AFIELD_d&, const AFIELD_d&, int&, double&);

  void invert_DdagD(AFIELD_d&, const AFIELD_d&, int&, double&);

  double flop_count();

  void reset_performance();

  void report_performance();

  void mult_performance(const std::string mode, const int Nrepeat);

 private:

  void init(const Parameters& params_fopr,
            const Parameters& params_solver);

  void init(const Parameters& params_fopr,
            const Parameters& params_solver,
            Director_Smear *dr_smrar);

  void tidyup();
};
#endif
