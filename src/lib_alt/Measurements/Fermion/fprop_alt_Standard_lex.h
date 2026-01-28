/*!
      @file    fprop_alt_Standard_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-07-05 20:48:53 #$
      @version $LastChangedRevision: 2530 $
*/

#ifndef FPROP_ALT_STANDARD_LEX_INCLUDED
#define FPROP_ALT_STANDARD_LEX_INCLUDED

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
#include "lib/Smear/director_Smear.h"

#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
#include "lib_alt/Smear/director_alt_Smear.h"
#endif

#include "lib_alt/Solver/asolver.h"

//! Get quark propagator for Fopr with lexical site index: alternative version.

/*!
    This is temporary implementation.
                                        [30 May 2017 H.Matsufuru]
 */

template<typename AFIELD>
class Fprop_alt_Standard_lex : public Fprop_alt<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;
  using Fprop_alt<AFIELD>::m_vl;
  using Fprop_alt<AFIELD>::m_mode;

 private:
  AFopr<AFIELD> *m_kernel;
  AFopr<AFIELD> *m_fopr;
  ASolver<AFIELD> *m_solver;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;

  Director_Smear *m_dr_smear;
#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
  Director_alt_Smear<AFIELD> *m_dr_smear_alt;
#else
  Director_Smear *m_dr_smear_alt;
#endif
  bool m_alt_director;

 public:
  Fprop_alt_Standard_lex(const Parameters& params_fopr,
                         const Parameters& params_solver)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver); }

  Fprop_alt_Standard_lex(const Parameters& params_fopr,
                         const Parameters& params_solver,
                         Director_Smear* dr_smear)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver, dr_smear); }

#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
  Fprop_alt_Standard_lex(const Parameters& params_fopr,
                         const Parameters& params_solver,
                         Director_alt_Smear<AFIELD>* dr_smear)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver, dr_smear); }
#endif

  ~Fprop_alt_Standard_lex()
  { tidyup(); }

  void set_config(Field *);

  void invert(Field&, const Field&, int&, double&);

  void invert_D(Field&, const Field&, int&, double&);

  void invert_DdagD(Field&, const Field&, int&, double&);

  void invert_D_prec(Field&, const Field&, int&, double&);

  void invert_DdagD_prec(Field&, const Field&, int&, double&);

  // inverter with AFIELD
  void invert(AFIELD&, const AFIELD&, int&, double&);

  void invert_D(AFIELD&, const AFIELD&, int&, double&);

  void invert_DdagD(AFIELD&, const AFIELD&, int&, double&);

  void invert_D_prec(AFIELD&, const AFIELD&, int&, double&);

  void invert_DdagD_prec(AFIELD&, const AFIELD&, int&, double&);


  double flop_count();

  void reset_performance();

  void get_performance(double& flop_count, double& elapsed_time);

  void report_performance();

  void mult_performance(const std::string mode, const int Nrepeat);

 private:
  void init(const Parameters& params_fopr,
            const Parameters& params_solver);

  void init(const Parameters& params_fopr,
            const Parameters& params_solver,
            Director_Smear *dr_smear);

#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
  void init(const Parameters& params_fopr,
            const Parameters& params_solver,
            Director_alt_Smear<AFIELD> *dr_smear);
#endif

  void tidyup();
};
#endif
