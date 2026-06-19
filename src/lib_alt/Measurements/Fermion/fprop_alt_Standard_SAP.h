/*!
      @file    fprop_alt_Standard_SAP.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef FPROP_ALT_STANDARD_LEX_SAP_INCLUDED
#define FPROP_ALT_STANDARD_LEX_SAP_INCLUDED

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt/Fopr/afopr_dd.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Tools/timer.h"
#include "lib/Smear/director_Smear.h"

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Solver/asolver_SAP.h"

//! Get quark propagator for domain-decomposed AFopr.

/*!
  This class determines the quark propagator for the
  domain-decomposed version of AFopr (alt-version).
                                      [10 Apr 2021 H.Matsufuru]
 */

template<typename AFIELD>
class Fprop_alt_Standard_SAP : public Fprop_alt<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;
  using Fprop_alt<AFIELD>::m_vl;
  using Fprop_alt<AFIELD>::m_mode;

 private:
  AFopr_dd<AFIELD> *m_kernel;
  AFopr_dd<AFIELD> *m_afopr;
  AFopr_dd<AFIELD> *m_afopr_dd;
  ASolver_SAP<AFIELD> *m_asolver;

  Director_Smear *m_dr_smear;

  std::vector<int> m_coarse_lattice;
  std::vector<int> m_fine_lattice;
  AIndex_block_lex<real_t, AFIELD::IMPL> *m_block_index;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;

 public:
  Fprop_alt_Standard_SAP(const Parameters& params_fopr,
                         const Parameters& params_solver)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver); }

  Fprop_alt_Standard_SAP(const Parameters& params_fopr,
                         const Parameters& params_solver,
                         Director_Smear *dr_smear)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver, dr_smear); }

  ~Fprop_alt_Standard_SAP()
  { tidyup(); }

  void set_config(Field *);

  void invert(Field&, const Field&, int&, double&);

  void invert_D(Field&, const Field&, int&, double&);

  void invert_DdagD(Field&, const Field&, int&, double&);

  void invert(AFIELD&, const AFIELD&, int&, double&);

  void invert_D(AFIELD&, const AFIELD&, int&, double&);

  void invert_DdagD(AFIELD&, const AFIELD&, int&, double&);

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

  void init_common(const Parameters& params_fopr,
                   const Parameters& params_solver);

  void tidyup();
};
#endif
