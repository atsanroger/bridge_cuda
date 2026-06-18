/*!
      @file    aprecond_MG_dw.h
      @brief   Multigrid preconditionor for Domainwall operator
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2024-04-01 12:13:09 +0900 (2024年04月01日 (月)) $
      @version $LastChangedRevision: 2595 $
*/

#ifndef APRECOND_MG_DW_H
#define APRECOND_MG_DW_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include  "lib/Parameters/commonParameters.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;
#include "lib/Tools/timer.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt/Solver/asolver.h"
#include "lib_alt/Solver/MultiGrid.h"
#include "lib_alt/Solver/aprecond.h"

/*!
    explanation to be added.
 */

template<typename AFIELD, typename AFIELD2>
class APrecond_MG_dw : public APrecond<AFIELD>
{
 public:
  typedef typename AFIELD::real_t         real_t0; // assumes double
  typedef typename AFIELD2::real_t        real_t;  // assumes float
  typedef MultiGrid<AFIELD2, AFIELD2>     MultiGrid_t;
  typedef typename MultiGrid_t::Index_t   block_index_t;

  static const Impl IMPL  = AFIELD::IMPL;
  static const Impl IMPL2 = AFIELD2::IMPL;

  static const std::string class_name;

 private:
  ASolver<AFIELD2> *m_coarse_solver;
  ASolver<AFIELD2> *m_smoother;
  ASolver<AFIELD2> *m_PV_solver;
  const MultiGrid_t *m_multigrid;

  AFopr<AFIELD> *m_afoprD;  // double
  AFopr<AFIELD2> *m_afoprF; // float

  AFIELD2 m_coarse_v, m_coarse_w;
  AFIELD2 m_fine_w, m_fine_v, m_fine_tmp, m_fine_r;
  AFIELD2 m_fine_tmp2;
  AFIELD2 m_v;
  AFIELD2 m_fine_v0, m_fine_v1, m_fine_z0, m_fine_z1;
  AFIELD m_fine_vecD, m_fine_debugD, m_fine_rD;
  //  const std::vector<AFIELD2> *m_testvectors;

  
  double m_accum_flop_coarse;
  double m_accum_flop_smoother;
  double m_accum_flop_other;
  double m_accum_flop_float;
  double m_accum_flop_double;
  double m_accum_flop_restrict;
  double m_accum_flop_prolong;

  int m_coarse_ncol;

  // elapsed time for each ingredient
  double m_time_f2d;
  double m_time_d2f;
  double m_time_residual;
  double m_time_restriction;
  double m_time_coarse_solver;
  double m_time_smoother;
  double m_time_prolongation;
  double m_time_mult_single_total;
  double m_time_mult_total;
  double m_time_double;
  double m_time_PV_solver;
  int m_num_mult_called;        //!< number of mult call
  int m_num_mult_single_called; //!< number of mult call

  //  unique_ptr<Timer> m_timer_f2d;  // float to double
  //  unique_ptr<Timer> m_timer_d2f;  // double to flaot
  //  unique_ptr<Timer> m_timer_residual;
  //  unique_ptr<Timer> m_timer_restrcition;
  //  unique_ptr<Timer> m_timer_coarse_solver;
  //  unique_ptr<Timer> m_timer_smoother;
  //  unique_ptr<Timer> m_timer_prolongation;
  //  unique_ptr<Timer> m_timer_mult_single_total;

 public:
  APrecond_MG_dw() { init(); }

  void mult(AFIELD&, const AFIELD&);
  void mult_as_setup(AFIELD2&, const AFIELD2&);
  //! dev8: public accessor to the single-RHS float V-cycle, so the AMG-local
  //! block-FGMRES driver can apply the (float) preconditioner column-by-column.
  void apply_vcycle_single(AFIELD2& out, const AFIELD2& in) { mult_single(out, in); }
  void reset_flop_count();
  void set_coarse_solver(ASolver<AFIELD2> *solver);
  void set_smoother(ASolver<AFIELD2> *solver);
  void set_PV_solver(ASolver<AFIELD2> *solver);

  void set_multigrid(const MultiGrid_t *multigrid)
  {
    m_multigrid = multigrid;
  }

  void set_fopr(AFopr<AFIELD> *foprD, AFopr<AFIELD2> *foprF) { m_afoprD = foprD; m_afoprF = foprF; }
  void report_timer();

  double flop_count() { return m_accum_flop_float + m_accum_flop_double; }
  double flop_count_coarse() const { return m_accum_flop_coarse; }
  double flop_count_smoother() const { return m_accum_flop_smoother; }
  double flop_count_other() const { return m_accum_flop_other; }
  double flop_count_double() const { return m_accum_flop_double; }

 private:
  void init();
  void tidyup();
  void update_flop_count();
  void mult_single(AFIELD2&, const AFIELD2&);
  void mult_single2(AFIELD2&, const AFIELD2&);
  void mult_single_gmres2(AFIELD2&, const AFIELD2&);
  void mult_single_smoother_only(AFIELD2&, const AFIELD2&);
};

#endif // include guard
//============================================================END=====
