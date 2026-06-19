/*!
      @file    apreccond_MG-tmpl.h
      @brief   multi grid peconditionor with afield
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate::  $
      @version $LastChangedRevision: 2499 $
 */

#include <cassert>
#include "lib_alt/Solver/aprecond_MG.h"
#include "lib/ResourceManager/threadManager.h"


template<typename AFIELD, typename AFIELD2>
const std::string APrecond_MG<AFIELD, AFIELD2>::class_name = "APrecond_MG";

//====================================================================
namespace {
#ifndef AFILED_HAS_SUB
  template<typename AFIELD>
  inline void sub(AFIELD& v, const AFIELD& w)
  {
    typedef typename AFIELD::real_t real_t;
    axpy(v, (real_t)-1.0, w);
  }
#endif

#ifndef AFILED_HAS_ADD
  template<typename AFIELD>
  inline void add(AFIELD& v, const AFIELD& w)
  {
    typedef typename AFIELD::real_t real_t;
    axpy(v, (real_t)1.0, w);
  }
#endif

  static inline int get_flop_restrict_site(int ncol) // ncol=2*nvec; spin proj + mult
  {
    return ncol * (2 * 6 + 8 * 12);
  }


  static inline int get_flop_prolong_site(int ncol)  // ncol=2*nvec; spin proj + mult
  {
    return ncol * (2 * 6 + 8 * 12);
  }


  template<typename INDEX, typename AFIELD>
  void convert(const INDEX&, AFIELD& w, const INDEX&, const AFIELD& v)
  {
    copy(w, v);
  }
}

//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::init()
{
  ThreadManager::assert_single_thread(class_name);

  reset_flop_count();
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::set_coarse_solver(ASolver<AFIELD2> *solver)
{
  ThreadManager::assert_single_thread(class_name);
  m_coarse_solver = solver;
  int nin  = solver->get_fopr()->field_nin();
  int nvol = solver->get_fopr()->field_nvol();
  int nex  = solver->get_fopr()->field_nex();
  m_coarse_v.reset(nin, nvol, nex);
  m_coarse_w.reset(nin, nvol, nex);
  m_coarse_ncol = nin / 2;
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::set_smoother(ASolver<AFIELD2> *solver)
{
  ThreadManager::assert_single_thread(class_name);
  m_smoother = solver;
  int nin  = solver->get_fopr()->field_nin();
  int nvol = solver->get_fopr()->field_nvol();
  int nex  = solver->get_fopr()->field_nex();
  m_fine_w.reset(nin, nvol, nex);
  m_fine_v.reset(nin, nvol, nex);
  m_fine_r.reset(nin, nvol, nex);
  m_fine_tmp.reset(nin, nvol, nex);
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::tidyup()
{
  // ThreadManager::assert_single_thread(class_name);
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::mult(AFIELD& v, const AFIELD& w)
{
  int nin  = v.nin();
  int nvol = v.nvol();
  int nex  = v.nex();
  assert(w.check_size(nin, nvol, nex));
  assert(m_fine_v.check_size(nin, nvol, nex));
  assert(m_fine_w.check_size(nin, nvol, nex));

  Timer timer;
  Timer timer_each;
  timer.start();
  //  want to solve |v> =Dinv |w>

#pragma omp master
  {
    m_coarse_solver->set_init_mode(ASolver<AFIELD2>::InitialGuess::ZERO);
  }

#pragma omp barrier
  timer_each.reset();
  timer_each.start();
  double time0 = 0.0;
  // use the normalized source
  double norm2_in = norm2(w);
  double scale_in = sqrt(norm2_in);
  copy(v, w);
  v.scal(1.0 / scale_in);
  timer_each.stop();
  double time1       = timer_each.elapsed_sec();
  double time_double = time1 - time0;
  time0 = time1;

  // float <- double on fine lattice assuming lexical indexing
  timer_each.start();
  AIndex_lex<real_t0, IMPL> index_d;
  AIndex_lex<real_t, IMPL2> index_f;      // assumes float
  convert(index_f, m_fine_w, index_d, v); // if AFILED==AFIELD2, calls the local version
  timer_each.stop();
  time1 = timer_each.elapsed_sec();
  double time_d2f = time1 - time0;
  time0 = time1;

  // single precision part
  //  timer_each.start();
  mult_single(m_fine_v, m_fine_w);
  //  timer_each.stop();
  //  time1=timer_each.elapsed_sec();
  //  double time_mult_single=time1-time0; time0=time1;// measured in mult_sinle

  // double <- float
  timer_each.start();
  convert(index_d, v, index_f, m_fine_v);  // if AFILED==AFIELD2, calls the local version
  timer_each.stop();
  time1 = timer_each.elapsed_sec();
  double time_f2d = time1 - time0;
  time0 = time1;

  timer_each.start();
  v.scal(scale_in);  // recover the normalization
  timer_each.stop();
  time1        = timer_each.elapsed_sec();
  time_double += (time1 - time0);
  time0        = time1;

#pragma omp master
  { // flop counting
    const size_t Lvol = CommonParameters::Lvol();
    const double flop_other_double = 4 * 24 * static_cast<double>(Lvol);
    m_accum_flop_double += flop_other_double;
  }
#pragma omp barrier

  timer.stop();

#pragma omp master
  {
    double elapsed_time = timer.elapsed_sec();
    //  vout.general("  time for mult in %s = %e sec\n",
    //               class_name.c_str(), elapsed_time);
    ++m_num_mult_called;
    m_time_mult_total += elapsed_time;
    m_time_d2f        += time_d2f;
    m_time_f2d        += time_f2d;
    m_time_double     += time_double;
  }
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::mult_single(AFIELD2& v,
                                               const AFIELD2& w)
{
  //  want to solve |v> =Dinv |w>
  //  float [pre]
  //     guess: |v> = |w>
  //            |r> = |w> - D|v>  = |v> -D|v>
  //  one more time
  //
  //
  //  float: coarse
  //            |ww> = Dinv[coarse]|r>
  //            |v> := |v> + |ww>
  //  float: calc res.
  //            |r> = |w> - D|v>
  //  float: smoother
  //            |ww> = Dinv[smoother]|rr>
  //            |v> := |v> + |ww>

  //initial barrier --- no need as m_afoprF has a barrier
  //#pragma omp barrier

  Timer  timer;
  double time0 = 0.0;

  //  |v> = |w>
  // -|r> = |w> - D|v>
  timer.start();
  m_afoprF->mult(m_fine_r, w);
  copy(v, w);
  sub(m_fine_r, w);
  timer.stop();
  double time1         = timer.elapsed_sec();
  double time_residual = time1 - time0;
  time0 = time1;

  // restriction: corase vector <- fine vector
  timer.start();
  m_multigrid->make_coarse_vector(m_coarse_v, m_fine_r);
#pragma omp barrier
  timer.stop();
  time1 = timer.elapsed_sec();
  double time_restriction = time1 - time0;
  time0 = time1;
  //  vout.general("  time for restriction = %e sec\n", time_restriction);

  // coarse solver: coarse_w = (Dinv) coarse_v
  timer.start();
  int    coarse_nconv = -1;
  real_t coarse_diff  = -1.0;

  m_coarse_solver->solve(m_coarse_w, m_coarse_v, coarse_nconv, coarse_diff);
  vout.detailed("%s: coarse:    nconv = %d  diff = %e\n",
                class_name.c_str(), coarse_nconv, coarse_diff);
  timer.stop();
  time1 = timer.elapsed_sec();
  double time_coarse_solver = time1 - time0;
  time0 = time1;
  //  vout.general("  time for coarse solver = %e sec\n",
  //               time_coarse_solver);

  // prolong the solution to the fine grid
  // coarse_w => fine_tmp
  timer.start();
  m_multigrid->make_fine_vector(m_fine_tmp, m_coarse_w);
  timer.stop();
  time1 = timer.elapsed_sec();
  double time_prolongation = time1 - time0;
  time0 = time1;
  //  vout.general("  time for prologation = %e sec\n", time_prolongation);

  // update solution and calculate residual vector
  // -|v> := |w> - D |v>
  timer.start();
  sub(v, m_fine_tmp);
  m_afoprF->mult(m_fine_r, v);
  sub(m_fine_r, w);
#ifdef DEBUG
  {
    double r2 = m_fine_r.norm2();
    vout.general("%s:  after the coarse solver, r2=%15.7e\n", class_name.c_str(), r2);
  }
#endif
  timer.stop();
  time1          = timer.elapsed_sec();
  time_residual += (time1 - time0);
  time0          = time1;

  // call smoother
  // |tmp> = Dinv[fine] |v>    smoother
  // |w> = |w> + |tmp>         update solution
  timer.start();
#pragma omp barrier

  int    smoother_nconv = -1;
  real_t smoother_diff  = -1.0;
  m_smoother->solve(m_fine_tmp, m_fine_r, smoother_nconv, smoother_diff);

  vout.detailed("%s: smoother:  nconv = %d  diff = %e\n",
                class_name.c_str(), smoother_nconv, smoother_diff);
  timer.stop();
  time1 = timer.elapsed_sec();
  double time_smoother = time1 - time0;
  time0 = time1;
  //  vout.general("  time for smoother = %e sec\n", time_smoother);

  // update the solution
  timer.start();
  sub(v, m_fine_tmp);
#ifdef DEBUG
  {
    m_afoprF->mult(m_fine_r, v);
    sub(m_fine_r, w);
    double r2 = m_fine_r.norm2();
    vout.general("%s:  after the smoother, r2=%15.7e\n", class_name.c_str(), r2);
  }
#endif
  timer.stop();
  time1          = timer.elapsed_sec();
  time_residual += (time1 - time0);
  time0          = time1;
  double time_single_total = time1;
  //  double elapsed_time = time1;
  //  vout.general("  time for mult_single (total) = %e sec\n",
  //               elapsed_time);

  update_flop_count();

#pragma omp master
  {
    m_time_restriction       += time_restriction;
    m_time_coarse_solver     += time_coarse_solver;
    m_time_prolongation      += time_prolongation;
    m_time_smoother          += time_smoother;
    m_time_residual          += time_residual;
    m_time_mult_single_total += time_single_total;
    ++m_num_mult_single_called;
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::mult_as_setup(AFIELD2& v,
                                                 const AFIELD2& w)
{
  //  m_coarse_solver->set_init_mode(ASolver<AFIELD2>::InitialGuess::ZERO);
  mult_single(v, w);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::update_flop_count()
{
#pragma omp master
  {
    const size_t Lvol = CommonParameters::Lvol();

    const double flop_restrict_site = get_flop_restrict_site(m_coarse_ncol);
    const double flop_prolong_site  = get_flop_prolong_site(m_coarse_ncol);
    const double flop_other_float
      = m_afoprF->flop_count()
        + (flop_restrict_site + flop_prolong_site + 4 * 24)
        * static_cast<double>(Lvol)
        + m_afoprF->flop_count();
    double tmp = 0.0;
    m_accum_flop_coarse   += m_coarse_solver->flop_count();
    m_accum_flop_smoother += m_smoother->flop_count();
    m_accum_flop_other    += flop_other_float;
    tmp += m_coarse_solver->flop_count();
    tmp += m_smoother->flop_count();
    tmp += flop_other_float;
    m_accum_flop_float    += tmp;
    m_accum_flop_restrict += flop_restrict_site * static_cast<double>(Lvol);
    m_accum_flop_prolong  += flop_prolong_site * static_cast<double>(Lvol);
    vout.general("update_flop_count: flop_other_float=%e, m_accum_flop_float=%e, tmp=%e\n",
                 flop_other_float, m_accum_flop_float, tmp);
    //#pragma omp flush(m_accum_flop_float, m_accum_flop_other, m_accum_flop_smoother, m_accum_flop_coarse)
  }
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::reset_flop_count()
{
  m_accum_flop_coarse   = 0.0;
  m_accum_flop_smoother = 0.0;
  m_accum_flop_other    = 0.0;
  m_accum_flop_float    = 0.0;
  m_accum_flop_double   = 0.0;

  m_time_f2d               = 0.0;
  m_time_d2f               = 0.0;
  m_time_double            = 0.0;
  m_time_residual          = 0.0;
  m_time_restriction       = 0.0;
  m_time_coarse_solver     = 0.0;
  m_time_smoother          = 0.0;
  m_time_prolongation      = 0.0;
  m_time_mult_total        = 0.0;
  m_time_mult_single_total = 0.0;

  m_num_mult_called        = 0;
  m_num_mult_single_called = 0;
}


//====================================================================
template<typename AFIELD, typename AFIELD2>
void APrecond_MG<AFIELD, AFIELD2>::report_timer()
{
  vout.general("%s: time budget\n", class_name.c_str());
  vout.general("Elapsed time:  restriction        : total %14.6e\n", m_time_restriction);
  vout.general("Elapsed time:  coarse solver      : total %14.6e\n", m_time_coarse_solver);
  vout.general("Elapsed time:  prolongation       : total %14.6e\n", m_time_prolongation);
  vout.general("Elapsed time:  smoother           : total %14.6e\n", m_time_smoother);
  vout.general("Elapsed time:  resudial(+lin.alg) : total %14.6e\n", m_time_residual);
  vout.general("Elapsed time:  convert(f2d)       : total %14.6e\n", m_time_f2d);
  vout.general("Elapsed time:  convert(d2f)       : total %14.6e\n", m_time_d2f);
  vout.general("Elapsed time:  double             : total %14.6e\n", m_time_double);
  vout.general("Elapsed time:  mult_single(total) : total %14.6e\n", m_time_mult_single_total);
  vout.general("Elapsed time:  mult (total)       : total %14.6e\n", m_time_mult_total);
  vout.general("  number of mult call : %d\n", m_num_mult_called);

  double flop_restrict = get_flop_restrict_site(m_coarse_ncol)
                         * m_num_mult_single_called * static_cast<double>(CommonParameters::Lvol());
  double flop_prolong = get_flop_prolong_site(m_coarse_ncol)
                        * m_num_mult_single_called * static_cast<double>(CommonParameters::Lvol());

  vout.general("  Flops: smoother (float)       : %14.6e GFlop/s\n", flop_count_smoother() / m_time_smoother * 1.0e-9);
  vout.general("  Flops: coarse solver (float)  : %14.6e GFlop/s\n", flop_count_coarse() / m_time_coarse_solver * 1.0e-9);
  vout.general("  Flops: float total            : %14.6e GFlop/s\n", m_accum_flop_float / m_time_mult_single_total * 1.0e-9);
  vout.general("  Flops: restrict (float)       : %14.6e GFlop/s\n", flop_restrict / m_time_restriction * 1.0e-9);
  vout.general("  Flops: prolong (float)        : %14.6e GFlop/s\n", flop_prolong / m_time_prolongation * 1.0e-9);
}


//============================================================END=====
