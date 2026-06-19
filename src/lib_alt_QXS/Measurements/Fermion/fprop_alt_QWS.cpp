/*!
      @file    fprop_alt_QWS.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/Measurements/Fermion/fprop_alt_QWS.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Fopr/afopr_Clover_QWS_dd.h"

#ifdef USE_QWSLIB
#include "lib_alt_QXS/extra/qws_lib.h"
#endif


template<typename AFIELD>
const std::string Fprop_alt_QWS<AFIELD>::class_name
  = "Fprop_alt_QWS";
//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::init(const Parameters& params_fopr,
                                 const Parameters& params_solver)
{   // this constructor assumes that the factories are available.
  // vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup (without link smearing).\n",
               class_name.c_str());

  //  typedef AFopr<AField<real_t> >  AltFopr;
  //  typedef ASolver<AField<real_t> >  AltSolver;
  typedef AFopr_dd<AFIELD> AltFopr;

  string fopr_type = params_fopr.get_string("fermion_type");
  //  m_fopr = AltFopr::New("AFopr_Clover_QWS_dd", params_fopr);
  m_fopr     = new AFopr_Clover_QWS_dd<AFIELD>(params_fopr);
  m_dr_smear = 0;

  m_kernel = 0;

  //  string solver_type = params_solver.get_string("solver_type");

  set_parameters(params_solver);
  reset_performance();

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::set_parameters(const Parameters& params_solver)
{
  //- fetch and check input parameters
  int    Niter, Nrestart, Niter_s;
  double Stop_cond;


  int err = 0;
  err += params_solver.fetch_int("maximum_number_of_iteration_single_prec", Niter_s);
  err += params_solver.fetch_int("maximum_number_of_iteration", Niter);
  err += params_solver.fetch_int("maximum_number_of_restart", Nrestart);
  err += params_solver.fetch_double("convergence_criterion_squared", Stop_cond);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Niter_d = Niter;
  m_Niter_s = Niter_s;
  m_nm      = 2;
  m_Nsap    = 4;
  if (Stop_cond > 0) {
    m_tol_d = sqrt(Stop_cond);
    m_tol_s = 1.0e-6;
  } else {
    m_tol_d = -1.0;
    m_tol_s = -1.0;
    vout.general(m_vl, "%s: negative convergence criterion is given, run as a fixed iteration sovler\n", class_name.c_str());
  }


  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter_d     = %d\n", m_Niter_d);
  vout.general(m_vl, "  Niter_s     = %d\n", m_Niter_s);
  vout.general(m_vl, "  nm          = %d\n", m_nm);
  vout.general(m_vl, "  nsap        = %d\n", m_Nsap);
  vout.general(m_vl, "  tol_d       = %e\n", m_tol_d);
  vout.general(m_vl, "  tol_s       = %e\n", m_tol_s);
}


/*
//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::init(const Parameters& params_fopr,
                                          const Parameters& params_solver,
                                          Director_Smear* dr_smear)
{  // this constructor assumes that the factories are available.

  // vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup (with link smearing).\n",
               class_name.c_str());

  typedef AFopr<AFIELD>  AltFopr;
  typedef ASolver<AFIELD>  AltSolver;

  m_dr_smear = dr_smear;

  string fopr_type = params_fopr.get_string("fermion_type");
  m_kernel = AltFopr::New(fopr_type, params_fopr);

  //  m_fopr = AltFopr::New("Smeared", m_kernel, m_dr_smear);
  m_fopr = new AFopr_Smeared<AFIELD>(m_kernel, m_dr_smear);

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());

}
*/

//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::tidyup()
{
  delete m_fopr;
  if (m_kernel != 0) delete m_kernel;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::invert(Field& xq, const Field& b,
                                   int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
    //  }else if(m_mode == "DdagD"){
    //    invert_DdagD(xq, b, nconv, diff);
    //  }else if(m_mode == "D_prec"){
    //    invert_D_prec(xq, b, nconv, diff);
    //  }else if(m_mode == "DdagD_prec"){
    //    invert_DdagD_prec(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::invert_D(Field& xq, const Field& b,
                                     int& nconv, double& diff)
{
  m_timer.reset();

  if (sizeof(real_t) != 8) {
    vout.crucial(m_vl, "%s: single prec is not supported\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  vout.paranoiac(m_vl, "nin = %d  nvol = %d nex = %d\n", nin, nvol, nex);

  AFIELD abq(nin, nvol, nex);
  AFIELD axq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

#pragma omp parallel
  {
    m_fopr->convert(abq, b);
    axq.set(0.0); // set initial guess
  }
  // copy(axq, abq);


#ifdef USE_QWSLIB
  int    nsap         = m_Nsap;
  int    nm           = m_nm;
  int    dd_maxiter_s = m_Niter_s;
  int    dd_maxiter   = m_Niter_d;
  double tol_s        = m_tol_s;
  double tol          = m_tol_d;
  int    iter         = -1;

  m_timer.start();
#pragma omp parallel
  {
    ((AFopr_Clover_QWS_dd<AFIELD> *)m_fopr)->mult_clv_inv(abq);
  }
  //  bicgstab_dd_mix_(x, b, &tol, &iter, &dd_maxiter, &tol_s, &dd_maxiter_s, &nsap, &nm);

  bicgstab_dd_mix_((scd_t *)axq.ptr(0), (scd_t *)abq.ptr(0), &tol, &iter, &dd_maxiter, &tol_s, &dd_maxiter_s, &nsap, &nm);

  m_timer.stop();
  nconv = iter;
  diff  = -1.0;

#pragma omp parallel
  {
    m_fopr->reverse(xq, axq);
  }
#endif


  m_elapsed_time += m_timer.elapsed_sec();
}


template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::invert_DdagD(Field& xq,
                                         const Field& b,
                                         int& nconv, double& diff)
{
  // dummy
}


/*
//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::invert(AFIELD& xq,
                                            const AFIELD& b,
                                            int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if(m_mode == "D"){
    invert_D(xq, b, nconv, diff);
    //  }else if(m_mode == "DdagD"){
    //    invert_DdagD(xq, b, nconv, diff);
    //  }else if(m_mode == "D_prec"){
    //    invert_D_prec(xq, b, nconv, diff);
    //  }else if(m_mode == "DdagD_prec"){
    ///    invert_DdagD_prec(xq, b, nconv, diff);
  }else{
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}
*/

//====================================================================
template<typename AFIELD>
double Fprop_alt_QWS<AFIELD>::flop_count()
{
  return 0.0;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::get_performance(double& flop_count,
                                            double& elapsed_time)
{
  flop_count   = 0.0;
  elapsed_time = m_elapsed_time;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::report_performance()
{
  //  double flops  = m_flop_count/m_elapsed_time;
  double flops  = 0;
  double gflops = flops * 1.0e-9;

  //  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: solver performance:\n", class_name.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", m_elapsed_time);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n", m_flop_count);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_QWS<AFIELD>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
  /*
  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  AFIELD axq(nin, nvol, nex), abq(nin, nvol, nex);
  abq.set(0.0);
  abq.set(0, 1.0);

  unique_ptr<Timer> timer(new Timer);

  std::string mode_prev = m_fopr->get_mode();

  m_fopr->set_mode(mode);

  timer->start();

#pragma omp parallel
  {
    for(int i = 0; i < Nrepeat; ++i){
      m_fopr->mult(axq, abq);
      m_fopr->mult(abq, axq);
    }
  }

  timer->stop();

  double flop_fopr = m_fopr->flop_count();
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops  = flop_total/elapsed_time;
  double gflops = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  mult mode = %s\n", mode.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(Fopr)   = %18.0f\n",flop_fopr);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n",flop_total);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);

  m_fopr->set_mode(mode_prev);
  */
}


//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Fprop_alt_QWS<AField<double, QXS> >::class_name
  = "Fprop_alt_QWS<Afield<double,QXS> >";


template class Fprop_alt_QWS<AField<double, QXS> >;

//============================================================END=====
