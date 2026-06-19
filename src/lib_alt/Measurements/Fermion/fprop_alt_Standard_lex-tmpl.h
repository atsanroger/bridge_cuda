/*!
      @file    fprop_alt_Standard_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-11-13 13:02:59 #$
      @version $LastChangedRevision: 2553 $
*/

template<typename AFIELD>
const std::string Fprop_alt_Standard_lex<AFIELD>::class_name
                                  = "Fprop_alt_Standard_lex<AFIELD>";
//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::init(const Parameters& params_fopr,
                                          const Parameters& params_solver)
{   // this constructor assumes that the factories are available.
  ThreadManager::assert_single_thread(class_name);

  // vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: construction (without link smearing).\n",
               class_name.c_str());
  vout.increase_indent();
  
  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  string fopr_type = params_fopr.get_string("fermion_type");
  m_fopr = AltFopr::New(fopr_type, params_fopr);

  m_dr_smear = 0;
  m_dr_smear_alt = 0;
  m_kernel = 0;
  m_alt_director = false;

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::init(const Parameters& params_fopr,
                                          const Parameters& params_solver,
                                          Director_Smear *dr_smear)
{  // this constructor assumes that the factories are available.
  // vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: construction (with link smearing: corelib).\n",
               class_name.c_str());
  vout.increase_indent();

  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  m_dr_smear = dr_smear;
  m_dr_smear_alt = 0;
  m_alt_director = false;

  string fopr_type = params_fopr.get_string("fermion_type");
  m_kernel = AltFopr::New(fopr_type, params_fopr);

  //  m_fopr = AltFopr::New("Smeared", m_kernel, m_dr_smear);
  m_fopr = new AFopr_Smeared<AFIELD>(m_kernel, m_dr_smear);

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());
}


//====================================================================
#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::init(
                               const Parameters& params_fopr,
                               const Parameters& params_solver,
                               Director_alt_Smear<AFIELD> *dr_smear)
{  // this constructor assumes that the factories are available.

  vout.general(m_vl, "%s: construction (with link smearing: alt).\n",
               class_name.c_str());
  vout.increase_indent();

  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  m_dr_smear = 0;
  m_dr_smear_alt = dr_smear;
  m_alt_director = true;

  string fopr_type = params_fopr.get_string("fermion_type");
  m_kernel = AltFopr::New(fopr_type, params_fopr);

  //  m_fopr = AltFopr::New("Smeared_alt", m_kernel, m_dr_smear_alt);
  m_fopr = new AFopr_Smeared_alt<AFIELD>(m_kernel, m_dr_smear_alt);

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());
}
#endif

//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::tidyup()
{
  delete m_solver;
  delete m_fopr;
  if (m_kernel != 0) delete m_kernel;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert(Field& xq, const Field& b,
                                            int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "D_prec") {
    invert_D_prec(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_prec") {
    invert_DdagD_prec(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_D(Field& xq, const Field& b,
                                              int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  vout.paranoiac(m_vl, "nin = %d  nvol = %d nex = %d\n", nin, nvol, nex);

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_fopr->convert(abq, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_alt, abq, b);
    }
  }

  vout.detailed(m_vl, "%s: convert finished.\n", class_name.c_str());

  m_fopr->set_mode("D");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  diff = double(diff2);

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->reverse(xq, axq);
    } else {
      reverse(index_alt, xq, axq);
    }
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_DdagD(Field& xq, const Field& b,
                                                  int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->convert(abq, b);
    } else {
      convert(index_alt, abq, b);
    }
  }

  real_t diff2;

  m_fopr->set_mode("DdagD");

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  diff = double(diff2);

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->reverse(xq, axq);
    } else {
      reverse(index_alt, xq, axq);
    }
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_D_prec(Field& xq,
                                                   const Field& b,
                                                   int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  if (m_fopr->needs_convert()) {
    m_fopr->convert(abq, b);
  } else {
    convert(index_alt, abq, b);
  }

  m_fopr->set_mode("D_prec");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  m_flop_count += m_solver->flop_count();

  m_fopr->mult(abq, axq, "Prec");

  if (m_fopr->needs_convert()) {
    m_fopr->reverse(xq, abq);
  } else {
    reverse(index_alt, xq, abq);
  }

  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  // 1 solve + 2 mult; flop for the solve has been already accumulated
  m_flop_count += 2 * m_fopr->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_DdagD_prec(Field& xq,
                                                       const Field& b,
                                                       int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "invert_DdagD_prec is called.\n");

  m_timer.reset();
  m_timer.start();

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  if (m_fopr->needs_convert()) {
    m_fopr->convert(axq, b);
  } else {
    convert(index_alt, axq, b);
  }

  m_fopr->mult(abq, axq, "Precdag");

  m_fopr->set_mode("DdagD_prec");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  m_flop_count += m_solver->flop_count();

  m_fopr->mult(abq, axq, "Prec");

  if (m_fopr->needs_convert()) {
    m_fopr->reverse(xq, abq);
  } else {
    reverse(index_alt, xq, abq);
  }

  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  // 1 solve + 2 mult; flop for the solve has been already accumulated
  m_flop_count += 2 * m_fopr->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert(AFIELD& xq,
                                            const AFIELD& b,
                                            int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "D_prec") {
    invert_D_prec(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_prec") {
    invert_DdagD_prec(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_D(AFIELD& axq,
                                              const AFIELD& abq,
                                              int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  m_fopr->set_mode("D");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_DdagD(AFIELD& axq,
                                                  const AFIELD& abq,
                                                  int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  real_t diff2;

  m_fopr->set_mode("DdagD");

#pragma omp parallel
  {
    m_solver->solve(axq, abq, nconv, diff2);
  }
  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_D_prec(AFIELD& axq,
                                                   const AFIELD& abq,
                                                   int& nconv, double& diff)
{
  m_timer.reset();
  m_timer.start();

  int    nin  = m_fopr->field_nin();
  int    nvol = m_fopr->field_nvol();
  int    nex  = m_fopr->field_nex();
  AFIELD atmp(nin, nvol, nex);

  m_fopr->set_mode("D_prec");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(atmp, abq, nconv, diff2);
  }
  m_flop_count += m_solver->flop_count();
  m_fopr->mult(axq, atmp, "Prec");

  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  // 1 solve + 2 mult; flop for the solve has been already accumulated
  m_flop_count += 2 * m_fopr->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::invert_DdagD_prec(AFIELD& axq,
                                                       const AFIELD& abq,
                                                       int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "invert_DdagD_prec is called.\n");

  m_timer.reset();
  m_timer.start();

  int    nin  = m_fopr->field_nin();
  int    nvol = m_fopr->field_nvol();
  int    nex  = m_fopr->field_nex();
  AFIELD atmp(nin, nvol, nex);

  m_fopr->mult(axq, abq, "Precdag");

  m_fopr->set_mode("DdagD_prec");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(atmp, axq, nconv, diff2);
  }
  m_flop_count += m_solver->flop_count();

  m_fopr->mult(axq, atmp, "Prec");

  diff = double(diff2);

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  // 1 solve + 2 mult; flop for the solve has been already accumulated
  m_flop_count += 2 * m_fopr->flop_count();
}


//====================================================================
template<typename AFIELD>
double Fprop_alt_Standard_lex<AFIELD>::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::get_performance(double& flop_count,
                                                     double& elapsed_time)
{
  flop_count   = m_flop_count;
  elapsed_time = m_elapsed_time;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::report_performance()
{
  double flops  = m_flop_count / m_elapsed_time;
  double gflops = flops * 1.0e-9;

  //  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: solver performance:\n", class_name.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", m_elapsed_time);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n", m_flop_count);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_lex<AFIELD>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
  ThreadManager::assert_single_thread(class_name);

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
    for (int i = 0; i < Nrepeat; ++i) {
      m_fopr->mult(axq, abq);
      m_fopr->mult(abq, axq);
    }
  }

  timer->stop();

  double flop_fopr  = m_fopr->flop_count();
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops        = flop_total / elapsed_time;
  double gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  mult mode = %s\n", mode.c_str());
  vout.general(m_vl, "  Number of mult = %18d\n", 2 * Nrepeat);
  vout.general(m_vl, "  Elapsed time   = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(Fopr)     = %18.0f\n", flop_fopr);
  vout.general(m_vl, "  Flop(total)    = %18.0f\n", flop_total);
  vout.general(m_vl, "  Performance    = %11.3f GFlops\n", gflops);

  m_fopr->set_mode(mode_prev);
}


//============================================================END=====
