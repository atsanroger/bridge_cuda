/*!
      @file    fprop_alt_Standard_eo_Mixedprec.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

template<typename AFIELD_d, typename AFIELD_f>
const std::string Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::class_name
  = "Fprop_alt_Standard_eo_Mixedprec";
//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::init(const Parameters& params_fopr,
                                                               const Parameters& params_solver)
{
  // this constructor assumes that the factories are available.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup.\n", class_name.c_str());

  string fopr_type = params_fopr.get_string("fermion_type");
  fopr_type += "_eo";

  m_kernel_d = 0;
  m_fopr_d   = AFopr<AFIELD_d>::New(fopr_type, params_fopr);

  m_kernel_f = 0;
  m_fopr_f   = AFopr<AFIELD_f>::New(fopr_type, params_fopr);

  // solver parameters for preconditioner solver
  Parameters params_solver2;
  params_solver2.set_int("maximum_number_of_iteration", 200);
  params_solver2.set_int("maximum_number_of_restart", 1);
  params_solver2.set_double("convergence_criterion_squared", 1.0e-8);
  //params_solver2.set_double("convergence_criterion_squared", 1.0e-14);
  string str_vlevel = params_solver.get_string("verbose_level");
  params_solver2.set_string("verbose_level", str_vlevel);

  string solver_type = params_solver.get_string("solver_type");
  m_solver_prec = ASolver<AFIELD_f>::New(solver_type, m_fopr_f);
  m_solver_prec->set_parameters(params_solver2);

  m_precond
    = new APrecond_Mixedprec<AFIELD_d, AFIELD_f>(m_solver_prec);

  m_solver
    = new ASolver_FBiCGStab<AFIELD_d>(m_fopr_d, m_precond);
  //      ASolver_BiCGStab_Precond<AFIELD_d>(m_fopr_d, m_precond);

  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::init(const Parameters& params_fopr,
                                                               const Parameters& params_solver,
                                                               Director_Smear *dr_smear)
{
  // this constructor assumes that the factories are available.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup.\n", class_name.c_str());

  m_dr_smear = dr_smear;

  string fopr_type = params_fopr.get_string("fermion_type");
  fopr_type += "_eo";

  m_kernel_d = AFopr<AFIELD_d>::New(fopr_type, params_fopr);
  m_fopr_d   = new AFopr_Smeared<AFIELD_d>(m_kernel_d, m_dr_smear);

  m_kernel_f = AFopr<AFIELD_f>::New(fopr_type, params_fopr);
  m_fopr_f   = new AFopr_Smeared<AFIELD_f>(m_kernel_f, m_dr_smear);

  // solver parameters for preconditioner solver
  Parameters params_solver2;
  params_solver2.set_int("maximum_number_of_iteration", 200);
  params_solver2.set_int("maximum_number_of_restart", 1);
  params_solver2.set_double("convergence_criterion_squared", 1.0e-8);
  string str_vlevel = params_solver.get_string("verbose_level");
  params_solver2.set_string("verbose_level", str_vlevel);

  string solver_type = params_solver.get_string("solver_type");
  m_solver_prec = ASolver<AFIELD_f>::New(solver_type, m_fopr_f);
  m_solver_prec->set_parameters(params_solver2);

  m_precond
    = new APrecond_Mixedprec<AFIELD_d, AFIELD_f>(m_solver_prec);

  m_solver = new
             ASolver_BiCGStab_Precond<AFIELD_d>(m_fopr_d, m_precond);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::tidyup()
{
  delete m_solver;
  delete m_precond;
  delete m_solver_prec;
  delete m_fopr_f;
  delete m_fopr_d;
  if (m_kernel_d != 0) delete m_kernel_d;
  if (m_kernel_f != 0) delete m_kernel_f;
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::set_config(Field *U)
{
  m_fopr_d->set_config(U);
  m_fopr_f->set_config(U);
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert(Field& xq, const Field& b,
                                                                 int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert(AFIELD_d& xq, const AFIELD_d& b,
                                                                 int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert_D(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac("invert_D start.\n");

  Timer timer;
  timer.reset();
  timer.start();

  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD_d axq(nin, nvol, nex);
  AFIELD_d abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr_d->needs_convert()) {
      m_fopr_d->convert(abq, b);
    } else {
      convert(index_alt, abq, b);
    }
  }

  timer.stop();
  invert_D(axq, abq, nconv, diff);
  timer.start();
#pragma omp parallel
  {
    if (m_fopr_d->needs_convert()) {
      m_fopr_d->reverse(xq, axq);
    } else {
      reverse(index_alt, xq, axq);
    }
  }


  timer.stop();

  m_elapsed_time += timer.elapsed_sec();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert_DdagD(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD_d axq(nin, nvol, nex);
  AFIELD_d abq(nin, nvol, nex);

  AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;

#pragma omp parallel
  {
    if (m_fopr_d->needs_convert()) {
      m_fopr_d->convert(abq, b);
    } else {
      convert(index_alt, abq, b);
    }
    convert(index_alt, abq, b);
  }


  invert_DdagD(axq, abq, nconv, diff);

#pragma omp parallel
  {
    if (m_fopr_d->needs_convert()) {
      m_fopr_d->reverse(xq, axq);
    } else {
      reverse(index_alt, xq, axq);
    }
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert_D(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac("invert_D start.\n");

  m_timer.reset();
  m_timer.start();

  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD_d be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD_d xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD_d y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
  AIndex_eo<real_t, AFIELD_d::IMPL>  index_eo;
#pragma omp parallel
  {
    index_eo.split(be, bo, b);
  }

  invert_De(xe, xo, be, bo, nconv, diff);

#pragma omp parallel
  {
    index_eo.merge(xq, xe, xo);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert_DdagD(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac("invert_DdagD start.\n");

  m_timer.reset();
  m_timer.start();

  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();

  AFIELD_d be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD_d xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD_d y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  AIndex_eo<real_t, AFIELD_f::IMPL> index_eo;

#pragma omp parallel
  {
    index_eo.split(be, bo, b);
#pragma omp barrier

    m_fopr_d->mult_gm5(y1, be);
    m_fopr_d->mult_gm5(y2, bo);
  }

  int    nconv1;
  double diff1;
  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv = nconv1;
  diff  = diff1;

#pragma omp parallel
  {
    m_fopr_d->mult_gm5(y1, xe);
    m_fopr_d->mult_gm5(y2, xo);
  }

  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv += nconv1;
  diff  += diff1;

#pragma omp parallel
  {
    index_eo.merge(xq, xe, xo);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::invert_De(
  AFIELD_d& xe, AFIELD_d& xo,
  AFIELD_d& be, AFIELD_d& bo,
  int& nconv, double& diff)
{
  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD_d y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  // set even source vector.
#pragma omp parallel
  {
    m_fopr_d->mult(y1, bo, "Doo_inv");
#pragma omp barrier

    m_fopr_d->mult(y2, y1, "Deo");
#pragma omp barrier

    axpy(be, real_t(-1.0), y2);
#pragma omp barrier

    m_fopr_d->mult(y1, be, "Dee_inv");
  }

  m_fopr_d->set_mode("D");
  m_fopr_f->set_mode("D");

#pragma omp parallel
  {
    m_solver->solve(xe, y1, nconv, diff);
#pragma omp barrier

    m_fopr_d->mult(y1, xe, "Doe");
#pragma omp barrier

    aypx(real_t(-1.0), y1, bo);
#pragma omp barrier

    m_fopr_d->mult(xo, y1, "Doo_inv");
#pragma omp barrier
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
double Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::report_performance()
{
  double flops  = m_flop_count / m_elapsed_time;
  double gflops = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: solver performance:\n", class_name.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", m_elapsed_time);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n", m_flop_count);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Mixedprec<AFIELD_d, AFIELD_f>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
  int nin  = m_fopr_d->field_nin();
  int nvol = m_fopr_d->field_nvol();
  int nex  = m_fopr_d->field_nex();

  unique_ptr<Timer> timer(new Timer);

  std::string mode_prev_d = m_fopr_d->get_mode();
  std::string mode_prev_f = m_fopr_f->get_mode();

  m_fopr_d->set_mode(mode);
  m_fopr_f->set_mode(mode);

  {
    AFIELD_d axq(nin, nvol, nex), abq(nin, nvol, nex);
    abq.set(0.0);
    abq.set(0, 1.0);

    timer->start();

#pragma omp parallel
    {
      for (int i = 0; i < Nrepeat; ++i) {
        m_fopr_d->mult(axq, abq);
        m_fopr_d->mult(abq, axq);
      }
    }

    timer->stop();
  }

  double flop_fopr  = m_fopr_d->flop_count();
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops        = flop_total / elapsed_time;
  double gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  mult mode = %s\n", mode.c_str());
  vout.general(m_vl, "  Double precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n", flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n", flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

  {
    AFIELD_f axq(nin, nvol, nex), abq(nin, nvol, nex);
    abq.set(0.0);
    abq.set(0, 1.0);

    timer->reset();
    timer->start();

#pragma omp parallel
    {
      for (int i = 0; i < Nrepeat; ++i) {
        m_fopr_f->mult(axq, abq);
        m_fopr_f->mult(abq, axq);
      }
    }

    timer->stop();
  }

  flop_fopr  = m_fopr_f->flop_count();
  flop_total = flop_fopr * double(2 * Nrepeat);

  elapsed_time = timer->elapsed_sec();
  flops        = flop_total / elapsed_time;
  gflops       = flops * 1.0e-9;

  vout.general(m_vl, "  Single precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n", flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n", flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

  m_fopr_d->set_mode(mode_prev_d);
  m_fopr_f->set_mode(mode_prev_f);
}


//============================================================END=====
