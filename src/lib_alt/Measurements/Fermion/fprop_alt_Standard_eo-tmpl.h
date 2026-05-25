/*!
      @file    fprop_alt_Standard_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-11-26 18:23:37 #$
      @version $LastChangedRevision: 2554 $
*/

template<typename AFIELD>
const std::string Fprop_alt_Standard_eo<AFIELD>::class_name
  = "Fprop_alt_Standard_eo";
//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::init(const Parameters& params_fopr,
                                         const Parameters& params_solver)
{
  // this constructor assumes that the factories are available.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup (without link smearing).\n",
               class_name.c_str());
  vout.increase_indent();

  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  string fopr_type = params_fopr.get_string("fermion_type");
  if(fopr_type.substr(fopr_type.size()-3 ,3) != "_eo")
    fopr_type += "_eo";

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
  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::init(const Parameters& params_fopr,
                                         const Parameters& params_solver,
                                         Director_Smear *dr_smear)
{
  ThreadManager::assert_single_thread(class_name);

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup  (with link smearing).\n",
               class_name.c_str());
  vout.increase_indent();

  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  m_dr_smear = dr_smear;
  m_dr_smear_alt = 0;
  m_alt_director = false;

  string fopr_type = params_fopr.get_string("fermion_type");
  if(fopr_type.substr(fopr_type.size()-3 ,3) != "_eo")
    fopr_type += "_eo";

  m_kernel   = AltFopr::New(fopr_type, params_fopr);

  //  m_fopr = AltFopr::New("Smeared", m_kernel, m_dr_smear);
  m_fopr = new AFopr_Smeared<AFIELD>(m_kernel, m_dr_smear);

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.decrease_indent();
  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}


//====================================================================
#ifdef DIRECTOR_ALT_SMEARED_IMPLEMENTED
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::init(
                                const Parameters& params_fopr,
                                const Parameters& params_solver,
                                Director_alt_Smear<AFIELD> *dr_smear)
{
  ThreadManager::assert_single_thread(class_name);

  vout.general(m_vl, "%s: being setup  (with link smearing).\n",
               class_name.c_str());
  vout.increase_indent();

  typedef AFopr<AFIELD>     AltFopr;
  typedef ASolver<AFIELD>   AltSolver;

  m_dr_smear = 0;
  m_dr_smear_alt = dr_smear;
  m_alt_director = true;

  string fopr_type = params_fopr.get_string("fermion_type");
  if(fopr_type.substr(fopr_type.size()-3 ,3) != "_eo")
    fopr_type += "_eo";

  m_kernel   = AltFopr::New(fopr_type, params_fopr);

  //  m_fopr = AltFopr::New("Smeared_alt", m_kernel, m_dr_smear_alt);
  m_fopr = new AFopr_Smeared_alt<AFIELD>(m_kernel, m_dr_smear_alt);

  string solver_type = params_solver.get_string("solver_type");
  m_solver = AltSolver::New(solver_type, m_fopr);
  m_solver->set_parameters(params_solver);

  reset_performance();

  vout.decrease_indent();
  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());
}
#endif

//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::tidyup()
{
  delete m_solver;
  delete m_fopr;
  if (m_kernel != 0) delete m_kernel;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert(Field& xq, const Field& b,
                                           int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if ((m_mode == "D") || (m_mode == "D_prec")) {
    invert_D(xq, b, nconv, diff);
  } else if ((m_mode == "DdagD") || (m_mode == "DdagD_prec")) {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "D_even") {
    m_fopr->set_mode("D");
    invert_De(xq, b, nconv, diff);
  } else if (m_mode == "Ddag_even") {
    m_fopr->set_mode("Ddag");
    invert_De(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_even") {
    m_fopr->set_mode("DdagD");
    invert_De(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert(AFIELD& xq, const AFIELD& b,
                                           int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "D_even") {
    m_fopr->set_mode("D");
    invert_De(xq, b, nconv, diff);
  } else if (m_mode == "Ddag_even") {
    m_fopr->set_mode("Ddag");
    invert_De(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_even") {
    m_fopr->set_mode("DdagD");
    invert_De(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert_D(Field& xq, const Field& b,
                                             int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.detailed(m_vl, "invert_D(Field) start.\n");

  m_timer.reset();
  m_timer.start();

  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;
  AIndex_eo<real_t, AFIELD::IMPL>  index_eo;

  bool is_qdw = (m_fopr->get_mw_mode() == AFopr<AFIELD>::MWMode::DW);

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_fopr->convert(abq, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_alt, abq, b);
    }
#pragma omp barrier
    index_eo.split(be, bo, abq, is_qdw);
  }

  invert_De(xe, xo, be, bo, nconv, diff);

  vout.detailed(m_vl, "%s: diff = %e\n", class_name.c_str(), diff);

#pragma omp parallel
  {
    index_eo.merge(axq, xe, xo, is_qdw);
#pragma omp barrier

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
void Fprop_alt_Standard_eo<AFIELD>::invert_DdagD(Field& xq, const Field& b,
                                                 int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.detailed(m_vl, "invert_DdagD start.\n");

  m_timer.reset();
  m_timer.start();

  //  xq.set(0.0);

  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD axq(nin, nvol, nex);
  AFIELD abq(nin, nvol, nex);

  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;
  AIndex_eo<real_t, AFIELD::IMPL>  index_eo;

  bool is_qdw = (m_fopr->get_mw_mode() == AFopr<AFIELD>::MWMode::DW);

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_fopr->convert(abq, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_alt, abq, b);
    }
#pragma omp barrier

    index_eo.split(be, bo, abq, is_qdw);
#pragma omp barrier

    m_fopr->mult_gm5(y1, be);
    m_fopr->mult_gm5(y2, bo);
  }

  int    nconv1;
  double diff1;
  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv = nconv1;
  diff  = diff1;

#pragma omp parallel
  {
    m_fopr->mult_gm5(y1, xe);
    m_fopr->mult_gm5(y2, xo);
  }

  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv += nconv1;
  diff  += diff1;

#pragma omp parallel
  {
    index_eo.merge(axq, xe, xo, is_qdw);
#pragma omp barrier

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
void Fprop_alt_Standard_eo<AFIELD>::invert_D(AFIELD& xq, const AFIELD& b,
                                             int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.detailed(m_vl, "invert_D start.\n");

  m_timer.reset();
  m_timer.start();

  //  xq.set(0.0);

  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();

  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);

  AIndex_eo<real_t, AFIELD::IMPL> index_eo;

  bool is_qdw = (m_fopr->get_mw_mode() == AFopr<AFIELD>::MWMode::DW);

#pragma omp parallel
  {
    index_eo.split(be, bo, b, is_qdw);
  }

  invert_De(xe, xo, be, bo, nconv, diff);

#pragma omp parallel
  {
    index_eo.merge(xq, xe, xo, is_qdw);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert_De(Field& xq, const Field& b,
                                              int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.detailed(m_vl, "invert_De(Field) start.\n");

  m_timer.reset();
  m_timer.start();

  //  xq.set(0.0);

  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();

  AFIELD abq(nin, nvol2, nex);
  AFIELD axq(nin, nvol2, nex);

  AIndex_eo<real_t, AFIELD::IMPL> index_eo;

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_fopr->convert(abq, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_eo, abq, b);
    }
  }

  invert_De(axq, abq, nconv, diff);

#pragma omp parallel
  {
    if (m_fopr->needs_convert()) {
      m_fopr->reverse(xq, axq);
    } else {
      reverse(index_eo, xq, axq);
    }
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert_DdagD(AFIELD& xq, const AFIELD& b,
                                                 int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.detailed(m_vl, "invert_DdagD start.\n");

  m_timer.reset();
  m_timer.start();

  //  xq.set(0.0);

  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();

  AFIELD be(nin, nvol2, nex), bo(nin, nvol2, nex);
  AFIELD xe(nin, nvol2, nex), xo(nin, nvol2, nex);
  AFIELD y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  AIndex_eo<real_t, AFIELD::IMPL> index_eo;

  bool is_qdw = (m_fopr->get_mw_mode() == AFopr<AFIELD>::MWMode::DW);

#pragma omp parallel
  {
    index_eo.split(be, bo, b, is_qdw);
#pragma omp barrier

    m_fopr->mult_gm5(y1, be);
    m_fopr->mult_gm5(y2, bo);
  }

  int    nconv1;
  double diff1;
  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv = nconv1;
  diff  = diff1;

#pragma omp parallel
  {
    m_fopr->mult_gm5(y1, xe);
    m_fopr->mult_gm5(y2, xo);
  }

  invert_De(xe, xo, y1, y2, nconv1, diff1);

  nconv += nconv1;
  diff  += diff1;

#pragma omp parallel
  {
    index_eo.merge(xq, xe, xo, is_qdw);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert_De(AFIELD& xe, AFIELD& xo,
                                              AFIELD& be, AFIELD& bo,
                                              int& nconv, double& diff)
{
  int nin   = m_fopr->field_nin();
  int nvol2 = m_fopr->field_nvol();
  int nex   = m_fopr->field_nex();
  int nvol  = 2 * nvol2;

  AFIELD y1(nin, nvol2, nex), y2(nin, nvol2, nex);

  vout.detailed(m_vl, "invert_De(AFIELD)(6arg) start.\n");

#pragma omp parallel
  {
    // set even source vector.
    m_fopr->mult(y1, bo, "Doo_inv");

#pragma omp barrier
    m_fopr->mult(y2, y1, "Deo");

#pragma omp barrier
    axpy(be, real_t(-1.0), y2);

#pragma omp barrier
    m_fopr->mult(y1, be, "Dee_inv");
  }

    m_fopr->set_mode("D");

  real_t diff2;

#pragma omp parallel
  {
    m_solver->solve(xe, y1, nconv, diff2);
#pragma omp barrier

    m_fopr->normalize_fprop(xe);

    m_fopr->mult(y1, xe, "Doe");
#pragma omp barrier

    aypx(real_t(-1.0), y1, bo);
#pragma omp barrier

    m_fopr->mult(xo, y1, "Doo_inv");

  }
  
  diff = double(diff2);
  vout.detailed(m_vl, "diff(invert_De) = %e\n", diff);

}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::invert_De(AFIELD& xe,
                                              const AFIELD& be,
                                              int& nconv, double& diff)
{
  real_t diff2;

  vout.detailed("invert_De(AFIELD)(4 arg) start.\n");

#pragma omp parallel
  {
    m_solver->solve(xe, be, nconv, diff2);
  }

  diff = double(diff2);

}


//====================================================================
template<typename AFIELD>
double Fprop_alt_Standard_eo<AFIELD>::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::report_performance()
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
template<typename AFIELD>
void Fprop_alt_Standard_eo<AFIELD>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
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
