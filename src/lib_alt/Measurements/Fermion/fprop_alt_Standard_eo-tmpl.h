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

  params_solver.fetch_int("refinement_iterations", m_n_refine);
  vout.general(m_vl, "  refinement_iterations = %d\n", m_n_refine);

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

  params_solver.fetch_int("refinement_iterations", m_n_refine);
  vout.general(m_vl, "  refinement_iterations = %d\n", m_n_refine);

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

  params_solver.fetch_int("refinement_iterations", m_n_refine);
  vout.general(m_vl, "  refinement_iterations = %d\n", m_n_refine);

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
  AFIELD be0(nin, nvol2, nex);  // saved even RHS (invert_De destroys be)

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;
  AIndex_eo<real_t, AFIELD::IMPL>  index_eo;

  int mw_mode = static_cast<int>(m_fopr->get_mw_mode());

  // Full-system refinement passes. The inner eo solves (initial + corrections)
  // must run WITHOUT their own preconditioned (kappa^2) refinement: on the tiny
  // correction RHS that nested loop divides by a near-zero residual norm and
  // produces NaN. Accuracy here is driven by the full-system loop below instead.
  const int n_full = m_n_refine;
  m_n_refine = 0;

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
    index_eo.split(be, bo, abq, mw_mode);
    if (n_full > 0) be0.copy(be);
  }

  invert_De(xe, xo, be, bo, nconv, diff);

  vout.detailed(m_vl, "%s: diff = %e\n", class_name.c_str(), diff);

  // --- full-system iterative refinement on M x = b (sees kappa, not kappa^2).
  // r = b - M x via the now-all-TW forward sub-ops (Dee/Deo/Doe/Doo), then solve
  // M dx = r with the same eo solver and accumulate. Residual is full multiword
  // precision so the float-triple path can drive the solution past FP64. -------
  if (n_full > 0) {
    AFIELD Mxe(nin, nvol2, nex), Mxo(nin, nvol2, nex), tt(nin, nvol2, nex);
    AFIELD re(nin, nvol2, nex),  ro(nin, nvol2, nex);
    AFIELD dxe(nin, nvol2, nex), dxo(nin, nvol2, nex);

    // Source norm, to test the relative residual. Once the base solve is at its
    // floor, r = b - Mx is pure underflow noise (~1e-37); refining it makes the
    // inner CG normalize by 1/||r|| on denormals -> NaN. Stop when ||r||^2/||b||^2
    // drops below this (i.e. nothing left to correct). When kappa is large (small
    // quark mass) the residual stays well above this and refinement proceeds.
    const real_t rr_rel_floor = real_t(1.0e-22);
    real_t bnorm = 0;
#pragma omp parallel
    {
      real_t bn = be0.norm2(mw_mode) + bo.norm2(mw_mode);
#pragma omp master
      bnorm = bn;
    }

    for (int ir = 0; ir < n_full; ++ir) {
      int    nc_r = 0;
      double df_r = 0.0;
      real_t rnorm = 0;
#pragma omp parallel
      {
        m_fopr->mult(Mxe, xe, "Dee");
#pragma omp barrier
        m_fopr->mult(tt,  xo, "Deo");
#pragma omp barrier
        Mxe.axpy(real_t(1.0), tt, mw_mode);     // Mxe = Dee xe + Deo xo
#pragma omp barrier
        m_fopr->mult(Mxo, xo, "Doo");
#pragma omp barrier
        m_fopr->mult(tt,  xe, "Doe");
#pragma omp barrier
        Mxo.axpy(real_t(1.0), tt, mw_mode);     // Mxo = Doe xe + Doo xo
#pragma omp barrier
        re.copy(be0); re.axpy(real_t(-1.0), Mxe, mw_mode);   // re = be - (Mx)_e
        ro.copy(bo);  ro.axpy(real_t(-1.0), Mxo, mw_mode);   // ro = bo - (Mx)_o
#pragma omp barrier
        real_t rn = re.norm2(mw_mode) + ro.norm2(mw_mode);
#pragma omp master
        rnorm = rn;
      }
      if (rnorm < rr_rel_floor * bnorm) {
#pragma omp master
        vout.general(m_vl, "  full-refine: ||r||^2/||b||^2=%.3e below floor, "
                     "solve already converged; stopping.\n",
                     double(bnorm > 0 ? rnorm / bnorm : 0.0));
        break;
      }
      invert_De(dxe, dxo, re, ro, nc_r, df_r);  // M dx = r
#pragma omp parallel
      {
        xe.axpy(real_t(1.0), dxe, mw_mode);
        xo.axpy(real_t(1.0), dxo, mw_mode);
      }
#pragma omp master
      vout.general(m_vl, "  full-refine[%d]: inner nconv=%d diff=%.3e\n",
                   ir, nc_r, df_r);
    }
  }
  m_n_refine = n_full;  // restore the configured refinement count

#pragma omp parallel
  {
    index_eo.merge(axq, xe, xo, mw_mode);
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

  int mw_mode = static_cast<int>(m_fopr->get_mw_mode());

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

    index_eo.split(be, bo, abq, mw_mode);
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
    index_eo.merge(axq, xe, xo, mw_mode);
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

  int mw_mode = static_cast<int>(m_fopr->get_mw_mode());

#pragma omp parallel
  {
    index_eo.split(be, bo, b, mw_mode);
  }

  invert_De(xe, xo, be, bo, nconv, diff);

#pragma omp parallel
  {
    index_eo.merge(xq, xe, xo, mw_mode);
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

  int mw_mode = static_cast<int>(m_fopr->get_mw_mode());

#pragma omp parallel
  {
    index_eo.split(be, bo, b, mw_mode);
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
    index_eo.merge(xq, xe, xo, mw_mode);
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
  AFIELD rr(nin, nvol2, nex), dd(nin, nvol2, nex);  // iterative-refinement temps

  // multiword mode (0=FP, 2=TW): the vector axpy/aypx below must run in this
  // mode so the eo RHS and the refinement residual are full triple-word, not
  // single-word (otherwise refinement is capped at ~double).
  const int mw_mode = static_cast<int>(m_fopr->get_mw_mode());

  vout.detailed(m_vl, "invert_De(AFIELD)(6arg) start.\n");

#pragma omp parallel
  {
    // set even source vector.
    m_fopr->mult(y1, bo, "Doo_inv");

#pragma omp barrier
    m_fopr->mult(y2, y1, "Deo");

#pragma omp barrier
    be.axpy(real_t(-1.0), y2, mw_mode);

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

    // --- iterative refinement on the eo-preconditioned system D_hat xe = y1 ---
    // Residual r = y1 - D_hat xe is recomputed in the field's own (multiword)
    // precision, so each pass restores accuracy the CGNR normal equations lose
    // to kappa^2. FP32-only for the QTW path (D_hat = D_eo + LU_inv, no double).
    for (int ir = 0; ir < m_n_refine; ++ir) {
#pragma omp barrier
      m_fopr->mult(rr, xe);                 // rr = D_hat xe  (mode "D")
#pragma omp barrier
      rr.aypx(real_t(-1.0), y1, mw_mode);   // rr = y1 - D_hat xe  (TW residual)
#pragma omp barrier
      real_t rnorm = rr.norm2(mw_mode);     // residual norm before correcting
      int    nc_r;
      real_t df_r;
      m_solver->solve(dd, rr, nc_r, df_r);  // D_hat dd = rr
#pragma omp barrier
      m_fopr->normalize_fprop(dd);
      xe.axpy(real_t(1.0), dd, mw_mode);    // xe += dd  (TW)
#pragma omp master
      vout.general(m_vl, "  refine[%d]: ||r||^2=%.3e  inner nconv=%d diff=%.3e\n",
                   ir, double(rnorm), nc_r, double(df_r));
    }
#pragma omp barrier

    m_fopr->mult(y1, xe, "Doe");
#pragma omp barrier

    y1.aypx(real_t(-1.0), bo, mw_mode);
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
