template<typename AFIELD>
const std::string ASolver_CG<AFIELD>::class_name = "ASolver_CG";
//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_x.reset(nin, nvol, nex);
  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);

  //  m_vl = Bridge::DETAILED;

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  InitialGuess init_guess_mode = InitialGuess::RHS;
  if (params.find_string("initial_guess_mode")) {
    const string initial_guess_mode = params.get_string("initial_guess_mode");
    vout.detailed(m_vl, "  initila_guess_mode %s\n", initial_guess_mode.c_str());
    if (initial_guess_mode == "RHS") {
      init_guess_mode = InitialGuess::RHS;
    } else if (initial_guess_mode == "GIVEN") {
      init_guess_mode = InitialGuess::GIVEN;
    } else if (initial_guess_mode == "ZERO") {
      init_guess_mode = InitialGuess::ZERO;
    } else {
      vout.crucial(m_vl, "Error at %s: unknown initial guess mode, %s\n", class_name.c_str(), initial_guess_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond, init_guess_mode);

  m_mw_mode = MWMode::FP;
  std::string mw_mode_str = "FP";
  params.fetch_string("multiword_mode", mw_mode_str);
  if (mw_mode_str == "DW") {
    m_mw_mode = MWMode::DW;
  } else if (mw_mode_str == "TW") {
    m_mw_mode = MWMode::TW;
  } else if (mw_mode_str != "FP") {
    vout.crucial(m_vl, "Error at %s: unknown multiword_mode '%s' (use FP/DW/TW)\n",
                 class_name.c_str(), mw_mode_str.c_str());
    exit(EXIT_FAILURE);
  }

  // Multiword is only implemented by even-odd fopr. Fopr that do not override
  // set_mw_mode (e.g. lexical Domainwall_5din) keep the base no-op and keep
  // reporting FP; driving them with enlarged QDW solver fields produces
  // illegal device access (CUDA 700). Probe support and fall back to FP.
  if (m_mw_mode != MWMode::FP && m_fopr) {
    m_fopr->set_mw_mode(m_mw_mode);
    const bool supported = (m_fopr->get_mw_mode() == m_mw_mode);
    m_fopr->set_mw_mode(MWMode::FP);  // restore so field_nin() reports FP size
    if (!supported) {
      vout.general(m_vl, "%s: fopr does not support multiword_mode=%s; "
                   "falling back to FP.\n", class_name.c_str(),
                   mw_mode_str.c_str());
      m_mw_mode = MWMode::FP;
    }
  }

  if (m_mw_mode == MWMode::DW) {
    int nin = 2 * m_fopr->field_nin();
    int nvol = m_fopr->field_nvol();
    int nex = m_fopr->field_nex();
    m_x.reset(nin, nvol, nex);
    m_r.reset(nin, nvol, nex);
    m_p.reset(nin, nvol, nex);
    m_s.reset(nin, nvol, nex);
  } else if (m_mw_mode == MWMode::TW) {
    int nin = 3 * m_fopr->field_nin();
    int nvol = m_fopr->field_nvol();
    int nex = m_fopr->field_nex();
    m_x.reset(nin, nvol, nex);
    m_r.reset(nin, nvol, nex);
    m_p.reset(nin, nvol, nex);
    m_s.reset(nin, nvol, nex);
  }

  if (m_fopr) {
    m_fopr->set_mw_mode(m_mw_mode);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond,
                                        const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter        = Niter;
  m_Stop_cond    = Stop_cond;
  m_initial_mode = init_guess_mode;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                               int& Nconv, real_t& diff)
{
  copy(m_s, b);

  real_t sr    = m_s.norm2(static_cast<int>(m_mw_mode));
  real_t snorm = 1.0 / sr;
  vout.detailed(m_vl, "  snorm = %22.15e\n", snorm);

  real_t rr, rrp;
  int    nconv = -1;

  const bool use_pair   = (m_mw_mode == MWMode::DW);
  const bool use_triple = (m_mw_mode == MWMode::TW);

  if (use_triple) {
    solve_CG_init_triple(rr);
  } else if (use_pair) {
    solve_CG_init_pair(rr);
  } else {
    solve_CG_init(rrp, rr);
  }
  vout.detailed(m_vl, "  init: %22.15e\n", rr * snorm);

  for (int iter = 0; iter < m_Niter; ++iter) {
    if (use_triple) {
      solve_CG_step_triple(rr);
    } else if (use_pair) {
      solve_CG_step_pair(rr);
    } else {
      solve_CG_step(rrp, rr);
    }
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr * snorm);

    if (rr * snorm < m_Stop_cond) {
      nconv = iter;
#pragma omp master
      {
        m_nconv = nconv + 1;
      }
      break;
    }
  }

  if (nconv == -1)  {
    vout.crucial(m_vl, "Error at %s: not converged\n",
                 class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                 m_Niter, rr * snorm);
#pragma omp master
    {
      m_nconv = nconv;
    }
#pragma omp barrier
    //exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "converged:\n");
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  copy(xq, m_x);

  AFIELD s_temp(b.nin(), b.nvol(), b.nex());
  m_fopr->mult(s_temp, xq);

  s_temp.axpy(real_t(-1.0), b, false);
  real_t diff2 = s_temp.norm2(false);

#pragma omp master
  {
    diff  = diff2;
    Nconv = m_nconv;
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init(real_t& rrp, real_t& rr)
{
  if (m_initial_mode == InitialGuess::RHS) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::RHS\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->set_mode("DdagD");
    vout.general(m_vl, "ASolver_CG solve_CG_init: m_x.nin() = %d, m_s.nin() = %d\n", m_x.nin(), m_s.nin());
    m_fopr->mult(m_s, m_x);
    m_s.normalize(static_cast<int>(m_mw_mode));
    m_r.axpy(real_t(-1.0), m_s, static_cast<int>(m_mw_mode));
    m_r.normalize(static_cast<int>(m_mw_mode));
    copy(m_p, m_r);
    rr  = m_r.norm2(static_cast<int>(m_mw_mode));
    rrp = rr;
  } else if (m_initial_mode == InitialGuess::GIVEN) {
    vout.crucial("%s: InitialGuess::GIVEN is not yet ready\n", class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (m_initial_mode == InitialGuess::ZERO) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::ZERO\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    m_s.set(0.0);
    m_x.set(0.0);
    copy(m_p, m_r);
    rr  = m_r.norm2(static_cast<int>(m_mw_mode));
    rrp = rr;
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step(real_t& rrp, real_t& rr)
{
  using complex_t = typename AFIELD::complex_t;

  m_fopr->mult(m_s, m_p);
  m_s.normalize(static_cast<int>(m_mw_mode));

  real_t pap = m_s.dot(m_p, static_cast<int>(m_mw_mode));
  //    m_fopr->mult_normA_dev(pap, m_s, m_p);
  real_t cr = rrp / pap;

  m_x.axpy(cr, m_p, static_cast<int>(m_mw_mode));
  m_x.normalize(static_cast<int>(m_mw_mode));

  m_r.axpy(-cr, m_s, static_cast<int>(m_mw_mode));
  m_r.normalize(static_cast<int>(m_mw_mode));
  rr = m_r.norm2(static_cast<int>(m_mw_mode));

  real_t bk = rr / rrp;

  aypx(bk, m_p, m_r, static_cast<int>(m_mw_mode));
  m_p.normalize(static_cast<int>(m_mw_mode));

  rrp = rr;
}


//====================================================================
// DD-pair CG: rrp, rr, cr, bk and pap are all (h, l) float-pairs. We use
// dw_add/dw_mul/dw_div on real_t (= float in the float build) so no FP64
// ALU is touched. snorm/rr*snorm at the host level is still computed in
// real_t since it only drives the convergence check.
namespace dd_priv {
template<typename T>
inline void dd_two_sum(T a, T b, T& s, T& e) {
    s = a + b;
    T v = s - a;
    e = (a - (s - v)) + (b - v);
}
template<typename T>
inline void dd_two_prod(T a, T b, T& p, T& e) {
    p = a * b;
    e = std::fma(a, b, -p);
}
// (a_h, a_l) / (b_h, b_l) via 2-iteration Newton refinement on the FP32 quotient.
template<typename T>
inline void dd_div(T a_h, T a_l, T b_h, T b_l, T& q_h, T& q_l) {
    T q1 = a_h / b_h;
    // r = a - q1*b   (DD)
    T p_h, p_e;
    dd_two_prod(q1, b_h, p_h, p_e);
    T t  = std::fma(q1, b_l, p_e);
    T r_h, r_e;
    dd_two_sum(a_h, -p_h, r_h, r_e);
    T r_l = a_l - t + r_e;
    T q2  = (r_h + r_l) / b_h;
    // q = q1 + q2 (renormalize)
    dd_two_sum(q1, q2, q_h, q_l);
}
} // namespace dd_priv

template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init_pair(real_t& rr)
{
  if (m_initial_mode == InitialGuess::RHS) {
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->set_mode("DdagD");
    m_fopr->mult(m_s, m_x);
    m_s.normalize(static_cast<int>(m_mw_mode));
    // m_r -= m_s : single-rep is fine here, b - A*x is bounded.
    m_r.axpy(real_t(-1.0), m_s, static_cast<int>(m_mw_mode));
    m_r.normalize(static_cast<int>(m_mw_mode));
    copy(m_p, m_r);
    real_t rr_h, rr_l;
    m_r.norm2_pair(rr_h, rr_l, static_cast<int>(m_mw_mode));
    m_rrp_h = rr_h;
    m_rrp_l = rr_l;
    rr      = rr_h + rr_l;
  } else if (m_initial_mode == InitialGuess::ZERO) {
    copy(m_r, m_s);
    m_s.set(0.0);
    m_x.set(0.0);
    copy(m_p, m_r);
    real_t rr_h, rr_l;
    m_r.norm2_pair(rr_h, rr_l, static_cast<int>(m_mw_mode));
    m_rrp_h = rr_h;
    m_rrp_l = rr_l;
    rr      = rr_h + rr_l;
  } else {
    vout.crucial("%s: unsupported initial guess mode in pair path\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step_pair(real_t& rr)
{
  m_fopr->mult(m_s, m_p);
  m_s.normalize(static_cast<int>(m_mw_mode));

  // pap = <s, p> as DD pair
  real_t pap_h, pap_l;
  m_s.dot_pair(pap_h, pap_l, m_p, static_cast<int>(m_mw_mode));

  // cr = rrp / pap (DD division, FP32-only)
  real_t cr_h, cr_l;
  dd_priv::dd_div(m_rrp_h, m_rrp_l, pap_h, pap_l, cr_h, cr_l);

  // m_x += cr * m_p  (DD axpy)
  m_x.axpy_pair(cr_h, cr_l, m_p, static_cast<int>(m_mw_mode));
  m_x.normalize(static_cast<int>(m_mw_mode));

  // m_r -= cr * m_s
  m_r.axpy_pair(-cr_h, -cr_l, m_s, static_cast<int>(m_mw_mode));
  m_r.normalize(static_cast<int>(m_mw_mode));

  // rr = <r, r> (DD)
  real_t rr_h, rr_l;
  m_r.norm2_pair(rr_h, rr_l, static_cast<int>(m_mw_mode));
  rr = rr_h + rr_l;

  // bk = rr / rrp
  real_t bk_h, bk_l;
  dd_priv::dd_div(rr_h, rr_l, m_rrp_h, m_rrp_l, bk_h, bk_l);

  // m_p = bk * m_p + m_r
  m_p.aypx_pair(bk_h, bk_l, m_r, static_cast<int>(m_mw_mode));
  m_p.normalize(static_cast<int>(m_mw_mode));

  m_rrp_h = rr_h;
  m_rrp_l = rr_l;
}


//====================================================================
// TW-triple CG: rrp, rr, cr, bk and pap are all (h, m, l) float-triples.
// Triple-word arithmetic via Hida QD-library style renormalize-after-each-op,
// FP32-only — never promotes to FP64. The CG inner-product BLAS layer
// (norm2_tw3 / dot_tw3) yields per-rank 3-word reductions; the device kernels
// preserve (h, m, l) at red[0..2] so the host-side coefficient algebra stays
// at ~72-bit precision.
namespace tw_priv {
// Error-free transforms (host). Same role as the device TwoSum/TwoProd.
template<typename T>
inline void tps(T a, T b, T& s, T& e) { s = a + b; T v = s - a; e = (a - (s - v)) + (b - v); }
template<typename T>
inline void tpp(T a, T b, T& p, T& e) { p = a * b; e = std::fma(a, b, -p); }

// Renormalize: enforce |m| <= ulp(h)/2, |l| <= ulp(m)/2.
template<typename T>
inline void tw_three_sum(T a, T b, T c, T& s_h, T& s_m, T& s_l) {
  // Two TwoSums + a fold: classic Hida 3-sum.
  T t1, e1;
  t1 = a + b; { T v = t1 - a; e1 = (a - (t1 - v)) + (b - v); }
  T t2, e2;
  t2 = t1 + c; { T v = t2 - t1; e2 = (t1 - (t2 - v)) + (c - v); }
  s_h = t2;
  T s2 = e1 + e2;        // no renorm needed — leading bits already separated.
  T sv = s2 - e1; T e3 = (e1 - (s2 - sv)) + (e2 - sv);
  s_m = s2;
  s_l = e3;
}

// (a_h, a_m, a_l) + (b_h, b_m, b_l) → (s_h, s_m, s_l), renormalized.
template<typename T>
inline void tw_add(T a_h, T a_m, T a_l, T b_h, T b_m, T b_l,
                   T& s_h, T& s_m, T& s_l) {
  T u_h, u_m; tps(a_h, b_h, u_h, u_m);
  T v_h, v_m; tps(a_m, b_m, v_h, v_m);
  T w_h = a_l + b_l;
  // Capture the eps^2 fold error of (u_m + v_h) instead of dropping it —
  // this is what makes the scalar genuinely triple-word, not dual.
  T mid, em; tps(u_m, v_h, mid, em);
  T lo = (v_m + w_h) + em;
  tw_three_sum(u_h, mid, lo, s_h, s_m, s_l);
}

// (a_h, a_m, a_l) * b (real, single-word) → (p_h, p_m, p_l) via TwoProd cascade.
template<typename T>
inline void tw_mul_real(T a_h, T a_m, T a_l, T b,
                        T& p_h, T& p_m, T& p_l) {
  T ph, pe; tpp(a_h, b, ph, pe);
  T mh, me; tpp(a_m, b, mh, me);
  T lh = a_l * b;
  T mid, em; tps(pe, mh, mid, em);     // capture fold error
  T lo = (me + lh) + em;
  tw_three_sum(ph, mid, lo, p_h, p_m, p_l);
}

// (a_h, a_m, a_l) * (b_h, b_m, b_l) → (p_h, p_m, p_l), genuine ~72-bit product.
template<typename T>
inline void tw_mul(T a_h, T a_m, T a_l, T b_h, T b_m, T b_l,
                   T& p_h, T& p_m, T& p_l) {
  T hh, he; tpp(a_h, b_h, hh, he);
  // Form the two mid cross-terms error-free so their eps^2 rounding lands in lo.
  T cm1_h, cm1_l; tpp(a_h, b_m, cm1_h, cm1_l);
  T cm2_h, cm2_l; tpp(a_m, b_h, cm2_h, cm2_l);
  T s1, e1; tps(cm1_h, cm2_h, s1, e1);
  T mid, e2; tps(he, s1, mid, e2);
  T lo = cm1_l + cm2_l + e1 + e2
       + std::fma(a_h, b_l, std::fma(a_m, b_m, std::fma(a_l, b_h, T(0))));
  tw_three_sum(hh, mid, lo, p_h, p_m, p_l);
}

// (a_h, a_m, a_l) / (b_h, b_m, b_l) via 2-step Newton refinement.
// FP32-only on real_t.
template<typename T>
inline void tw_div(T a_h, T a_m, T a_l, T b_h, T b_m, T b_l,
                   T& q_h, T& q_m, T& q_l) {
  T q1 = a_h / b_h;
  // r = a - q1*b
  T qb_h, qb_m, qb_l;
  tw_mul_real(b_h, b_m, b_l, q1, qb_h, qb_m, qb_l);
  T r_h, r_m, r_l;
  tw_add(a_h, a_m, a_l, -qb_h, -qb_m, -qb_l, r_h, r_m, r_l);

  T q2 = r_h / b_h;
  T qb2_h, qb2_m, qb2_l;
  tw_mul_real(b_h, b_m, b_l, q2, qb2_h, qb2_m, qb2_l);
  T r2_h, r2_m, r2_l;
  tw_add(r_h, r_m, r_l, -qb2_h, -qb2_m, -qb2_l, r2_h, r2_m, r2_l);

  T q3 = r2_h / b_h;

  // q = q1 + q2 + q3 (renormalized to (h, m, l)).
  tw_three_sum(q1, q2, q3, q_h, q_m, q_l);
}
} // namespace tw_priv

template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init_triple(real_t& rr)
{
  if (m_initial_mode == InitialGuess::RHS) {
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->set_mode("DdagD");
    m_fopr->mult(m_s, m_x);
    m_s.normalize(static_cast<int>(m_mw_mode));
    m_r.axpy(real_t(-1.0), m_s, static_cast<int>(m_mw_mode));
    m_r.normalize(static_cast<int>(m_mw_mode));
    copy(m_p, m_r);
    real_t rr_h, rr_m, rr_l;
    m_r.norm2_triple(rr_h, rr_m, rr_l, static_cast<int>(m_mw_mode));
    m_rrp_th = rr_h;
    m_rrp_tm = rr_m;
    m_rrp_tl = rr_l;
    rr       = rr_h + rr_m + rr_l;
  } else if (m_initial_mode == InitialGuess::ZERO) {
    copy(m_r, m_s);
    m_s.set(0.0);
    m_x.set(0.0);
    copy(m_p, m_r);
    real_t rr_h, rr_m, rr_l;
    m_r.norm2_triple(rr_h, rr_m, rr_l, static_cast<int>(m_mw_mode));
    m_rrp_th = rr_h;
    m_rrp_tm = rr_m;
    m_rrp_tl = rr_l;
    rr       = rr_h + rr_m + rr_l;
  } else {
    vout.crucial("%s: unsupported initial guess mode in triple path\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step_triple(real_t& rr)
{
  m_fopr->mult(m_s, m_p);
  m_s.normalize(static_cast<int>(m_mw_mode));

  // pap = <s, p> as TW triple
  real_t pap_h, pap_m, pap_l;
  m_s.dot_triple(pap_h, pap_m, pap_l, m_p, static_cast<int>(m_mw_mode));

  // cr = rrp / pap (TW division, FP32-only)
  real_t cr_h, cr_m, cr_l;
  tw_priv::tw_div(m_rrp_th, m_rrp_tm, m_rrp_tl,
                  pap_h, pap_m, pap_l,
                  cr_h, cr_m, cr_l);

  // m_x += cr * m_p
  m_x.axpy_triple(cr_h, cr_m, cr_l, m_p, static_cast<int>(m_mw_mode));
  m_x.normalize(static_cast<int>(m_mw_mode));

  // m_r -= cr * m_s
  m_r.axpy_triple(-cr_h, -cr_m, -cr_l, m_s, static_cast<int>(m_mw_mode));
  m_r.normalize(static_cast<int>(m_mw_mode));

  // rr = <r, r> (TW)
  real_t rr_h, rr_m, rr_l;
  m_r.norm2_triple(rr_h, rr_m, rr_l, static_cast<int>(m_mw_mode));
  rr = rr_h + rr_m + rr_l;

  // bk = rr / rrp
  real_t bk_h, bk_m, bk_l;
  tw_priv::tw_div(rr_h, rr_m, rr_l,
                  m_rrp_th, m_rrp_tm, m_rrp_tl,
                  bk_h, bk_m, bk_l);

  // m_p = bk * m_p + m_r
  m_p.aypx_triple(bk_h, bk_m, bk_l, m_r, static_cast<int>(m_mw_mode));
  m_p.normalize(static_cast<int>(m_mw_mode));

  m_rrp_th = rr_h;
  m_rrp_tm = rr_m;
  m_rrp_tl = rr_l;
}


//====================================================================
template<typename AFIELD>
double ASolver_CG<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field  = static_cast<double>(Nin * Nvol * Nex) * NPE;
  double flop_vector = (6 + ninit * 4 + m_nconv * 11) * flop_field;
  double flop_fopr   = (1 + ninit + m_nconv) * m_fopr->flop_count();

  double flop = flop_vector + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
