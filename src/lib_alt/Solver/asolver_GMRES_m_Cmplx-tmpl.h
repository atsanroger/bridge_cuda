template<typename AFIELD>
const std::string ASolver_GMRES_m_Cmplx<AFIELD>::class_name = "ASolver_GMRES_m_Cmplx";

namespace {
  template<typename AFIELD>
    std::string InitialGuess_str(typename ASolver<AFIELD>::InitialGuess init_guess){
    using InitialGuess = typename ASolver<AFIELD>::InitialGuess;
    switch(init_guess) {
    case InitialGuess::RHS:
      return "RHS";
      break;
    case InitialGuess::GIVEN:
      return "GIVEN";
      break;
    case InitialGuess::ZERO:
      return "ZERO";
      break;
    default:
      break;  // nothing to do
    }
  return "unknown";
  }

  template<typename AFIELD>
    int fetch_initial_guess(const Parameters &params,
                            const std::string& key,
                            typename ASolver<AFIELD>::InitialGuess &init_guess_mode){
    if (!params.find_string(key)) {
      return EXIT_FAILURE;
    }
    using InitialGuess = typename ASolver<AFIELD>::InitialGuess;
    const string initial_guess_mode = params.get_string(key);
    vout.detailed("  initila_guess_mode %s\n", initial_guess_mode.c_str());
    if (initial_guess_mode == "RHS") {
      init_guess_mode = InitialGuess::RHS;
    } else if (initial_guess_mode == "GIVEN") {
      init_guess_mode = InitialGuess::GIVEN;
    } else if (initial_guess_mode == "ZERO") {
      init_guess_mode = InitialGuess::ZERO;
    } else {
      vout.crucial("Error: unsupported initial guess mode, %s\n", initial_guess_mode.c_str());
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

}

//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_initial_mode = InitialGuess::RHS;

  m_s.reset(nin, nvol, nex);
  m_r.reset(nin, nvol, nex);
  m_x.reset(nin, nvol, nex);
  m_v_tmp.reset(nin, nvol, nex);

}

//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::set_parameters(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;  // fetching as double
  int    N_M;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_int("number_of_orthonormal_vectors", N_M);

  InitialGuess init_guess_mode;
  int err1 = fetch_initial_guess<AFIELD>(params, "initial_guess_mode", init_guess_mode);
  if (err1 ){ // not a mandatory parameter
    init_guess_mode = InitialGuess::RHS;
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond, init_guess_mode);
  set_parameters_GMRES_m(N_M);
}


/*
int ASolver_GMRES_m_Cmplx::fetch_initial_guess(const Parameters &param,
                                              const std::string& key,
                                              InitialGuess &init_guess_mode){
  if (!params.find_string(key)) {
    return EXIT_FAILURE;
  }

  const string initial_guess_mode = params.get_string(key);
  vout.detailed(m_vl, "  initila_guess_mode %s\n", initial_guess_mode.c_str());
  if (initial_guess_mode == "RHS") {
    init_guess_mode = InitialGuess::RHS;
  } else if (initial_guess_mode == "GIVEN") {
    init_guess_mode = InitialGuess::GIVEN;
  } else if (initial_guess_mode == "ZERO") {
    init_guess_mode = InitialGuess::ZERO;
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported initial guess mode, %s\n", class_name.c_str(), initial_guess_mode.c_str());
    return EXIT_FAILURE;
  }(

  return EXIT_SUCCESS;

}
*/

//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_int("maximum_number_of_restart", m_Nrestart);
  params.set_double("convergence_criterion_squared", m_Stop_cond);
  params.set_int("number_of_orthonormal_vectors", m_N_M);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  std::string init_guess_str = InitialGuess_str<AFIELD>(m_initial_mode);
  params.set_string("initial_guess_mode", init_guess_str);

}

//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::set_parameters(const int Niter, const int Nrestart, const real_t Stop_cond)
{
  ThreadManager::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Nrestart  = %d\n", Nrestart);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_negative(Nrestart);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter     = Niter;
  m_Nrestart  = Nrestart;
  m_Stop_cond = Stop_cond;
}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::set_parameters(const int Niter, const int Nrestart, const real_t Stop_cond, const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  Niter          = %d\n", Niter);
  vout.general(m_vl, "  Nrestart       = %d\n", Nrestart);
  vout.general(m_vl, "  Stop_cond      = %8.2e\n", Stop_cond);
  //  vout.general(m_vl, "  init_guess_mode = %d (0:RHS, 1:GIVEN, 2:ZERO)\n", init_guess_mode);
  vout.general(m_vl, "  init_guess_mode = %s\n", InitialGuess_str<AFIELD>(init_guess_mode).c_str());
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_negative(Nrestart);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter          = Niter;
  m_Nrestart       = Nrestart;
  m_Stop_cond      = Stop_cond;
  m_initial_mode   = init_guess_mode;
}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::set_parameters_GMRES_m(const int N_M)
{
  ThreadManager::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "  N_M   = %d\n", N_M);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(N_M);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_N_M = N_M;

  m_v.resize(m_N_M + 1);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();
  for (int i = 0; i < m_N_M + 1; ++i) {
    m_v[i].reset(nin, nvol, nex);
  }

}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::set_parameters(const int Niter,
                                          const int Nrestart,
                                          const real_t Stop_cond,
                                          const InitialGuess init_guess_mode,
                                          const int N_M)
{
  set_parameters(Niter, Nrestart, Stop_cond, init_guess_mode);
  set_parameters_GMRES_m(N_M);
}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                                          int& Nconv, real_t& diff)
{
  const real_t bnorm2 = b.norm2();
  const int     bsize = b.size();

  vout.paranoiac(m_vl, "%s: starts\n", class_name.c_str());
  vout.paranoiac(m_vl, "  norm of b = %16.8e\n", bnorm2);
  vout.paranoiac(m_vl, "  size of b = %d\n", bsize);

  bool   is_converged = false;
  int    Nconv2       = 0;
  real_t diff2        = 1.0;  // superficial initialization
  real_t rr;
  coeff_t coeff;

  int Nconv_unit = 1;
  // if (m_fopr->get_mode() == "DdagD" || m_fopr->get_mode() == "DDdag") {
  //   Nconv_unit = 2;
  // }

  solve_init(b, xq, rr, coeff, m_initial_mode);
  Nconv2 += Nconv_unit;

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);


  for (int i_restart = 0; i_restart < m_Nrestart; i_restart++) {
    for (int iter = 0; iter < m_Niter; iter++) {
      if (rr / bnorm2 < m_Stop_cond) break;

      solve_step(b, rr, coeff);
      Nconv2 += Nconv_unit * m_N_M;

      vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);
    }

    //- calculate true residual
    m_fopr->mult(m_s, m_x);  // s  = m_fopr->mult(x);
    axpy(m_s, -1.0, b);      // s -= b;
    diff2 = m_s.norm2();

    if (diff2 < m_Stop_cond * bnorm2) {
      vout.detailed(m_vl, "%s: converged.\n", class_name.c_str());
      vout.detailed(m_vl, "  iter(final): %8d  %22.15e\n", Nconv2, diff2 / bnorm2);

      is_converged = true;

      m_Nrestart_count = i_restart;
      m_Nconv_count    = Nconv2;

      break;
    } else {
      //- restart with new approximate solution
      solve_init(b, m_x, rr, coeff, InitialGuess::GIVEN);

      vout.detailed(m_vl, "%s: restarted.\n", class_name.c_str());
    }
  }


  if (!is_converged) {
    vout.crucial(m_vl, "Error at %s: not converged.\n", class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n", Nconv2, diff2 / bnorm2);
    //exit(EXIT_FAILURE);
  }


  copy(xq, m_x);  // xq = x;

#pragma omp barrier
#pragma omp master
  {
    diff  = sqrt(diff2 / bnorm2);
    Nconv = Nconv2;
  }
#pragma omp barrier
}



//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::solve_init(const AFIELD& b, const AFIELD& xq,
                                              real_t& rr, coeff_t &coeff,
                                              const InitialGuess initial_guess)
{

  coeff.h.resize((m_N_M + 1) * m_N_M);
  coeff.y.resize(m_N_M);
  coeff.g.resize(m_N_M+1);

  if(initial_guess == InitialGuess::GIVEN){
    copy(m_s, xq);  // s = xq;
  } else if (initial_guess == InitialGuess::ZERO){
    m_s.set(0.0);
  } else { // default: RHS
    copy(m_s, b);  // s = xq;
  }

  copy(m_x, m_s);  // x  = s;

  for (int i = 0; i < m_N_M + 1; ++i) { // redundant
    m_v[i].set(0.0);           // m_v[i] = 0.0;
  }

  // r = b - A x_0
  if (initial_guess == InitialGuess::ZERO){
    copy(m_r, b);                // r  = b;
  } else {
    m_fopr->mult(m_v_tmp, m_s);  // v_tmp = m_fopr->mult(s);
    copy(m_r, b);                // r  = b;
    axpy(m_r, -1.0, m_v_tmp);    // r -= v_tmp;
  }

  rr = m_r.norm2();            // rr = r * r;

  coeff.beta_prev=sqrt(rr);

  //- v[0] = (1.0 / beta_prev) * r;
  copy(m_v[0], m_r);                      // v[0] = r;
  scal(m_v[0], (1.0 / coeff.beta_prev));  // v[0] = (1.0 / beta_p) * v[0];
}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::solve_step(const AFIELD& b, real_t& rr,
                                              coeff_t &coeff)
{
  std::vector<dcomplex> &h = coeff.h;
  std::vector<dcomplex> &y = coeff.y;
  for(int ij=0; ij<h.size(); ++ij) { h[ij] = dcomplex(0.0, 0.0); }
  for(int i=0; i<y.size(); ++i) { y[i] = dcomplex(0.0, 0.0); }

  for (int j = 0; j < m_N_M; ++j) {
    m_fopr->mult(m_v_tmp, m_v[j]);  // v_tmp = m_fopr->mult(v[j]);

    for (int i = 0; i < j + 1; ++i) {
      int ij = index_ij(i, j);
      h[ij] = dotc(m_v[i], m_v_tmp);  // h[ij] = (v[i], A v[j]);
    }

    //- v[j+1] = A v[j] - \Sum_{i=0}^{j-1} h[i,j] * v[i]
    copy(m_v[j + 1],  m_v_tmp);

    for (int i = 0; i < j + 1; ++i) {
      int ij = index_ij(i, j);
      complex_t a = complex_t( real(h[ij]), imag(h[ij]));
      axpy(m_v[j + 1], -a, m_v[i]);  // v[j+1] -= h[ij] * v[i];
    }

    double v_norm2 = m_v[j + 1].norm2();

    int j1j = index_ij(j + 1, j);
    h[j1j] = dcomplex(sqrt(v_norm2), 0.0);

    scal(m_v[j + 1], real_t(1.0 / sqrt(v_norm2)));  // v[j+1] /= sqrt(v_norm2);
  }


  // Compute y, which minimizes J := |r_new| = |beta_p - h * y|
  min_J(coeff);


  // x += Sum_{i=0}^{N_M-1} y[i] * v[i];
  for (int i = 0; i < m_N_M; ++i) {
    complex_t a = complex_t(real(y[i]), imag(y[i]));
    axpy(m_x, a, m_v[i]);  // x += y[i] * v[i];
  }


  // r = b - m_fopr->mult(x);
  //  copy(m_s, m_x);  // s = x;
  //  solve_init(b, rr);

  copy(m_s, m_x);  // s  = x;

  //  for (int i = 0; i < m_N_M + 1; ++i) { // redundant
  //    m_v[i].set(0.0);           // m_v[i] = 0.0;
  //  }

  // r = b - A x_0
  m_fopr->mult(m_v_tmp, m_s);  // v_tmp = m_fopr->mult(s);
  copy(m_r, b);                // r  = b;
  axpy(m_r, -1.0, m_v_tmp);    // r -= v_tmp;

  rr = m_r.norm2();            // rr = r * r;

  coeff.beta_prev=sqrt(rr);

  //- v[0] = (1.0 / beta_prev) * r;
  copy(m_v[0], m_r);                      // v[0] = r;
  scal(m_v[0], (1.0 / coeff.beta_prev));  // v[0] = (1.0 / beta_p) * v[0];

}


//====================================================================
template<typename AFIELD>
void ASolver_GMRES_m_Cmplx<AFIELD>::min_J(coeff_t& coeff)
{
  // Compute y, which minimizes J := |r_new| = |beta_p - h * y|

  std::vector<dcomplex> &h = coeff.h;
  std::vector<dcomplex> &y = coeff.y;
  std::vector<dcomplex> &g = coeff.g;

  for(int i=0; i<g.size(); ++i) { g[i] = dcomplex(0.0, 0.0); }
  g[0] = dcomplex(coeff.beta_prev, 0.0);

  for (int i = 0; i < m_N_M; ++i) {
    int    ii    = index_ij(i, i);
    double h_1_r2 = real(h[ii])*real(h[ii]) + imag(h[ii])*imag(h[ii]);

    int    i1i   = index_ij(i + 1, i);
    double h_2_r2 = real(h[i1i])*real(h[i1i]) + imag(h[i1i])*imag(h[i1i]);

    double denomi = sqrt(h_1_r2 + h_2_r2);

    dcomplex cs = h[ii] / denomi;
    dcomplex sn = h[i1i] / denomi;

    for (int j = i; j < m_N_M; ++j) {
      int ij  = index_ij(i, j);
      int i1j = index_ij(i + 1, j);

      dcomplex const_1_c = conj(cs) * h[ij] + sn * h[i1j];
      dcomplex const_2_c = -sn * h[ij] + cs * h[i1j];

      h[ij]  = const_1_c;
      h[i1j] = const_2_c;
    }

    dcomplex const_1_c = conj(cs) * g[i] + sn * g[i + 1];
    dcomplex const_2_c = -sn * g[i] + cs * g[i + 1];

    g[i]     = const_1_c;
    g[i + 1] = const_2_c;
  }


  for (int i = m_N_M - 1; i > -1; --i) {
    for (int j = i + 1; j < m_N_M; ++j) {
      int ij = index_ij(i, j);
      g[i] -= h[ij] * y[j];
    }

    int ii = index_ij(i, i);
    y[i] = g[i] / h[ii];
  }
}


//====================================================================
template<typename AFIELD>
double ASolver_GMRES_m_Cmplx<AFIELD>::flop_count()
{
  const int NPE = CommonParameters::NPE();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol = CommonParameters::Nvol()/2 for eo
  const int Nin  = m_x.nin();
  const int Nvol = m_x.nvol();
  const int Nex  = m_x.nex();

  const double gflop_fopr = m_fopr->flop_count();

  if (gflop_fopr < CommonParameters::epsilon_criterion()) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0\n", class_name.c_str());
    return 0.0;
  }

  const double gflop_axpy = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_dotc = (Nin * Nex * 4) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_norm = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_scal = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);

  int N_M_part = 0;
  for (int j = 0; j < m_N_M; ++j) {
    for (int i = 0; i < j + 1; ++i) {
      N_M_part += 1;
    }
  }

  const double gflop_init = gflop_fopr + gflop_axpy + gflop_norm + gflop_scal;
  const double gflop_step = m_N_M * gflop_fopr + N_M_part * gflop_dotc
                            + (N_M_part + m_N_M) * gflop_axpy
                            + m_N_M * gflop_scal
                            + gflop_init;
  const double gflop_true_residual = gflop_fopr + gflop_axpy + gflop_norm;

  const int    N_iter = (m_Nconv_count - 1) / m_N_M;
  const double gflop  = gflop_norm + gflop_init
                        + gflop_step * N_iter
                        + gflop_true_residual * (m_Nrestart_count + 1)
                        + gflop_init * m_Nrestart_count;

  return gflop;
}


//====================================================================
//============================================================END=====
