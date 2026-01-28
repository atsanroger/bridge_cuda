template<typename AFIELD>
const std::string ASolver_FGMRES<AFIELD>::class_name = "ASolver_FGMRES";

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
void ASolver_FGMRES<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_initial_mode = InitialGuess::ZERO;

  m_r.reset(nin, nvol, nex);
  m_x.reset(nin, nvol, nex);
  m_w.reset(nin, nvol, nex);

}

//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::tidyup(void)
{
  // nothing to do
}

//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::set_parameters(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Niter;
  double Stop_cond;  // fetching as double
  int    N_M;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  //  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_int("number_of_orthonormal_vectors", N_M);

  InitialGuess init_guess_mode;
  int err1 = fetch_initial_guess<AFIELD>(params, "initial_guess_mode", init_guess_mode);
  if (err1 ){ // not a mandatory parameter
    init_guess_mode = InitialGuess::ZERO;
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Stop_cond, init_guess_mode);
  set_parameters_GMRES_m(N_M);
}


/*
int ASolver_FGMRES::fetch_initial_guess(const Parameters &param,
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
void ASolver_FGMRES<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_int("maximum_number_of_iteration", m_Niter);
  //  params.set_int("maximum_number_of_restart", m_Nrestart);
  params.set_double("convergence_criterion_squared", m_Stop_cond);
  params.set_int("number_of_orthonormal_vectors", m_N_M);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  std::string init_guess_str = InitialGuess_str<AFIELD>(m_initial_mode);
  params.set_string("initial_guess_mode", init_guess_str);

}

//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::set_parameters(const int Niter, const real_t Stop_cond)
{
  //  set_parameters(Niter, Stop_cond, InitialGuess::ZERO);
  set_parameters(Niter, Stop_cond, m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter          = %d\n", Niter);
  //  vout.general(m_vl, "  Nrestart       = %d\n", Nrestart);
  vout.general(m_vl, "  Stop_cond      = %8.2e\n", Stop_cond);
  //  vout.general(m_vl, "  init_guess_mode = %d (0:RHS, 1:GIVEN, 2:ZERO)\n", init_guess_mode);
  vout.general(m_vl, "  init_guess_mode = %s\n", InitialGuess_str<AFIELD>(init_guess_mode).c_str());
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  //  err += ParameterCheck::non_negative(Nrestart);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter          = Niter;
  //  m_Nrestart       = Nrestart;
  m_Stop_cond      = Stop_cond;
  m_initial_mode   = init_guess_mode;
}


//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::set_parameters_GMRES_m(const int N_M)
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

  m_z.resize(m_N_M + 1);
  m_v.resize(m_N_M + 1);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();
  for (int i = 0; i < m_v.size(); ++i) {
    m_v[i].reset(nin, nvol, nex);
  }
  for (int i = 0; i < m_z.size(); ++i) {
    m_z[i].reset(nin, nvol, nex);
  }

}


//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::set_parameters(const int Niter,
                                            const real_t Stop_cond,
                                            const InitialGuess init_guess_mode,
                                          const int N_M)
{
  set_parameters(Niter,Stop_cond, init_guess_mode);
  set_parameters_GMRES_m(N_M);
}


//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::init_coeff(coeff_t &coeff)
{
  assert(m_N_M>0);

  int m=m_N_M;

  coeff.ht.resize(m);  // H(i+1,i)
  coeff.h.resize(m*m); // H(i,j)
  coeff.r.resize(m*m);
  coeff.y.resize(m);
  coeff.c.resize(m);
  coeff.s.resize(m);
  coeff.b.resize(m+1);

}
//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
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

  //  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);


  int i = 0;
  int niter = -1;
  for (int iter = 0; iter < m_Niter; iter++) {

    if( i == 0 ){
      //      double r2 = m_r.norm2();
      double r=sqrt(rr);
      coeff.b[0] = dcomplex(r, 0.0);
      copy(m_v[0], m_r);
      scal(m_v[0], 1.0/r);
    }

    solve_step1(rr, coeff, i);

    vout.detailed(m_vl, "FGMRES    iter: %8d  %22.15e\n", iter, rr / bnorm2);

    if (rr < m_Stop_cond * bnorm2) { // congerged, maybe?
      vout.detailed(m_vl, "  seems converged, constructing the solution\n");
      solve_step2(rr, coeff, i, b);  // construct solution
      if (rr  < m_Stop_cond * bnorm2) {
        niter = iter;
        break;
      }
      i = 0;
      continue;
    }

    ++i;
    if( i == m_N_M){
      int nr = m_N_M - 1;
      solve_step2(rr, coeff, nr, b);
      if (rr < m_Stop_cond * bnorm2) {
        niter = iter;
        break;
      }
      vout.detailed(m_vl, "  i = %d, rr=%23.15e, restarting\n", i, rr);
      i = 0;
    }
  } // iter


  Nconv2 = niter;
  diff2 = rr;
  if (niter < 0) {
    vout.crucial(m_vl, "Error at %s: not converged.\n", class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n", m_Niter, diff2 / bnorm2);
    //exit(EXIT_FAILURE);
  } else {
    vout.detailed(m_vl, "FGMRES converged: %d  %22.15e\n", Nconv2, diff2/bnorm2);
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
void ASolver_FGMRES<AFIELD>::solve_init(const AFIELD& b, const AFIELD& xq,
                                        real_t& rr, coeff_t &coeff,
                                        const InitialGuess initial_guess)
{

  init_coeff(coeff);


  if(initial_guess == InitialGuess::GIVEN){
    copy(m_x, xq);  // s = xq;
  } else if (initial_guess == InitialGuess::ZERO){
    m_x.set(0.0);
  } else { // default: RHS
    copy(m_x, b);  // s = xq;
  }


  // r = b - A x_0
  if (initial_guess == InitialGuess::ZERO){
    copy(m_r, b);                // r  = b;
  } else {
    m_fopr->mult(m_w, m_x);  // v_tmp = m_fopr->mult(s);
    copy(m_r, b);                // r  = b;
    axpy(m_r, -1.0, m_w);    // r -= v_tmp;
  }

  rr = m_r.norm2();          // rr = r * r;

}


//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::solve_step1(real_t& rr, coeff_t &coeff,
                                         const int i)
{

  // apply preconditinor
  m_prec->mult(m_z[i], m_v[i]);

  // mult
  m_fopr->mult(m_v[i+1], m_z[i]);
#ifdef DEBUG
  {
    copy(m_r, m_v[i+1]);
    axpy(m_r, -1.0, m_v[i]);
    double b2 = m_v[i].norm2();
    double r2 = m_r.norm2();
    vout.general("after the prec: i=%d, b2=%e, r2=%e (realtvie: %e)\n",
                 i, b2, r2, r2/b2);
  }
#endif

  // G-S
  for(int k=0; k<i+1; ++k){
    dcomplex h_ki = dotc(m_v[k], m_v[i+1]);
    coeff.h[index_ij(k,i)] = h_ki;
    axpy(m_v[i+1], -complex_t( real(h_ki), imag(h_ki)), m_v[k]);

    dcomplex dh_ki = dotc(m_v[k], m_v[i+1]);
    coeff.h[index_ij(k,i)] += dh_ki;
    axpy(m_v[i+1], -complex_t( real(dh_ki), imag(dh_ki)), m_v[k]);
  }
  double w2 = m_v[i+1].norm2();
  coeff.ht[i] = sqrt(w2);
  m_v[i+1].scal(1.0/coeff.ht[i]);

  coeff.r[index_ij(0,i)] = coeff.h[index_ij(0,i)];

  for(int k=1; k<i+1; ++k){
    double   c = coeff.c[k-1];
    dcomplex s = coeff.s[k-1];
    dcomplex r = coeff.r[index_ij(k-1,i)];
    dcomplex h = coeff.h[index_ij(k,i)];

    coeff.r[index_ij(k-1,i)] =  c*r + std::conj(s)*h;
    coeff.r[index_ij(k,i)]   = -s*r + c*h;
  }
  double rii2 = norm( coeff.r[index_ij(i,i)] );
  double rii  = sqrt(rii2);
  double hti  = coeff.ht[i]; // h(i+1, i)
  double delta_inv = 1.0/sqrt( rii2 + hti*hti );
  dcomplex mu;
  if( rii < hti){
    mu = std::conj(rii/hti);
  } else {
    mu = hti/rii;
  }
  double   c = rii * delta_inv;
  dcomplex s = hti * delta_inv * mu/abs(mu);
  coeff.c[i] = c;
  coeff.s[i] = s;
  coeff.r[index_ij(i,i)] = c*coeff.r[index_ij(i,i)] + std::conj(s)*hti;
  coeff.b[i+1] = - s * coeff.b[i];
  coeff.b[i]   = c*coeff.b[i];

  rr = norm(coeff.b[i+1]);

}

//====================================================================
template<typename AFIELD>
void ASolver_FGMRES<AFIELD>::solve_step2(real_t& rr, coeff_t &coeff,
                                         const int nr,
                                         const AFIELD &b)
{

  vout.detailed(m_vl, "%s: constructing solution, nr = %d\n",
                class_name.c_str(), nr);

  coeff.y[nr] = coeff.b[nr]/coeff.r[index_ij(nr,nr)];
  for(int k = nr-1; k > -1; --k){
    dcomplex tmp = coeff.b[k];
    for(int i=k+1; i<nr+1; ++i){
      tmp -= coeff.r[index_ij(k,i)] * coeff.y[i];
    }
    coeff.y[k] = tmp/coeff.r[index_ij(k,k)];
  }

  for(int i=0; i<nr+1; ++i){
    dcomplex yi =coeff.y[i];
    axpy(m_x, complex_t(real(yi), imag(yi)), m_z[i]);
  }

  m_fopr->mult(m_w, m_x);
  copy(m_r, b);
  axpy(m_r, -1.0, m_w);

  rr = m_r.norm2();

}

//====================================================================
template<typename AFIELD>
double ASolver_FGMRES<AFIELD>::flop_count()
{

  vout.crucial(m_vl, "%s: flop_count is not yet ready\n",
               class_name.c_str());
  return 0.0;


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

  //  const int    N_iter = (m_Nconv_count - 1) / m_N_M;
  //  const double gflop  = gflop_norm + gflop_init
  //                        + gflop_step * N_iter
  //                        + gflop_true_residual * (m_Nrestart_count + 1)
  //                        + gflop_init * m_Nrestart_count;
  const double gflop=0;

  return gflop;
}


//====================================================================
//============================================================END=====
