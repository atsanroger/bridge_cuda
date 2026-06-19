/*!
      @file    afopr_Domainwall_General_5din_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Field/afield-inc.h"

template<typename AFIELD>
const std::string AFopr_Domainwall_General_5din_eo<AFIELD>::class_name
                                     = "AFopr_Domainwall_General_5din_eo";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  // switches
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  vout.detailed(m_vl, "%s: initalization starts.\n",
                class_name.c_str());

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  m_Nvcd = 2 * Nc * Nd;
  m_Ndf  = 2 * Nc * Nc;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nx = CommonParameters::Nx();
  m_Ny = CommonParameters::Ny();
  m_Nz = CommonParameters::Nz();
  m_Nt = CommonParameters::Nt();

  check_Nc(Nc);
  check_setup();

  // condition check
  if(m_Nx % 2 != 0 || m_Ny % 2 != 0){
    vout.crucial(m_vl, "%s: Nx and Ny must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if(m_Nx % (2*VLEN2) != 0){
    vout.crucial(m_vl, "%s: Nx must be mulriple of 2*VLEN2.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nx2 = m_Nx/2;
  m_Nst2 = m_Nvol/2;

  m_Nx2v  = m_Nx2/VLEN2;
  m_Nst2v = m_Nst2/VLEN2;

  vout.general(m_vl, "  Nx2v = %d  Nst2v = %d\n", m_Nx2v, m_Nst2v);

  m_Leo.resize(m_Ny * m_Nz * m_Nt);

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  for (int t = 0; t < m_Nt; ++t) {
    for (int z = 0; z < m_Nz; ++z) {
      for (int y = 0; y < m_Ny; ++y) {
        int t2 = ipe3 * m_Nt + t;
        int z2 = ipe2 * m_Nz + z;
	int y2 = ipe1 * m_Ny + y;
        m_Leo[y + m_Ny * (z + m_Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

  // setup verbose level
  string vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(vlevel);

  vout.general(m_vl, "Initialization of %s:\n", class_name.c_str());
  int err = 0;

  // setup kernel operator
  std::string kernel_type;
  err += params.fetch_string("kernel_type", kernel_type);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
		 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel = params;
  double M0;
  err += params.fetch_double("domain_wall_height", M0);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
		 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  m_M0 = real_t(M0);

  double kappa = 1.0/(8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  // Factory is assumed to work.
  m_foprw = AFopr<AFIELD>::New(kernel_type, params_kernel);
  m_foprw->set_mode("D");

  m_Ns = 0; // temporary set

  set_parameters(params);

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  int Nbdin = (m_Nvcd/2) * m_Ns;
  m_Nbdsize[0] = Nbdin * ((m_Ny * m_Nz * m_Nt + 1)/2);
  m_Nbdsize[1] = Nbdin * m_Nx2 * m_Nz * m_Nt;
  m_Nbdsize[2] = Nbdin * m_Nx2 * m_Ny * m_Nt;
  m_Nbdsize[3] = Nbdin * m_Nx2 * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nvol, m_Ndim);

  vout.detailed(m_vl, "%s: initalization finished.\n",
                class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::tidyup()
{
  delete m_foprw;
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu] = Communicator::send_init(Nvsize, mu, -1);
    chsend_up[mu] = Communicator::send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1);
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu]->ptr();
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu]->ptr();
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if(do_comm[mu] == 1){
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_parameters(
                                             const Parameters& params)
{

  string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string gmset_type;
  double mq, M0;
  int Ns;
  std::vector<int> bc;
  double b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c);

  if(real_t(M0) != m_M0) set_kernel_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_parameters(
                                const double mq,
                                const double M0,
                                const int Ns,
                                const vector<int> bc,
                                const double b,
                                const double c)
{
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());

  m_M0 = real_t(M0);
  m_mq = real_t(mq);
  m_Ns = Ns;

  m_NinF = m_Nvcd * m_Ns;

  assert(bc.size() == m_Ndim);
  if(m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
  }

  if(m_b.size() != m_Ns){
    m_b.resize(m_Ns);
    m_c.resize(m_Ns);
  }
  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = real_t(b);
    m_c[is] = real_t(c);
  }

  vout.general(m_vl, "  mq   = %8.4f\n",m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n",m_M0);
  vout.general(m_vl, "  Ns   = %4d\n",m_Ns);
  for(int mu = 0; mu < m_Ndim; ++mu){
    vout.general(m_vl, "  boundary[%d] = %2d\n",mu,m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
                                                  m_b[0], m_c[0]);

  set_precond_parameters();

  // working 5d vectors.
  if(m_y1.nin() != m_NinF){
    m_y1.reset(m_NinF, m_Nst2, 1);
    m_v1.reset(m_NinF, m_Nst2, 1);
    m_v2.reset(m_NinF, m_Nst2, 1);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_kernel_parameters(
                                             const Parameters& params)
{
  Parameters params_kernel = params;

  double M0;
  params.fetch_double("domain_wall_height", M0);

  double kappa = 1.0/(8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  m_foprw->set_parameters(params_kernel);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_precond_parameters()
{

  if(m_dp.size() != m_Ns){
    m_dp.resize(m_Ns);
    m_dm.resize(m_Ns);
    m_e.resize(m_Ns-1);
    m_f.resize(m_Ns-1);
  }

  for(int is=0; is < m_Ns; ++is){
    m_dp[is] = 1.0 + m_b[is]*(4.0 - m_M0);
    m_dm[is] = 1.0 - m_c[is]*(4.0 - m_M0);
  }

  m_e[0] = m_mq * m_dm[m_Ns-1]/m_dp[0];
  m_f[0] = m_mq * m_dm[0];
  for(int is=1; is < m_Ns-1; ++is){
    m_e[is] = m_e[is-1]*m_dm[is-1]/m_dp[is];
    m_f[is] = m_f[is-1]*m_dm[is]/m_dp[is-1];
  }

  m_g = m_e[m_Ns-2]*m_dm[m_Ns-2];

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_coefficients(
                                  const std::vector<double> vec_b,
                                  const std::vector<double> vec_c)
{
  if(vec_b.size() != m_Ns || vec_c.size() != m_Ns){
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
                 class_name.c_str());

  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = real_t(vec_b[is]);
    m_c[is] = real_t(vec_c[is]);
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                       is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_config(Field* u)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  AFIELD Ulex(m_Ndf, m_Nvol, m_Ndim);

#pragma omp parallel
 {
  Index_lex_alt<real_t,SIMD2> index_lex;
  convert_gauge(index_lex, Ulex, *u);

#pragma omp barrier

  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1) set_boundary_config(Ulex, mu);
  }

#pragma omp barrier

  Index_eo_alt<real_t,SIMD2> index_eo;

  //convert(index_eo, m_Ueo, index_lex, Ulex);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nvol);

  for(int ex = 0; ex < m_Ndim; ++ex){
    for(int site = is; site < ns; ++site){
      for(int in = 0; in < m_Ndf; ++in){
        int iv1 = index_lex.idx( in, m_Ndf, site, ex);
        int iv2 = index_eo.idx(in, m_Ndf, site, ex);
        m_Ueo.set(iv2, Ulex.cmp(iv1));
      }
    }
  }

 } // #pragma omp parallel

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::set_boundary_config(
                                                 AFIELD& U,
						 const int mu)
{
  int ipe[4], Nsize[4], Lsize[4];
  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

  real_t bc = real_t(m_boundary[mu]);

  for(int i = 0; i < m_Ndim; ++i){
    ipe[i]   = Communicator::ipe(i);
    Nsize[i] = CommonParameters::Nsize(i);
    Lsize[i] = CommonParameters::Lsize(i);
  }

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nvol);

  int j[4];
  for(int site = is; site < ns; ++site){
    int x = site % m_Nx;
    int yzt = site/m_Nx;
    int y = yzt % m_Ny;
    int z = (yzt/m_Ny) % m_Nz;
    int t = yzt/(m_Ny*m_Nz);
    j[0] = ipe[0] * Nsize[0] + x;
    j[1] = ipe[1] * Nsize[1] + y;
    j[2] = ipe[2] * Nsize[2] + z;
    j[3] = ipe[3] * Nsize[3] + t;

    if(j[mu] == Lsize[mu]-1){
      for(int in = 0; in < m_Ndf; ++in){
        int i = index_lex.idx_G(in, site, mu);
        real_t uv = bc * U.cmp(i);
        U.set(i, uv);
      }
    }

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::convert(AFIELD& v,
                                                  const Field& w)
{
  //  ThreadManager_OpenMP::assert_single_thread(class_name);

  vout.paranoiac(m_vl, " convert start.\n");

  //#pragma omp parallel
 {
  Index_lex_alt<real_t,AFIELD::IMPL> index;

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nvol);

  for(int site = site0; site < site1; ++site){
    for(int is = 0; is < m_Ns; ++is){
      for(int ivcd = 0; ivcd < m_Nvcd; ++ivcd){
        int in = ivcd + m_Nvcd * is;
        v.set(index.idx(in, m_NinF, site, 0),
              real_t(w.cmp(ivcd, site, is)));
      }
    }
  }

 }

  vout.paranoiac(m_vl, " convert finished.\n");

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::reverse(Field& v,
                                                  const AFIELD& w)
{
  //  ThreadManager_OpenMP::assert_single_thread(class_name);

  vout.paranoiac(m_vl, " reverse start.\n");

  //#pragma omp parallel
 {
  Index_lex_alt<real_t,AFIELD::IMPL> index;

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nvol);

  for(int site = site0; site < site1; ++site){
    for(int is = 0; is < m_Ns; ++is){
      for(int ivcd = 0; ivcd < m_Nvcd; ++ivcd){
        int in = ivcd + m_Nvcd * is;
        v.set(ivcd, site, is,
              double(w.cmp(index.idx(in, m_NinF, site, 0))) );
      }
    }
  }

 }

  vout.paranoiac(m_vl, " reverse finished.\n");

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::mult(AFIELD& v,
                                               const AFIELD& w)
{
  if(m_mode=="D"){
    D(v,w);
  }else if(m_mode=="Ddag"){
    Ddag(v, w);
  }else if(m_mode=="DdagD"){
    DdagD(v, w);
  }else{
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::mult_dag(AFIELD& v,
                                                   const AFIELD& w)
{
  if(m_mode=="D"){
    Ddag(v, w);
  }else if(m_mode=="Ddag"){
    D(v, w);
  }else if(m_mode=="DdagD"){
    DdagD(v, w);
  }else{
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::mult(AFIELD& v,
                                                    const AFIELD& w,
                                                    std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

#pragma omp barrier

  if(mode=="Deo"){
    D_eo(v, w, 0);
  }else if(mode=="Doe"){
    D_eo(v, w, 1);
  }else if(mode=="Dee"){
    D_ee(v, w, 0);
  }else if(mode=="Doo"){
    D_ee(v, w, 1);
  }else if(mode=="Dee_inv"){
    L_inv(m_v1, w);
    U_inv(v, m_v1);
  }else if(mode=="Doo_inv"){
#pragma omp barrier
    L_inv(m_v1, w);
#pragma omp barrier
    U_inv(v, m_v1);
#pragma omp barrier
  }else{
    vout.crucial(m_vl, "Error at %s: mode undefined in mult_dag.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::mult_dag(AFIELD& v,
                                                     const AFIELD& w,
                                                     std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

#pragma omp barrier

  if(mode=="Deo"){
    Ddag_eo(v, w, 0);
  }else if(mode=="Doe"){
    Ddag_eo(v, w, 1);
  }else if(mode=="Dee"){
    Ddag_ee(v, w, 0);
  }else if(mode=="Doo"){
    Ddag_ee(v, w, 1);
  }else if(mode=="Dee_inv"){
    Udag_inv(m_v1, w);
    Ldag_inv(v, m_v1);
  }else if(mode=="Doo_inv"){
    Udag_inv(m_v1, w);
    Ldag_inv(v, m_v1);
  }else{
    vout.crucial(m_vl, "Error at %s: mode undefined in mult_dag.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::DdagD(AFIELD& v,
                                                  const AFIELD& w)
{
#pragma omp barrier

  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#pragma omp barrier

  Udag_inv(v, m_v1);
  Ldag_inv(m_v2, v);
  Ddag_eo(v, m_v2, 1);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
  Ddag_eo(m_v2, v, 0);

  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  Udag_inv(m_v1, w);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::D_eo(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even

  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], wL[NCD], yL[NCD];
    clear_vec(vL, NCD);

    for(int is = 0; is < m_Ns; ++is){

      load_vec(wL, &wp2[Nin4*is], NCD);
      set_vec(yL, m_b[is], wL, NCD);

      int is_up = (is + 1) % m_Ns;
      real_t Fup = 0.5 * m_c[is];
      if(is == m_Ns-1) Fup *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_up], NCD);
      add_aPm5_dirac_vec(yL, Fup, wL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn = 0.5 * m_c[is];
      if(is == 0) Fdn *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_dn], NCD);
      add_aPp5_dirac_vec(yL, Fdn, wL, NC);

      scal_vec(yL, real_t(-0.5), NCD);
      save_vec(&yp2[Nin4*is], yL, NCD);

      save_vec(&vp2[Nin4*is], vL, NCD);

    }

  }
#pragma omp barrier

  Dhop(vp, yp, ieo);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Ddag_eo(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even

  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  // v <- w  here gamma_5 is multiplied.
  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], yL[NCD];

    clear_vec(yL, NCD);
    for(int is = 0; is < m_Ns; ++is){
      load_mult_gm5_dirac_vec(vL, &wp2[Nin4*is], NC);
      save_vec(&vp2[Nin4*is], vL, NCD);
      save_vec(&yp2[Nin4*is], yL, NCD);
    }

  }

#pragma omp barrier

  Dhop(yp, vp, ieo);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], yL[NCD];

    for(int is = 0; is < m_Ns; ++is){

      load_mult_gm5_dirac_vec(yL, &yp2[Nin4*is], NCD);
      set_vec(vL, real_t(-0.5)*m_b[is], yL, NCD);

      int is_up = (is + 1) % m_Ns;
      real_t Fup = 0.5 * (-0.5) * m_c[is_up];
      if(is == m_Ns-1) Fup *= - m_mq;
      load_mult_gm5_dirac_vec(yL, &yp2[Nin4*is_up], NCD);
      add_aPp5_dirac_vec(vL, Fup, yL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn = 0.5 * (-0.5) * m_c[is_dn];
      if(is == 0) Fdn *= - m_mq;
      load_mult_gm5_dirac_vec(yL, &yp2[Nin4*is_dn], NCD);
      add_aPm5_dirac_vec(vL, Fdn, yL, NC);

      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::D_ee(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even (no difference)

  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NCD], wL[NCD];

    for(int is = 0; is < m_Ns; ++is){

      real_t FF1 = m_b[is] *(4.0 - m_M0) + 1.0;
      real_t FF2 = m_c[is] *(4.0 - m_M0) - 1.0;

      load_vec(wL, &wp2[Nin4*is], NCD);
      set_vec(vL, FF1, wL, NCD);

      int is_up = (is + 1) % m_Ns;
      real_t Fup = 0.5 * FF2;
      if(is == m_Ns-1) Fup *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_up], NCD);
      add_aPm5_dirac_vec(vL, Fup, wL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn = 0.5 * FF2;
      if(is == 0) Fdn *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_dn], NCD);
      add_aPp5_dirac_vec(vL, Fdn, wL, NC);

      save_vec(&vp2[Nin4*is], vL, NCD);

    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Ddag_ee(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even

  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NCD], wL[NCD];

    for(int is = 0; is < m_Ns; ++is){

      real_t FF1 = m_b[is] *(4.0 - m_M0) + 1.0;
      load_vec(vL, &wp2[Nin4*is], NCD);
      scal_vec(vL, FF1, NCD);

      int is_up = (is + 1) % m_Ns;
      real_t Fup = 0.5 * (m_c[is_up] *(4.0 - m_M0) - 1.0);
      if(is == m_Ns-1) Fup *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_up], NCD);
      add_aPp5_dirac_vec(vL, Fup, wL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn = 0.5 * (m_c[is_dn] *(4.0 - m_M0) - 1.0);
      if(is == 0) Fdn *= - m_mq;
      load_vec(wL, &wp2[Nin4*is_dn], NCD);
      add_aPm5_dirac_vec(vL, Fdn, wL, NC);

      save_vec(&vp2[Nin4*is], vL, NCD);

    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Dhop(real_t* vp,
                                                    real_t* wp,
                                                    const int ieo)
{

  Dhop_1(vp, wp, ieo);

#pragma omp master
 {
  chset_send.start();
  chset_recv.start();
 }

  Dhop_b(vp, wp, ieo);

#pragma omp master
 {
  chset_send.wait();
  chset_recv.wait();
 }

#pragma omp barrier

  Dhop_2(vp, wp, ieo);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Dhop_1(real_t* vp,
                                                      real_t* wp,
                                                      const int ieo)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU2  = m_Ndf * m_Nst2;

  int Nxy2 = m_Nx2v * m_Ny;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  real_t bufL[Nin5];

  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    int idir;

    real_t *wp2 = &wp[Nin5 * site];

    idir = 0; // x-direction boundary
    if(do_comm[idir] == 1){

      if(ix == 0 && Leo == 1){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * (iyzt/2);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_xp1(&buf1[ibf + NVC*ND2*is], &wp2[Nin4*is], NC);
        }
      }

      if(ix == m_Nx2v-1 && Leo == 0){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * (iyzt/2);
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_xm1(&buf1[ibf + NVC*ND2*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 1; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iy == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nx2v * izt);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_yp1(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(iy == m_Ny-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nx2v * izt);
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_ym1(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 2; // z-direction boundary
    if(do_comm[idir] == 1){

      if(iz == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy2 * it);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_zp1(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(iz == m_Nz-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy2 * it);
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_zm1(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 3; // t-direction boundary
    if(do_comm[idir] == 1){

      if(it == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * ixyz;
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_tp1_dirac(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(it == m_Nt-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * ixyz;
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_tm1_dirac(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

  } // site loop (boundary)

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Dhop_b(real_t* vp,
                                                      real_t* wp,
                                                      const int ieo)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int NvU2  = m_Ndf * m_Nst2;

  int Nxy2  = m_Nx2v * m_Ny;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  real_t bufL[Nin5];

  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    int idir;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN*NCD];
    Vsimd_t vL[NCD];

    idir = 0;

    if(Leo == 0){
      real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_xpb(vL, up, &wp2[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }else{
      if(ix < m_Nx2v-1 || do_comm[idir] == 0){
        int ix2 = (ix + 1) % m_Nx2v;
        int nei = ix2 + m_Nx2v * iyzt;
        real_t *wpn = &wp[Nin5 * nei];
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo + 2*idir));
        for(int is = 0; is < m_Ns; ++is){
          shift_vec2_bw(z4, &wp2[Nin4*is], &wpn[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xpb(vL, up, z4, NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }
    }

    if(Leo == 1){
      real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_xmb(vL, up, &wp2[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }else{
      if(ix > 0 || do_comm[idir] == 0){
        int ix2 = (ix - 1 + m_Nx2v) % m_Nx2v;
        int nei = ix2 + m_Nx2v * iyzt;
        real_t *wpn = &wp[Nin5 * nei];
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo + 2*idir));
        real_t *un = m_Ueo.ptr(VLEN*NDF2*nei  + NvU2*(1-ieo + 2*idir));
        real_t uL[VLEN*NDF2];
        shift_vec2_fw(uL, up, un, NDF2);
        for(int is = 0; is < m_Ns; ++is){
          shift_vec2_fw(z4, &wp2[Nin4*is], &wpn[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xmb(vL, uL, z4, NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }
    }

    idir = 1; // y-direction bulk

    if(iy < m_Ny-1 || do_comm[idir] == 0){
      int iy2 = (iy + 1) % m_Ny;
      int nei = ix + m_Nx2v * (iy2 + m_Ny * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_ypb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(iy > 0 || do_comm[idir] == 0){
      int iy2 = (iy - 1 + m_Ny) % m_Ny;
      int nei = ix + m_Nx2v * (iy2 + m_Ny * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*nei + NvU2*(1-ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_ymb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    idir = 2; // z-direction bulk

    if(iz < m_Nz-1 || do_comm[idir] == 0){
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxy2 * (iz2 + m_Nz*it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_zpb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(iz > 0 || do_comm[idir] == 0){
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy2 * (iz2 + m_Nz*it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*nei + NvU2*(1-ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_zmb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    idir = 3; // t-direction bulk

    if(it < m_Nt-1 || do_comm[idir] == 0){
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyz2 * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_tpb_dirac(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(it > 0 || do_comm[idir] == 0){
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz2 * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_Ueo.ptr(VLEN*NDF2*nei + NvU2*(1-ieo + 2*idir));
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_tmb_dirac(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Dhop_2(real_t* vp,
                                                      real_t* wp,
                                                      const int ieo)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU2  = m_Ndf * m_Nst2;

  int Nxy2  = m_Nx2v * m_Ny;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  real_t bufL[Nin5];

  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN*NCD];
    Vsimd_t vL[NCD];

    idir = 0;
    if(do_comm[idir] == 1){

      if(Leo == 1 && ix == m_Nx2v-1){
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo+2*idir));
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * (iyzt/2);
        for(int is = 0; is < m_Ns; ++is){
          shift_vec0_bw(z4, &wp2[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xpb(vL, up, z4, NC);
          mult_wilson_xp2(vL, up, &buf2[ibf + NVC*ND2*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(Leo == 0 && ix == 0){
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(1-ieo+2*idir));
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * (iyzt/2);
        real_t uL[VLEN*NDF2];
        shift_vec0_fw(uL, up, NDF2);
        for(int is = 0; is < m_Ns; ++is){
          shift_vec0_fw(z4, &wp2[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xmb(vL, uL, z4, NC);
          mult_wilson_xm2(vL, &buf2[ibf + NVC*ND2*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 1; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iy == m_Ny-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nx2v * izt);
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo+2*idir));
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_yp2(vL, up, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(iy == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nx2v * izt);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_ym2(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 2; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iz == m_Nz-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy2 * it);
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo+2*idir));
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_zp2(vL, up, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(iz == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy2 * it);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_zm2(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 3; // y-direction boundary
    if(do_comm[idir] == 1){

      if(it == m_Nt-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * ixyz;
        real_t *up = m_Ueo.ptr(VLEN*NDF2*site + NvU2*(ieo+2*idir));
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_tp2_dirac(vL, up, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(it == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * ixyz;
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_tm2_dirac(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

  } // site loop (boundary)

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::L_inv(AFIELD& v,
                                                  const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NCD], xL[NCD], yL[NCD];

    load_vec(vL, &wp2[0], NCD);
    save_vec(&vp2[0], vL, NCD);
    set_vec(yL, m_e[0], vL, NCD);

    for(int is = 1; is < m_Ns-1; ++is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is]/m_dp[is-1], xL, NC);
      save_vec(&vp2[Nin4*is], vL, NCD);
      axpy_vec(yL, m_e[is], vL, NCD);
    }

    int is = m_Ns-1;
    copy_vec(xL, vL, NCD);
    load_vec(vL, &wp2[Nin4*is], NCD);
    add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is]/m_dp[is-1], xL, NC);
    add_aPm5_dirac_vec(vL, real_t(-0.5), yL, NC);
    save_vec(&vp2[Nin4*is], vL, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::U_inv(AFIELD& v,
                                                  const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    int is = m_Ns-1;
    load_vec(vL, &wp2[Nin4*is], NCD);
    scal_vec(vL, real_t(1.0)/(m_dp[is]+m_g), NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);

    set_aPp5_dirac_vec(yL, real_t(0.5), vL, NC);

    for(int is = m_Ns-2; is >= 0; --is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is], xL, NC);
      axpy_vec(vL, -m_f[is], yL, NCD);
      scal_vec(vL, real_t(1.0)/m_dp[is], NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Udag_inv(AFIELD& v,
                                                     const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    load_vec(vL, &wp2[0], NCD);
    scal_vec(vL, real_t(1.0)/m_dp[0], NCD);
    save_vec(&vp2[0], vL, NCD);
    set_vec(yL, m_f[0], vL, NCD);

    for(int is = 1; is < m_Ns-1; ++is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is-1], xL, NC);
      scal_vec(vL, real_t(1.0)/m_dp[is], NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
      axpy_vec(yL, m_f[is], vL, NCD);
    }

    int is = m_Ns - 1;
    copy_vec(xL, vL, NCD);
    load_vec(vL, &wp2[Nin4*is], NCD);
    add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is-1], xL, NC);
    add_aPp5_dirac_vec(vL, real_t(-0.5), yL, NC);
    scal_vec(vL, real_t(1.0)/(m_dp[is] + m_g), NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din_eo<AFIELD>::Ldag_inv(AFIELD& v,
                                                     const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nst2v);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    int is = m_Ns-1;
    load_vec(vL, &wp2[Nin4*is], NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);
    set_aPm5_dirac_vec(yL, real_t(0.5), vL, NC);

    for(int is = m_Ns-2; is >= 0; --is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is+1]/m_dp[is], xL, NC);
      axpy_vec(vL, -m_e[is], yL, NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_General_5din_eo<AFIELD>::flop_count(std::string mode)
{
  int Lvol = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);

  // double flop_Wilson = m_foprw->flop_count("D");
  double flop_Wilson = m_foprw->flop_count();

  double axpy1 = static_cast<double>(2 * m_NinF);
  double scal1 = static_cast<double>(1 * m_NinF);

  double flop_DW = vNs * (flop_Wilson + vsite*(6*axpy1 + 2*scal1));
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW. 

  double flop_LU_inv = 2.0 * vsite * 
            ((3.0*axpy1 + scal1)*(vNs-1.0) + axpy1 + 2.0*scal1);

  double flop = 0.0;
  if(mode == "Prec"){
    flop = flop_LU_inv;
  } else if ((mode == "D") || (mode == "Ddag")){
    flop = flop_DW;
  } else if (mode == "DdagD"){
    flop = 2.0 * flop_DW;
  } else if ((mode == "D_prec") || (mode == "Ddag_prec")){
    flop = flop_LU_inv + flop_DW;
  } else if (mode == "DdagD_prec"){
    flop = 2.0 * (flop_LU_inv + flop_DW);
  }else{
    vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;

}

//============================================================END=====
