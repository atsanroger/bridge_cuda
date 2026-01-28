/*!
      @file    afopr_Wilson_dd-tmpl.h
      @brief     
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#include "lib_alt_Accel/inline/afield_th-inc.h"
#include "lib_alt_Accel/Field/afield-inc.h"

template<typename AFIELD>
const std::string AFopr_Wilson_dd<AFIELD>::class_name
                                         = "AFopr_Wilson_dd<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches
  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 1 if communication forced any time

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nc = CommonParameters::Nc();
  m_Nd = CommonParameters::Nd();
  m_Ndim = CommonParameters::Ndim();
  m_Nvc  = m_Nc * 2;
  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Ndm2 = m_Nd * m_Nd / 2,

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();

  m_Nsize[0] = m_Nx;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  //  check_Nc(m_Nc);
  //  check_setup();

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  int Nd2 = m_Nd/2;
  m_Nbdsize[0] = m_Nvc * Nd2 * ceil_nwp(m_Ny * m_Nz * m_Nt);
  m_Nbdsize[1] = m_Nvc * Nd2 * ceil_nwp(m_Nx * m_Nz * m_Nt);
  m_Nbdsize[2] = m_Nvc * Nd2 * ceil_nwp(m_Nx * m_Ny * m_Nt);
  m_Nbdsize[3] = m_Nvc * Nd2 * ceil_nwp(m_Nx * m_Ny * m_Nz);

  setup_channels();

  vout.increase_indent();

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nst, m_Ndim);

  // gauge configuration with block condition.
  m_Ublock.reset(NDF, m_Nst, m_Ndim);

  m_Uhop.reset(NDF, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v2.reset(NinF, m_Nst, 1);

  set_parameters(params);

  vout.decrease_indent();
  vout.detailed(m_vl, "%s: initalization finished.\n",
  		      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    int Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu,  1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if(do_comm[mu] == 1){
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

    // openacc device memory allocation
#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_init(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up2, m_Nbdsize[mu]);
#else
    BridgeACC::afield_init((real_t*)buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_init((real_t*)buf_dn, m_Nbdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  // openacc device memory clean up
  for(int mu = 0; mu < m_Ndim; ++mu){

#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up2, m_Nbdsize[mu]);
#else
    real_t* buf_up = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn, m_Nbdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double kappa;
  std::vector<int> bc;
  std::vector<int> block_size;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_int_vector("block_size", block_size);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                                                class_name.c_str());
  exit(EXIT_FAILURE);
 }

  //- setting gamma matrix representation
  std::string repr;
  err = params.fetch_string("gamma_matrix_type", repr);
  if(err){
    vout.general(m_vl, "  gamma_matrix_type is not given - set to Dirac\n");
    m_repr = DIRAC;
  }else if(repr == "Dirac"){
    m_repr = DIRAC;
    vout.general(m_vl, "  gamma_matrix_type: Dirac\n");
  }else if(repr == "Chiral"){  // currently not available.
    m_repr = CHIRAL;
    vout.crucial(m_vl, "%s: chiral repr. is not implemented\n",
		 class_name.c_str());
  }else{
    vout.crucial(m_vl, "Error in %s: irrelevant gamma_matrix_type: %s\n",
		 class_name.c_str(), repr.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  set_parameters is called.\n");
  set_parameters(real_t(kappa), bc, block_size);
  vout.general(m_vl, "  set_parameters finished.\n");

  vout.general(m_vl, "  all set_parameters finished.\n");
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_parameters(
                                   const real_t CKs,
                                   const std::vector<int> bc,
                                   const std::vector<int> block_size)
{
  ThreadManager::assert_single_thread(class_name);

  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
    m_boundary.resize(m_Ndim);
    m_block_size.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu]   = bc[mu];
      m_block_size[mu] = block_size[mu];
    }

    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_bc[mu]  = 1;
      if(do_comm[mu] > 0){ // do communication
        if(Communicator::ipe(mu) == 0) m_bc[mu]  = m_boundary[mu];
        m_bc2[mu] = 0;
      }else{  // no communication
        m_bc[mu]  = 0;    // for boundar part (dummy)
        m_bc2[mu] = m_boundary[mu];  // for bulk part
      }
    }
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  if(m_repr == DIRAC){
    vout.general(m_vl, "  gamma-matrix type = Dirac\n");
  }else{
    vout.general(m_vl, "  gamma-matrix type = Chiral\n");
  }
  vout.general(m_vl, "  kappa = %8.4f\n", m_CKs);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  block_size[%d] = %2d\n", mu, m_block_size[mu]);
  }

  int rem =  m_Nx % m_block_size[0]
           + m_Ny % m_block_size[1]
           + m_Nz % m_block_size[2]
    + m_Nt % m_block_size[3];
  if(rem != 0){
    vout.crucial(m_vl, "%s: block_size is irelevant.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int NBx = m_Nx/m_block_size[0];
  int NBy = m_Ny/m_block_size[1];
  int NBz = m_Nz/m_block_size[2];
  int NBt = m_Nt/m_block_size[3];
  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);
  int ipet = Communicator::ipe(3);
  m_Ieo = (NBx * ipex + NBy * ipey + NBz * ipez + NBt * ipet) % 2;

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_int_vector("block_size", m_block_size);

  std::string repr;
  if(m_repr == DIRAC) repr  = "Dirac";
  if(m_repr == CHIRAL) repr = "Chiral";
  params.set_string("gamma_matrix_type", repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_config(Field* u)
{
  int nth = ThreadManager::get_num_threads();
  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_config_omp(Field* u)
{
#pragma omp parallel
 {
  set_config_impl(u);
 }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_config_impl(Field* u)
{
#pragma omp barrier

  m_timer.reset();
  m_timer.start();

  vout.general(m_vl, "%s: set_config start\n", class_name.c_str());

  m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  convert_gauge(index_lex, m_U, *u);

  copy(m_Ublock, m_U);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    BridgeACC::set_block_config(m_Ublock.ptr(0), m_Nsize, &m_block_size[0]);
  }

  copy(m_Uhop, m_U);
  m_Uhop.axpy(real_t(-1.0), m_Ublock);

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();
  vout.general(m_vl, "%s: set_config finished in %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_block_config(AFIELD& U)
{
#pragma omp barrier

  AIndex_lex<real_t,AFIELD::IMPL> index;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  for(int site = is; site < ns; ++site){
    int ix = site % m_Nx;
    int iy = (site/m_Nx) % m_Ny;
    int iz = (site/(m_Nx * m_Ny)) % m_Nz;
    int it = site/(m_Nx * m_Ny * m_Nz);

    int mu = 0;
    if((ix + 1) % m_block_size[mu] == 0){
      for(int in = 0; in < NDF; ++in){
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 1;
    if((iy + 1) % m_block_size[mu] == 0){
      for(int in = 0; in < NDF; ++in){
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 2;
    if((iz + 1) % m_block_size[mu] == 0){
      for(int in = 0; in < NDF; ++in){
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

    mu = 3;
    if((it + 1) % m_block_size[mu] == 0){
      for(int in = 0; in < NDF; ++in){
        int i = index.idx_G(in, site, mu);
        U.set(i, real_t(0.0));
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0) m_mode = mode;

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
std::string AFopr_Wilson_dd<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult(AFIELD &v, const AFIELD &w)
{
  if(m_mode == "D") {
    //D_alt(v, w);
    D(v, w);
    // mult_sap(v, w, 0);
    // mult_sap(v, w, 1);
  }else if(m_mode == "DdagD") {
    DdagD(v, w);
  }else if(m_mode == "Ddag") {
    Ddag(v, w);
  }else if(m_mode == "H") {
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w)
{
  if(m_mode == "D"){
    Ddag(v, w);
  }else if(m_mode == "DdagD"){
    DdagD(v, w);
  }else if(m_mode == "Ddag"){
    D(v, w);
  }else if(m_mode == "H"){
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult(AFIELD &v, const AFIELD &w,
                                   const std::string mode)
{
  if(mode == "D"){
    D(v, w);
  }else if(mode == "Ddag"){
    Ddag(v, w);
  }else if(mode == "DdagD"){
    DdagD(v, w);
  }else if(mode == "H"){
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w,
                                       const std::string mode)
{
  if(mode == "D"){
    Ddag(v, w);
  }else if(mode == "Ddag"){
    D(v, w);
  }else if(mode == "DdagD"){
    DdagD(v, w);
  }else if(mode == "H"){
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_sap(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{
#pragma omp barrier

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *up = m_Ublock.ptr(0);

  int jeo = (ieo + m_Ieo) % 2;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::mult_wilson_dd_dirac(v2, up, v1,
                                    m_CKs, m_bc2,
                                    m_Nsize, &m_block_size[0], jeo);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_dd(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *up = m_Ublock.ptr(0);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::mult_wilson_D_dirac(v2, up, v1,
                                   m_Nsize, m_bc2, m_CKs);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_gm5(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm5(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_up(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(mu == 0){
    mult_xp(vp, wp, 0);
  }else if(mu == 1){
    mult_yp(vp, wp, 0);
  }else if(mu == 2){
    mult_zp(vp, wp, 0);
  }else if(mu == 3){
    mult_tp(vp, wp, 0);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_dn(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(mu == 0){
    mult_xm(vp, wp, 0);
  }else if(mu == 1){
    mult_ym(vp, wp, 0);
  }else if(mu == 2){
    mult_zm(vp, wp, 0);
  }else if(mu == 3){
    mult_tm(vp, wp, 0);
  }else{
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_dup(AFIELD &v, const AFIELD &w,
                                       const int mu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  //clear(vp);
  v.set(real_t(0.0));

  if(mu == 0){
    mult_xp(vp, wp, 1);
  }else if(mu == 1){
    mult_yp(vp, wp, 1);
  }else if(mu == 2){
    mult_zp(vp, wp, 1);
  }else if(mu == 3){
    mult_tp(vp, wp, 1);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

  scal_local(vp, -m_CKs);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_ddn(AFIELD &v, const AFIELD &w,
                                       const int mu)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  v.set(real_t(0.0));

  if(mu == 0){
    mult_xm(vp, wp, 1);
  }else if(mu == 1){
    mult_ym(vp, wp, 1);
  }else if(mu == 2){
    mult_zm(vp, wp, 1);
  }else if(mu == 3){
    mult_tm(vp, wp, 1);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

  scal_local(vp, -m_CKs);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::DdagD(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::Ddag(AFIELD &v, const AFIELD &w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::scal_local(real_t *v, real_t a)
{
#pragma omp master
  {
    int Nin = 2 * m_Nc * m_Nd;
    BridgeACC::scal(v, 0, a, Nin, m_Nst);
  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::D_alt(AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *xp = m_v1.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  m_v1.set(real_t(0.0));
  mult_xp(xp, wp, 0);
  mult_xm(xp, wp, 0);
  mult_yp(xp, wp, 0);
  mult_ym(xp, wp, 0);
  mult_zp(xp, wp, 0);
  mult_zm(xp, wp, 0);
  mult_tp(xp, wp, 0);
  mult_tm(xp, wp, 0);

  axpy(v, -m_CKs, m_v1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::D(AFIELD &v, const AFIELD &w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u  = m_U.ptr(0);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    if(do_comm_any > 0){

      real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
      real_t *buf1xm = (real_t*)chsend_up[0].ptr();

      real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
      real_t *buf1ym = (real_t*)chsend_up[1].ptr();

      real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
      real_t *buf1zm = (real_t*)chsend_up[2].ptr();

      real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
      real_t *buf1tm = (real_t*)chsend_up[3].ptr();

      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_1_dirac(
                            buf1xp, buf1xm, buf1yp, buf1ym,
                            buf1zp, buf1zm, buf1tp, buf1tm,
                            u, v1, m_Nsize, m_bc, do_comm, NC);
      }else{
        BridgeACC::mult_wilson_1_chiral(
                            buf1xp, buf1xm, buf1yp, buf1ym,
                            buf1zp, buf1zm, buf1tp, buf1tm,
                            u, v1, m_Nsize, m_bc, do_comm, NC);
      }

    }

    chset_send.start();
    chset_recv.start();

    // bulk part
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_D_dirac(v2, u, v1, m_Nsize, m_bc2, m_CKs);
    }else{
      BridgeACC::mult_wilson_D_chiral(v2, u, v1, m_Nsize, m_bc2, m_CKs);
    }

    chset_send.wait();
    chset_recv.wait();

    if(do_comm_any > 0){

      real_t *buf2xp = (real_t*)chrecv_up[0].ptr();
      real_t *buf2xm = (real_t*)chrecv_dn[0].ptr();

      real_t *buf2yp = (real_t*)chrecv_up[1].ptr();
      real_t *buf2ym = (real_t*)chrecv_dn[1].ptr();

      real_t *buf2zp = (real_t*)chrecv_up[2].ptr();
      real_t *buf2zm = (real_t*)chrecv_dn[2].ptr();

      real_t *buf2tp = (real_t*)chrecv_up[3].ptr();
      real_t *buf2tm = (real_t*)chrecv_dn[3].ptr();

      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_2_dirac(v2, u,
                          buf2xp, buf2xm, buf2yp, buf2ym,
                          buf2zp, buf2zm, buf2tp, buf2tm,
                          m_CKs, m_Nsize, m_bc, do_comm, NC);
      }else{
        BridgeACC::mult_wilson_2_chiral(v2, u,
                          buf2xp, buf2xm, buf2yp, buf2ym,
                          buf2zp, buf2zm, buf2tp, buf2tm,
                          m_CKs, m_Nsize, m_bc, do_comm, NC);
      }

    }

  } // if(ith == 0)

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::H(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::H_alt(AFIELD &v, const AFIELD &w)
{
  D_alt(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_gm5_dirac(v, w, m_Nsize, NC);
    }else{
      BridgeACC::mult_wilson_gm5_chiral(v, w, m_Nsize, NC);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_xp(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 0;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_xp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_wilson_xpb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_xp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_xm(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 0;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_xm1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_wilson_xmb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_xm2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_yp(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 1;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_yp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_wilson_ypb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_yp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_ym(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 1;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_ym1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_wilson_ymb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_ym2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_zp(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 2;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_zp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_wilson_zpb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_zp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_zm(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 2;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_zm1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_wilson_zmb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_zm2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_tp(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 3;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){

      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_tp1_dirac(buf1, v1, m_Nsize, m_bc, NC);
      }else{
        BridgeACC::mult_wilson_tp1_chiral(buf1, v1, m_Nsize, m_bc, NC);
      }
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tpb_dirac(v2, u, v1, m_Nsize, m_bc2, NC);
    }else{
      BridgeACC::mult_wilson_tpb_chiral(v2, u, v1, m_Nsize, m_bc2, NC);
    }

    if(do_comm[idir] > 0){

      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_tp2_dirac(v2, u, buf2, m_Nsize, m_bc, NC);
      }else{
        BridgeACC::mult_wilson_tp2_chiral(v2, u, buf2, m_Nsize, m_bc, NC);
      }
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_dd<AFIELD>::mult_tm(real_t *v2, real_t *v1, int isap)
{
#pragma omp barrier

  int idir = 3;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));
    if(isap != 0) u = m_Uhop.ptr(m_Ndf * m_Nst * idir);

    if(do_comm[idir] > 0){

      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_tm1_dirac(buf1, u, v1, m_Nsize, m_bc, NC);
      }else{
        BridgeACC::mult_wilson_tm1_chiral(buf1, u, v1, m_Nsize, m_bc, NC);
      }
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tmb_dirac(v2, u, v1, m_Nsize, m_bc2, NC);
    }else{
      BridgeACC::mult_wilson_tmb_chiral(v2, u, v1, m_Nsize, m_bc2, NC);
    }

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      if(m_repr == DIRAC){
        BridgeACC::mult_wilson_tm2_dirac(v2, buf2, m_Nsize, m_bc, NC);
      }else{
        BridgeACC::mult_wilson_tm2_chiral(v2, buf2, m_Nsize, m_bc, NC);
      }
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Wilson_dd<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == DIRAC) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1)));
  } else if (m_repr == CHIRAL) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}

//====================================================================
template<typename AFIELD>
double AFopr_Wilson_dd<AFIELD>::flop_count_sap()
{
  // The following implentastion was copied from SIMD implementation
  // by Kanamori-san. To be confirmed. [08 May 2021 H.Matsufuru]

  int NPE  = CommonParameters::NPE();
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();

  int NBx = m_Nx/m_block_size[0];
  int NBy = m_Ny/m_block_size[1];
  int NBz = m_Nz/m_block_size[2];
  int NBt = m_Nt/m_block_size[3];
  size_t nvol2 = NBx * NBy * NBz * NBt /2;

  int block_x = m_block_size[0];
  int block_y = m_block_size[1];
  int block_z = m_block_size[2];
  int block_t = m_block_size[3];

  int hop_x_site = Nc * Nd * (4 * Nc + 2);
  int hop_y_site = Nc * Nd * (4 * Nc + 2);
  int hop_z_site = Nc * Nd * (4 * Nc + 2);
  int hop_t_site = Nc * Nd * (4 * Nc + 1);
  int accum_site = 4 * Nc * Nd;

  // Dirac representation is assumed.
  double flop_x = 2.0 * static_cast<double>(hop_x_site
                      * (block_x-1) * block_y * block_z * block_t);
  double flop_y = 2.0 * static_cast<double>(hop_y_site
                      * block_x * (block_y-1) * block_z * block_t);
  double flop_z = 2.0 * static_cast<double>(hop_z_site
                      * block_x * block_y * (block_z-1) * block_t);
  double flop_t = 2.0 * static_cast<double>(hop_t_site
                      * block_x * block_y * block_z * (block_t-1));
  double flop_sap_mult = flop_x + flop_y + flop_z + flop_t
                   + static_cast<double>(accum_site);

  return flop_sap_mult * static_cast<double>(nvol2)
                       * static_cast<double>(NPE);

}

//============================================================END=====
