/*!
      @file    afopr_Clover-tmpl.h
      @brief     
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/


template<typename AFIELD>
const std::string AFopr_Clover<AFIELD>::class_name
                                           = "AFopr_Clover<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches for coomunication
  // int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  m_memory_saved = false;  // memory not saved
  //m_memory_saved = true;  // memory saved
  // Cation: set_config() must be called outside OMP parallel region

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  //  m_vl = CommonParameters::Vlevel();

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

  if(!m_memory_saved){
    m_fopr_ct = new AFopr_CloverTerm<AFIELD>(params);
  }
  m_T.reset(m_Ndf, m_Nst, m_Ndm2);

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nst, m_Ndim);


  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  if(!m_memory_saved) m_v1.reset(NinF, m_Nst, 1);
  m_v2.reset(NinF, m_Nst, 1);

  set_parameters(params);

  vout.decrease_indent();
  vout.detailed(m_vl, "%s: initalization finished.\n",
  		      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::setup_channels()
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

    // device memory allocation
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
void AFopr_Clover<AFIELD>::tidyup()
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

  if(!m_memory_saved)  delete m_fopr_ct;

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);

  err += params.fetch_int_vector("boundary_condition", bc);
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
  }else if(repr == "Chiral"){
    m_repr = CHIRAL;
  }else{
    vout.crucial(m_vl, "Error in %s: irrelevant gamma_matrix_type: %s\n",
		 class_name.c_str(), repr.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(kappa), real_t(cSW), bc);

  if(!m_memory_saved) m_fopr_ct->set_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_parameters(const real_t CKs,
                                          const real_t cSW,
                                          const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
    m_cSW = cSW;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
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
  vout.general(m_vl, "  cSW   = %8.4f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_double("clover_coefficient", double(m_cSW));
  params.set_int_vector("boundary_condition", m_boundary);

  std::string repr;
  if(m_repr == DIRAC) repr  = "Dirac";
  if(m_repr == CHIRAL) repr = "Chiral";
  params.set_string("gamma_matrix_type", repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_config(Field* u)
{
  if(m_memory_saved){
    set_config_saved(u);
  }else{
    int nth = ThreadManager::get_num_threads();
    if (nth > 1) {
      set_config_impl(u);
    } else {
      set_config_omp(u);
    }
  }

}
//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_config_omp(Field *u)
{
#pragma omp parallel
 {
  set_config_impl(u);
 }
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_config_impl(Field* u)
{
#pragma omp barrier

  m_timer.reset();
  m_timer.start();

  vout.general(m_vl, "%s: set_config start\n", class_name.c_str());

  m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  convert_gauge(index_lex, m_U, *u);

  m_timer.stop();
  double elapsed_time1 = m_timer.elapsed_sec();
  vout.general(m_vl, "  convert_gauge: %10.6f [sec]\n", elapsed_time1);

  m_timer.reset();
  m_timer.start();

  vout.increase_indent();
  m_fopr_ct->set_config(u);
  vout.decrease_indent();

  m_timer.stop();
  double elapsed_time2 = m_timer.elapsed_sec();
  vout.general(m_vl, "  ct set_config: %10.6f [sec]\n", elapsed_time2);

  m_timer.reset();
  m_timer.start();

  m_fopr_ct->get_csw(m_T);

  m_timer.stop();
  double elapsed_time3 = m_timer.elapsed_sec();
  vout.general(m_vl, "  set_csw:       %10.6f [sec]\n", elapsed_time3);

  double elapsed_time = elapsed_time1 + elapsed_time2 + elapsed_time3;
  vout.general(m_vl, "%s: set_config finished in %10.6f [sec]\n",
               class_name.c_str(), elapsed_time);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_config_saved(Field* u)
{
  ThreadManager::assert_single_thread(class_name);
  vout.crucial(m_vl, "%s: set_config_saved start\n", class_name.c_str());

  m_timer.reset();
  m_timer.start();

  vout.general(m_vl, "%s: set_config start\n", class_name.c_str());

  m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
 {
  convert_gauge(index_lex, m_U, *u);
 }

  m_timer.stop();
  double elapsed_time1 = m_timer.elapsed_sec();
  vout.general(m_vl, "  convert_gauge: %10.6f [sec]\n", elapsed_time1);

  m_timer.reset();
  m_timer.start();

  vout.increase_indent();

  Parameters params;
  get_parameters(params);
  m_fopr_ct = new AFopr_CloverTerm<AFIELD>(params);

  m_fopr_ct->set_config(u);

  vout.decrease_indent();

  m_timer.stop();
  double elapsed_time2 = m_timer.elapsed_sec();
  vout.general(m_vl, "  ct set_config: %10.6f [sec]\n", elapsed_time2);

  m_timer.reset();
  m_timer.start();

  m_fopr_ct->get_csw(m_T);

  delete m_fopr_ct;

  m_timer.stop();
  double elapsed_time3 = m_timer.elapsed_sec();
  vout.general(m_vl, "  set_csw:       %10.6f [sec]\n", elapsed_time3);

  double elapsed_time = elapsed_time1 + elapsed_time2 + elapsed_time3;
  vout.general(m_vl, "%s: set_config finished in %10.6f [sec]\n",
               class_name.c_str(), elapsed_time);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
std::string AFopr_Clover<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult(AFIELD &v, const AFIELD &w)
{
  if(m_mode == "D") {
    D(v, w);
    //D_alt(v, w);
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
void AFopr_Clover<AFIELD>::mult_dag(AFIELD &v,
                                  const AFIELD &w)
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
void AFopr_Clover<AFIELD>::mult(AFIELD &v, const AFIELD &w,
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
void AFopr_Clover<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w,
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
void AFopr_Clover<AFIELD>::mult_gm5(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm5(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_up(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(mu == 0){
    mult_xp(vp, wp);
  }else if(mu == 1){
    mult_yp(vp, wp);
  }else if(mu == 2){
    mult_zp(vp, wp);
  }else if(mu == 3){
    mult_tp(vp, wp);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_dn(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(mu == 0){
    mult_xm(vp, wp);
  }else if(mu == 1){
    mult_ym(vp, wp);
  }else if(mu == 2){
    mult_zm(vp, wp);
  }else if(mu == 3){
    mult_tm(vp, wp);
  }else{
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::DdagD(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::Ddag(AFIELD &v, const AFIELD &w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::D_alt(AFIELD &v, const AFIELD &w)
{
  if(m_memory_saved){

    vout.detailed("%s: D_alt is not available in memory saved mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

  }else{

#pragma omp barrier
    real_t *vp = v.ptr(0);
    real_t *xp = m_v1.ptr(0);
    real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

    m_v1.set(real_t(0.0));
    mult_xp(xp, wp);
    mult_xm(xp, wp);
    mult_yp(xp, wp);
    mult_ym(xp, wp);
    mult_zp(xp, wp);
    mult_zm(xp, wp);
    mult_tp(xp, wp);
    mult_tm(xp, wp);

    m_fopr_ct->mult_csw(v, w);

    axpy(v, -m_CKs, m_v1);
#pragma omp barrier

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::D(AFIELD &v, const AFIELD &w)
{
#pragma omp barrier

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u  = m_U.ptr(0);
  real_t *csw = m_T.ptr(0);

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

      chset_recv.start();
      chset_send.start();
    }

    // bulk part
    if(m_repr == DIRAC){
      BridgeACC::mult_clover_D_dirac(v2, u, csw, v1, m_Nsize, m_bc2, m_CKs);
    }else{
      BridgeACC::mult_clover_D_chiral(v2, u, csw, v1, m_Nsize, m_bc2, m_CKs);
    }

    if(do_comm_any > 0){

      chset_send.wait();
      chset_recv.wait();

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

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::H(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::H_alt(AFIELD &v, const AFIELD &w)
{
  D_alt(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
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
void AFopr_Clover<AFIELD>::multadd_csw(real_t *v2, real_t *v1)
{
#pragma omp barrier

  vout.crucial(m_vl, "%s: multadd_csw is not implemented.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::multadd_csw_dirac(real_t *v2,
                                               real_t *v1,
                                               int site)
{
#pragma omp barrier

  vout.crucial(m_vl, "%s: multadd_csw is not implemented.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
  int idir = 1;
 
  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
  int idir = 2;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
  int idir = 3;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Clover<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir = 3;

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
double AFopr_Clover<AFIELD>::flop_count()
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == DIRAC) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd); // <- clover term
  } else if (m_repr == CHIRAL) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd);   // <- clover term
  } else {
    vout.crucial(m_vl, "%s: input repr is undefined.\n",
                 class_name.c_str());
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

  return flop;
}

//============================================================END=====
