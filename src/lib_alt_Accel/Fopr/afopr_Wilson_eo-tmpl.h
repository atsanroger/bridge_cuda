/*!
      @file    afopr_Wilson_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2542 $
*/

template<typename AFIELD>
const std::string AFopr_Wilson_eo<AFIELD>::class_name
                                       = "AFopr_Wilson_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches for coomunication
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nc = CommonParameters::Nc();
  m_Nd = CommonParameters::Nd();
  m_vl = CommonParameters::Vlevel();

  m_Ncol = m_Nc;
  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Nvc  = m_Nc * 2;
  m_Nd   = m_Nd;
  m_Ndim = CommonParameters::Ndim();
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();

  m_Nx2  = m_Nx/2;
  m_Nst2 = m_Nst/2;

  m_Nsize[0] = m_Nx2;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  //  check_Nc(m_Nc);
  //  check_setup();

  // condition check
  if(m_Nx % 2 != 0){
    vout.crucial(m_vl, "%s: Nx must be even.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  m_Ieo_origin = (ipe1 * m_Ny + ipe2 * m_Nz + ipe3 * m_Nt) % 2;

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  // setup of communication buffers
  m_Nbdsize.resize(m_Ndim);
  int Nd2 = m_Nd/2;
  m_Nbdsize[0] = m_Nvc * Nd2 * ceil_nwp((m_Ny * m_Nz * m_Nt + 1)/2);
  m_Nbdsize[1] = m_Nvc * Nd2 * ceil_nwp(m_Nx2 * m_Nz * m_Nt);
  m_Nbdsize[2] = m_Nvc * Nd2 * ceil_nwp(m_Nx2 * m_Ny * m_Nt);
  m_Nbdsize[3] = m_Nvc * Nd2 * ceil_nwp(m_Nx2 * m_Ny * m_Nz);

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst2, 1);
  m_v2.reset(NinF, m_Nst2, 1);

  set_parameters(params);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::setup_channels()
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
void AFopr_Wilson_eo<AFIELD>::tidyup()
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
void AFopr_Wilson_eo<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
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

  set_parameters(real_t(kappa), bc);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::set_parameters(const real_t CKs,
                                            const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
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
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", double(m_CKs));
  params.set_int_vector("boundary_condition", m_boundary);

  std::string repr;
  if(m_repr == DIRAC) repr  = "Dirac";
  if(m_repr == CHIRAL) repr = "Chiral";
  params.set_string("gamma_matrix_type", repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::set_config(Field *u)
{
  int nth = ThreadManager::get_num_threads();

  vout.paranoiac(m_vl, "%s: set_config is called: num_threads = %d\n",
                 class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.paranoiac(m_vl, "%s: set_config finished\n", class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::set_config_omp(Field *u)
{
  vout.paranoiac(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
 {
  set_config_impl(u);
 }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::set_config_impl(Field* u)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) m_conf = u;

  convert(m_index_eo, m_Ueo, *u);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult(AFIELD &v,
                                  const AFIELD &w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_dag(AFIELD &v,
                                  const AFIELD &w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult(AFIELD &v, const AFIELD &w,
                                   const std::string mode)
{
  if(mode == "Dee"){
    copy(v, w);
  }else if(mode == "Doo"){
    copy(v, w);
  }else if(mode == "Dee_inv"){
    copy(v, w);
  }else if(mode == "Doo_inv"){
    copy(v, w);
  }else if(mode == "Deo_test"){
    Meo(v, w, 0);
  }else if(mode == "Deo"){
    Meo(v, w, 0);
  }else if(mode == "Doe"){
    Meo(v, w, 1);
  }else if(mode == "convert"){
    //mult_gm4(v, w);
    copy(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_gm5(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm5(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::DdagD(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::Ddag(AFIELD &v, const AFIELD &w)
{
#pragma omp barrier

  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::D(AFIELD &v, const AFIELD &w)
{
  Meo(m_v1, w, 1);
  Meo(v, m_v1, w, 0, 1);
  //Meo_alt(m_v1, w, 1);
  //Meo_alt(v, m_v1, w, 0, 1);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::aypx(real_t a, AFIELD &v,
                                   const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  aypx(a, vp, wp);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::Meo_alt(AFIELD &v, const AFIELD &w,
                                      const int ieo)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  clear(vp);
  mult_xp(vp, wp, ieo);
  mult_xm(vp, wp, ieo);
  mult_yp(vp, wp, ieo);
  mult_ym(vp, wp, ieo);
  mult_zp(vp, wp, ieo);
  mult_zm(vp, wp, ieo);
  mult_tp(vp, wp, ieo);
  mult_tm(vp, wp, ieo);
  scal(vp, -m_CKs);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::Meo_alt(AFIELD &v,
                                      const AFIELD &w, const AFIELD &x,
                                      const int ieo, const int iflag)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *xp = const_cast<AFIELD*>(&x)->ptr(0);

  clear(vp);
  mult_xp(vp, wp, ieo);
  mult_xm(vp, wp, ieo);
  mult_yp(vp, wp, ieo);
  mult_ym(vp, wp, ieo);
  mult_zp(vp, wp, ieo);
  mult_zm(vp, wp, ieo);
  mult_tp(vp, wp, ieo);
  mult_tm(vp, wp, ieo);

  if(iflag == 1){
    aypx(m_CKs, vp, xp);
  }else{
    scal(vp, -m_CKs);
  }
  
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::Meo(AFIELD &v, const AFIELD &w,
                                  const int ieo)
{
  Meo(v, w, w, ieo, 0);
  //Meo_alt(v, w, ieo);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::Meo(AFIELD &v,
                                  const AFIELD &w, const AFIELD &x,
                                  const int ieo, const int iflag)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
    real_t *xp = const_cast<AFIELD*>(&x)->ptr(0);
    real_t *up = m_Ueo.ptr(0);
    int jeo    = (m_Ieo_origin + ieo) % 2;

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
        BridgeACC::mult_wilson_eo_1_dirac(
                             buf1xp, buf1xm, buf1yp, buf1ym,
                             buf1zp, buf1zm, buf1tp, buf1tm,
                             up, wp, ieo, jeo,
                             m_Nsize, m_bc, do_comm, NC);
      }else{
        BridgeACC::mult_wilson_eo_1_chiral(
                              buf1xp, buf1xm, buf1yp, buf1ym,
                              buf1zp, buf1zm, buf1tp, buf1tm,
                              up, wp, ieo, jeo,
                              m_Nsize, m_bc, do_comm, NC);
      }

      chset_send.start();
      chset_recv.start();
    }

    // bulk part
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_Meo_dirac(vp, up, wp, xp,
                          ieo, jeo, iflag, m_Nsize, m_bc2, m_CKs);
    }else{
      BridgeACC::mult_wilson_Meo_chiral(vp, up, wp, xp,
                           ieo, jeo, iflag, m_Nsize, m_bc2, m_CKs);
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
        BridgeACC::mult_wilson_eo_2_dirac(vp, up,
                             buf2xp, buf2xm, buf2yp, buf2ym,
                             buf2zp, buf2zm, buf2tp, buf2tm,
                             m_CKs, ieo, jeo, iflag,
                             m_Nsize, m_bc, do_comm, NC);
      }else{
        BridgeACC::mult_wilson_eo_2_chiral(vp, up,
                              buf2xp, buf2xm, buf2yp, buf2ym,
                              buf2zp, buf2zm, buf2tp, buf2tm,
                              m_CKs, ieo, jeo, iflag,
                              m_Nsize, m_bc, do_comm, NC);
      }

    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

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
void AFopr_Wilson_eo<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    BridgeACC::mult_wilson_aypx(a, v, w, m_Nsize, NC);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::clear(real_t *v)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    BridgeACC::mult_wilson_clear(v, m_Nsize, NC);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::scal(real_t *v, const real_t a)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    BridgeACC::mult_wilson_scal(v, a, m_Nsize, NC);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_xp(real_t *v2, real_t *v1,
                                      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 0;
    int jeo = (m_Ieo_origin + ieo) % 2;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2*idir) );

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_xp1_eo(buf1, v1, m_Nsize, m_bc, jeo, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_wilson_xpb_eo(v2, u, v1, m_Nsize, m_bc2, jeo, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_xp2_eo(v2, u, buf2, m_Nsize, m_bc, jeo, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_xm(real_t *v2, real_t *v1,
                                        const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 0;
    int jeo = (m_Ieo_origin + ieo) % 2;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2*idir) );

    if(do_comm[idir] > 0){
      BridgeACC::mult_wilson_xm1_eo(buf1, u, v1, m_Nsize, m_bc, jeo, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_wilson_xmb_eo(v2, u, v1, m_Nsize, m_bc2, jeo, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_wilson_xm2_eo(v2, buf2, m_Nsize, m_bc, jeo, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_eo<AFIELD>::mult_yp(real_t *v2, real_t *v1,
				      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 1;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

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
void AFopr_Wilson_eo<AFIELD>::mult_ym(real_t *v2, real_t *v1,
				      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 1;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

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
void AFopr_Wilson_eo<AFIELD>::mult_zp(real_t *v2, real_t *v1,
                                        const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 2;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

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
void AFopr_Wilson_eo<AFIELD>::mult_zm(real_t *v2, real_t *v1,
				      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 2;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

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
void AFopr_Wilson_eo<AFIELD>::mult_tp(real_t *v2, real_t *v1,
				      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 3;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

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
void AFopr_Wilson_eo<AFIELD>::mult_tm(real_t *v2, real_t *v1,
				      const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){

    int idir = 3;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

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
double AFopr_Wilson_eo<AFIELD>::flop_count()
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == DIRAC) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1)));
  } else if (m_repr == CHIRAL) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));
  } else {
    vout.crucial(m_vl, "%s: gamma matrix repr is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

  if (m_mode == "mult_test") flop = flop/2.0;

  return flop;
}

//============================================================END=====
