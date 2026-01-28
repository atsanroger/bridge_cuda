/*!
      @file    afopr_Wilson_CLE-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

template<typename AFIELD>
const std::string AFopr_Wilson_CLE<AFIELD>::class_name
                                           = "AFopr_Wilson_CLE<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  // switch for communication
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication only when necessary

  vout.general(m_vl, "%s: being setup:\n", class_name.c_str());

  m_vl = CommonParameters::Vlevel();

  m_Nc = CommonParameters::Nc();
  m_Nd = CommonParameters::Nd();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  int Nvc = m_Nc * 2;
  int Ndf = 2 * m_Nc * m_Nc;
  int Nx  = CommonParameters::Nx();
  int Ny  = CommonParameters::Ny();
  int Nz  = CommonParameters::Nz();
  int Nt  = CommonParameters::Nt();

  m_Nsize[0] = Nx;
  m_Nsize[1] = Ny;
  m_Nsize[2] = Nz;
  m_Nsize[3] = Nt;

  // Nc check: temporary commented out.
  // check_Nc(m_Nc);
  // check_setup();

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general(m_vl, "  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  int Nd2 = m_Nd/2;
  m_Nbdsize[0] = Nvc * Nd2 * Ny * Nz * Nt;
  m_Nbdsize[1] = Nvc * Nd2 * Nx * Nz * Nt;
  m_Nbdsize[2] = Nvc * Nd2 * Nx * Ny * Nt;
  m_Nbdsize[3] = Nvc * Nd2 * Nx * Ny * Nz;

  setup_channels();

  // gauge configuration.
  m_U.reset(Ndf, m_Nst, m_Ndim);
  m_Uinvd.reset(Ndf, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v2.reset(NinF, m_Nst, 1);
  m_w1.reset(NinF, m_Nst, 1);

  vout.detailed(m_vl, "%s: initalization finished.\n",
  		      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::setup_channels()
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
void AFopr_Wilson_CLE<AFIELD>::tidyup()
{

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
void AFopr_Wilson_CLE<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double kappa, chpot;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("chemical_potential", chpot);
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
    vout.general(m_vl, "  gamma_matrix_type: Dirac\n");
  }else if(repr == "Chiral"){
    m_repr = CHIRAL;
    vout.general(m_vl, "  gamma_matrix_type: Chiral\n");
  }else{
    vout.crucial(m_vl, "Error in %s: irrelevant gamma_matrix_type: %s\n",
 		class_name.c_str(), repr.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(kappa), bc, real_t(chpot));

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::set_parameters(const real_t CKs,
                                          const std::vector<int> bc,
                                          const real_t chpot)
{
  ThreadManager::assert_single_thread(class_name);

  assert(bc.size() == m_Ndim);

  m_CKs = CKs;
  m_chpot = chpot;

  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_bc[mu]  = 1;
    if(do_comm[mu] > 0){ // do communication
      if(Communicator::ipe(mu) == 0) m_bc[mu] = m_boundary[mu];
      m_bc2[mu] = 0;
    }else{  // no communication
      m_bc[mu]  = 0;    // for boundar part (dummy)
      m_bc2[mu] = m_boundary[mu];  // for bulk part
    }
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  CKs   = %8.4f\n", m_CKs);
  vout.general(m_vl, "  chpot = %8.4f\n", m_chpot);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::set_config(Field* u)
{
  ThreadManager::assert_single_thread(class_name);

  m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
 {
  convert_gauge(index_lex, m_U, *u);
 }

  inverse_dag(m_Uinvd, m_U);

  real_t fug_up = exp( m_chpot);
  real_t fug_dn = exp(-m_chpot);

  m_U.scal(    fug_up, m_Ndim-1);
  m_Uinvd.scal(fug_dn, m_Ndim-1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::inverse_dag(AFIELD& Uinvd, AFIELD& U)
{  // take (U^inv)^\dag of gauge field

  real_t e_crit;
  if(sizeof(real_t) == 4){
    e_crit = 1.e-14;
  }else if(sizeof(real_t) == 8){
    e_crit = 1.e-28;
  }else{
    vout.crucial(m_vl, "%s: unsupported precision,\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  real_t *up = Uinvd.ptr(0);
  real_t *vp = U.ptr(0);

  for(int idir = 0; idir < m_Ndim; ++idir){
    BridgeACC::inverse_dag(up, idir, vp, idir, m_Nst);
  }

  // check of inverse_dag
  AFIELD wt(NDF, m_Nst, 1);
  real_t *wp = wt.ptr(0);

  vout.detailed(m_vl, "%s: check of iverse:\n", class_name.c_str());
  for(int idir = 0; idir < m_Ndim; ++idir){
    BridgeACC::mult_Gnd(wp, 0, vp, idir, up, idir, m_Nst);
    BridgeACC::add_unit(wp, 0, real_t(-1.0), m_Nst);
    real_t ww2 = wt.norm2()/real_t(NDF * m_Nst);
    vout.detailed(m_vl, "  idir = %d  norm2 = %e\n", idir, ww2);
    if(ww2 > e_crit){
      vout.detailed(m_vl, "%s: too large deviation: %e\n",
                    class_name.c_str(), idir, ww2);

    }
  }


}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::set_mode(std::string mode)
{
  m_mode = mode;

  vout.paranoiac(m_vl, "%s: mode is set to %s\n",
                 class_name.c_str(), m_mode.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult(AFIELD &v, const AFIELD &w)
{
  if (m_mode == "D") {
    mult_D(v, w);
  } else if (m_mode == "DdagD") {
    mult_DdagD(v, w);
  } else if (m_mode == "Ddag") {
    mult_Ddag(v, w);
  } else if (m_mode == "H") {
    mult_H(v, w);
  } else {
    vout.crucial(m_vl, "%s: undefined mode = %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w)
{
  if (m_mode == "D") {
    mult_Ddag(v, w);
  } else if (m_mode == "DdagD") {
    mult_DdagD(v, w);
  } else if (m_mode == "Ddag") {
    mult_D(v, w);
  } else if (m_mode == "H") {
    mult_H(v, w);
  } else {
    vout.crucial(m_vl, "%s: undefined mode = %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult(AFIELD &v, const AFIELD &w,
                                const std::string mode)
{
  if(mode == "D"){
    mult_D(v, w);
  }else if(mode == "Ddag"){
    mult_Ddag(v, w);
  }else if(mode == "DdagD"){
    mult_DdagD(v, w);
  }else if(mode == "H"){
    mult_H(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w,
                                    const std::string mode)
{
  if(mode == "D"){
    mult_Ddag(v, w);
  }else if(mode == "Ddag"){
    mult_D(v, w);
  }else if(mode == "DdagD"){
    mult_DdagD(v, w);
  }else if(mode == "H"){
    mult_H(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_gm5(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm5(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_DdagD(AFIELD &v, const AFIELD &w)
{
  mult_D(m_v2, w);
  mult_Ddag(v, m_v2);

  //  mult_D_alt(m_v2, w);
  //  mult_Ddag_alt(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_D(AFIELD &v, const AFIELD &w)
{
  if(m_repr != DIRAC){
    mult_D_alt(v, w);
    return;
  }

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u_up = m_U.ptr(0);
  real_t *u_dn = m_Uinvd.ptr(0);

  if(do_comm_any > 0){

    real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
    real_t *buf1xm = (real_t*)chsend_up[0].ptr();

    real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
    real_t *buf1ym = (real_t*)chsend_up[1].ptr();

    real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
    real_t *buf1zm = (real_t*)chsend_up[2].ptr();

    real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
    real_t *buf1tm = (real_t*)chsend_up[3].ptr();

    BridgeACC::mult_wilson_cle_1_dirac(
                             buf1xp, buf1xm, buf1yp, buf1ym,
                             buf1zp, buf1zm, buf1tp, buf1tm,
                             u_dn, v1, m_Nsize, m_bc, do_comm, NC);

  }

  chset_send.start();
  chset_recv.start();

  // bulk part
  BridgeACC::mult_wilson_cle_D_dirac(v2, u_up, u_dn, v1,
                                       m_Nsize, m_bc2, m_CKs);

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

      BridgeACC::mult_wilson_cle_2_dirac(v2, u_up,
                                 buf2xp, buf2xm, buf2yp, buf2ym,
                                 buf2zp, buf2zm, buf2tp, buf2tm,
                                 m_CKs, m_Nsize, m_bc, do_comm, NC);

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_Ddag(AFIELD &v, const AFIELD &w)
{
  if(m_repr != DIRAC){
    mult_Ddag_alt(v, w);
    return;
  }

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);

  real_t *u_up = m_Uinvd.ptr(0);
  real_t *u_dn = m_U.ptr(0);

  mult_gm5(vp, wp);

  if(do_comm_any > 0){

    real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
    real_t *buf1xm = (real_t*)chsend_up[0].ptr();

    real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
    real_t *buf1ym = (real_t*)chsend_up[1].ptr();

    real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
    real_t *buf1zm = (real_t*)chsend_up[2].ptr();

    real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
    real_t *buf1tm = (real_t*)chsend_up[3].ptr();

    BridgeACC::mult_wilson_cle_1_dirac(
                             buf1xp, buf1xm, buf1yp, buf1ym,
                             buf1zp, buf1zm, buf1tp, buf1tm,
                             u_dn, vp, m_Nsize, m_bc, do_comm, NC);

  }

  chset_send.start();
  chset_recv.start();

  // bulk part
  BridgeACC::mult_wilson_cle_D_dirac(yp, u_up, u_dn, vp,
                                     m_Nsize, m_bc2, m_CKs);

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

    BridgeACC::mult_wilson_cle_2_dirac(yp, u_up,
                                  buf2xp, buf2xm, buf2yp, buf2ym,
                                  buf2zp, buf2zm, buf2tp, buf2tm,
                                  m_CKs, m_Nsize, m_bc, do_comm, NC);

  }

  mult_gm5(vp, yp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_D_alt(AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int jd = 1;
  clear(vp);
  mult_xp(vp, wp, jd);
  mult_xm(vp, wp, jd);
  mult_yp(vp, wp, jd);
  mult_ym(vp, wp, jd);
  mult_zp(vp, wp, jd);
  mult_zm(vp, wp, jd);
  mult_tp(vp, wp, jd);
  mult_tm(vp, wp, jd);
  aypx(-m_CKs, vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_Ddag_alt(AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);

  mult_gm5(vp, wp);

  int jd = -1;
  clear(yp);
  mult_xp(yp, vp, jd);
  mult_xm(yp, vp, jd);
  mult_yp(yp, vp, jd);
  mult_ym(yp, vp, jd);
  mult_zp(yp, vp, jd);
  mult_zm(yp, vp, jd);
  mult_tp(yp, vp, jd);
  mult_tm(yp, vp, jd);
  aypx(-m_CKs, yp, vp);

  mult_gm5(vp, yp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_up(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int jd = 1;
  real_t fug_up = exp(m_chpot);

  if(mu == 0){
    mult_xp(vp, wp, jd);
  }else if(mu == 1){
    mult_yp(vp, wp, jd);
  }else if(mu == 2){
    mult_zp(vp, wp, jd);
  }else if(mu == 3){
    mult_tp(vp, wp, jd);
    scal(vp, fug_up);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_dn(int mu, AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int jd = 1;
  real_t fug_dn = exp(-m_chpot);

  if(mu == 0){
    mult_xm(vp, wp, jd);
  }else if(mu == 1){
    mult_ym(vp, wp, jd);
  }else if(mu == 2){
    mult_zm(vp, wp, jd);
  }else if(mu == 3){
    mult_tm(vp, wp, jd);
    scal(vp, fug_dn);
  }else{
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_H(AFIELD &v, const AFIELD &w)
{
  mult_D(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
  if(m_repr == DIRAC){
    BridgeACC::mult_wilson_gm5_dirac(v, w, m_Nsize, NC);
  }else{
    BridgeACC::mult_wilson_gm5_chiral(v, w, m_Nsize, NC);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::gm5_aypx(real_t a, real_t *v, real_t *w)
{
  if(m_repr == DIRAC){
    BridgeACC::mult_wilson_gm5_aypx_dirac(a, v, w, m_Nsize, NC);
  }else{
    BridgeACC::mult_wilson_gm5_aypx_chiral(a, v, w, m_Nsize, NC);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::scal(real_t *v, real_t a)
{
  BridgeACC::mult_wilson_scal(v, a, m_Nsize, NC);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
  BridgeACC::mult_wilson_aypx(a, v, w, m_Nsize, NC);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::axpy(real_t *v, real_t a, real_t *w)
{
  BridgeACC::mult_wilson_axpy(v, a, w, m_Nsize, NC);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::clear(real_t *v)
{
  BridgeACC::mult_wilson_clear(v, m_Nsize, NC);
}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_xp(real_t *v2, real_t *v1, int jd)
{
  int idir = 0;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_xp1(buf1, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeACC::mult_wilson_xpb(v2, u, v1, m_Nsize, m_bc2, NC);

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_xp2(v2, u, buf2, m_Nsize, m_bc, NC);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_xm(real_t *v2, real_t *v1, int jd)
{
  int idir = 0;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_xm1(buf1, u, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeACC::mult_wilson_xmb(v2, u, v1, m_Nsize, m_bc2, NC);

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_xm2(v2, buf2, m_Nsize, m_bc, NC);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_yp(real_t *v2, real_t *v1, int jd)
{
  int idir = 1;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_yp1(buf1, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeACC::mult_wilson_ypb(v2, u, v1, m_Nsize, m_bc2, NC);

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_yp2(v2, u, buf2, m_Nsize, m_bc, NC);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_ym(real_t *v2, real_t *v1, int jd)
{
  int idir = 1;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_ym1(buf1, u, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

   BridgeACC::mult_wilson_ymb(v2, u, v1, m_Nsize, m_bc2, NC);

   if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_ym2(v2, buf2, m_Nsize, m_bc, NC);
   }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_zp(real_t *v2, real_t *v1, int jd)
{
  int idir = 2;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_zp1(buf1, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeACC::mult_wilson_zpb(v2, u, v1, m_Nsize, m_bc2, NC);

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_zp2(v2, u, buf2, m_Nsize, m_bc, NC);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_zm(real_t *v2, real_t *v1, int jd)
{
  int idir = 2;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    BridgeACC::mult_wilson_zm1(buf1, u, v1, m_Nsize, m_bc, NC);
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeACC::mult_wilson_zmb(v2, u, v1, m_Nsize, m_bc2, NC);

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeACC::mult_wilson_zm2(v2, buf2, m_Nsize, m_bc, NC);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_tp(real_t *v2, real_t *v1, int jd)
{
  int idir = 3;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tp1_dirac(buf1, v1, m_Nsize, m_bc, NC);
    }else{
      BridgeACC::mult_wilson_tp1_chiral(buf1, v1, m_Nsize, m_bc, NC);
    }
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  if(m_repr == DIRAC){
    BridgeACC::mult_wilson_tpb_dirac(v2, u, v1, m_Nsize, m_bc2, NC);
  }else{
    BridgeACC::mult_wilson_tpb_chiral(v2, u, v1, m_Nsize, m_bc2, NC);
  }

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tp2_dirac(v2, u, buf2, m_Nsize, m_bc, NC);
    }else{
      BridgeACC::mult_wilson_tp2_chiral(v2, u, buf2, m_Nsize, m_bc, NC);
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Wilson_CLE<AFIELD>::mult_tm(real_t *v2, real_t *v1, int jd)
{
  int idir = 3;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  real_t *u;
  if(jd == 1){
    u = m_Uinvd.ptr(index_lex.idx_G(0, 0, idir));
  }else{
    u = m_U.ptr(index_lex.idx_G(0, 0, idir));
  }

  if(do_comm[idir] > 0){

    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tm1_dirac(buf1, u, v1, m_Nsize, m_bc, NC);
    }else{
      BridgeACC::mult_wilson_tm1_chiral(buf1, u, v1, m_Nsize, m_bc, NC);
    }
    BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tmb_dirac(v2, u, v1, m_Nsize, m_bc2, NC);
    }else{
      BridgeACC::mult_wilson_tmb_chiral(v2, u, v1, m_Nsize, m_bc2, NC);
    }

  if(do_comm[idir] > 0){
    BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_tm2_dirac(v2, buf2, m_Nsize, m_bc, NC);
    }else{
      BridgeACC::mult_wilson_tm2_chiral(v2, buf2, m_Nsize, m_bc, NC);
    }
  }

}

//====================================================================
template<typename AFIELD>
double AFopr_Wilson_CLE<AFIELD>::flop_count()
{
  return flop_count(m_mode);
}

//====================================================================
template<typename AFIELD>
double AFopr_Wilson_CLE<AFIELD>::flop_count(const std::string mode)
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
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}

//====================================================================
//============================================================END=====
