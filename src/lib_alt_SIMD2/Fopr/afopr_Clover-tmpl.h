/*!
      @file    afopr_Clover-tmpl.h
      @brief     
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/
// afopr_Clover-tmpl.h :  template file for AFopr_Clover class

template<typename AFIELD>
const std::string AFopr_Clover<AFIELD>::class_name
                                           = "AFopr_Clover<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::init()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  // switches
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  m_vl = AFopr<AFIELD>::m_vl;

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

  check_Nc(m_Nc);
  check_setup();

  m_Nxv  = m_Nx/VLEN2;
  m_Nstv = m_Nst/VLEN2;
  vout.general(m_vl, "  Nxv = %d  Nstv = %d\n",
               m_Nxv, m_Nstv);

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  int Nd2 = m_Nd/2;
  m_Nbdsize[0] = m_Nvc * Nd2 * m_Ny * m_Nz * m_Nt;
  m_Nbdsize[1] = m_Nvc * Nd2 * m_Nx * m_Nz * m_Nt;
  m_Nbdsize[2] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nt;
  m_Nbdsize[3] = m_Nvc * Nd2 * m_Nx * m_Ny * m_Nz;

  setup_channels();

  m_fopr_ct = new AFopr_CloverTerm<AFIELD>;

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nst, m_Ndim);

  // clover term.
  m_T.reset(m_Ndf, m_Nst, m_Ndm2);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst, 1);
  m_v2.reset(NinF, m_Nst, 1);

  vout.detailed(m_vl, "%s: initalization finished.\n",
  		      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::tidyup()
{
  //  ThreadManager_OpenMP::assert_single_thread(class_name);

  delete m_fopr_ct;

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
void AFopr_Clover<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

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

  set_parameters(real_t(kappa), real_t(cSW), bc);

  m_fopr_ct->set_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_parameters(const real_t CKs,
                                          const real_t cSW,
                                          const std::vector<int> bc)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  assert(bc.size() == m_Ndim);

  m_CKs = CKs;
  m_cSW = cSW;
  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  CKs  = %8.4f\n", m_CKs);
  vout.general(m_vl, "  cSW  = %8.4f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_config(Field* u)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_timer.reset();
  m_timer.start();

  vout.paranoiac(m_vl, "%s: set_config start\n", class_name.c_str());

  m_conf = u;

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
 {
  convert_gauge(index_lex, m_U, *u);
 }

  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1) set_boundary_config(m_U, mu);
  }

  m_fopr_ct->set_config(u);
  m_fopr_ct->get_csw(m_T);

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();
  vout.paranoiac(m_vl, "%s: set_config finished in %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::set_boundary_config(AFIELD& U,
                                                    const int mu)
{
  int ipe[4], Nsize[4], Lsize[4];
  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

  real_t bc = (real_t)(m_boundary[mu]);

  for(int i = 0; i < m_Ndim; ++i){
    ipe[i]   = Communicator::ipe(i);
    Nsize[i] = CommonParameters::Nsize(i);
    Lsize[i] = CommonParameters::Lsize(i);
  }

#pragma omp parallel
  {
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nst);

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

}

//====================================================================
template<typename AFIELD>
std::string AFopr_Clover<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult(AFIELD &v,
                                  const AFIELD &w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_dag(AFIELD &v,
                                  const AFIELD &w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
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
    return DdagD(v, w);
  }else if(mode == "H"){
    return H(v, w);
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
void AFopr_Clover<AFIELD>::mult_gm4(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm4(vp, wp);

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
  real_t *vp = m_v1.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  clear(vp);
  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);

  m_fopr_ct->mult_csw(v, w);
  axpy(v, -m_CKs, m_v1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::D(AFIELD &v, const AFIELD &w)
{
#pragma omp barrier

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  int Nxy = m_Nxv * m_Ny;
  int Nxyz = m_Nxv * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nstv);

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int iy  = iyzt % m_Ny;
    int izt = site/Nxy;
    int iz  = izt % m_Nz;
    int it  = izt/m_Nz;
    int ixy = ix + m_Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    if(do_comm[0] == 1){
      if(ix == 0){
        real_t *buf1 = (real_t*)chsend_dn[0]->ptr();
        mult_wilson_xp1(&buf1[NVC*ND2*iyzt], &v1[VLEN * NCD*site], NC);
      }
      if(ix == m_Nxv-1 ){
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 0);
        real_t *buf1 = (real_t*)chsend_up[0]->ptr();
        mult_wilson_xm1(&buf1[NVC*ND2*iyzt], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == 0){
        real_t *buf1 = (real_t*)chsend_dn[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_yp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iy == m_Ny-1){
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 1);
        real_t *buf1 = (real_t*)chsend_up[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_ym1(&buf1[ibf], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == 0){
        real_t *buf1 = (real_t*)chsend_dn[2]->ptr();
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iz == m_Nz-1){
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 2);
        real_t *buf1 = (real_t*)chsend_up[2]->ptr();
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zm1(&buf1[ibf], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[3] == 1){
      if(it == 0){
        real_t *buf1 = (real_t*)chsend_dn[3]->ptr();
        mult_wilson_tp1_dirac(&buf1[VLEN*NC*ND2*ixyz],
                              &v1[VLEN*NCD*site], NC);
      }
      if(it == m_Nt-1){
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 3);
        real_t *buf1 = (real_t*)chsend_up[3]->ptr();
        mult_wilson_tm1_dirac(&buf1[VLEN*NC*ND2*ixyz],
                              &u[VLEN*NDF2*site], &v1[VLEN*NCD*site], NC);
      }
    }

  }

#pragma omp barrier

#pragma omp master
  {
    chset_send.start();
    chset_recv.start();
  }

  // bulk
  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int iy  = iyzt % m_Ny;
    int izt = site/Nxy;
    int iz  = izt % m_Nz;
    int it  = izt/m_Nz;
    int ixy = ix + m_Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    clear_vec(v2v, NCD);

    if(ix < m_Nxv-1 || do_comm[0] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 0);
      int nei = ix+1 + m_Nxv * iyzt;
      if(ix == m_Nxv-1) nei = 0 + m_Nxv * iyzt;
      shift_vec2_bw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
      mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
    }

    if(ix > 0 ||  do_comm[0] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 0);
      int nei = ix-1 + m_Nxv * iyzt;
      if(ix == 0) nei = m_Nxv-1 + m_Nxv * iyzt;
      shift_vec2_fw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
      shift_vec2_fw(uL, &u[VLEN*NDF2*site], &u[VLEN*NDF2*nei], NDF2);
      mult_wilson_xmb(v2v, uL, zL, NC);
    }

    if(iy < m_Ny-1 ||  do_comm[1] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 1);
      int iy2 = (iy + 1) % m_Ny;
      int nei = ix + m_Nxv * (iy2 + m_Ny * izt);
      mult_wilson_ypb(v2v, &u[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }

    if(iy > 0 ||  do_comm[1] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 1);
      real_t *buf2 = (real_t*)chrecv_dn[1]->ptr();
      int iy2 = (iy - 1 + m_Ny) % m_Ny; 
      int nei = ix + m_Nxv * (iy2 + m_Ny * izt);
      mult_wilson_ymb(v2v, &u[VLEN*NDF2*nei], &v1[VLEN*NCD*nei], NC);
    }

    if(iz < m_Nz-1 ||  do_comm[2] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 2);
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz*it);
      mult_wilson_zpb(v2v, &u[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }

    if(iz > 0 ||  do_comm[2] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 2);
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz*it);
      mult_wilson_zmb(v2v, &u[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }

    if(it < m_Nt-1 ||  do_comm[3] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 3);
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tpb_dirac(v2v, &u[VLEN*NDF2*site],
                            &v1[VLEN*NCD*nei], NC);
    }

    if(it > 0 ||  do_comm[3] == 0){
      real_t *u = m_U.ptr(m_Ndf * m_Nst * 3);
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tmb_dirac(v2v, &u[VLEN*NDF2*nei],
                            &v1[VLEN*NCD*nei], NC);
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp master
  {
    chset_send.wait();
    chset_recv.wait();
  }

#pragma omp barrier

  // boundary
  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int iy  = iyzt % m_Ny;
    int izt = site/Nxy;
    int iz  = izt % m_Nz;
    int it  = izt/m_Nz;
    int ixy = ix + m_Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    load_vec(v2v, &v2[VLEN*NCD*site], NCD);

    if(do_comm[0] == 1){
      if(ix == m_Nxv-1){
	real_t *u = m_U.ptr(m_Ndf * m_Nst * 0);
        real_t *buf2 = (real_t*)chrecv_up[0]->ptr();
	shift_vec0_bw(zL, &v1[VLEN*NCD*site], NCD);
        mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
        mult_wilson_xp2(v2v, &u[VLEN*NDF2*site],
			&buf2[NVC*ND2*iyzt], NC);
      }
      if(ix == 0){
	real_t *buf2 = (real_t*)chrecv_dn[0]->ptr();
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 0);
	shift_vec0_fw(zL, &v1[VLEN*NCD*site], NCD);
        shift_vec0_fw(uL, &u[VLEN*NDF2*site], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
        mult_wilson_xm2(v2v, &buf2[NVC*ND2*iyzt], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == m_Ny-1){
        real_t *buf2 = (real_t*)chrecv_up[1]->ptr();
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 1);
        int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_yp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
      }
      if(iy == 0){
        real_t *buf2 = (real_t*)chrecv_dn[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
        mult_wilson_ym2(v2v, &buf2[ibf], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == m_Nz-1){
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 2);
        real_t *buf2 = (real_t*)chrecv_up[2]->ptr();
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
      }
      if(iz == 0){
        real_t *buf2 = (real_t*)chrecv_dn[2]->ptr();
        int ibf = VLEN * NC *ND2 * (ixy + Nxy * it);
        mult_wilson_zm2(v2v, &buf2[ibf], NC);
      }
    }

    if(do_comm[3] == 1){
      if(it == m_Nt-1){
        real_t *buf2 = (real_t*)chrecv_up[3]->ptr();
        real_t *u = m_U.ptr(m_Ndf * m_Nst * 3);
        mult_wilson_tp2_dirac(v2v, &u[VLEN*NDF2*site],
                              &buf2[VLEN*NC*ND2*ixyz], NC);
      }
      if(it == 0){
        real_t *buf2 = (real_t*)chrecv_dn[3]->ptr();
        mult_wilson_tm2_dirac(v2v, &buf2[VLEN*NC*ND2*ixyz], NC);
      }
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

  for(int site = is; site < ns; ++site){

    load_vec(v2v, &v2[VLEN*NCD*site], NCD);

    scal_vec(v2v, -m_CKs, NCD);

    Vsimd_t v1v[NCD];
    load_vec(v1v, &v1[VLEN*NCD*site], NCD);

    // this is in Dirac representation
    Vsimd_t ut[NDF2], wt1, wt2;
    real_t *u = m_T.ptr(0);
    for(int jd = 0; jd < ND2; ++jd){
      for(int id = 0; id < ND; ++id){
        int ig = VLEN * NDF2 * (site + m_Nstv * (id + ND * jd));
        load_vec(ut, &u[ig], NDF2);
        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          int id2 = (id + ND2) % ND;
          mult_uv(wt1, &ut[ic2], &v1v[NC*id],  NC);
          mult_uv(wt2, &ut[ic2], &v1v[NC*id2], NC);
          int icd1 = ic + NC * jd;
          int icd2 = ic + NC * (jd + ND2);
          axpy_vec(&v2v[icd1], real_t(1.0), &wt1, 1);
          axpy_vec(&v2v[icd2], real_t(1.0), &wt2, 1);
	}
      }
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

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
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  clear(vp);

  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);

  gm5_aypx(-m_CKs, vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  for(int site = is; site < ns; ++site){
    for(int ic = 0; ic < NC; ++ic){
      v[index.idx_SPr(ic,0,site,0)] = w[index.idx_SPr(ic,2,site,0)];
      v[index.idx_SPi(ic,0,site,0)] = w[index.idx_SPi(ic,2,site,0)];
      v[index.idx_SPr(ic,1,site,0)] = w[index.idx_SPr(ic,3,site,0)];
      v[index.idx_SPi(ic,1,site,0)] = w[index.idx_SPi(ic,3,site,0)];
      v[index.idx_SPr(ic,2,site,0)] = w[index.idx_SPr(ic,0,site,0)];
      v[index.idx_SPi(ic,2,site,0)] = w[index.idx_SPi(ic,0,site,0)];
      v[index.idx_SPr(ic,3,site,0)] = w[index.idx_SPr(ic,1,site,0)];
      v[index.idx_SPi(ic,3,site,0)] = w[index.idx_SPi(ic,1,site,0)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_gm4(real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  for(int site = is; site < ns; ++site){
    for(int ic = 0; ic < NC; ++ic){
      v[index.idx_SPr(ic,0,site,0)] =  w[index.idx_SPr(ic,0,site,0)];
      v[index.idx_SPi(ic,0,site,0)] =  w[index.idx_SPi(ic,0,site,0)];
      v[index.idx_SPr(ic,1,site,0)] =  w[index.idx_SPr(ic,1,site,0)];
      v[index.idx_SPi(ic,1,site,0)] =  w[index.idx_SPi(ic,1,site,0)];
      v[index.idx_SPr(ic,2,site,0)] = -w[index.idx_SPr(ic,2,site,0)];
      v[index.idx_SPi(ic,2,site,0)] = -w[index.idx_SPi(ic,2,site,0)];
      v[index.idx_SPr(ic,3,site,0)] = -w[index.idx_SPr(ic,3,site,0)];
      v[index.idx_SPi(ic,3,site,0)] = -w[index.idx_SPi(ic,3,site,0)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::gm5_aypx(real_t a, real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  int id1 = 0;
  int id2 = m_Nvc;
  int id3 = m_Nvc*2;
  int id4 = m_Nvc*3;
  int Ncd = m_Nvc * m_Nd;

  for(int site = is; site < ns; ++site){
    int i = Ncd * site;
    for(int ivc = 0; ivc < m_Nvc; ++ivc){
      real_t v3 = v[ivc + id3 + i];
      real_t v4 = v[ivc + id4 + i];
      v[ivc + id3 + i] = w[ivc + id1 + i] + a * v[ivc + id1 + i];
      v[ivc + id4 + i] = w[ivc + id2 + i] + a * v[ivc + id2 + i];
      v[ivc + id1 + i] = w[ivc + id3 + i] + a * v3;
      v[ivc + id2 + i] = w[ivc + id4 + i] + a * v4;
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::gm5_aypx(real_t a, real_t *v,
                                      real_t *w,
                                      const int is, const int ns)
{
  int id1 = 0;
  int id2 = m_Nvc;
  int id3 = m_Nvc*2;
  int id4 = m_Nvc*3;
  int Ncd = m_Nvc * m_Nd;

  for(int site = is; site < ns; ++site){
    int i = Ncd * site;
    for(int ivc = 0; ivc < m_Nvc; ++ivc){
      real_t v3 = v[ivc + id3 + i];
      real_t v4 = v[ivc + id4 + i];
      v[ivc + id3 + i] = w[ivc + id1 + i] + a * v[ivc + id1 + i];
      v[ivc + id4 + i] = w[ivc + id2 + i] + a * v[ivc + id2 + i];
      v[ivc + id1 + i] = w[ivc + id3 + i] + a * v3;
      v[ivc + id2 + i] = w[ivc + id4 + i] + a * v4;
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      int i2 = icd + Ncd * i;
      v[i2] = a * v[i2] + w[i2];
    }
  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::aypx(real_t a, real_t *v,
                                  real_t *w,
                                  const int is, const int ns)
{
  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      int i2 = icd + Ncd * i;
      v[i2] = a * v[i2] + w[i2];
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  real_t zero = 0.0;
  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      v[icd + Ncd * i] = zero;
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::clear(real_t *v,
                                   const int is, const int ns)
{
  real_t zero = 0.0;
  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      v[icd + Ncd * i] = zero;
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::multadd_csw(real_t *v2, real_t *v1)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  for(int site = is; site < ns; ++site){
    multadd_csw_dirac(v2, v1, site);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::multadd_csw_dirac(real_t *v2,
                                               real_t *v1,
                                               int site)
{
  int Nd2  = m_Nd / 2;

  real_t* u = m_T.ptr(0);

  real_t vt[VLEN*NCD], wt[VLEN*NCD];
  real_t ut[VLEN*NDF2];
  real_t wt1[VLEN], wt2[VLEN];

  int iv = VLEN * NCD  * site;
  load_vec(vt, &v1[iv], NCD);
  load_vec(wt, &v2[iv], NCD);

  for(int icd = 0; icd < NC * m_Nd; ++icd){
    for(int k = 0; k < VLEN; ++k){
      wt[k+VLEN*icd] *= -m_CKs;
    }
  }

  for(int jd = 0; jd < Nd2; ++jd){
    for(int id = 0; id < m_Nd; ++id){
      int ig = VLEN * NDF2 * (site + m_Nstv * (id + m_Nd * jd));
      load_vec(ut, &u[ig], NDF2);
      for(int ic = 0; ic < NC; ++ic){
        int ic2 = VLEN * NC * ic;
        int id2 = (id + Nd2) % m_Nd;
        mult_uv(wt1, &ut[ic2], &vt[VLEN*NC*id],  NC);
        mult_uv(wt2, &ut[ic2], &vt[VLEN*NC*id2], NC);
        int icd1 = ic + NC * jd;
        int icd2 = ic + NC * (jd + Nd2);
        for(int k = 0; k < VLEN; ++k){
          wt[k+VLEN*icd1] += wt1[k];
          wt[k+VLEN*icd2] += wt2[k];
        }
      }
    }
  }

  save_vec(&v2[iv], wt, NCD);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[0]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[0]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  real_t zL[VLEN*NCD];

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    if(ix == 0){
      mult_wilson_xp1(&buf1[NVC*ND2*iyzt], &v1[VLEN * NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_dn[0]->start();
  chrecv_up[0]->start();
  chsend_dn[0]->wait();
  chrecv_up[0]->wait();
 }
#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;

    clear_vec(v2v, NCD);

    if(ix != m_Nxv-1){
      int nei = ix+1 + m_Nxv * iyzt;
      // if(ix==m_Nxv-1) nei = ix+1 + m_Nxv * iyzt;
      shift_vec2_bw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NC*ND);
      mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
    }else{
      shift_vec0_bw(zL, &v1[VLEN*NCD*site], NCD);
      mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
      mult_wilson_xp2(v2v, &u[VLEN*NDF2*site],
                      &buf2[NVC*ND2*iyzt], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[0]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[0]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  real_t zL[VLEN*NCD];
  real_t uL[VLEN*NDF2];

  if(do_comm[0] > 0){

    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nxv;
      if(ix == m_Nxv-1 && do_comm[0] > 0){
        int iyzt = site/m_Nxv;
        mult_wilson_xm1(&buf1[NVC*ND2*iyzt], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

#pragma omp barrier

#pragma omp master
   {
    chsend_up[0]->start();
    chrecv_dn[0]->start();
    chsend_up[0]->wait();
    chrecv_dn[0]->wait();
   }
#pragma omp barrier

  } // end of if(do_comm[0] > 0)

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;

    clear_vec(v2v, NCD);

    if(ix > 0){
      int nei = ix-1 + m_Nxv * iyzt;
      // if(ix == 0) nei = m_Nxv-1 + m_Nxv * iyzt;
      shift_vec2_fw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
      shift_vec2_fw(uL, &u[VLEN*NDF2*site], &u[VLEN*NDF2*nei], NDF2);
      mult_wilson_xmb(v2v, uL, zL, NC);
    }else{
      shift_vec0_fw(zL, &v1[VLEN*NCD*site], NCD);
      shift_vec0_fw(uL, &u[VLEN*NDF2*site], NDF2);
      mult_wilson_xmb(v2v, uL, zL, NC);
      mult_wilson_xm2(v2v, &buf2[NVC*ND2*iyzt], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy = m_Nxv * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[1]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[1]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nxv;
    int iy  = (site/m_Nxv) % m_Ny;
    int izt = site/Nxy;
    if(iy == 0){
      int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
      mult_wilson_yp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_dn[1]->start();
  chrecv_up[1]->start();

  chsend_dn[1]->wait();
  chrecv_up[1]->wait();
 }
#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nxv;
    int iy  = (site/m_Nxv) % m_Ny;
    int izt = site/Nxy;

    clear_vec(v2v, NCD);

    if(iy != m_Ny-1){
      int nei = ix + m_Nxv * (iy+1 + m_Ny * izt);
      mult_wilson_ypb(v2v, &u[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
      mult_wilson_yp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy = m_Nxv * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[1]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[1]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nxv;
    int iy  = (site/m_Nxv) % m_Ny;
    int izt = site/Nxy;
    if(iy == m_Ny-1){
      int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
      mult_wilson_ym1(&buf1[ibf], &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_up[1]->start();
  chrecv_dn[1]->start();

  chsend_up[1]->wait();
  chrecv_dn[1]->wait();
 }

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nxv;
    int iy  = (site/m_Nxv) % m_Ny;
    int izt = site/Nxy;

    clear_vec(v2v, NCD);

    if(iy != 0){
      int nei = ix + m_Nxv * (iy-1 + m_Ny * izt);
      mult_wilson_ymb(v2v, &u[VLEN*NDF2*nei], &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC * ND2 * (ix + m_Nxv * izt);
      mult_wilson_ym2(v2v, &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy = m_Nxv * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[2]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[2]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);

    if(iz == 0){
      int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
      mult_wilson_zp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_dn[2]->start();
  chrecv_up[2]->start();

  chsend_dn[2]->wait();
  chrecv_up[2]->wait();
 }

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);

    clear_vec(v2v, NCD);

    if(iz != m_Nz-1){
      int nei = ixy + Nxy * (iz+1 + m_Nz*it);
      mult_wilson_zpb(v2v, &u[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
      mult_wilson_zp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy = m_Nxv * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[2]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[2]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);
    if(iz == m_Nz-1){
      int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
      mult_wilson_zm1(&buf1[ibf], &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_up[2]->start();
  chrecv_dn[2]->start();

  chsend_up[2]->wait();
  chrecv_dn[2]->wait();
 }

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);

    clear_vec(v2v, NCD);

    if(iz != 0){
      int nei = ixy + Nxy * (iz-1 + m_Nz*it);
      mult_wilson_zmb(v2v, &u[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC *ND2 * (ixy + Nxy * it);
      mult_wilson_zm2(v2v, &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[3]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[3]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if(it == 0){
      mult_wilson_tp1_dirac(&buf1[VLEN*NC*ND2*ixyz],
                            &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_dn[3]->start();
  chrecv_up[3]->start();

  chsend_dn[3]->wait();
  chrecv_up[3]->wait();
}

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    clear_vec(v2v, NCD);

    if(it != m_Nt-1){
      int nei = ixyz + Nxyz*(it+1);
      mult_wilson_tpb_dirac(v2v, &u[VLEN*NDF2*site],
                            &v1[VLEN*NCD*nei], NC);
    }else{
      mult_wilson_tp2_dirac(v2v, &u[VLEN*NDF2*site],
                            &buf2[VLEN*NC*ND2*ixyz], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nxv * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[3]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[3]->ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if(it == m_Nt-1){
       mult_wilson_tm1_dirac(&buf1[VLEN*NC*ND2*ixyz], 
                     &u[VLEN*NDF2*site], &v1[VLEN*NCD*site], NC);
    }
  }

#pragma omp barrier

#pragma omp master
 {
  chsend_up[3]->start();
  chrecv_dn[3]->start();

  chsend_up[3]->wait();
  chrecv_dn[3]->wait();
 }
#pragma omp barrier

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    clear_vec(v2v, NCD);

    if(it != 0){
      int nei = ixyz + Nxyz*(it-1);
      mult_wilson_tmb_dirac(v2v, &u[VLEN*NDF2*nei],
                            &v1[VLEN*NCD*nei], NC);
    }else{
      mult_wilson_tm2_dirac(v2v, &buf2[VLEN*NC*ND2*ixyz], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

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

  if (m_repr == "Dirac") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd); // <- clover term
  } else if (m_repr == "Chiral") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd);   // <- clover term
  } else {
    //    vout.crucial(m_vl, "%s: input repr is undefined.\n",
    //                 class_name.c_str());
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

  return flop;
}

//============================================================END=====
