/*!
      @file    afopr_Clover_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

template<typename AFIELD>
const std::string AFopr_Clover_eo<AFIELD>::class_name
                                       = "AFopr_Clover_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::init()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  // switches
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  m_vl = AFopr<AFIELD>::m_vl;

  vout.detailed(m_vl, "%s: initalization starts.\n",
		class_name.c_str());

  m_Nc = CommonParameters::Nc();
  if(m_Nc != 3){
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nd = CommonParameters::Nd();
  m_vl = CommonParameters::Vlevel();

  // the follwoing member variables(without m_) are used in
  // implementation of mult.
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
  m_Nx2  = m_Nx/2;
  m_Nst2 = m_Nst/2;
  m_Nx2v  = m_Nx2/VLEN2;
  m_Nst2v = m_Nst2/VLEN2;
  vout.general(m_vl, "  Nx2v = %d  Nst2v = %d\n",
                m_Nx2v, m_Nst2v);

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

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 0;
    if(Communicator::npe(mu) > 1) do_comm[mu] = 1;
    // if(req_comm == 1) do_comm[mu] = 1;
    if(req_comm == 1 && mu > 0) do_comm[mu] = 1;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  // setup of communication buffers
  m_Nbdsize.resize(m_Ndim);
  int Nd2 = m_Nd/2;
  m_Nbdsize[0] = m_Nvc * Nd2 * (m_Ny/2) * m_Nz * m_Nt;
  m_Nbdsize[1] = m_Nvc * Nd2 * m_Nx2 * m_Nz * m_Nt;
  m_Nbdsize[2] = m_Nvc * Nd2 * m_Nx2 * m_Ny * m_Nt;
  m_Nbdsize[3] = m_Nvc * Nd2 * m_Nx2 * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nst, m_Ndim);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst2, 1);
  m_v2.reset(NinF, m_Nst2, 1);

  m_fopr_ct = new AFopr_CloverTerm<AFIELD>;

  m_fee_inv.reset(m_Ndf, m_Nst2, m_Nd*m_Nd);
  m_foo_inv.reset(m_Ndf, m_Nst2, m_Nd*m_Nd);

  vout.detailed(m_vl, "%s: initalization finished.\n",
		class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::setup_channels()
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
void AFopr_Clover_eo<AFIELD>::tidyup()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  delete m_fopr_ct;

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_parameters(const Parameters& params)
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
void AFopr_Clover_eo<AFIELD>::set_parameters(const real_t CKs,
                                               const real_t cSW,
                                               const std::vector<int> bc)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  assert(bc.size() == m_Ndim);

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  CKs  = %8.4f\n", CKs);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  m_CKs = CKs;
  m_cSW = cSW;
  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_config(Field* u)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_timer.reset();
  m_timer.start();

  vout.detailed(m_vl, "%s: set_config start\n", class_name.c_str());

  m_conf = u;

  AFIELD Ulex(m_Ndf, m_Nst, m_Ndim);

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;
  Index_eo_alt<real_t,AFIELD::IMPL> index;

#pragma omp parallel
 {
  convert_gauge(index_lex, Ulex, *u);

#pragma omp barrier

  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1){
      set_boundary_config(Ulex, mu);
    }
  }

  convert(index, m_Ueo, index_lex, Ulex);

 } // omp parallel

  m_fopr_ct->set_config(u);

  int Ncsw = m_Nd * m_Nd/2;

  AFIELD Tinv(m_Ndf, m_Nst, 1);
  AFIELD Tinv_e(m_Ndf, m_Nst, 1), Tinv_o(m_Ndf, m_Nst, 1);

  for(int j = 0; j < Ncsw; ++j){

    m_fopr_ct->get_csw_inv(Tinv, j);

#pragma omp parallel
    {
      index.split(Tinv_e, Tinv_o, Tinv);
#pragma omp barrier

      copy(m_fee_inv, j, Tinv_e, 0);
      copy(m_foo_inv, j, Tinv_o, 0);
     }

  }

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();
  vout.detailed(m_vl, "%s: set_config finished in %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::set_boundary_config(AFIELD& U,
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

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult(AFIELD &v,
                                  const AFIELD &w)
{
  if(m_mode == "D"){
    return D(v, w);
  }else if(m_mode == "DdagD"){
    return DdagD(v, w);
  }else if(m_mode == "Ddag"){
    return Ddag(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_dag(AFIELD &v,
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
void AFopr_Clover_eo<AFIELD>::mult_gm5(AFIELD &v,
                                  const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

#pragma omp barrier

  mult_gm5(vp, wp);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_gm4(AFIELD &v,
                                  const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

#pragma omp barrier

  mult_gm4(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult(AFIELD &v,
                                     const AFIELD &w,
                                     const std::string mode)
{
  if(mode == "Dee_inv"){
    mult_cswinv(v, w, 0);
  }else if(mode == "Doo_inv"){
    mult_cswinv(v, w, 1);
  }else if(mode == "Deo"){
    Meo(v, w, 0);
  }else if(mode == "Doe"){
    Meo(v, w, 1);
  }else if(mode == "convert"){
    mult_gm4(v, w);
    // copy(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::DdagD(AFIELD &v,
                                   const AFIELD &w)
{
  D(m_v2, w);
  Ddag(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Ddag(AFIELD &v,
                                   const AFIELD &w)
{
#pragma omp barrier
  mult_gm5(m_v1, w);
#pragma omp barrier

  mult_cswinv(v, m_v1, 0);
#pragma omp barrier
  Meo(m_v1, v, 1);
#pragma omp barrier
  mult_cswinv(v, m_v1, 1);
#pragma omp barrier
  Meo(m_v1, v, 0);
#pragma omp barrier

  mult_gm5(v, m_v1);
#pragma omp barrier

  aypx(real_t(-1.0), v, w);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::D(AFIELD &v,
                                  const AFIELD &w)
{
#pragma omp barrier
  Meo(m_v1, w, 1);
#pragma omp barrier
  mult_cswinv(v, m_v1, 1);
#pragma omp barrier
  Meo(m_v1, v, 0);
#pragma omp barrier
  mult_cswinv(v, m_v1, 0);
#pragma omp barrier

  aypx(real_t(-1.0), v, w);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::aypx(real_t a,
                                     AFIELD &v,
                                     const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  aypx(a, vp, wp);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Clover_eo<AFIELD>::flop_count()
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == "Dirac") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd);
  } else if (m_repr == "Chiral") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))
      + 4 * m_Nc * m_Nc * m_Nd * m_Nd);
  } else {
    //    vout.crucial(m_vl, "%s: input repr is undefined.\n",
    //                 class_name.c_str());
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    exit(EXIT_FAILURE);
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

  return flop;
}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Meo_alt(AFIELD &v, const AFIELD &w,
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
void AFopr_Clover_eo<AFIELD>::mult_gm5(real_t *v,
                                      real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2);

  Index_eo_alt<real_t,AFIELD::IMPL> index;

  for(int site = is; site < ns; ++site){
    for(int ic = 0; ic < m_Nc; ++ic){
      v[index.idxh_SPr(ic,0,site,0)] = w[index.idxh_SPr(ic,2,site,0)];
      v[index.idxh_SPi(ic,0,site,0)] = w[index.idxh_SPi(ic,2,site,0)];
      v[index.idxh_SPr(ic,1,site,0)] = w[index.idxh_SPr(ic,3,site,0)];
      v[index.idxh_SPi(ic,1,site,0)] = w[index.idxh_SPi(ic,3,site,0)];
      v[index.idxh_SPr(ic,2,site,0)] = w[index.idxh_SPr(ic,0,site,0)];
      v[index.idxh_SPi(ic,2,site,0)] = w[index.idxh_SPi(ic,0,site,0)];
      v[index.idxh_SPr(ic,3,site,0)] = w[index.idxh_SPr(ic,1,site,0)];
      v[index.idxh_SPi(ic,3,site,0)] = w[index.idxh_SPi(ic,1,site,0)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_gm4(real_t *v,
                                         real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2);

  Index_eo_alt<real_t,AFIELD::IMPL> index;

  for(int site = is; site < ns; ++site){
    for(int ic = 0; ic < m_Nc; ++ic){
      v[index.idxh_SPr(ic,0,site,0)] =  w[index.idxh_SPr(ic,0,site,0)];
      v[index.idxh_SPi(ic,0,site,0)] =  w[index.idxh_SPi(ic,0,site,0)];
      v[index.idxh_SPr(ic,1,site,0)] =  w[index.idxh_SPr(ic,1,site,0)];
      v[index.idxh_SPi(ic,1,site,0)] =  w[index.idxh_SPi(ic,1,site,0)];
      v[index.idxh_SPr(ic,2,site,0)] = -w[index.idxh_SPr(ic,2,site,0)];
      v[index.idxh_SPi(ic,2,site,0)] = -w[index.idxh_SPi(ic,2,site,0)];
      v[index.idxh_SPr(ic,3,site,0)] = -w[index.idxh_SPr(ic,3,site,0)];
      v[index.idxh_SPi(ic,3,site,0)] = -w[index.idxh_SPi(ic,3,site,0)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::aypx(real_t a, real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2);

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
void AFopr_Clover_eo<AFIELD>::aypx(real_t a, real_t *v,
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
void AFopr_Clover_eo<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2);

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
void AFopr_Clover_eo<AFIELD>::scal(real_t *v, const real_t a)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2);

  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      v[icd + Ncd * i] *= a;
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::clear(real_t *v,
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
void AFopr_Clover_eo<AFIELD>::scal(real_t *v, const real_t a,
                                   const int is, const int ns)
{
  real_t zero = 0.0;
  int Ncd = m_Nvc * m_Nd;

  for(int i = is; i < ns; ++i){
    for(int icd = 0; icd < Ncd; ++icd){
      v[icd + Ncd * i] *= a;
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_cswinv(AFIELD& v,
                                      const AFIELD& w, int ieo)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

#pragma omp barrier

  real_t *csw_inv = 0;
  if(ieo == 0){
    csw_inv = m_fee_inv.ptr(0);
  }else if(ieo == 1){
    csw_inv = m_foo_inv.ptr(0);
  }

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t vt[NCD], wt[NCD];
  Vsimd_t ut[NDF2];
  Vsimd_t wt1, wt2;

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NCD  * site;

    load_vec(vt, &v1[iv], NCD);
    clear_vec(wt, NCD);

    for(int jd = 0; jd < ND2; ++jd){
      for(int id = 0; id < ND; ++id){
        int id2 = (id + ND2) % m_Nd;
        int ig  = VLEN * NDF2 * (site + m_Nst2v * (id + ND * jd));
        load_vec(ut, &csw_inv[ig], NDF2);
        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          mult_uv(wt1, &ut[ic2], &vt[NC*id],  m_Nc);
          mult_uv(wt2, &ut[ic2], &vt[NC*id2], m_Nc);
          int icd1 = ic + m_Nc * jd;
          int icd2 = ic + m_Nc * (jd + ND2);
          axpy_vec(&wt[icd1], 1.0, &wt1, 1);
          axpy_vec(&wt[icd2], 1.0, &wt2, 1);
        }
      }
    }

    save_vec(&v2[iv], wt, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::Meo(AFIELD &v,
                                  const AFIELD &w, const int ieo)
{
#pragma omp barrier

  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

  int Nxy2 = m_Nx2v * m_Ny;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t zL[VLEN*NCD];
  real_t uL[VLEN*NDF2];

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    if(do_comm[0] == 1){
      if(ix == 0 && Leo == 1){
        int in = VLEN * NCD * (ix + m_Nx2v * iyzt);
        real_t *buf1 = (real_t*)chsend_dn[0]->ptr();
        int ibf_up = iyzt/2;
        mult_wilson_xp1(&buf1[NVC*ND2*ibf_up], &v1[in], NC);
      }
      if(ix == m_Nx2v-1 && Leo == 0){ 
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2*0));
        real_t *buf1 = (real_t*)chsend_up[0]->ptr();
        int ibf_dn = iyzt/2;
        mult_wilson_xm1(&buf1[NVC*ND2*ibf_dn], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == 0){
        real_t *buf1 = (real_t*)chsend_dn[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
        mult_wilson_yp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iy == m_Ny-1){
        real_t *buf1 = (real_t*)chsend_up[1]->ptr();
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * 1));
        int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
        mult_wilson_ym1(&buf1[ibf], &u[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == 0){
        real_t *buf1 = (real_t*)chsend_dn[2]->ptr();
        int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
        mult_wilson_zp1(&buf1[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iz == m_Nz-1){
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * 2));
        real_t *buf1 = (real_t*)chsend_up[2]->ptr();
        int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
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
        real_t *u = m_Ueo.ptr(NDF * m_Nst2 *(1-ieo + 2 * 3));
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

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    clear_vec(v2v, NCD);

    if(Leo == 0){
      real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 0));
      mult_wilson_xpb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }else{
      if(ix < m_Nx2v-1 || do_comm[0] == 0){
        int nei = ix+1 + m_Nx2v * iyzt;
	if(ix == m_Nx2v-1) nei = 0 + m_Nx2v * iyzt;
        real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 0));
        shift_vec2_bw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NC*ND);
        mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
      }
    }

    if(Leo == 1){
      real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2*0));
      mult_wilson_xmb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }else{
      if(ix > 0 || do_comm[0] == 0){
        int nei = ix-1 + m_Nx2v * iyzt;
	if(ix == 0) nei = m_Nx2v-1 + m_Nx2v * iyzt;
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2*0));
        shift_vec2_fw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
        shift_vec2_fw(uL, &u[VLEN*NDF2*site], &u[VLEN*NDF2*nei], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
      }
    }

    if(iy != m_Ny-1 || do_comm[1] == 0){
      int iy2 = (iy + 1) % m_Ny;
      int nei = ix + m_Nx2v * (iy2 + m_Ny * izt);
      real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 1));
      mult_wilson_ypb(v2v, &u[VLEN*NDF2*site],
                    &v1[VLEN*NCD*nei], NC);
    }

    if(iy != 0 || do_comm[1] == 0){
      int iy2 = (iy - 1 + m_Ny) % m_Ny;
      int nei = ix + m_Nx2v * (iy2 + m_Ny * izt);
      real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * 1));
      mult_wilson_ymb(v2v, &u[VLEN*NDF2*nei],
                           &v1[VLEN*NCD*nei], NC);
    }

    if(iz != m_Nz-1 || do_comm[2] == 0){
      real_t* u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * 2));
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxy2 * (iz2 + m_Nz*it);
      mult_wilson_zpb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*nei], NC);
    }

    if(iz != 0 || do_comm[2] == 0){
      real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * 2));
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy2 * (iz2 + m_Nz*it);
      mult_wilson_zmb(v2v, &u[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }

    if(it != m_Nt-1 || do_comm[3] == 0){
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyz2 * it2;
      real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * 3));
      mult_wilson_tpb_dirac(v2v, &u[VLEN*NDF2*site],
                            &v1[VLEN*NCD*nei], NC);
    }

    if(it != 0 || do_comm[3] == 0){
      real_t *u = m_Ueo.ptr(NDF * m_Nst2 *(1-ieo + 2 * 3));
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz2 * it2;
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

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy2;
    int iz   = izt % m_Nz;
    int it   = izt/m_Nz;
    int ixy = ix + m_Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    load_vec(v2v, &v2[VLEN*NCD*site], NCD);

    if(do_comm[0] == 1){
      if(Leo == 1 && ix == m_Nx2v-1){
        real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 0));
        real_t *buf2 = (real_t*)chrecv_up[0]->ptr();
        int ibf_up = iyzt/2;
        shift_vec0_bw(zL, &v1[VLEN*NCD*site], NCD);
        mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
        mult_wilson_xp2(v2v, &u[VLEN*NDF2*site],
                        &buf2[NVC*ND2*ibf_up], NC);
      }
      if(Leo == 0 && ix == 0){
        real_t *buf2 = (real_t*)chrecv_dn[0]->ptr();
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2*0));
        int ibf_dn = iyzt/2;
        shift_vec0_fw(zL, &v1[VLEN*NCD*site], NCD);
        shift_vec0_fw(uL, &u[VLEN*NDF2*site], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
        mult_wilson_xm2(v2v, &buf2[NVC*ND2*ibf_dn], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == m_Ny-1){
        real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 1));
        real_t *buf2 = (real_t*)chrecv_up[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
        mult_wilson_yp2(v2v, &u[VLEN*NDF2*site],
                        &buf2[ibf], NC);
      }
      if(iy == 0){
        real_t *buf2 = (real_t*)chrecv_dn[1]->ptr();
        int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
        mult_wilson_ym2(v2v, &buf2[ibf], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == m_Nz-1){
        real_t *buf2 = (real_t*)chrecv_up[2]->ptr();
        real_t* u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * 2));
        int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
        mult_wilson_zp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
      }
      if(iz == 0){
        real_t *buf2 = (real_t*)chrecv_dn[2]->ptr();
        int ibf = VLEN * NC *ND2 * (ixy + Nxy2 * it);
        mult_wilson_zm2(v2v, &buf2[ibf], NC);
      }
    }

    if(do_comm[3] == 1){
      if(it == m_Nt-1){
        real_t *buf2 = (real_t*)chrecv_up[3]->ptr();
        real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * 3));
        mult_wilson_tp2_dirac(v2v, &u[VLEN*NDF2*site],
                              &buf2[VLEN*NC*ND2*ixyz], NC);
      }
      if(it == 0){
        real_t *buf2 = (real_t*)chrecv_dn[3]->ptr();
        mult_wilson_tm2_dirac(v2v, &buf2[VLEN*NC*ND2*ixyz], NC);
      }
    }

    scal_vec(v2v, -m_CKs, NCD);
    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

 }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_xp(real_t *v2, real_t *v1,
                                        const int ieo)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[0]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[0]->ptr();
  real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

  real_t zL[VLEN*NCD];

  if(do_comm[0] == 1){
    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nx2v;
      int iyzt = site/m_Nx2v;
      int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];
      if(ix == 0 && Leo == 1){
        int in = VLEN * NCD * (ix + m_Nx2v * iyzt);
        int ibf_up = iyzt/2;
        mult_wilson_xp1(&buf1[NVC*ND2*ibf_up], &v1[in], NC);
      }
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
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt];

    clear_vec(v2v, NCD);

    if(Leo == 0){
      mult_wilson_xpb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }else{
      if(ix < m_Nx2v-1 || do_comm[0] == 0){
        int nei = ix+1 + m_Nx2v * iyzt;
        if(ix == m_Nx2v-1) nei = 0 + m_Nx2v * iyzt;
        real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * 0));
        shift_vec2_bw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NC*ND);
        mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
      }else{
        int ibf_up = iyzt/2;
        shift_vec0_bw(zL, &v1[VLEN*NCD*site], NCD);
        mult_wilson_xpb(v2v, &u[VLEN*NDF2*site], zL, NC);
        mult_wilson_xp2(v2v, &u[VLEN*NDF2*site],
                        &buf2[NVC*ND2*ibf_up], NC);
      }
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_xm(real_t *v2, real_t *v1,
                                        const int ieo)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[0]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[0]->ptr();
  real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2*idir));

  real_t zL[VLEN*NCD];
  real_t uL[VLEN*NDF2];

  real_t yt1[NVC], yt2[NVC];

  if(do_comm[0] == 1){
    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nx2v;
      int iyzt = site/m_Nx2v;
      int Leo = ieo + (1 - 2*ieo) * m_Leo[iyzt];
      if(Leo == 0){ 
        if(ix == m_Nx2v-1){
          int ibf_dn = iyzt/2;
          mult_wilson_xm1(&buf1[NVC*ND2*ibf_dn], &u[VLEN*NDF2*site],
                          &v1[VLEN*NCD*site], NC);
        }
      }
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

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx2v;
    int iyzt = site/m_Nx2v;
    int Leo = ieo + (1 - 2*ieo) * m_Leo[iyzt];

    clear_vec(v2v, NCD);

    if(Leo == 1){
      mult_wilson_xmb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*site], NC);
    }else{
      if(ix > 0 || do_comm[0] == 0){
        int nei = ix-1 + m_Nx2v * iyzt;
        if(ix == 0) nei = m_Nx2v-1 + m_Nx2v * iyzt;
        shift_vec2_fw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
        shift_vec2_fw(uL, &u[VLEN*NDF2*site], &u[VLEN*NDF2*nei], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
      }else{
        int ibf_dn = iyzt/2;
        shift_vec0_fw(zL, &v1[VLEN*NCD*site], NCD);
        shift_vec0_fw(uL, &u[VLEN*NDF2*site], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
        mult_wilson_xm2(v2v, &buf2[NVC*ND2*ibf_dn], NC);
      }
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_yp(real_t *v2, real_t *v1,
				      const int ieo)
{
  int idir = 1;
  int Nxy2 = m_Nx2v * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[1]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[1]->ptr();
  real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nx2v;
    int iy  = (site/m_Nx2v) % m_Ny;
    int izt = site/Nxy2;
    if(iy == 0){
      int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
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
    int ix  = site % m_Nx2v;
    int iy  = (site/m_Nx2v) % m_Ny;
    int izt = site/Nxy2;

    clear_vec(v2v, NCD);

    if(iy != m_Ny-1){
      int nei = ix + m_Nx2v * (iy+1 + m_Ny * izt);
      mult_wilson_ypb(v2v, &u[VLEN*NDF2*site],
                    &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
      mult_wilson_yp2(v2v, &u[VLEN*NDF2*site],
		      &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_ym(real_t *v2, real_t *v1,
				      const int ieo)
{
  int idir = 1;
  int Nxy2 = m_Nx2v * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[1]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[1]->ptr();
  real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nx2v;
    int iy  = (site/m_Nx2v) % m_Ny;
    int izt = site/Nxy2;
    if(iy == m_Ny-1){
      int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
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
    int ix  = site % m_Nx2v;
    int iy  = (site/m_Nx2v) % m_Ny;
    int izt = site/Nxy2;

    clear_vec(v2v, NCD);

    if(iy != 0){
      int nei = ix + m_Nx2v * (iy-1 + m_Ny * izt);
      mult_wilson_ymb(v2v, &u[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }else{
      int ibf = VLEN * NC * ND2 * (ix + m_Nx2v * izt);
      mult_wilson_ym2(v2v, &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_zp(real_t *v2, real_t *v1,
                                        const int ieo)
{
  int idir = 2;
  int Nxy2 = m_Nx2v * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[2]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[2]->ptr();
  real_t* u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy2;
    int iz  = (site/Nxy2) % m_Nz;
    int it  = site/(Nxy2 * m_Nz);
    if(iz == 0){
      int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
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
    int ixy = site % Nxy2;
    int iz  = (site/Nxy2) % m_Nz;
    int it  = site/(Nxy2 * m_Nz);

    clear_vec(v2v, NCD);

    if(iz != m_Nz-1){
      int nei = ixy + Nxy2 * (iz+1 + m_Nz*it);
      mult_wilson_zpb(v2v, &u[VLEN*NDF2*site],
                      &v1[VLEN*NCD*nei], NC);
    }
    if(iz == m_Nz-1){
      int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
      mult_wilson_zp2(v2v, &u[VLEN*NDF2*site], &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_zm(real_t *v2, real_t *v1,
				      const int ieo)
{
  int idir = 2;
  int Nxy2 = m_Nx2v * m_Ny;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[2]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[2]->ptr();
  real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (1-ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy2;
    int iz  = (site/Nxy2) % m_Nz;
    int it  = site/(Nxy2 * m_Nz);
    int ibf = ixy + Nxy2 * it;
    if(iz == m_Nz-1){
      int ibf = VLEN * NC * ND2 * (ixy + Nxy2 * it);
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
    int ixy = site % Nxy2;
    int iz  = (site/Nxy2) % m_Nz;
    int it  = site/(Nxy2 * m_Nz);

    clear_vec(v2v, NCD);

    if(iz != 0){
      int nei = ixy + Nxy2 * (iz-1 + m_Nz*it);
      mult_wilson_zmb(v2v, &u[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }
    if(iz == 0){
      int ibf = VLEN * NC *ND2 * (ixy + Nxy2 * it);
      mult_wilson_zm2(v2v, &buf2[ibf], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_tp(real_t *v2, real_t *v1,
				      const int ieo)
{
  int idir = 3;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_dn[3]->ptr();
  real_t *buf2 = (real_t*)chrecv_up[3]->ptr();
  real_t *u = m_Ueo.ptr(m_Ndf * m_Nst2 * (ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz2;
    int it   = site / Nxyz2;
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
    int ixyz = site % Nxyz2;
    int it   = site / Nxyz2;

    clear_vec(v2v, NCD);

    if(it != m_Nt-1){
      int nei = ixyz + Nxyz2*(it+1);
      mult_wilson_tpb_dirac(v2v, &u[VLEN*NDF2*site],
                            &v1[VLEN*NCD*nei], NC);
    }
    if(it == m_Nt-1){
      mult_wilson_tp2_dirac(v2v, &u[VLEN*NDF2*site],
                            &buf2[VLEN*NC*ND2*ixyz], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Clover_eo<AFIELD>::mult_tm(real_t *v2, real_t *v1,
				      const int ieo)
{
  int idir = 3;
  int Nxyz2 = m_Nx2v * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

  Vsimd_t v2v[NCD];

  real_t *buf1 = (real_t*)chsend_up[3]->ptr();
  real_t *buf2 = (real_t*)chrecv_dn[3]->ptr();
  real_t *u = m_Ueo.ptr(NDF * m_Nst2 *(1-ieo + 2 * idir));

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz2;
    int it   = site / Nxyz2;
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
    int ixyz = site % Nxyz2;
    int it   = site / Nxyz2;

    clear_vec(v2v, NCD);

    if(it != 0){
      int nei = ixyz + Nxyz2*(it-1);
      mult_wilson_tmb_dirac(v2v, &u[VLEN*NDF2*nei],
                                 &v1[VLEN*NCD*nei], NC);
    }
    if(it == 0){
      mult_wilson_tm2_dirac(v2v, &buf2[VLEN*NC*ND2*ixyz], NC);
    }

    add_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//============================================================END=====
