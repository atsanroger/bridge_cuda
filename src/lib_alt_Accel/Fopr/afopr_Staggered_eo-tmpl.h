/*!
      @file    afopr_Staggered_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2554 $
*/

#include "lib_alt_Accel/Fopr/afopr_Staggered_eo.h"

template<typename AFIELD>
const std::string AFopr_Staggered_eo<AFIELD>::class_name
                                      = "AFopr_Staggered_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches
  //  int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication only when necessary

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nc = CommonParameters::Nc();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nvc = m_Nc * 2;
  int Ndf = 2 * m_Nc * m_Nc;
  int Nx  = CommonParameters::Nx();
  int Ny  = CommonParameters::Ny();
  int Nz  = CommonParameters::Nz();
  int Nt  = CommonParameters::Nt();
  int Nx2 = Nx/2;
  m_Nst2  = m_Nst/2;

  m_Nsize[0] = Nx2;
  m_Nsize[1] = Ny;
  m_Nsize[2] = Nz;
  m_Nsize[3] = Nt;

  // Nc check: temporary commented out.
  // check_Nc(m_Nc);
  // check_setup();

  if(Nx % 2 != 0){
    vout.crucial(m_vl, "%s: Nx must be even.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ipe1 = Communicator::ipe(1);
  int ipe2 = Communicator::ipe(2);
  int ipe3 = Communicator::ipe(3);
  m_Ieo_origin = (ipe1*Ny + ipe2*Nz + ipe3*Nt) % 2;

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general(m_vl, "  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  m_Nbdsize[0] = m_Nvc * ceil_nwp((Ny * Nz * Nt + 1)/2);
  m_Nbdsize[1] = m_Nvc * ceil_nwp(Nx2 * Nz * Nt);
  m_Nbdsize[2] = m_Nvc * ceil_nwp(Nx2 * Ny * Nt);
  m_Nbdsize[3] = m_Nvc * ceil_nwp(Nx2 * Ny * Nz);

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(Ndf, m_Nst2, 2*m_Ndim);

  // working vectors.
  m_v1.reset(m_Nvc, m_Nst2, 1);
  m_v2.reset(m_Nvc, m_Nst2, 1);

  m_staggered_phase.reset(1, m_Nst, m_Ndim);
  m_parity.reset(1, m_Nst, 1);

  set_staggered_phase();

  vout.increase_indent();

  set_parameters(params);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

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
void AFopr_Staggered_eo<AFIELD>::tidyup()
{
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
void AFopr_Staggered_eo<AFIELD>::set_parameters(const Parameters& params)
{

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  double mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(real_t(mq), bc);

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_parameters(const real_t mq,
                                        const std::vector<int> bc)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_mq = mq;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
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

  //- print input parameters
  vout.general(m_vl, "%s: set parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", double(m_mq));
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_staggered_phase()
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {

    int Nt   = CommonParameters::Nt();
    int Nz   = CommonParameters::Nz();
    int Ny   = CommonParameters::Ny();
    int Nx   = CommonParameters::Nx();
    int Nvol = CommonParameters::Nvol();

    int ipex = Communicator::ipe(0);
    int ipey = Communicator::ipe(1);
    int ipez = Communicator::ipe(2);
    int ipet = Communicator::ipe(3);

    real_t* stgph = m_staggered_phase.ptr(0);
    real_t* prty  = m_parity.ptr(0);

    AIndex_lex<real_t,AFIELD::IMPL> index;

    for (int t = 0; t < Nt; ++t) {
      int t2 = t + ipet * Nt;
      for (int z = 0; z < Nz; ++z) {
        int z2 = z + ipez * Nz;
        for (int y = 0; y < Ny; ++y) {
          int y2 = y + ipey * Ny;
          for (int x = 0; x < Nx; ++x) {
            int x2 = x + ipex * Nx;
            int is = index.site(x, y, z, t);

            stgph[index.idx(0, 1, is, 0)] = 1.0;
            stgph[index.idx(0, 1, is, 1)] = 1.0;
            stgph[index.idx(0, 1, is, 2)] = 1.0;
            stgph[index.idx(0, 1, is, 3)] = 1.0;
            prty[ index.idx(0, 1, is, 0)] = 1.0;

            if (( x2 % 2) == 1){
              stgph[index.idx(0, 1, is, 1)] = -1.0;
            }
            if (((x2 + y2) % 2) == 1){
              stgph[index.idx(0, 1, is, 2)] = -1.0;
            }
            if (((x2 + y2 + z2) % 2) == 1){
              stgph[index.idx(0, 1, is, 3)] = -1.0;
            }
            if (((x2 + y2 + z2 + t2) % 2) == 1){
              prty[index.idx(0, 1, is, 0)] = -1.0;
            }

          }
        }
      }
    }

    update_device(m_staggered_phase);
    update_device(m_parity);
  }

#pragma omp barrier

}

//====================================================================
/*
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_staggered_phase(AFIELD& v, int mu)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {

    int Nvol  = CommonParameters::Nvol();
    int Nin = v.nin();
    int Nex = v.nex();

    assert(v.nvol() == Nvol);

    real_t* stgph = m_staggered_phase.ptr(0);
    real_t* vp = v.ptr(0);

    AIndex_lex<real_t,AFIELD::IMPL> index;

    BridgeACC::copy_from_device(vp, Nin*Nvol*Nex);

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = 0; site < Nvol; ++site) {
        real_t sph = stgph[index.idx(0, 1, site, mu)];
        for (int in = 0; in < Nin; ++in) {
          real_t vt = vp[index.idx(in, Nin, site, ex)];
          vp[index.idx(in, Nin, site, ex)] = sph * vt;
        }
      }
    }

    BridgeACC::copy_to_device(vp, Nin*Nvol*Nex);
  }

#pragma omp barrier

}
*/
//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_config(Field *U)
{
  ThreadManager::assert_single_thread(class_name);

  int Ndf = m_Nvc * m_Nc;
  AFIELD Ulex(Ndf, m_Nst, m_Ndim);

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  convert(index_lex, Ulex, *U);

  int Nsize2[NDIM];
  Nsize2[0] = m_Nsize[0] * 2;
  for(int mu = 1; mu < NDIM; ++mu) Nsize2[mu] = m_Nsize[mu];

  real_t* stgph = m_staggered_phase.ptr(0);
  real_t* up = Ulex.ptr(0);

  BridgeACC::mult_staggered_phase_dev(up, stgph, Nsize2, NC);

  AIndex_eo<real_t,AFIELD::IMPL> index_eo;
  index_eo.split_gauge(m_Ueo, Ulex);

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if(m_mode == "D") {
    D(v, w);
  }else if(m_mode == "Ddag"){
    Ddag(v, w);
  }else if(m_mode == "DdagD"){
    DdagD(v, w);
    //  }else if(m_mode == "H"){
    //    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if(m_mode == "D"){
    Ddag(v, w);
  }else if(m_mode == "Ddag"){
    D(v, w);
  }else if(m_mode == "DdagD"){
    DdagD(v, w);
    //  } else if(m_mode == "H"){
    //    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                      const std::string mode)
{
  if(mode == "Dee"){
    copy(v, w);
    scal(v, m_mq);
  }else if(mode == "Doo"){
    copy(v, w);
    scal(v, m_mq);
  }else if(mode == "Dee_inv"){
    copy(v, w);
    scal(v, real_t(1.0/m_mq));
  }else if(mode == "Doo_inv"){
    copy(v, w);
    scal(v, real_t(1.0/m_mq));
  }else if(mode == "Deo"){
    Meo(v, w, 0);
    //Meo_alt(v, w, 0);
  }else if(mode == "Doe"){
    Meo(v, w, 1);
    //Meo_alt(v, w, 1);
  }else if(mode == "convert"){
    copy(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  Ddag(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  Meo(m_v1, w, 1);
  Meo(v, m_v1, w, 0, 1);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  //  is the following implementation right ?
  D(v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Meo_alt(AFIELD &v, const AFIELD &w,
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
  // scal(vp, 0.5/m_mq);  // hopping normalization
  scal(vp, 0.5); // mass normalizaton

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Meo(AFIELD &v, const AFIELD& w,
                                     const int ieo)
{
  Meo(v, w, w, ieo, 0);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::Meo(AFIELD &v, const AFIELD& w,
                                     const AFIELD& x,
                                     const int ieo, const int iflag)
{  // iflag = 0: v = Meo * w,  others: v = mq2 * x - Meo * w.

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
    real_t *xp = const_cast<AFIELD*>(&x)->ptr(0);
    real_t *up = m_Ueo.ptr(0);
    int jeo = (m_Ieo_origin + ieo) % 2;

    if(do_comm_any > 0){

      real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
      real_t *buf1xm = (real_t*)chsend_up[0].ptr();
      real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
      real_t *buf1ym = (real_t*)chsend_up[1].ptr();
      real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
      real_t *buf1zm = (real_t*)chsend_up[2].ptr();
      real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
      real_t *buf1tm = (real_t*)chsend_up[3].ptr();

      BridgeACC::mult_staggered_1eo(buf1xp, buf1xm, buf1yp, buf1ym,
                                    buf1zp, buf1zm, buf1tp, buf1tm,
                                    up, wp, m_bc, m_Nsize, do_comm,
                                    ieo, jeo);

      chset_recv.start();
      chset_send.start();
    }

    real_t mq2 = m_mq * m_mq;
    BridgeACC::mult_staggered_Meo(vp, up, up, wp, xp,
                                  mq2, m_bc2, m_Nsize, ieo, jeo, iflag);

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

      BridgeACC::mult_staggered_2eo(vp, up,
                                    buf2xp, buf2xm, buf2yp, buf2ym,
                                    buf2zp, buf2zm, buf2tp, buf2tm,
                                    m_mq, m_bc, m_Nsize, do_comm,
                                    ieo, jeo, iflag);
    }

  } // if(ith == 0)

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_gm5 is not defined.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_gm5(AFIELD& v)
{
  vout.crucial(m_vl, "%s: mult_gm5 is not defined.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::scal(AFIELD& v, real_t a)
{
  real_t *vp = v.ptr(0);
  scal(vp, a);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::clear(real_t *v2)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){
    BridgeACC::mult_staggered_clear(v2, m_Nsize, NC);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::aypx(real_t a, real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){
    BridgeACC::mult_staggered_aypx(a, v2, v1, m_Nsize, NC);
  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::scal(real_t *v, real_t a)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){
    BridgeACC::mult_staggered_scal(v, a, m_Nsize, NC);
  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::normalize_fprop(AFIELD& v)
{
  // D is defined with mass normalization:
  // even-odd solver assumes diagonal element is unity, thus
  // solution is to be multiplied by diagonal term
  #pragma omp barrier
  scal(v, m_mq * m_mq);
  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::normalize_fopr(AFIELD& v)
{
  // this function is define the inverse of normalize_fprop().
  #pragma omp barrier
  scal(v, real_t(1.0) / (m_mq * m_mq));
  #pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_xp(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 0;
    int jeo = (m_Ieo_origin + ieo) % 2;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

    if(do_comm[idir] > 0){

      BridgeACC::mult_staggered_xp1_eo(buf1, v1, m_Nsize, m_bc, jeo, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_staggered_xpb_eo(v2, u, v1, m_Nsize, m_bc2, jeo, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_xp2_eo(v2, u, buf2, m_Nsize, m_bc, jeo, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_xm(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 0;
    int jeo = (m_Ieo_origin + ieo) % 2;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_xm1_eo(buf1, u, v1, m_Nsize, m_bc, jeo, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_staggered_xmb_eo(v2, u, v1, m_Nsize, m_bc2, jeo, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_xm2_eo(v2, buf2, m_Nsize, m_bc, jeo, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_yp(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 1;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_yp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_staggered_ypb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_yp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_ym(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 1;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_ym1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_staggered_ymb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_ym2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_zp(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 2;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_zp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_staggered_zpb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_zp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_zm(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 2;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_zm1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_staggered_zmb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_zm2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_tp(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 3;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_tp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_staggered_tpb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_tp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered_eo<AFIELD>::mult_tm(real_t *v2, real_t *v1,
                                         const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 3;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_Ueo.ptr(NDF * m_Nst2 * (1-ieo + 2 * idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_tm1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_staggered_tmb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){ 
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_tm2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
double AFopr_Staggered_eo<AFIELD>::flop_count()
{
  return flop_count(m_mode);
}

//====================================================================
template<typename AFIELD>
double AFopr_Staggered_eo<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The following is based on rev.1976.   [21 Jul 2019 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  flop_site = static_cast<double>(m_Nvc * (2  +  8 * 2 * m_Nvc));
                              //  #comp   aypx  dir FMA  #comp

  flop = flop_site * static_cast<double>(Lvol);

  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;

}

//============================================================END=====
