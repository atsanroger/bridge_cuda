/*!
      @file    afopr_Staggered-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#include "lib_alt_Accel/Fopr/afopr_Staggered.h"

template<typename AFIELD>
const std::string AFopr_Staggered<AFIELD>::class_name
                                         = "AFopr_Staggered<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  // switches
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication only when necessary

  vout.general(m_vl, "%s: being setup:\n", class_name.c_str());

  m_vl = CommonParameters::Vlevel();

  m_Nc = CommonParameters::Nc();
  m_Nst  = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nvc = m_Nc * 2;
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
  m_Nbdsize[0] = m_Nvc * Ny * Nz * Nt;
  m_Nbdsize[1] = m_Nvc * Nx * Nz * Nt;
  m_Nbdsize[2] = m_Nvc * Nx * Ny * Nt;
  m_Nbdsize[3] = m_Nvc * Nx * Ny * Nz;

  setup_channels();

  // gauge configuration.
  m_U.reset(Ndf, m_Nst, m_Ndim);

  // working vectors.
  m_v2.reset(m_Nvc, m_Nst, 1);

  m_staggered_phase.reset(1, m_Nst, m_Ndim);
  m_parity.reset(1, m_Nst, 1);

  set_staggered_phase();

  set_parameters(params);

  //vout.detailed(m_vl, "%s: initalization finished.\n",
  vout.general(m_vl, "%s: initalization finished.\n",
		class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::setup_channels()
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
void AFopr_Staggered<AFIELD>::tidyup()
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
void AFopr_Staggered<AFIELD>::set_parameters(const Parameters& params)
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
void AFopr_Staggered<AFIELD>::set_parameters(const real_t mq,
                                        const std::vector<int> bc)
{

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == m_Ndim);

  //- store values
  m_mq = mq;
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
      m_bc[mu]  = 0;    // for boundary part (dummy)
      m_bc2[mu] = m_boundary[mu];  // for bulk part
    }
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }


}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", double(m_mq));
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::set_staggered_phase()
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

    BridgeACC::copy_to_device(stgph, m_Nst*m_Ndim);
    BridgeACC::copy_to_device(prty,  m_Nst);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_staggered_phase(AFIELD& v, int mu)
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

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::set_config(Field *U)
{
  ThreadManager::assert_single_thread(class_name);

  AIndex_lex<real_t,AFIELD::IMPL> index;

  convert(index, m_U, *U);

  real_t* stgph = m_staggered_phase.ptr(0);
  real_t* up = m_U.ptr(0);

  BridgeACC::mult_staggered_phase_dev(up, stgph, m_Nsize, NC);

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undeifined.\n", class_name.c_str());
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(v, w);
  mult_gm5(v);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  Ddag(v, m_v2);
  //H(m_v2, w);
  //H(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::D_alt(AFIELD& v, const AFIELD& w)
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

  // hopping normalization
  // aypx(0.5/m_mq, vp, wp);

  // mass normalization
  scal(vp, 0.5);
  axpy(vp, m_mq, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::Ddag_alt(AFIELD& v, const AFIELD& w)
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

  // hopping normalization
  // aypx(-0.5/m_mq, vp, wp);

  // mass normalization
  scal(vp, -0.5);
  axpy(vp, m_mq, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
    real_t *up = m_U.ptr(0);

    int jd = 1;

    if(do_comm_any > 0){
      real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
      real_t *buf1xm = (real_t*)chsend_up[0].ptr();

      real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
      real_t *buf1ym = (real_t*)chsend_up[1].ptr();

      real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
      real_t *buf1zm = (real_t*)chsend_up[2].ptr();

      real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
      real_t *buf1tm = (real_t*)chsend_up[3].ptr();

      BridgeACC::mult_staggered_1(buf1xp, buf1xm, buf1yp, buf1ym,
                                  buf1zp, buf1zm, buf1tp, buf1tm,
                                  up, wp, m_bc, m_Nsize, do_comm);
 
      chset_send.start();
      chset_recv.start();
    }

    BridgeACC::mult_staggered_D(vp, up, up, wp, m_mq, m_bc2, m_Nsize, jd);

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

      BridgeACC::mult_staggered_2(vp, up,
                                  buf2xp, buf2xm, buf2yp, buf2ym,
                                  buf2zp, buf2zm, buf2tp, buf2tm,
                                  m_mq, m_bc, m_Nsize, do_comm, jd);
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
    real_t *up = m_U.ptr(0);

    int jd = -1;

    if(do_comm_any > 0){

      real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
      real_t *buf1xm = (real_t*)chsend_up[0].ptr();

      real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
      real_t *buf1ym = (real_t*)chsend_up[1].ptr();

      real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
      real_t *buf1zm = (real_t*)chsend_up[2].ptr();

      real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
      real_t *buf1tm = (real_t*)chsend_up[3].ptr();

      BridgeACC::mult_staggered_1(buf1xp, buf1xm, buf1yp, buf1ym,
                                  buf1zp, buf1zm, buf1tp, buf1tm,
                                  up, wp, m_bc, m_Nsize, do_comm);

      chset_send.start();
      chset_recv.start();
    }

    BridgeACC::mult_staggered_D(vp, up, up, wp, m_mq, m_bc2, m_Nsize, jd);

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

      BridgeACC::mult_staggered_2(vp, up,
                                  buf2xp, buf2xm, buf2yp, buf2ym,
                                  buf2zp, buf2zm, buf2tp, buf2tm,
                                  m_mq, m_bc, m_Nsize, do_comm, jd);

    }
  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int Nin  = w.nin();
    int Nvol = w.nvol();
    int Nex  = w.nex();

    AIndex_lex<real_t,AFIELD::IMPL> index;

    real_t* prty = m_parity.ptr(0);
    real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);
    real_t* vp = v.ptr(0);

    BridgeACC::copy_from_device(wp, Nin*Nvol*Nex);

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = 0; site < Nvol; ++site) {
        real_t ph = prty[index.idx(0, 1, site, 0)];
        for (int in = 0; in < Nin; ++in) {
          real_t wt = wp[index.idx(in, Nin, site, ex)];
          vp[index.idx(in, Nin, site, ex)] = ph * wt;
        }
      }
    }

    BridgeACC::copy_to_device(vp, Nin*Nvol*Nex);

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_gm5(AFIELD& v)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int Nin  = v.nin();
    int Nvol = v.nvol();
    int Nex  = v.nex();

    AIndex_lex<real_t,AFIELD::IMPL> index;

    real_t* prty = m_parity.ptr(0);
    real_t* vp = v.ptr(0);

    BridgeACC::copy_from_device(vp, Nin*Nvol*Nex);

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = 0; site < Nvol; ++site) {
        real_t ph = prty[index.idx(0, 1, site, 0)];
        for (int in = 0; in < Nin; ++in) {
          real_t vt = vp[index.idx(in, Nin, site, ex)];
          vp[index.idx(in, Nin, site, ex)] = ph * vt;
        }
      }
    }

    BridgeACC::copy_to_device(vp, Nin*Nvol*Nex);

  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::scal(AFIELD& v, real_t a)
{
  real_t *vp = v.ptr(0);
  scal(vp, a);
}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::clear(real_t *v2)
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
void AFopr_Staggered<AFIELD>::aypx(real_t a, real_t *v2, real_t *v1)
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
void AFopr_Staggered<AFIELD>::axpy(real_t *v2, real_t a, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){
    BridgeACC::mult_staggered_axpy(v2, a, v1, m_Nsize, NC);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::scal(real_t *v, real_t a)
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
void AFopr_Staggered<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 0;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_xp1(buf1, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_up[idir].start();
      chsend_dn[idir].start();
    }

    BridgeACC::mult_staggered_xpb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_dn[idir].wait();
      chrecv_up[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_xp2(v2, u, buf2, m_Nsize, m_bc, NC);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 0;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

    if(do_comm[idir] > 0){
      BridgeACC::mult_staggered_xm1(buf1, u, v1, m_Nsize, m_bc, NC);
      BridgeACC::copy_from_device(buf1, m_Nbdsize[idir]);

      chrecv_dn[idir].start();
      chsend_up[idir].start();
    }

    BridgeACC::mult_staggered_xmb(v2, u, v1, m_Nsize, m_bc2, NC);

    if(do_comm[idir] > 0){
      chsend_up[idir].wait();
      chrecv_dn[idir].wait();

      BridgeACC::copy_to_device(buf2, m_Nbdsize[idir]);
      BridgeACC::mult_staggered_xm2(v2, buf2, m_Nsize, m_bc, NC);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Staggered<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 1;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Staggered<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 1;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Staggered<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 2;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Staggered<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 2;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_up[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Staggered<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 3;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
void AFopr_Staggered<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0){

    int idir = 3;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
    real_t *buf2 = (real_t*)chrecv_up[idir].ptr();
    real_t *u = m_U.ptr(index_lex.idx_G(0, 0, idir));

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
double AFopr_Staggered<AFIELD>::flop_count()
{
  return flop_count(m_mode);
}

//====================================================================
template<typename AFIELD>
double AFopr_Staggered<AFIELD>::flop_count(const std::string mode)
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

//====================================================================
//============================================================END=====
