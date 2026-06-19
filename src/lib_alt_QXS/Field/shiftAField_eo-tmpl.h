/*!
      @file    shiftAField_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


template<typename AFIELD>
const std::string ShiftAField_eo<AFIELD>::class_name =
  "ShiftAField_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::init(int Nin)
{
  int              Ndim = CommonParameters::Ndim();
  std::vector<int> bc(Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    bc[mu] = 1;
  }

  init(Nin, bc);
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::init(int Nin, std::vector<int>& bc)
{
  m_vl = CommonParameters::Vlevel();

  int req_comm = 0;
  //int req_comm = 1;  // set 1 if communication forced any time

  vout.general(m_vl, "%s: being constructed.\n", class_name.c_str());

  m_Nin = Nin;
  vout.general(m_vl, "  Nin = %d\n", m_Nin);

  m_Nx    = CommonParameters::Nx();
  m_Ny    = CommonParameters::Ny();
  m_Nz    = CommonParameters::Nz();
  m_Nt    = CommonParameters::Nt();
  m_Nvol  = m_Nx * m_Ny * m_Nz * m_Nt;
  m_Nx2   = m_Nx / 2;
  m_Nvol2 = m_Nvol / 2;

  m_Ndim = CommonParameters::Ndim();

  m_Nx2v  = m_Nx2 / VLENX;
  m_Nyv   = m_Ny / VLENY;
  m_Nst2v = m_Nvol2 / VLEN;

  vout.general(m_vl, "  VLENX = %d  VLENY = %d\n", VLENX, VLENY);
  vout.general(m_vl, "  Nx2  = %d  Nvol2 = %d\n", m_Nx2, m_Nvol2);
  vout.general(m_vl, "  Nx2v = %d  Nyv  = %d\n", m_Nx2v, m_Nyv);

  if (bc.size() != m_Ndim) {
    vout.crucial(m_vl, "%s: incorrect size of boundary condition\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_boundary.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  m_Nbdsize[0] = m_Nin * m_Ny * m_Nz * m_Nt;
  m_Nbdsize[1] = m_Nin * m_Nx * m_Nz * m_Nt;
  m_Nbdsize[2] = m_Nin * m_Nx * m_Ny * m_Nt;
  m_Nbdsize[3] = m_Nin * m_Nx * m_Ny * m_Nz;

  setup_channels();

  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::tidyup()
{
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    int Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif
    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::backward(AFIELD& v, const AFIELD& w,
                                      const int mu, const int ieo)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_eo<real_t, AFIELD::IMPL> index;

  for (int ex = 0; ex < Nex; ++ex) {
    real_t *vp = v.ptr(index.idxh(0, m_Nin, 0, ex));
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idxh(0, m_Nin, 0, ex));

    if (mu == 0) {
      //up_xh_simd(vp, wp, ieo);
      up_xh_naive(vp, wp, ieo);
    } else if (mu == 1) {
      up_yh_nv(vp, wp, ieo);
    } else if (mu == 2) {
      up_zh_nv(vp, wp, ieo);
    } else if (mu == 3) {
      up_th_nv(vp, wp, ieo);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::backward(AFIELD& v, const int ex1,
                                      const AFIELD& w, const int ex2,
                                      const int mu, const int ieo)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_eo<real_t, AFIELD::IMPL> index;

  real_t *vp = v.ptr(index.idxh(0, m_Nin, 0, ex1));
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idxh(0, m_Nin, 0, ex2));

  if (mu == 0) {
    //up_xh_simd(vp, wp, ieo);
    up_xh_naive(vp, wp, ieo);
  } else if (mu == 1) {
    up_yh_nv(vp, wp, ieo);
  } else if (mu == 2) {
    up_zh_nv(vp, wp, ieo);
  } else if (mu == 3) {
    up_th_nv(vp, wp, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::forward(AFIELD& v, const AFIELD& w,
                                     const int mu, const int ieo)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_eo<real_t, AFIELD::IMPL> index;

  for (int ex = 0; ex < Nex; ++ex) {
    real_t *vp = v.ptr(index.idxh(0, m_Nin, 0, ex));
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idxh(0, m_Nin, 0, ex));

    if (mu == 0) {
      //dn_xh_simd(vp, wp, ieo);
      dn_xh_naive(vp, wp, ieo);
    } else if (mu == 1) {
      dn_yh_nv(vp, wp, ieo);
    } else if (mu == 2) {
      dn_zh_nv(vp, wp, ieo);
    } else if (mu == 3) {
      dn_th_nv(vp, wp, ieo);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::forward(AFIELD& v, const int ex1,
                                     const AFIELD& w, const int ex2,
                                     const int mu, const int ieo)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_eo<real_t, AFIELD::IMPL> index;

  real_t *vp = v.ptr(index.idxh(0, m_Nin, 0, ex1));
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idxh(0, m_Nin, 0, ex2));

  if (mu == 0) {
    //dn_xh_simd(vp, wp, ieo);
    dn_xh_naive(vp, wp, ieo);
  } else if (mu == 1) {
    dn_yh_nv(vp, wp, ieo);
  } else if (mu == 2) {
    dn_zh_nv(vp, wp, ieo);
  } else if (mu == 3) {
    dn_th_nv(vp, wp, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
#ifdef USE_QXS_ACLE

template<typename AFIELD>
void ShiftAField_eo<AFIELD>::up_xh_simd(real_t *vp, real_t *wp,
                                        const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = real_t(m_boundary[0]);

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  // qqqqqq

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

#pragma omp barrier

  if (do_comm[0] == 1) {
#pragma omp barrier

#pragma omp master
    {
      chrecv_up[0].start();
      chsend_dn[0].start();
      // communication is not overlapped in x-direction.
      chrecv_up[0].wait();
      chsend_dn[0].wait();
    }

#pragma omp barrier
  }

  Vsimd_t xt;
  for (int site = is; site < ns; ++site) {
    int ix2  = site % m_Nx2v;
    int iyzt = site / m_Nx2v;
    int leo  = index_alt.leo(VLENY * iyzt);
    int jeo  = (ieo + leo) % 2;
    if ((ix2 < m_Nx2v - 1) || (do_comm[0] == 0)) {
      int nei = site + 1;
      if (ix2 == m_Nx2v - 1) nei = 0 + m_Nx2v * iyzt;
      for (int in = 0; in < m_Nin; ++in) {
        int iw1 = VLEN * (in + m_Nin * site);
        int iw2 = VLEN * (in + m_Nin * nei);
        shift_vec2_xbw_eo(&xt, &wp[iw1], &wp[iw2], jeo, 1);
        save_vec(&vp[iw1], &xt, 1);
      }
    }
  }


#pragma omp barrier
}


#endif

//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::up_xh_naive(real_t *vp, real_t *wp,
                                         const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = real_t(m_boundary[0]);

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  if (do_comm[0] == 1) {
    for (int site = is; site < ns; ++site) {
      int ix2   = site % m_Nx2;
      int iyzt  = site / m_Nx2;
      int iyzt2 = iyzt / 2;
      int leo   = index_alt.leo(iyzt);
      if ((ix2 == 0) && (leo == (1 - ieo))) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idxh(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt2] = bc2 * wp[index];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_up[0].start();
      chsend_dn[0].start();
      // communication is not overlapped in x-direction.
      chrecv_up[0].wait();
      chsend_dn[0].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ix2  = site % m_Nx2;
    int iyzt = site / m_Nx2;
    int leo  = index_alt.leo(iyzt);
    if (leo == ieo) {
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        vp[iv] = wp[iv];
      }
    } else {
      if ((ix2 < m_Nx2 - 1) || (do_comm[0] == 0)) {
        int nei = ix2 + 1 + m_Nx2 * iyzt;
        int bc3 = 1.0;
        if (ix2 == m_Nx2 - 1) {
          nei = 0 + m_Nx2 * iyzt;
          bc3 = bc2;
        }
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          int iw = index_alt.idxh(in, m_Nin, nei, ex);
          vp[iv] = bc3 * wp[iw];
        }
      } else {
        for (int in = 0; in < m_Nin; ++in) {
          int iv    = index_alt.idxh(in, m_Nin, site, ex);
          int iyzt2 = iyzt / 2;
          vp[iv] = buf2[in + m_Nin * iyzt2];
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
#ifdef USE_QXS_ACLE

template<typename AFIELD>
void ShiftAField_eo<AFIELD>::dn_xh_simd(real_t *vp, real_t *wp,
                                        const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == Communicator::npe(0) - 1) {
    bc2 = real_t(m_boundary[0]);
  }

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst2v);

#pragma omp barrier

  if (do_comm[0] == 1) {
#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[0].start();
      chsend_up[0].start();
      // communication is not overlapped in x-direction.
      chrecv_dn[0].wait();
      chsend_up[0].wait();
    }

#pragma omp barrier
  }

  Vsimd_t xt;
  for (int site = is; site < ns; ++site) {
    int ix2  = site % m_Nx2v;
    int iyzt = site / m_Nx2v;
    int leo  = index_alt.leo(VLENY * iyzt);
    int jeo  = (ieo + leo) % 2;
    if ((ix2 > 0) || (do_comm[0] == 0)) {
      int nei = site - 1;
      if (ix2 == 0) nei = m_Nx2v - 1 + m_Nx2v * iyzt;
      for (int in = 0; in < m_Nin; ++in) {
        int iw1 = VLEN * (in + m_Nin * site);
        int iw2 = VLEN * (in + m_Nin * nei);
        shift_vec2_xfw_eo(&xt, &wp[iw1], &wp[iw2], jeo, 1);
        save_vec(&vp[iw1], &xt, 1);
      }
    }
  }



#pragma omp barrier
}


#endif

//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::dn_xh_naive(real_t *vp, real_t *wp,
                                         const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == Communicator::npe(0) - 1) {
    bc2 = real_t(m_boundary[0]);
  }

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  if (do_comm[0] == 1) {
    for (int site = is; site < ns; ++site) {
      int ix2   = site % m_Nx2;
      int iyzt  = site / m_Nx2;
      int iyzt2 = iyzt / 2;
      int leo   = index_alt.leo(iyzt);
      if ((ix2 == m_Nx2 - 1) && (leo == ieo)) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idxh(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt2] = bc2 * wp[index];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chrecv_dn[0].start();
      chsend_up[0].start();
      // communication is not overlapped in x-direction.
      chrecv_dn[0].wait();
      chsend_up[0].wait();
    }

#pragma omp barrier
  }

  for (int site = is; site < ns; ++site) {
    int ix2  = site % m_Nx2;
    int iyzt = site / m_Nx2;
    int leo  = index_alt.leo(iyzt);
    if (leo == (1 - ieo)) {
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        vp[iv] = wp[iv];
      }
    } else {
      if ((ix2 > 0) || (do_comm[0] == 0)) {
        int nei = ix2 - 1 + m_Nx2 * iyzt;
        int bc3 = 1.0;
        if (ix2 == 0) {
          nei = m_Nx2 - 1 + m_Nx2 * iyzt;
          bc3 = bc2;
        }
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          int iw = index_alt.idxh(in, m_Nin, nei, ex);
          vp[iv] = bc3 * wp[iw];
        }
      } else {
        for (int in = 0; in < m_Nin; ++in) {
          int iv    = index_alt.idxh(in, m_Nin, site, ex);
          int iyzt2 = iyzt / 2;
          vp[iv] = buf2[in + m_Nin * iyzt2];
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::up_yh_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_up[1].start();
    }
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nx2;
      int iy  = (site / m_Nx2) % m_Ny;
      int izt = site / (m_Nx2 * m_Ny);
      if (iy == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idxh(in, m_Nin, site, ex);
          int ixzt = ix + m_Nx2 * izt;
          buf1[in + m_Nin * ixzt] = bc2 * wp[iw];
        }
      }
    }

#pragma omp barrier
#pragma omp master
    {
      chsend_dn[1].start();
    }
  } // if(do_comm[1] == 1)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nx2;
    int iy  = (site / m_Nx2) % m_Ny;
    int izt = site / (m_Nx2 * m_Ny);
    if ((iy < m_Ny - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % m_Ny;
      int    nei = ix + m_Nx2 * (iy2 + m_Ny * izt);
      real_t bc3 = 1.0;
      if (iy == m_Ny - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_up[1].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iy   = (site / m_Nx2) % m_Ny;
      int izt  = site / (m_Nx2 * m_Ny);
      int ixzt = ix + m_Nx2 * izt;
      if (iy == m_Ny - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          vp[iv] = buf2[in + m_Nin * ixzt];
        }
      }
    }
#pragma omp master
    {
      chsend_dn[1].wait();
    }
  } // if(do_comm[1] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::dn_yh_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_dn[1].start();
    }
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nx2;
      int iy  = (site / m_Nx2) % m_Ny;
      int izt = site / (m_Nx2 * m_Ny);
      if (iy == m_Ny - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idxh(in, m_Nin, site, ex);
          int ixzt = ix + m_Nx2 * izt;
          buf1[in + m_Nin * ixzt] = wp[iw];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_up[1].start();
    }
  } // if(do_comm[1] == 1)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nx2;
    int iy  = (site / m_Nx2) % m_Ny;
    int izt = site / (m_Nx2 * m_Ny);
    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + m_Ny) % m_Ny;
      int    nei = ix + m_Nx2 * (iy2 + m_Ny * izt);
      real_t bc3 = 1.0;
      if (iy == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_dn[1].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx2;
      int iy   = (site / m_Nx2) % m_Ny;
      int izt  = site / (m_Nx2 * m_Ny);
      int ixzt = ix + m_Nx2 * izt;
      if (iy == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          vp[iv] = bc2 * buf2[in + m_Nin * ixzt];
        }
      }
    }

#pragma omp master
    {
      chsend_up[1].wait();
    }
  } // if(do_comm[1] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::up_zh_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

  if (do_comm[2] == 1) {
#pragma omp master
    {
      chrecv_up[2].start();
    }
    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idxh(in, m_Nin, site, ex);
          int ixyt = ixy + Nxy * it;
          buf1[in + m_Nin * ixyt] = bc2 * wp[iw];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_dn[2].start();
    }
  } // if(do_comm[2] == 1)

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);
    if ((iz < m_Nz - 1) || (do_comm[2] == 0)) {
      int    iz2 = (iz + 1) % m_Nz;
      int    nei = ixy + Nxy * (iz2 + m_Nz * it);
      real_t bc3 = 1.0;
      if (iz == m_Nz - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[2] == 1) {
#pragma omp master
    {
      chrecv_up[2].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == m_Nz - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv   = index_alt.idxh(in, m_Nin, site, ex);
          int ixyt = ixy + Nxy * it;
          vp[iv] = buf2[in + m_Nin * ixyt];
        }
      }
    }
#pragma omp master
    {
      chsend_dn[2].wait();
    }
  } // if(do_comm[2] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::dn_zh_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);



  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

  if (do_comm[2] == 1) {
#pragma omp master
    {
      chrecv_dn[2].start();
    }

    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == m_Nz - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idxh(in, m_Nin, site, ex);
          int ixyt = ixy + Nxy * it;
          buf1[in + m_Nin * ixyt] = wp[iw];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_up[2].start();
    }
  } // if(do_comm[2] == 1)

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxy;
    int iz  = (site / Nxy) % m_Nz;
    int it  = site / (Nxy * m_Nz);
    if ((iz > 0) || (do_comm[2] == 0)) {
      int    iz2 = (iz - 1 + m_Nz) % m_Nz;
      int    nei = ixy + Nxy * (iz2 + m_Nz * it);
      real_t bc3 = 1.0;
      if (iz == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[2] == 1) {
#pragma omp master
    {
      chrecv_dn[2].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxy;
      int iz  = (site / Nxy) % m_Nz;
      int it  = site / (Nxy * m_Nz);
      if (iz == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv   = index_alt.idxh(in, m_Nin, site, ex);
          int ixyt = ixy + Nxy * it;
          vp[iv] = bc2 * buf2[in + m_Nin * ixyt];
        }
      }
    }

#pragma omp master
    {
      chsend_up[2].wait();
    }
  } // if(do_comm[2] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::up_th_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

  if (do_comm[3] == 1) {
#pragma omp master
    {
      chrecv_up[3].start();
    }
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw = index_alt.idxh(in, m_Nin, site, ex);
          buf1[in + m_Nin * ixyz] = bc2 * wp[iw];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_dn[3].start();
    }
  } // if(do_comm[3] == 1)


  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int    it2 = (it + 1) % m_Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t bc3 = 1.0;
      if (it == m_Nt - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[3] == 1) {
#pragma omp master
    {
      chrecv_up[3].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          vp[iv] = buf2[in + m_Nin * ixyz];
        }
      }
    }
#pragma omp master
    {
      chsend_dn[3].wait();
    }
  } // if(do_comm[3] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_eo<AFIELD>::dn_th_nv(real_t *vp, real_t *wp,
                                      const int ieo)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  AIndex_eo<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

  if (do_comm[3] == 1) {
#pragma omp master
    {
      chrecv_dn[3].start();
    }
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw = index_alt.idxh(in, m_Nin, site, ex);
          buf1[in + m_Nin * ixyz] = wp[iw];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_up[3].start();
    }
  } // if(do_comm[3] == 1)

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;
    if ((it > 0) || (do_comm[3] == 0)) {
      int    it2 = (it - 1 + m_Nt) % m_Nt;
      int    nei = ixyz + Nxyz * it2;
      real_t bc3 = 1.0;
      if (it == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idxh(in, m_Nin, site, ex);
        int iw = index_alt.idxh(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[3] == 1) {
#pragma omp master
    {
      chrecv_dn[3].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idxh(in, m_Nin, site, ex);
          vp[iv] = bc2 * buf2[in + m_Nin * ixyz];
        }
      }
    }

#pragma omp master
    {
      chsend_up[3].wait();
    }
  } // if(do_comm[3] == 1)

#pragma omp barrier
}


//============================================================END=====
