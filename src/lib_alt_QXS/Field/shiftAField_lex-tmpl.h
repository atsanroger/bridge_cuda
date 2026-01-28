/*!
      @file    shiftAField_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


template<typename AFIELD>
const std::string ShiftAField_lex<AFIELD>::class_name =
  "ShiftAField_lex<AFIELD>";
//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::init(int Nin)
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
void ShiftAField_lex<AFIELD>::init(int Nin, std::vector<int>& bc)
{
  m_vl = CommonParameters::Vlevel();

  int req_comm = 0;
  //int req_comm = 1;  // set 1 if communication forced any time

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_Nin = Nin;
  vout.general(m_vl, "  Nin = %d\n", m_Nin);

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;

  m_Ndim = CommonParameters::Ndim();

  m_Nxv  = m_Nx / VLENX;
  m_Nyv  = m_Ny / VLENY;
  m_Nstv = m_Nvol / VLEN;
  vout.general(m_vl, "  VLENX = %d  Nxv  = %d\n", VLENX, m_Nxv);
  vout.general(m_vl, "  VLENY = %d  Nyv  = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %d  Nstv = %d\n", VLEN, m_Nstv);

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
void ShiftAField_lex<AFIELD>::tidyup()
{
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::setup_channels()
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
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const AFIELD& w,
                                       const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int ex = 0; ex < Nex; ++ex) {
    real_t *vp = v.ptr(index.idx(0, m_Nin, 0, ex));
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idx(0, m_Nin, 0, ex));

    if (mu == 0) {
      //up_x(vp, wp);
      up_x_nv(vp, wp);
    } else if (mu == 1) {
      //up_y(vp, wp);
      up_y_nv(vp, wp);
    } else if (mu == 2) {
      //up_z(vp, wp);
      up_z_nv(vp, wp);
    } else if (mu == 3) {
      //up_t(vp, wp);
      up_t_nv(vp, wp);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const int ex1,
                                       const AFIELD& w, const int ex2,
                                       const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  real_t *vp = v.ptr(index.idx(0, m_Nin, 0, ex1));
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idx(0, m_Nin, 0, ex2));

  if (mu == 0) {
    //up_x(vp, wp);
    up_x_nv(vp, wp);
  } else if (mu == 1) {
    //up_y(vp, wp);
    up_y_nv(vp, wp);
  } else if (mu == 2) {
    //    up_z(vp, wp);
    up_z_nv(vp, wp);
  } else if (mu == 3) {
    //up_t(vp, wp);
    up_t_nv(vp, wp);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const AFIELD& w,
                                      const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int ex = 0; ex < Nex; ++ex) {
    real_t *vp = v.ptr(index.idx(0, m_Nin, 0, ex));
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idx(0, m_Nin, 0, ex));

    if (mu == 0) {
      // dn_x(vp, wp);
      dn_x_nv(vp, wp);
    } else if (mu == 1) {
      // dn_y(vp, wp);
      dn_y_nv(vp, wp);
    } else if (mu == 2) {
      // dn_z(vp, wp);
      dn_z_nv(vp, wp);
    } else if (mu == 3) {
      // dn_t(vp, wp);
      dn_t_nv(vp, wp);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const int ex1,
                                      const AFIELD& w, const int ex2,
                                      const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  real_t *vp = v.ptr(index.idx(0, m_Nin, 0, ex1));
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(index.idx(0, m_Nin, 0, ex2));

  if (mu == 0) {
    //dn_x(vp, wp);
    dn_x_nv(vp, wp);
  } else if (mu == 1) {
    //dn_y(vp, wp);
    dn_y_nv(vp, wp);
  } else if (mu == 2) {
    //dn_z(vp, wp);
    dn_z_nv(vp, wp);
  } else if (mu == 3) {
    // dn_t(vp, wp);
    dn_t_nv(vp, wp);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_x(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = real_t(m_boundary[0]);

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[0] == 1) {
    /*
    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nxv;
      int iyzt = site/m_Nxv;
      if(ix == 0){
        real_t buf[m_Nin];
        load_vec1(buf, &wp[VLEN*Nin2*site], 0, Nin2);
        scal_th(buf, bc2, m_Nin);
        copy_th(&buf1[m_Nin*iyzt], buf, m_Nin);
      }
    }
    */

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idx(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt] = bc2 * wp[index];
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
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if ((ix < m_Nx - 1) || (do_comm[0] == 0)) {
      int    ix2 = (ix + 1) % m_Nx;
      int    nei = ix2 + m_Nx * iyzt;
      real_t bc3 = 1.0;
      if (ix == m_Nx - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    } else {
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        vp[iv] = buf2[in + m_Nin * iyzt];
      }
    }
  }

  /*
for(int site = is; site < ns; ++site){
  int ix   = site % m_Nxv;
  int iyzt = site/m_Nxv;
  int iv = VLEN * Nin2 * site;
  if(ix < m_Nxv-1){
    int nei = VLEN * Nin2 * (site + 1);
    shift_vec2_bw(&vp[iv], &wp[iv], &wp[nei], Nin2);
  }else if(do_comm[0] == 0){
    int nei = VLEN * Nin2 * (0 + m_Nxv * iyzt);
    real_t buf[m_Nin];
    load_vec1(buf, &wp[nei], 0, Nin2);
    scal_th(buf, bc2, m_Nin);
    shift_vec1_bw(&vp[iv], &wp[iv], buf, Nin2);
  }else{
    shift_vec1_bw(&vp[iv], &wp[iv], &buf2[m_Nin*iyzt], Nin2);
  }
}
  */

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_x_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = real_t(m_boundary[0]);

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[0] == 1) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idx(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt] = bc2 * wp[index];
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
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if ((ix < m_Nx - 1) || (do_comm[0] == 0)) {
      int    ix2 = (ix + 1) % m_Nx;
      int    nei = ix2 + m_Nx * iyzt;
      real_t bc3 = 1.0;
      if (ix == m_Nx - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    } else {
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        vp[iv] = buf2[in + m_Nin * iyzt];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_x(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == Communicator::npe(0) - 1) {
    bc2 = real_t(m_boundary[0]);
  }

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[0] == 1) {
    /*
    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nxv;
      int iyzt = site/m_Nxv;
      if(ix == m_Nxv-1){
        real_t buf[m_Nin];
        load_vec1(buf, &wp[VLEN*Nin2*site], VLEN2-1, Nin2);
        scal_th(buf, bc2, m_Nin);
        copy_th(&buf1[m_Nin*iyzt], buf, m_Nin);
      }
    }
    */

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == m_Nx - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idx(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt] = bc2 * wp[index];
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
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if ((ix > 0) || (do_comm[0] == 0)) {
      int    ix2 = (ix - 1 + m_Nx) % m_Nx;
      int    nei = ix2 + m_Nx * iyzt;
      real_t bc3 = 1.0;
      if (ix == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    } else {
      for (int in = 0; in < m_Nin; ++in) {
        int index = index_alt.idx(in, m_Nin, site, ex);
        vp[index] = buf2[in + m_Nin * iyzt];
      }
    }
  }

  /*
  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int iv = VLEN * Nin2 * site;
    if(ix > 0){
      int nei = VLEN * Nin2 * (site - 1);
      shift_vec2_fw(&vp[iv], &wp[iv], &wp[nei], Nin2);
    }else if(do_comm[0] == 0){
      int nei = VLEN * Nin2 * (m_Nxv-1 + m_Nxv * iyzt);
      real_t buf[m_Nin];
      load_vec1(buf, &wp[nei], VLEN2-1, Nin2);
      scal_th(buf, bc2, m_Nin);
      shift_vec1_fw(&vp[iv], &wp[iv], buf, Nin2);
    }else{
      shift_vec1_fw(&vp[iv], &wp[iv], &buf2[m_Nin*iyzt], Nin2);
    }
  }
  */

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_x_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(0) == Communicator::npe(0) - 1) {
    bc2 = real_t(m_boundary[0]);
  }

  int ex = 0; // ex loop is outside this method

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[0] == 1) {
    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == m_Nx - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int index = index_alt.idx(in, m_Nin, site, ex);
          buf1[in + m_Nin * iyzt] = bc2 * wp[index];
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
    int ix   = site % m_Nx;
    int iyzt = site / m_Nx;
    if ((ix > 0) || (do_comm[0] == 0)) {
      int    ix2 = (ix - 1 + m_Nx) % m_Nx;
      int    nei = ix2 + m_Nx * iyzt;
      real_t bc3 = 1.0;
      if (ix == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    } else {
      for (int in = 0; in < m_Nin; ++in) {
        int index = index_alt.idx(in, m_Nin, site, ex);
        vp[index] = buf2[in + m_Nin * iyzt];
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_y(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_up[1].start();
    }
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Ny;
      int izt = site / (m_Nxv * m_Ny);
      if (iy == 0) {
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ix + m_Nxv * izt);
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
      }
    }
#pragma omp barrier

#pragma omp master
    {
      chsend_dn[1].start();
    }
  } // if(do_comm[1] == 1)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Ny;
    int izt = site / (m_Nxv * m_Ny);
    int iv  = VLEN * Nin2 * site;
    if ((iy < m_Ny - 1) || (do_comm[1] == 0)) {
      int     iyn = (iy + 1) % m_Ny;
      int     nei = VLEN * Nin2 * (ix + m_Nxv * (iyn + m_Ny * izt));
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (iy == m_Ny - 1) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
    }
  }

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_up[1].wait();
    }
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Ny;
      int izt = site / (m_Nxv * m_Ny);
      if (iy == m_Ny - 1) {
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ix + m_Nxv * izt);
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        save_vec(&vp[iv], vt, Nin2);
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
void ShiftAField_lex<AFIELD>::up_y_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_up[1].start();
    }

    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nx;
      int iy  = (site / m_Nx) % m_Ny;
      int izt = site / (m_Nx * m_Ny);
      if (iy == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idx(in, m_Nin, site, ex);
          int ixzt = ix + m_Nx * izt;
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
    int ix  = site % m_Nx;
    int iy  = (site / m_Nx) % m_Ny;
    int izt = site / (m_Nx * m_Ny);
    if ((iy < m_Ny - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % m_Ny;
      int    nei = ix + m_Nx * (iy2 + m_Ny * izt);
      real_t bc3 = 1.0;
      if (iy == m_Ny - 1) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
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
      int ix   = site % m_Nx;
      int iy   = (site / m_Nx) % m_Ny;
      int izt  = site / (m_Nx * m_Ny);
      int ixzt = ix + m_Nx * izt;
      if (iy == m_Ny - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idx(in, m_Nin, site, ex);
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
void ShiftAField_lex<AFIELD>::dn_y(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_dn[1].start();
    }
    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Ny;
      int izt = site / (m_Nxv * m_Ny);
      if (iy == m_Ny - 1) {
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ix + m_Nxv * izt);
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      chsend_up[1].start();
    }
  } // if(do_comm[1] == 1)

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Ny;
    int izt = site / (m_Nxv * m_Ny);
    int iv  = VLEN * Nin2 * site;
    if ((iy > 0) || (do_comm[1] == 0)) {
      int     iyn = (iy - 1 + m_Ny) % m_Ny;
      int     nei = VLEN * Nin2 * (ix + m_Nxv * (iyn + m_Ny * izt));
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (iy == 0) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
    }
  }

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_dn[1].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Ny;
      int izt = site / (m_Nxv * m_Ny);
      if (iy == 0) {
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ix + m_Nxv * izt);
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&vp[iv], vt, Nin2);
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
void ShiftAField_lex<AFIELD>::dn_y_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = real_t(m_boundary[1]);

  int Nin2 = m_Nin / 2;

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chrecv_dn[1].start();
    }

    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nx;
      int iy  = (site / m_Nx) % m_Ny;
      int izt = site / (m_Nx * m_Ny);
      if (iy == m_Ny - 1) {
        for (int in = 0; in < m_Nin; ++in) {
          int iw   = index_alt.idx(in, m_Nin, site, ex);
          int ixzt = ix + m_Nx * izt;
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
    int ix  = site % m_Nx;
    int iy  = (site / m_Nx) % m_Ny;
    int izt = site / (m_Nx * m_Ny);
    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + m_Ny) % m_Ny;
      int    nei = ix + m_Nx * (iy2 + m_Ny * izt);
      real_t bc3 = 1.0;
      if (iy == 0) bc3 = bc2;
      for (int in = 0; in < m_Nin; ++in) {
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
        vp[iv] = bc3 * wp[iw];
      }
    }
  }

  if (do_comm[1] == 1) {
#pragma omp master
    {
      chsend_up[1].wait();
      chrecv_dn[1].wait();
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iy   = (site / m_Nx) % m_Ny;
      int izt  = site / (m_Nx * m_Ny);
      int ixzt = ix + m_Nx * izt;
      if (iy == 0) {
        for (int in = 0; in < m_Nin; ++in) {
          int iv = index_alt.idx(in, m_Nin, site, ex);
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
void ShiftAField_lex<AFIELD>::up_z(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  int Nxy = m_Nxv * m_Ny;
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ixy + Nxy * it);
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
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
    int iv  = VLEN * Nin2 * site;
    if ((iz < m_Nz - 1) || (do_comm[2] == 0)) {
      int     izn = (iz + 1) % m_Nz;
      int     nei = VLEN * Nin2 * (ixy + Nxy * (izn + m_Nz * it));
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (iz == m_Nz - 1) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ixy + Nxy * it);
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        save_vec(&vp[iv], vt, Nin2);
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
void ShiftAField_lex<AFIELD>::up_z_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxy = m_Nx * m_Ny;

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
          int iw   = index_alt.idx(in, m_Nin, site, ex);
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
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
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
          int iv   = index_alt.idx(in, m_Nin, site, ex);
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
void ShiftAField_lex<AFIELD>::dn_z(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  int Nxy = m_Nxv * m_Ny;
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ixy + Nxy * it);
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
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
    int iv  = VLEN * Nin2 * site;
    if ((iz > 0) || (do_comm[2] == 0)) {
      int     izn = (iz - 1 + m_Nz) % m_Nz;
      int     nei = VLEN * Nin2 * (ixy + Nxy * (izn + m_Nz * it));
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (iz == 0) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * (ixy + Nxy * it);
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&vp[iv], vt, Nin2);
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
void ShiftAField_lex<AFIELD>::dn_z_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = real_t(m_boundary[2]);



  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxy = m_Nx * m_Ny;

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
          int iw   = index_alt.idx(in, m_Nin, site, ex);
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
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
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
          int iv   = index_alt.idx(in, m_Nin, site, ex);
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
void ShiftAField_lex<AFIELD>::up_t(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  int Nxyz = m_Nxv * m_Ny * m_Nz;

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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * ixyz;
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
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
    int iv   = VLEN * Nin2 * site;
    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int     itn = (it + 1) % m_Nt;
      int     nei = VLEN * Nin2 * (ixyz + Nxyz * itn);
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (it == m_Nt - 1) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * ixyz;
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        save_vec(&vp[iv], vt, Nin2);
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
void ShiftAField_lex<AFIELD>::up_t_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxyz = m_Nx * m_Ny * m_Nz;

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
          int iw = index_alt.idx(in, m_Nin, site, ex);
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
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
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
          int iv = index_alt.idx(in, m_Nin, site, ex);
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
void ShiftAField_lex<AFIELD>::dn_t(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int Nin2 = m_Nin / 2;

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  int Nxyz = m_Nxv * m_Ny * m_Nz;
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * ixyz;
        Vsimd_t vt[Nin2];
        load_vec(vt, &wp[iv], Nin2);
        save_vec(&buf1[ibf], vt, Nin2);
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
    int iv   = VLEN * Nin2 * site;
    if ((it > 0) || (do_comm[3] == 0)) {
      int     itn = (it - 1 + m_Nt) % m_Nt;
      int     nei = VLEN * Nin2 * (ixyz + Nxyz * itn);
      Vsimd_t vt[Nin2];
      load_vec(vt, &wp[nei], Nin2);
      if (it == 0) scal_vec(vt, bc2, Nin2);
      save_vec(&vp[iv], vt, Nin2);
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
        int     iv  = VLEN * Nin2 * site;
        int     ibf = VLEN * Nin2 * ixyz;
        Vsimd_t vt[Nin2];
        load_vec(vt, &buf2[ibf], Nin2);
        scal_vec(vt, bc2, Nin2);
        save_vec(&vp[iv], vt, Nin2);
      }
    }
#pragma omp master
    {
      chsend_up[3].wait();
    }
  } // if(do_comm[3] == 1)

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_t_nv(real_t *vp, real_t *wp)
{
  real_t bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = real_t(m_boundary[3]);

  int ex = 0;

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  AIndex_lex<real_t, AFIELD::IMPL> index_alt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  int Nxyz = m_Nx * m_Ny * m_Nz;

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
          int iw = index_alt.idx(in, m_Nin, site, ex);
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
        int iv = index_alt.idx(in, m_Nin, site, ex);
        int iw = index_alt.idx(in, m_Nin, nei, ex);
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
          int iv = index_alt.idx(in, m_Nin, site, ex);
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
