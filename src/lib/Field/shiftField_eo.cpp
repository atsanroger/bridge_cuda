/*!
        @file    shiftField_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Field/shiftField_eo.h"
#include "Field/field_thread-inc.h"
#include "ResourceManager/threadManager.h"

const std::string ShiftField_eo::class_name = "ShiftField_eo";

//====================================================================
void ShiftField_eo::init(const int Nin)
{
  m_vl = CommonParameters::Vlevel();

  vout.paranoiac(m_vl, "%s: construction\n", class_name.c_str());

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;

  if ((m_Nx % 2) != 0) {
    vout.crucial("error in %s: Nx must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nx2   = m_Nx / 2;
  m_Nvol2 = m_Nvol / 2;

  m_Nin = Nin;

  m_yzt_eo.resize(m_Ny * m_Nz * m_Nt);
  for (int it = 0; it < m_Nt; ++it) {
    for (int iz = 0; iz < m_Nz; ++iz) {
      for (int iy = 0; iy < m_Ny; ++iy) {
        int it_global = it + Communicator::ipe(3) * m_Nt;
        int iz_global = iz + Communicator::ipe(2) * m_Nz;
        int iy_global = iy + Communicator::ipe(1) * m_Ny;
        m_yzt_eo[iy + m_Ny * (iz + m_Nz * it)]
          = (iy_global + iz_global + it_global) % 2;
      }
    }
  }

  if (m_Nin > 0) {
    vout.paranoiac(m_vl, "  Nin = %d: resetting working vectors.\n",
                   m_Nin);

    int Nvbufx = (m_Ny * m_Nz * m_Nt + 1) / 2;
    m_wt_x.reset(m_Nin, Nvbufx, 1);
    m_vt_x.reset(m_Nin, Nvbufx, 1);
    m_wt_y.reset(m_Nin, m_Nvol2 / m_Ny, 1);
    m_vt_y.reset(m_Nin, m_Nvol2 / m_Ny, 1);
    m_wt_z.reset(m_Nin, m_Nvol2 / m_Nz, 1);
    m_vt_z.reset(m_Nin, m_Nvol2 / m_Nz, 1);
    m_wt_t.reset(m_Nin, m_Nvol2 / m_Nt, 1);
    m_vt_t.reset(m_Nin, m_Nvol2 / m_Nt, 1);

    m_we.reset(m_Nin, m_Nvol2, 1);
    m_wo.reset(m_Nin, m_Nvol2, 1);
    m_ve.reset(m_Nin, m_Nvol2, 1);
    m_vo.reset(m_Nin, m_Nvol2, 1);

    m_w1.reset(m_Nin, m_Nvol, 1);
  } else {
    vout.paranoiac(m_vl, "  Nin = %d: working vectors are not set.\n",
                   m_Nin);
  }

  vout.paranoiac(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
}


//====================================================================
void ShiftField_eo::backward_h(Field& v, const Field& w,
                               const int mu, const int ieo)
{
  const int boundary_condition = 1;

  if (mu == 0) {        // x-direction
    up_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    up_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    up_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    up_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::forward_h(Field& v, const Field& w,
                              const int mu, const int ieo)
{
  const int boundary_condition = 1;

  if (mu == 0) {        // x-direction
    dn_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    dn_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    dn_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    dn_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::backward_h(Field& v, const Field& w,
                               const int boundary_condition, const int mu, const int ieo)
{
  if (mu == 0) {        // x-direction
    up_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    up_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    up_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    up_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::forward_h(Field& v, const Field& w,
                              const int boundary_condition, const int mu, const int ieo)
{
  if (mu == 0) {        // x-direction
    dn_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    dn_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    dn_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    dn_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::backward(Field& v, const Field& w, const int mu)
{
  backward(v, w, 1, mu);
}


//====================================================================
void ShiftField_eo::forward(Field& v, const Field& w, const int mu)
{
  forward(v, w, 1, mu);
}


//====================================================================
void ShiftField_eo::backward(Field& v, const Field& w,
                             const int boundary_condition, const int mu)
{
  const int Nex = w.nex();

  if (m_Nin > 0) {
    assert(w.check_size(m_Nin, m_Nvol, Nex));
    assert(v.check_size(m_Nin, m_Nvol, Nex));

    for (int ex = 0; ex < Nex; ++ex) {
      copy(m_w1, 0, w, ex);
#pragma omp barrier

      m_index_eo.splitField(m_we, m_wo, m_w1);

      backward_h(m_ve, m_wo, boundary_condition, mu, 0);
      backward_h(m_vo, m_we, boundary_condition, mu, 1);

      m_index_eo.mergeField(m_w1, m_ve, m_vo);

      copy(v, ex, m_w1, 0);
#pragma omp barrier
    }
  } else {
    ThreadManager::assert_single_thread(class_name);

    const int Nin   = w.nin();
    const int Nvol2 = w.nvol() / 2;

    Field field_e(Nin, Nvol2, Nex);
    Field field_o(Nin, Nvol2, Nex);

    Field field_se(Nin, Nvol2, Nex);
    Field field_so(Nin, Nvol2, Nex);

    m_index_eo.splitField(field_e, field_o, w);

    backward_h(field_se, field_o, boundary_condition, mu, 0);
    backward_h(field_so, field_e, boundary_condition, mu, 1);

    m_index_eo.mergeField(v, field_se, field_so);
  }
}


//====================================================================
void ShiftField_eo::forward(Field& v, const Field& w,
                            const int boundary_condition, const int mu)
{
  const int Nex = w.nex();

  if (m_Nin > 0) {
    assert(w.check_size(m_Nin, m_Nvol, Nex));
    assert(v.check_size(m_Nin, m_Nvol, Nex));

    for (int ex = 0; ex < Nex; ++ex) {
      copy(m_w1, 0, w, ex);
#pragma omp barrier

      m_index_eo.splitField(m_we, m_wo, m_w1);

      forward_h(m_ve, m_wo, boundary_condition, mu, 0);
      forward_h(m_vo, m_we, boundary_condition, mu, 1);

      m_index_eo.mergeField(m_w1, m_ve, m_vo);

      copy(v, ex, m_w1, 0);
#pragma omp barrier
    }
  } else {
    ThreadManager::assert_single_thread(class_name);

    const int Nin   = w.nin();
    const int Nvol2 = w.nvol() / 2;
    const int Nex   = w.nex();
    assert(v.nvol() == m_Nvol2);
    assert(v.check_size(Nin, m_Nvol, Nex));

    Field wt_e(Nin, Nvol2, Nex);
    Field wt_o(Nin, Nvol2, Nex);

    Field vt_e(Nin, Nvol2, Nex);
    Field vt_o(Nin, Nvol2, Nex);

    m_index_eo.splitField(wt_e, wt_o, w);

    forward_h(vt_e, wt_o, mu, 0);
    forward_h(vt_o, wt_e, mu, 1);

    m_index_eo.mergeField(v, vt_e, vt_o);
  }
}


//====================================================================
void ShiftField_eo::up_xh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = bc;

  int Nin;
  int Nex    = w.nex();
  int Nvbufx = (m_Ny * m_Nz * m_Nt + 1) / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_x.ptr(0);
    vt  = m_vt_x.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * Nvbufx];
      vt = new double[Nin * Nvbufx];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      if ((ix2 == 0) && (keo == 1)) {
        int iyzt2 = iyzt / 2;
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * iyzt2] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      int size = Nin * Nvbufx;
      Communicator::exchange(size, &vt[0], &wt[0], 0, 1, 0);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      int ix2n = ix2 + keo;
      int nei  = ix2n + m_Nx2 * iyzt;
      if (ix2n < m_Nx2) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        int iyzt2 = iyzt / 2;
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = vt[in + Nin * iyzt2];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::dn_xh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = bc;

  int Nin;
  int Nex    = w.nex();
  int Nvbufx = (m_Ny * m_Nz * m_Nt + 1) / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_x.ptr(0);
    vt  = m_vt_x.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * Nvbufx];
      vt = new double[Nin * Nvbufx];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      if ((ix2 == m_Nx2 - 1) && (keo == 0)) {
        int iyzt2 = iyzt / 2;
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * iyzt2] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      int size = Nin * Nvbufx;
      Communicator::exchange(size, &vt[0], &wt[0], 0, -1, 4);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int keo  = (ieo + m_yzt_eo[iyzt]) % 2;
      int ix2n = ix2 - (1 - keo);
      int nei  = ix2n + m_Nx2 * iyzt;
      if (ix2n >= 0) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        int iyzt2 = iyzt / 2;
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = bc2 * vt[in + Nin * iyzt2];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::up_yh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_y.ptr(0);
    vt  = m_vt_y.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Ny)];
      vt = new double[Nin * (m_Nvol2 / m_Ny)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix2 + m_Nx2 * izt;
      if (iy == 0) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixzt] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Ny);
      Communicator::exchange(size, &vt[0], &wt[0], 1, 1, 1);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix2 + m_Nx2 * izt;
      int nei  = ix2 + m_Nx2 * (iy + 1 + m_Ny * izt);

      if (iy < m_Ny - 1) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = vt[in + Nin * ixzt];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::dn_yh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_y.ptr(0);
    vt  = m_vt_y.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Ny)];
      vt = new double[Nin * (m_Nvol2 / m_Ny)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix2 + m_Nx2 * izt;
      if (iy == m_Ny - 1) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixzt] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Ny);
      Communicator::exchange(size, &vt[0], &wt[0], 1, -1, 5);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix2  = site % m_Nx2;
      int iyzt = site / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix2 + m_Nx2 * izt;
      int nei  = ix2 + m_Nx2 * (iy - 1 + m_Ny * izt);

      if (iy > 0) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = bc2 * vt[in + Nin * ixzt];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::up_zh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

  int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_z.ptr(0);
    vt  = m_vt_z.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Nz)];
      vt = new double[Nin * (m_Nvol2 / m_Nz)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      if (iz == 0) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixyt] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Nz);
      Communicator::exchange(size, &vt[0], &wt[0], 2, 1, 2);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      int nei  = ixy + Nxy * (iz + 1 + m_Nz * it);

      if (iz < m_Nz - 1) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = vt[in + Nin * ixyt];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::dn_zh(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

  int Nxy = m_Nx2 * m_Ny;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_z.ptr(0);
    vt  = m_vt_z.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Nz)];
      vt = new double[Nin * (m_Nvol2 / m_Nz)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      if (iz == m_Nz - 1) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixyt] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Nz);
      Communicator::exchange(size, &vt[0], &wt[0], 2, -1, 6);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxy;
      int izt  = site / Nxy;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int ixyt = ixy + Nxy * it;
      int nei  = ixy + Nxy * (iz - 1 + m_Nz * it);

      if (iz > 0) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = bc2 * vt[in + Nin * ixyt];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::up_th(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

  int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_t.ptr(0);
    vt  = m_vt_t.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Nt)];
      vt = new double[Nin * (m_Nvol2 / m_Nt)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == 0) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixyz] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Nt);
      Communicator::exchange(size, &vt[0], &wt[0], 3, 1, 3);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it + 1);

      if (it < m_Nt - 1) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = vt[in + Nin * ixyz];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//====================================================================
void ShiftField_eo::dn_th(Field& v, const Field& w,
                          const int bc, const int ieo)
{
  double bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  double *wt, *vt;

  int Nxyz = m_Nx2 * m_Ny * m_Nz;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_t.ptr(0);
    vt  = m_vt_t.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol2 / m_Nt)];
      vt = new double[Nin * (m_Nvol2 / m_Nt)];
      ns = m_Nvol2;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol2 * ex);
    const double *wp = w.ptr(Nin * m_Nvol2 * ex);

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if (it == m_Nt - 1) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixyz] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol2 / m_Nt);
      Communicator::exchange(size, &vt[0], &wt[0], 3, -1, 7);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      int nei  = ixyz + Nxyz * (it - 1);

      if (it > 0) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = bc2 * vt[in + Nin * ixyz];
        }
      }
    }
  }

  if ((m_Nin == 0) && (ith == 0)) {
    delete[] wt;
    delete[] vt;
  }

#pragma omp barrier
}


//============================================================END=====
