/*!
        @file    shiftField_lex.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Field/shiftField_lex.h"
#include "Field/field_thread-inc.h"
#include "ResourceManager/threadManager.h"

const std::string ShiftField_lex::class_name = "ShiftField_lex";

//====================================================================
void ShiftField_lex::init(const int Nin)
{
  m_vl = CommonParameters::Vlevel();

  vout.paranoiac(m_vl, "%s: construction\n", class_name.c_str());

  m_Nx = CommonParameters::Nx();
  m_Ny = CommonParameters::Ny();
  m_Nz = CommonParameters::Nz();
  m_Nt = CommonParameters::Nt();

  m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;

  m_Nin = Nin;

  if (m_Nin > 0) {
    vout.paranoiac(m_vl, "  Nin = %d: resetting working vectors.\n",
                   m_Nin);

    m_wt_x.reset(m_Nin, m_Nvol / m_Nx, 1);
    m_vt_x.reset(m_Nin, m_Nvol / m_Nx, 1);
    m_wt_y.reset(m_Nin, m_Nvol / m_Ny, 1);
    m_vt_y.reset(m_Nin, m_Nvol / m_Ny, 1);
    m_wt_z.reset(m_Nin, m_Nvol / m_Nz, 1);
    m_vt_z.reset(m_Nin, m_Nvol / m_Nz, 1);
    m_wt_t.reset(m_Nin, m_Nvol / m_Nt, 1);
    m_vt_t.reset(m_Nin, m_Nvol / m_Nt, 1);
  } else {
    vout.paranoiac(m_vl, "  Nin = %d: working vectors are not set.\n",
                   m_Nin);
  }

  vout.paranoiac(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
}


//====================================================================
void ShiftField_lex::backward(Field& v, const Field& w, const int mu)
{
  const int boundary_condition = 1;

  if (mu == 0) {  // x-direction
    up_x(v, w, boundary_condition);
  } else if (mu == 1) {
    up_y(v, w, boundary_condition);
  } else if (mu == 2) {
    up_z(v, w, boundary_condition);
  } else if (mu == 3) {
    up_t(v, w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::forward(Field& v, const Field& w, const int mu)
{
  const int boundary_condition = 1;

  if (mu == 0) {
    dn_x(v, w, boundary_condition);
  } else if (mu == 1) {
    dn_y(v, w, boundary_condition);
  } else if (mu == 2) {
    dn_z(v, w, boundary_condition);
  } else if (mu == 3) {
    dn_t(v, w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::backward(Field& v, const Field& w,
                              const int boundary_condition, const int mu)
{
  if (mu == 0) {  // x-direction
    up_x(v, w, boundary_condition);
  } else if (mu == 1) {
    up_y(v, w, boundary_condition);
  } else if (mu == 2) {
    up_z(v, w, boundary_condition);
  } else if (mu == 3) {
    up_t(v, w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::forward(Field& v, const Field& w,
                             const int boundary_condition, const int mu)
{
  if (mu == 0) {
    dn_x(v, w, boundary_condition);
  } else if (mu == 1) {
    dn_y(v, w, boundary_condition);
  } else if (mu == 2) {
    dn_z(v, w, boundary_condition);
  } else if (mu == 3) {
    dn_t(v, w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::up_x(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

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
      wt = new double[Nin * (m_Nvol / m_Nx)];
      vt = new double[Nin * (m_Nvol / m_Nx)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == 0) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * iyzt] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol / m_Nx);
      Communicator::exchange(size, &vt[0], &wt[0], 0, 1, 0);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int nei  = ix + 1 + m_Nx * iyzt;
      if (ix < m_Nx - 1) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = vt[in + Nin * iyzt];
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
void ShiftField_lex::dn_x(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(0) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

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
      wt = new double[Nin * (m_Nvol / m_Nx)];
      vt = new double[Nin * (m_Nvol / m_Nx)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      if (ix == m_Nx - 1) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * iyzt] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol / m_Nx);
      Communicator::exchange(size, &vt[0], &wt[0], 0, -1, 4);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int nei  = ix - 1 + m_Nx * iyzt;
      if (ix > 0) {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = wp[in + Nin * nei];
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          vp[in + Nin * site] = bc2 * vt[in + Nin * iyzt];
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
void ShiftField_lex::up_y(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

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
      wt = new double[Nin * (m_Nvol / m_Ny)];
      vt = new double[Nin * (m_Nvol / m_Ny)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      if (iy == 0) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixzt] = bc2 * wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol / m_Ny);
      Communicator::exchange(size, &vt[0], &wt[0], 1, 1, 1);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      int nei  = ix + m_Nx * (iy + 1 + m_Ny * izt);

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
void ShiftField_lex::dn_y(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(1) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

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
      wt = new double[Nin * (m_Nvol / m_Ny)];
      vt = new double[Nin * (m_Nvol / m_Ny)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      if (iy == m_Ny - 1) {
        for (int in = 0; in < Nin; ++in) {
          wt[in + Nin * ixzt] = wp[in + Nin * site];
        }
      }
    }

#pragma omp barrier

#pragma omp master
    {
      const int size = Nin * (m_Nvol / m_Ny);
      Communicator::exchange(size, &vt[0], &wt[0], 1, -1, 5);
    }

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nx;
      int iyzt = site / m_Nx;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int ixzt = ix + m_Nx * izt;
      int nei  = ix + m_Nx * (iy - 1 + m_Ny * izt);

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
void ShiftField_lex::up_z(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double *wt, *vt;

  int Nxy = m_Nx * m_Ny;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_z.ptr(0);
    vt  = m_vt_z.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol / m_Nz)];
      vt = new double[Nin * (m_Nvol / m_Nz)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

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
      const int size = Nin * (m_Nvol / m_Nz);
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
void ShiftField_lex::dn_z(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(2) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double *wt, *vt;

  int Nxy = m_Nx * m_Ny;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_z.ptr(0);
    vt  = m_vt_z.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol / m_Nz)];
      vt = new double[Nin * (m_Nvol / m_Nz)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

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
      const int size = Nin * (m_Nvol / m_Nz);
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
void ShiftField_lex::up_t(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double *wt, *vt;

  int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_t.ptr(0);
    vt  = m_vt_t.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol / m_Nt)];
      vt = new double[Nin * (m_Nvol / m_Nt)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

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
      const int size = Nin * (m_Nvol / m_Nt);
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
void ShiftField_lex::dn_t(Field& v, const Field& w, const int bc)
{
  double bc2 = 1.0;
  if (Communicator::ipe(3) == 0) bc2 = bc;

  int Nin;
  int Nex = w.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double *wt, *vt;

  int Nxyz = m_Nx * m_Ny * m_Nz;

#pragma omp barrier

  if (m_Nin > 0) {
    wt  = m_wt_t.ptr(0);
    vt  = m_vt_t.ptr(0);
    Nin = m_Nin;
  } else {
    is  = 0;
    Nin = w.nin();
    if (ith == 0) {
      wt = new double[Nin * (m_Nvol / m_Nt)];
      vt = new double[Nin * (m_Nvol / m_Nt)];
      ns = m_Nvol;
    } else {
      wt = 0;
      vt = 0;
      ns = 0;
    }
  }

  for (int ex = 0; ex < Nex; ++ex) {
    double       *vp = v.ptr(Nin * m_Nvol * ex);
    const double *wp = w.ptr(Nin * m_Nvol * ex);

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
      const int size = Nin * (m_Nvol / m_Nt);
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
