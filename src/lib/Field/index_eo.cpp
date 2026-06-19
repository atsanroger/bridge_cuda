/*!
        @file    index_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include <assert.h>

#include "Field/index_eo.h"
#include "Field/field_thread-inc.h"
#include "ResourceManager/threadManager.h"

const std::string Index_eo::class_name = "Index_eo";

//====================================================================
void Index_eo::init()
{
  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nvol = CommonParameters::Nvol();

  m_Nx2   = m_Nx / 2;
  m_Nvol2 = m_Nvol / 2;

  m_vl = CommonParameters::Vlevel();

  if ((m_Nx % 2) == 1) {
    vout.crucial("Error at %s: Nx = %d, which must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // grid coordinate
  const int ipex = Communicator::ipe(0);
  const int ipey = Communicator::ipe(1);
  const int ipez = Communicator::ipe(2);
  const int ipet = Communicator::ipe(3);

  // node even/odd
  m_node_eo
    = (ipex * m_Nx + ipey * m_Ny + ipez * m_Nz + ipet * m_Nt) % 2;

  m_yzt_eo.resize(m_Ny * m_Nz * m_Nt);
  m_site_up.resize(m_Nx2 * m_Ny * m_Nz * m_Nt * 2);
  m_site_dn.resize(m_Nx2 * m_Ny * m_Nz * m_Nt * 2);

  for (int it = 0; it < m_Nt; ++it) {
    for (int iz = 0; iz < m_Nz; ++iz) {
      for (int iy = 0; iy < m_Ny; ++iy) {
        int it_global = it + ipet * m_Nt;
        int iz_global = iz + ipez * m_Nz;
        int iy_global = iy + ipey * m_Ny;
        m_yzt_eo[iy + m_Ny * (iz + m_Nz * it)]
          = (iy_global + iz_global + it_global) % 2;
      }
    }
  }

  for (int it = 0; it < m_Nt; ++it) {
    for (int iz = 0; iz < m_Nz; ++iz) {
      for (int iy = 0; iy < m_Ny; ++iy) {
        int iyzt = iy + m_Ny * (iz + m_Nz * it);
        for (int ix2 = 0; ix2 < m_Nx2; ++ix2) {
          int site2 = ix2 + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
          m_site_up[site2]
            = ((ix2 + m_yzt_eo[iyzt]) % m_Nx2) + m_Nx2 * iyzt;
          m_site_up[site2 + m_Nvol2]
            = ((ix2 + 1 - m_yzt_eo[iyzt]) % m_Nx2) + m_Nx2 * iyzt;
          m_site_dn[site2]
            = ((ix2 - 1 + m_yzt_eo[iyzt] + m_Nx2) % m_Nx2) + m_Nx2 * iyzt;
          m_site_dn[site2 + m_Nvol2]
            = ((ix2 - m_yzt_eo[iyzt] + m_Nx2) % m_Nx2) + m_Nx2 * iyzt;
        }
      }
    }
  }
}


//====================================================================
void Index_eo::convertField(Field& field_eo, const Field& field_lex)
{
  const int Nin = field_lex.nin();
  const int Nex = field_lex.nex();
  assert(field_lex.nvol() == m_Nvol);
  assert(field_eo.check_size(Nin, m_Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int ix    = site % m_Nx;
      int iyzt  = site / m_Nx;
      int iy    = iyzt % m_Ny;
      int izt   = iyzt / m_Ny;
      int iz    = izt % m_Nz;
      int it    = izt / m_Nz;
      int ix2   = ix / 2;
      int ieo   = (ix + iy + iz + it + m_node_eo) % 2;
      int site2 = this->site(ix2, iy, iz, it, ieo);
      for (int in = 0; in < Nin; ++in) {
        double vt = field_lex.cmp(in, site, ex);
        field_eo.set(in, site2, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Index_eo::convertField(Field& field_eo, const Field& field_lex,
                            const int ieo)
{
  const int Nin = field_lex.nin();
  const int Nex = field_lex.nex();
  assert(field_lex.nvol() == m_Nvol);
  assert(field_eo.check_size(Nin, m_Nvol2, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site2 = is; site2 < ns; ++site2) {
      int ix2  = site2 % m_Nx2;
      int iyzt = site2 / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int keo  = (iy + iz + it + m_node_eo + ieo) % 2;
      int ix   = 2 * ix2 + keo;
      int site = m_index_lex.site(ix, iy, iz, it);
      for (int in = 0; in < Nin; ++in) {
        double vt = field_lex.cmp(in, site, ex);
        field_eo.set(in, site2, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Index_eo::reverseField(Field& field_lex, const Field& field_eo,
                            const int ieo)
{
  const int Nin = field_lex.nin();
  const int Nex = field_lex.nex();
  assert(field_lex.nvol() == m_Nvol);
  assert(field_eo.check_size(Nin, m_Nvol2, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site2 = is; site2 < ns; ++site2) {
      int ix2  = site2 % m_Nx2;
      int iyzt = site2 / m_Nx2;
      int iy   = iyzt % m_Ny;
      int izt  = iyzt / m_Ny;
      int iz   = izt % m_Nz;
      int it   = izt / m_Nz;
      int keo  = (iy + iz + it + m_node_eo + ieo) % 2;
      int ix   = 2 * ix2 + keo;
      int site = m_index_lex.site(ix, iy, iz, it);
      for (int in = 0; in < Nin; ++in) {
        double vt = field_eo.cmp(in, site2, ex);
        field_lex.set(in, site, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Index_eo::reverseField(Field& field_lex, const Field& field_eo)
{
  const int Nin = field_lex.nin();
  const int Nex = field_lex.nex();
  assert(field_lex.nvol() == m_Nvol);
  assert(field_eo.check_size(Nin, m_Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int ix    = site % m_Nx;
      int iyzt  = site / m_Nx;
      int iy    = iyzt % m_Ny;
      int izt   = iyzt / m_Ny;
      int iz    = izt % m_Nz;
      int it    = izt / m_Nz;
      int ix2   = ix / 2;
      int ieo   = (ix + iy + iz + it + m_node_eo) % 2;
      int site2 = this->site(ix2, iy, iz, it, ieo);
      for (int in = 0; in < Nin; ++in) {
        double vt = field_eo.cmp(in, site2, ex);
        field_lex.set(in, site, ex, vt);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Index_eo::splitField(Field& field_e, Field& field_o,
                          const Field& field_eo)
{
  const int Nin = field_eo.nin();
  const int Nex = field_eo.nex();
  assert(field_eo.nvol() == m_Nvol);
  assert(field_e.check_size(Nin, m_Nvol2, Nex));
  assert(field_o.check_size(Nin, m_Nvol2, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site2 = is; site2 < ns; ++site2) {
      for (int in = 0; in < Nin; ++in) {
        double ve = field_eo.cmp(in, site2, ex);
        field_e.set(in, site2, ex, ve);
        double vo = field_eo.cmp(in, site2 + m_Nvol2, ex);
        field_o.set(in, site2, ex, vo);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Index_eo::mergeField(Field& field_eo,
                          const Field& field_e, const Field& field_o)
{
  const int Nin = field_eo.nin();
  const int Nex = field_eo.nex();
  assert(field_eo.nvol() == m_Nvol);
  assert(field_e.check_size(Nin, m_Nvol2, Nex));
  assert(field_o.check_size(Nin, m_Nvol2, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site2 = is; site2 < ns; ++site2) {
      for (int in = 0; in < Nin; ++in) {
        double ve = field_e.cmp(in, site2, ex);
        field_eo.set(in, site2, ex, ve);
        double vo = field_o.cmp(in, site2, ex);
        field_eo.set(in, site2 + m_Nvol2, ex, vo);
      }
    }
  }

#pragma omp barrier
}


//============================================================END=====
