/*!
        @file    aindex_coares_lex.h
        @brief
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2499 $
 */

#ifndef AINDEX_COARSE_LEX_QXS_H
#define AINDEX_COARSE_LEX_QXS_H

#include "lib_alt/Field/aindex_coarse_lex_base.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// uses the same indexing as the fine lattice
#include "lib_alt_QXS/Field/aindex_lex.h"

template<typename REALTYPE>
class AIndex_coarse_lex<REALTYPE, QXS> {
 protected:
  int m_Nc, m_Nd, m_Nin, m_Nin2, m_Ncol;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;
  int m_Nxv, m_Nyv, m_Nzt, m_Nstv;

 private:
  AIndex_coarse_lex(); // coare lattice size is always needed

 public:
  AIndex_coarse_lex(const int Nx, const int Ny, const int Nz, const int Nt, const int Nc, const int Nd)
  {
    m_Nx   = Nx;
    m_Ny   = Ny;
    m_Nz   = Nz;
    m_Nt   = Nt;
    m_Nc   = Nc;          // "color"
    m_Nd   = Nd;          // "spin"
    m_Ncol = Nc * Nd;
    m_Nin  = 2 * Nc * Nd; // "color" x "spin" x "complex"
    m_Nin2 = 2 * m_Ncol * m_Ncol;
    m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;
    // m_Nxv = Index_lex_alt_qxs::set_Nxv<REALTYPE>(m_Nx);
    // m_Nyv = Index_lex_alt_qxs::set_Nyv<REALTYPE>(m_Ny);
    m_Nxv  = AIndex_lex_qxs::set_Nxv<REALTYPE>(m_Nx);
    m_Nyv  = AIndex_lex_qxs::set_Nyv<REALTYPE>(m_Ny);
    m_Nzt  = m_Nz * m_Nt;
    m_Nstv = m_Nxv * m_Nyv * m_Nzt;
  }

  int site(const int x, const int y, const int z, const int t) const
  { return x + m_Nx * (y + m_Ny * (z + m_Nz * t)); }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    //return Index_lex_alt_qxs::idx<REALTYPE>(in, Nin, ist, m_Nx, m_Ny,
    return AIndex_lex_qxs::idx<REALTYPE>(in, Nin, ist, m_Nx, m_Ny,
                                         m_Nxv, m_Nyv, m_Nstv, ex);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  { return idx(idf, m_Nin2, ist, ex); }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 2 * (ic2 + m_Ncol * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 1 + 2 * (ic2 + m_Ncol * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_SP(const int in, const int ist, const int ex) const
  { return idx(in, m_Nin, ist, ex); }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 2 * (id + m_Nd * ic);
    return idx_SP(in, ist, ex);
  }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 1 + 2 * (id + m_Nd * ic);
    return idx_SP(in, ist, ex);
  }

  /*
  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    return in+ Nin*ist + Nin*m_Nvol*ex;
  }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  {
    int in=2*(id+m_Nd*ic);
    return in + m_Nin*ist + m_Nin*m_Nvol*ex;
  }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  {
    int in=2*(id+m_Nd*ic)+1;
    return in + m_Nin*ist + m_Nin*m_Nvol*ex;
  }
  */
};

#endif
