/*!
      @file    aindex_coares_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate::  $
      @version $LastChangedRevision: 2229 $
 */

#ifndef ACCEL_AINDEX_COARAE_LEX_INCLUDED
#define ACCEL_AINDEX_COARAE_LEX_INCLUDED

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Field/aindex_coarse_lex_base.h"

// uses the same indexing as the fine lattice
#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

template<typename REALTYPE>
class AIndex_coarse_lex<REALTYPE, ACCEL> {
 protected:
  int m_Nc, m_Nd, m_Nin, m_Nin2, m_Ncol;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;
  int m_Nvol_pad;

 private:
  AIndex_coarse_lex(); // coare lattice size is always needed

 public:
  AIndex_coarse_lex(const int Nx, const int Ny, const int Nz,
                    const int Nt, const int Nc, const int Nd)
    {
     m_Nx = Nx;
     m_Ny = Ny;
     m_Nz = Nz;
     m_Nt = Nt;
     m_Nc = Nc;
     m_Nd = Nd;
     m_Ncol = Nc * Nd;
     m_Nin = 2 * Nc * Nd;
     m_Nin2 = 2 * m_Ncol * m_Ncol;
     m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;
     m_Nvol_pad = ceil_nwp(m_Nvol);
    }

  int site(const int x, const int y, const int z, const int t) const
  {  return m_Nx * (m_Ny * (m_Nz * t + z) + y) + x;  }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    int istx = ist + m_Nvol_pad * ex;
    return IDX2(Nin, in, istx);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  {  return idx(idf, m_Nin2, ist, ex);  }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 2*(ic1 + m_Ncol *ic2);
    return idx_G(idf, ist, ex); }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 1 + 2*(ic1 + m_Ncol *ic2);
    return idx_G(idf, ist, ex); }

  int idx_SP(const int in, const int ist, const int ex) const
  {  return idx(in, m_Nin, ist, ex);  }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  { int in = 2*(id + m_Nd *ic);
    return idx_SP(in, ist, ex); }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  { int in = 1 + 2*(id + m_Nd *ic);
    return idx_SP(in, ist, ex); }

};

#endif  // ACCEL_AINDEX_COARAE_LEX_INCLUDED
