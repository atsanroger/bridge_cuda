/*!
        @file    index_lex_alt.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/


#ifndef INDEX_LEX_ALT_VECTOR_INCLUDED
#define INDEX_LEX_ALT_VECTOR_INCLUDED

#include <string>

#include "lib_alt/Field/aindex_lex_base.h"

#include "lib/Parameters/commonParameters.h"

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

//! Lexical site index.

/*!
  This class defines lexicographical site index for alternative
  code set.
  For the Vector implementation, the index is actually defined
  in define_index.h as the macro definition IDX2.
                                       [13 Aug 2019 H.Matsufuru]
*/
template<typename REALTYPE>
class AIndex_lex<REALTYPE, VECTOR> {
 protected:
  int m_Nc, m_Nd, m_Ndf, m_Nvcd;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol;        //!< lattice volume (number of sites).
  int m_Nvol_pad;    //!< Nvol after padding.

 public:
  AIndex_lex() {
    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    m_Ndf  = 2 * m_Nc * m_Nc;
    m_Nvcd = 2 * m_Nc * m_Nd;
    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();
    m_Nvol = CommonParameters::Nvol();
    m_Nvol_pad = ceil_nwp(m_Nvol);
   }

  AIndex_lex(int Nx, int Ny, int Nz, int Nt)
  {
    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    m_Ndf  = 2 * m_Nc * m_Nc;
    m_Nvcd = 2 * m_Nc * m_Nd;
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
    m_Nvol = Nx * Ny * Nz * Nt;
    m_Nvol_pad = ceil_nwp(m_Nvol);
  }

  int site(const int x, const int y, const int z, const int t) const
  {  return m_Nx * (m_Ny * (m_Nz * t + z) + y) + x;  }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    return IDXV(Nin, m_Nvol, in, ist, ex);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  {  return idx(idf, NDF, ist, ex);  }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 2*(ic1 + NC *ic2);
    return idx(idf, NDF, ist, ex); }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 1 + 2*(ic1 + NC *ic2);
    return idx(idf, NDF, ist, ex); }

  int idx_SP(const int in, const int ist, const int ex) const
  {  return idx(in, NVCD, ist, ex);  }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  { int in = 2*(ic + NC *id);
    return idx(in, NVCD, ist, ex); }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  { int in = 1 + 2*(ic + NC *id);
    return idx(in, NVCD, ist, ex); }

};

#endif
