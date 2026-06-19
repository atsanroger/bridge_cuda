/*!
        @file    $Id:: index_lex_alt.h#$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$

        @version $LastChangedRevision: 2160 $
*/


#ifndef INDEX_LEX_ALT_SIMD2_INCLUDED
#define INDEX_LEX_ALT_SIMD2_INCLUDED

#include <string>

#include  "lib_alt/Field/index_lex_alt_base.h"

#include  "lib_alt_SIMD2/inline/define_vlen.h"

#include  "lib/Parameters/commonParameters.h"


namespace Index_lex_alt_simd2 {

template<typename REALTYPE>
inline int idx(const int in, const int Nin,
                     const int ist, const int Nvol, const int ex)
 {  return (in % 2) + 2*(ist % VLEND2) +
           VLEND*( (in/2) + (Nin/2)*((ist + Nvol*ex)/VLEND2) ); }



template<>
inline int idx<float>(const int in, const int Nin,
                     const int ist, const int Nvol, const int ex)
 {  return (in % 2) + 2*(ist % VLENS2) +
           VLENS*( (in/2) + (Nin/2)*((ist + Nvol*ex)/VLENS2) ); }

}

//! Lexical site index.

/*!
  This class defines lexicographycal site index for alternative
  code set.
                                       [06 Nov 2016 H.Matsufuru]
*/
template<typename REALTYPE>
class Index_lex_alt<REALTYPE, SIMD2> {
 protected:
  int m_Nc, m_Nd, m_Ndf, m_Nvcd;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;

 public:
  Index_lex_alt() {
    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    m_Ndf  = 2 * m_Nc * m_Nc;
    m_Nvcd = 2 * m_Nc * m_Nd;
    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();
    m_Nvol = CommonParameters::Nvol();
   }

  Index_lex_alt(int Nx, int Ny, int Nz, int Nt)
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
  }

  int site(const int x, const int y, const int z, const int t) const
  {  return m_Nx * (m_Ny * (m_Nz * t + z) + y) + x;  }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {  return Index_lex_alt_simd2::idx<REALTYPE>(in, Nin, ist, m_Nvol, ex); }

  int idx_G(const int idf, const int ist, const int ex) const
  {  return idx(idf, m_Ndf, ist, ex);  }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 2*(ic1 + m_Nc *ic2);
    return idx(idf, m_Ndf, ist, ex); }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 1 + 2*(ic1 + m_Nc *ic2);
    return idx(idf, m_Ndf, ist, ex); }

  int idx_SP(const int in, const int ist, const int ex) const
  {  return idx(in, m_Nvcd, ist, ex);  }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  { int in = 2*(ic + m_Nc *id);
    return idx_SP(in, ist, ex); }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  { int in = 1 + 2*(ic + m_Nc *id);
    return idx_SP(in, ist, ex); }

};
#endif
