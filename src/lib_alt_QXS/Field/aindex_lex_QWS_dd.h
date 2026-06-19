/*!
      @file    aindex_lex_QWS.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AINDEX_LEX_QWS_DD_INCLUDED
#define QXS_AINDEX_LEX_QWS_DD_INCLUDED

#include <string>

#include "lib_alt/alt_impl.h"

#include "lib/Parameters/commonParameters.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

namespace AIndex_lex_qxs_qws_dd {
  template<typename REALTYPE>
  inline int idx(const int in, const int Nin, const int ist,
                 const int Nx, const int Ny, const int Nxv,
                 const int Nyv_, const int Nstv, const int ex)
  {
    int ix   = ist % Nx;
    int Nxh  = Nx / 2;
    int deo  = ix / Nxh;
    int ixh  = ix % Nxh;
    int kx   = ixh % VLEND;
    int ixv  = ixh / VLEND;
    int iyzt = ist / Nx;
    int istv = ixv + Nxv * iyzt / 2;
    return kx + VLEND * (in + Nin * (istv + deo * Nstv / 2 + Nstv * ex));
  }


  template<>
  inline int idx<float>(const int in, const int Nin, const int ist,
                        const int Nx, const int Ny, const int Nxv,
                        const int Nyv_, const int Nstv, const int ex)
  {
    int ix   = ist % Nx;
    int Nxh  = Nx / 2;
    int deo  = ix / Nxh;
    int ixh  = ix % Nxh;
    int kx   = ixh % VLENS;
    int ixv  = ixh / VLENS;
    int iyzt = ist / Nx;
    int istv = ixv + Nxv * iyzt / 2;
    return kx + VLENS * (in + Nin * (istv + deo * Nstv / 2 + Nstv * ex));
  }


  template<typename REALTYPE>
  inline int set_Nxv(const int Nx) { return Nx / VLEND; }

  template<>
  inline int set_Nxv<float>(const int Nx) { return Nx / VLENS; }

  //template<typename REALTYPE>
  //  inline int set_Nyv(const int Ny){ return Ny/VLENYD; }

  //template<>
  //  inline int set_Nyv<float>(const int Ny){ return Ny/VLENYS; }
}

//! Lexical site index.

/*!
  This class defines lexicographycal site index for alternative
  code set: QXS2 version.
                                       [23 Feb 2021 H.Matsufuru]
  Special case with VLENYD=VLENYS=1 to use qws
                                       [20 Jul 2021 I.Kanamori]

*/

template<typename REALTYPE, Impl impl>
class AIndex_lex_QWS_dd;

template<typename REALTYPE>
class AIndex_lex_QWS_dd<REALTYPE, QXS>
{
 protected:
  int m_Nc, m_Nd, m_Ndf, m_Nvcd;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nzt, m_Nstv;

 public:
  AIndex_lex_QWS_dd()
  {
    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();
    init_common();
  }

  AIndex_lex_QWS_dd(int Nx, int Ny, int Nz, int Nt)
  {
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
    init_common();
  }

 private:
  void init_common()
  {
    m_Nc   = CommonParameters::Nc();
    m_Nd   = CommonParameters::Nd();
    m_Ndf  = 2 * m_Nc * m_Nc;
    m_Nvcd = 2 * m_Nc * m_Nd;
    m_Nst  = m_Nx * m_Ny * m_Nz * m_Nt;
    m_Nxv  = AIndex_lex_qxs_qws_dd::set_Nxv<REALTYPE>(m_Nx);
    m_Nyv  = m_Ny;
    m_Nzt  = m_Nz * m_Nt;
    m_Nstv = m_Nxv * m_Ny * m_Nzt;
    vout.paranoiac("AIndex_lex_QWS:  VLEND=%d, VLNES=%d\n", VLEND, VLENS);
  }

 public:

  int site(const int x, const int y, const int z, const int t) const
  { return x + m_Nxv * (y + m_Ny * (z + m_Nz * t)); }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    return AIndex_lex_qxs_qws_dd::idx<REALTYPE>(in, Nin, ist, m_Nx, m_Ny,
                                                m_Nxv, m_Nyv, m_Nstv, ex);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  { return idx(idf, m_Ndf, ist, ex); }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 2 * (ic2 + m_Nc * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 1 + 2 * (ic2 + m_Nc * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_SP(const int in, const int ist, const int ex) const
  { return idx(in, m_Nvcd, ist, ex); }

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
};

#endif
