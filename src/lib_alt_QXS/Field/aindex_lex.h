/*!
      @file    aindex_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
      @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_AINDEX_LEX_INCLUDED
#define QXS_AINDEX_LEX_INCLUDED

#include <string>

#include "lib_alt/Field/aindex_lex_base.h"

#include "lib/Parameters/commonParameters.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

namespace AIndex_lex_qxs {
  template<typename REALTYPE>
  inline int idx(const int in, const int Nin, const int ist,
                 const int Nx, const int Ny, const int Nxv,
                 const int Nyv, const int Nstv, const int ex)
  {
    int ix   = ist % Nx;
    int iy   = (ist / Nx) % Ny;
    int kx   = ix % VLENXD;
    int ky   = iy % VLENYD;
    int k    = kx + VLENXD * ky;
    int ixv  = ix / VLENXD;
    int iyv  = iy / VLENYD;
    int izt  = ist / (Nx * Ny);
    int istv = ixv + Nxv * (iyv + Nyv * izt);
    return k + VLEND * (in + Nin * (istv + Nstv * ex));
  }

  template<>
  inline int idx<float>(const int in, const int Nin, const int ist,
                        const int Nx, const int Ny, const int Nxv,
                        const int Nyv, const int Nstv, const int ex)
  {
    int ix   = ist % Nx;
    int iy   = (ist / Nx) % Ny;
    int kx   = ix % VLENXS;
    int ky   = iy % VLENYS;
    int k    = kx + VLENXS * ky;
    int ixv  = ix / VLENXS;
    int iyv  = iy / VLENYS;
    int izt  = ist / (Nx * Ny);
    int istv = ixv + Nxv * (iyv + Nyv * izt);
    return k + VLENS * (in + Nin * (istv + Nstv * ex));
  }

  template<typename REALTYPE>
  inline int set_Nxv(const int Nx) { return Nx / VLENXD; }

  template<>
  inline int set_Nxv<float>(const int Nx) { return Nx / VLENXS; }

  template<typename REALTYPE>
  inline int set_Nyv(const int Ny) { return Ny / VLENYD; }

  template<>
  inline int set_Nyv<float>(const int Ny) { return Ny / VLENYS; }

#ifdef USE_QXS_FP16
  template<>
  inline int idx<half>(const int in, const int Nin, const int ist,
                       const int Nx, const int Ny, const int Nxv,
                       const int Nyv, const int Nstv, const int ex)
  {
    int ix   = ist % Nx;
    int iy   = (ist / Nx) % Ny;
    int kx   = ix % VLENXH;
    int ky   = iy % VLENYH;
    int k    = kx + VLENXH * ky;
    int ixv  = ix / VLENXH;
    int iyv  = iy / VLENYH;
    int izt  = ist / (Nx * Ny);
    int istv = ixv + Nxv * (iyv + Nyv * izt);
    return k + VLENH * (in + Nin * (istv + Nstv * ex));
  }

  template<>
  inline int set_Nxv<half>(const int Nx) { return Nx / VLENXH; }

  template<>
  inline int set_Nyv<half>(const int Ny) { return Ny / VLENYH; }
#endif
}

//! Lexical site index.

/*!
  This class defines lexicographycal site index for alternative
  code set: QXS version.
                                       [23 Feb 2021 H.Matsufuru]
*/
template<typename REALTYPE>
class AIndex_lex<REALTYPE, QXS>
{
 protected:
  int m_Nc, m_Nd, m_Ndf, m_Nvcd;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nzt, m_Nstv;

 public:
  AIndex_lex()
  {
    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();
    init_common();
  }

  AIndex_lex(int Nx, int Ny, int Nz, int Nt)
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
    m_Nxv  = AIndex_lex_qxs::set_Nxv<REALTYPE>(m_Nx);
    m_Nyv  = AIndex_lex_qxs::set_Nyv<REALTYPE>(m_Ny);
    m_Nzt  = m_Nz * m_Nt;
    m_Nstv = m_Nxv * m_Nyv * m_Nzt;
  }

 public:

  int site(const int x, const int y, const int z, const int t) const
  { return x + m_Nx * (y + m_Ny * (z + m_Nz * t)); }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    return AIndex_lex_qxs::idx<REALTYPE>(in, Nin, ist, m_Nx, m_Ny,
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
