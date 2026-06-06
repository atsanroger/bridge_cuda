/*!
        @file    aindex_eo.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_AINDEX_EO_INCLUDED
#define QXS_AINDEX_EO_INCLUDED

#include <vector>

#include "lib_alt/Field/aindex_eo_base.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"


namespace AIndex_eo_qxs {
  template<typename REALTYPE>
  inline int idx(const int in, const int Nin, const int ist,
                 const int Nx2, const int Ny,
                 const int leo, const int Nvol2, const int ex)
  {
    int ist2   = ist / 2;
    int ieo    = (ist + leo) % 2;
    int offset = (ieo + 2 * ex) * Nin * Nvol2;
    int ix2    = ist2 % Nx2;
    int iy     = (ist2 / Nx2) % Ny;
    int kx2    = ix2 % VLENXD;
    int ky     = iy % VLENYD;
    int k      = kx2 + VLENXD * ky;
    int ix2v   = ix2 / VLENXD;
    int iyv    = iy / VLENYD;
    int izt    = ist2 / (Nx2 * Ny);
    int ist2v  = ix2v + (Nx2 / VLENXD) * (iyv + (Ny / VLENYD) * izt);
    return k + VLEND * (in + Nin * ist2v) + offset;
  }


  template<typename REALTYPE>
  inline int idxh(const int in, const int Nin, const int ist2,
                  const int Nx2, const int Ny,
                  const int Nvol2, const int ex)
  {
    int ix2   = ist2 % Nx2;
    int iy    = (ist2 / Nx2) % Ny;
    int kx2   = ix2 % VLENXD;
    int ky    = iy % VLENYD;
    int k     = kx2 + VLENXD * ky;
    int ix2v  = ix2 / VLENXD;
    int iyv   = iy / VLENYD;
    int izt   = ist2 / (Nx2 * Ny);
    int ist2v = ix2v + (Nx2 / VLENXD) * (iyv + (Ny / VLENYD) * izt);
    return k + VLEND * (in + Nin * ist2v) + Nin * Nvol2 * ex;
  }


  template<>
  inline int idx<float>(const int in, const int Nin, const int ist,
                        const int Nx2, const int Ny,
                        const int leo, const int Nvol2, const int ex)
  {
    int ist2   = ist / 2;
    int ieo    = (ist + leo) % 2;
    int offset = (ieo + 2 * ex) * Nin * Nvol2;
    int ix2    = ist2 % Nx2;
    int iy     = (ist2 / Nx2) % Ny;
    int kx2    = ix2 % VLENXS;
    int ky     = iy % VLENYS;
    int k      = kx2 + VLENXS * ky;
    int ix2v   = ix2 / VLENXS;
    int iyv    = iy / VLENYS;
    int izt    = ist2 / (Nx2 * Ny);
    int ist2v  = ix2v + (Nx2 / VLENXS) * (iyv + (Ny / VLENYS) * izt);
    return k + VLENS * (in + Nin * ist2v) + offset;
  }


  template<>
  inline int idxh<float>(const int in, const int Nin, const int ist2,
                         const int Nx2, const int Ny,
                         const int Nvol2, const int ex)
  {
    int ix2   = ist2 % Nx2;
    int iy    = (ist2 / Nx2) % Ny;
    int kx2   = ix2 % VLENXS;
    int ky    = iy % VLENYS;
    int k     = kx2 + VLENXS * ky;
    int ix2v  = ix2 / VLENXS;
    int iyv   = iy / VLENYS;
    int izt   = ist2 / (Nx2 * Ny);
    int ist2v = ix2v + (Nx2 / VLENXS) * (iyv + (Ny / VLENYS) * izt);
    return k + VLENS * (in + Nin * ist2v) + Nin * Nvol2 * ex;
  }

#ifdef USE_QXS_FP16
  template<>
  inline int idx<half>(const int in, const int Nin, const int ist,
                       const int Nx2, const int Ny,
                       const int leo, const int Nvol2, const int ex)
  {
    int ist2   = ist / 2;
    int ieo    = (ist + leo) % 2;
    int offset = (ieo + 2 * ex) * Nin * Nvol2;
    int ix2    = ist2 % Nx2;
    int iy     = (ist2 / Nx2) % Ny;
    int kx2    = ix2 % VLENXH;
    int ky     = iy % VLENYH;
    int k      = kx2 + VLENXH * ky;
    int ix2v   = ix2 / VLENXH;
    int iyv    = iy / VLENYH;
    int izt    = ist2 / (Nx2 * Ny);
    int ist2v  = ix2v + (Nx2 / VLENXH) * (iyv + (Ny / VLENYH) * izt);
    return k + VLENH * (in + Nin * ist2v) + offset;
  }

  template<>
  inline int idxh<half>(const int in, const int Nin, const int ist2,
                        const int Nx2, const int Ny,
                        const int Nvol2, const int ex)
  {
    int ix2   = ist2 % Nx2;
    int iy    = (ist2 / Nx2) % Ny;
    int kx2   = ix2 % VLENXH;
    int ky    = iy % VLENYH;
    int k     = kx2 + VLENXH * ky;
    int ix2v  = ix2 / VLENXH;
    int iyv   = iy / VLENYH;
    int izt   = ist2 / (Nx2 * Ny);
    int ist2v = ix2v + (Nx2 / VLENXH) * (iyv + (Ny / VLENYH) * izt);
    return k + VLENH * (in + Nin * ist2v) + Nin * Nvol2 * ex;
  }
#endif

}

//! Even-odd site index.

/*!
    This class defines even-odd site index for alternative
    implementation.
    For use of QXS version.
                                      [23 Feb 2021 H.Matsufuru]
*/
template<typename REALTYPE>
class AIndex_eo<REALTYPE, QXS>
{
 private:
  int Nx, Ny, Nz, Nt, Nvol;
  int Nx2, Nvol2;
  int Nc, Nd, Ndf, Nvcd;
  std::vector<int> Leo;
  Bridge::VerboseLevel m_vl;

  //! initial setup.
  void init();

 public:
  //! constructor.
  AIndex_eo() { init(); }

  int site(const int x, const int y, const int z, const int t) const
  {
    int ieo = (x + leo(y, z, t)) % 2;
    return (x / 2) + Nx2 * (y + Ny * (z + Nz * t)) + ieo * Nvol2;
  }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  {
    int ist2 = ist / 2;
    int leot = Leo[ist2 / Nx2];
    return AIndex_eo_qxs::idx<REALTYPE>(in, Nin, ist,
                                        Nx2, Ny, leot, Nvol2, ex);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  { return idx(idf, Ndf, ist, ex); }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 2 * (ic2 + Nc * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  {
    int idf = 1 + 2 * (ic2 + Nc * ic1);
    return idx_G(idf, ist, ex);
  }

  int idx_SP(const int in, const int ist, const int ex) const
  { return idx(in, Nvcd, ist, ex); }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 2 * (id + Nd * ic);
    return idx_SP(in, ist, ex);
  }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 1 + 2 * (id + Nd * ic);
    return idx_SP(in, ist, ex);
  }

  int idxh(const int in, const int Nin, const int ist2, const int ex) const
  {
    return AIndex_eo_qxs::idxh<REALTYPE>(in, Nin, ist2,
                                         Nx2, Ny, Nvol2, ex);
  }

  int idxh_SP(const int in, const int ist2, const int ex) const
  { return idxh(in, Nvcd, ist2, ex); }

  int idxh_SPr(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 2 * (id + Nd * ic);
    return idxh_SP(in, ist, ex);
  }

  int idxh_SPi(const int ic, const int id, const int ist, const int ex) const
  {
    int in = 1 + 2 * (id + Nd * ic);
    return idxh_SP(in, ist, ex);
  }

  int site(const int x2, const int y, const int z, const int t,
           const int ieo) const
  { return x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo; }

  int site(const int is, const int ieo) const
  { return is + Nvol2 * ieo; }

  int siteh(const int x2, const int y, const int z, const int t)
  const
  { return x2 + Nx2 * (y + Ny * (z + Nz * t)); }

  int leo(const int y, const int z, const int t) const
  { return Leo[y + Ny * (z + Nz * t)]; }

  int leo(const int iyzt) const { return Leo[iyzt]; }

  template<typename AFIELD>
  void split(AFIELD& v_e, AFIELD& v_o, const AFIELD& v, bool qdw = false);

  template<typename AFIELD>
  void merge(AFIELD& v, const AFIELD& v_e, const AFIELD& v_o, bool qdw = false);
};

#endif
