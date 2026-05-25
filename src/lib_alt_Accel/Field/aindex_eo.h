/*!
        @file    aindex_eo.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-09-10 17:39:07 #$
        @version $LastChangedRevision: 2537 $
*/

#ifndef ACCEL_AINDEX_EO_INCLUDED
#define ACCEL_AINDEX_EO_INCLUDED

#include "lib_alt/Field/aindex_eo_base.h"

#include <vector>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/aindex_lex.h"

namespace AIndex_eo_accel {

  template<typename REALTYPE>
  inline int idx(const int in,  const int Nin,   const int ist,
		 const int leo, const int Nvol2, const int ex) {
    int ist2 = ist/2;
    int ieo = (ist + leo) % 2;
    int offset = (ieo + 2*ex)* Nvol2;
    return (ist2 % NWP) + NWP*( in + Nin*( (ist2 + offset)/NWP) );
  }

  template<typename REALTYPE>
  inline int idxh(const int in, const int Nin,
                  const int ist2, const int Nvol2, const int ex) {
    return (ist2 % NWP) + NWP*( in + Nin*( (ist2 + Nvol2*ex)/NWP) );
  }

  template<>
  inline int idx<float>(const int in,  const int Nin,   const int ist,
		 const int leo, const int Nvol2, const int ex) {
    int ist2 = ist/2;
    int ieo = (ist + leo) % 2;
    int offset = (ieo + 2*ex)* Nvol2;
    return (ist2 % NWP) + NWP*( in + Nin*( (ist2 + offset)/NWP) );
  }

  template<>
  inline int idxh<float>(const int in, const int Nin,
                  const int ist2, const int Nvol2, const int ex) {
    return (ist2 % NWP) + NWP*( in + Nin*( (ist2 + Nvol2*ex)/NWP) );
  }

}

//! Even-odd site index.

/*!
    This class defines even-odd site index for alternative
    implementation.
    Accel version.
                                      [19 Jun 2019 H.Matsufuru]
*/
template<typename REALTYPE>
class AIndex_eo<REALTYPE, ACCEL>{

 public:
  static const std::string class_name;

 private:
  int Nx, Ny, Nz, Nt, Nvol;
  int Nx2, Nvol2;
  int Nc, Nd, Ndf, Nvcd;
  int m_ieo_origin;    //!< parity of local origin site
  int m_Nsize[4];
  int m_Nsize2[4];
  std::vector<int> Leo;
  Bridge::VerboseLevel m_vl;

  //! initial setup.
  void init();

 public:
  //! constructor.
  AIndex_eo(){ init(); }

  int site(const int x, const int y, const int z, const int t) const
  { int ieo = (x + leo(y,z,t)) % 2;
    return (x/2) + Nx2 * (y + Ny * (z + Nz * t)) + ieo * Nvol2; }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  { int ist2 = ist/2;
    int leo = Leo[ist2/Nx2];
    return AIndex_eo_accel::idx<REALTYPE>(in, Nin, ist, leo, Nvol2, ex);
  }

  int idx_G(const int idf, const int ist, const int ex) const
  { return idx(idf, Ndf, ist, ex);  }

  int idx_Gr(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 2*(ic1 + Nc * ic2);
    return idx(idf, Ndf, ist, ex);   }

  int idx_Gi(const int ic1, const int ic2, const int ist, const int ex) const
  { int idf = 1 + 2*(ic1 + Nc * ic2);
    return idx(idf, Ndf, ist, ex);   }

  int idx_SP(const int in, const int ist, const int ex) const
  { return idx(in, Nvcd, ist, ex);  }

  int idxh(const int in, const int Nin, const int ist2, const int ex) const
  { return AIndex_eo_accel::idxh<REALTYPE>(in, Nin, ist2, Nvol2, ex); }

  int idx_SPr(const int ic, const int id, const int ist, const int ex) const
    { int in = 2*(ic + Nc * id);
      return idx_SP(in, ist, ex); }

  int idx_SPi(const int ic, const int id, const int ist, const int ex) const
    { int in = 1 + 2*(ic + Nc * id);
      return idx_SP(in, ist, ex); }

  int idxh_SP(const int in, const int ist2, const int ex) const
  { return idxh(in, Nvcd, ist2, ex); }

  int idxh_SPr(const int ic, const int id, const int ist, const int ex) const
    { int in = 2*(ic + Nc * id);
      return idxh_SP(in, ist, ex); }

  int idxh_SPi(const int ic, const int id, const int ist, const int ex) const
    { int in = 1 + 2*(ic + Nc * id);
      return idxh_SP(in, ist, ex); }

  int idxh_Gr(const int ic1, const int ic2, const int ist, const int ex) const
    { int in = 2*(ic1 + Nc * ic2);
      return idxh(in, Ndf, ist, ex); }

  int idxh_Gi(const int ic1, const int ic2, const int ist, const int ex) const
    { int in = 1 + 2*(ic1 + Nc * ic2);
      return idxh(in, Ndf, ist, ex); }

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

  int ieo_origin() const
  { return m_ieo_origin; }

  template <typename AFIELD>
  void split(AFIELD& v_e, AFIELD& v_o, const AFIELD& v, bool qdw = false);

  template <typename AFIELD>
  void split(AFIELD& v_e, const int ex_e,
             AFIELD& v_o, const int ex_o,
             const AFIELD& v, const int ex);

  template <typename AFIELD>
  void split_gauge(AFIELD& ueo, const AFIELD& ulex);

  template <typename AFIELD>
  void merge(AFIELD& v, const AFIELD& v_e, const AFIELD& v_o, bool qdw = false);

};

#endif
