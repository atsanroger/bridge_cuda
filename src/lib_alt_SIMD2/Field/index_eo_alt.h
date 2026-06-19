/*!
        @file    index_eo_alt.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef INDEX_EO_ALT_INCLUDED
#define INDEX_EO_ALT_INCLUDED

#include <vector>

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Field/index_eo_alt_base.h"

#include "lib_alt_SIMD2/inline/define_vlen.h"


namespace Index_eo_alt_simd2 {

  template<typename REALTYPE>
  inline int idx(const int in,  const int Nin,   const int ist,
		 const int leo, const int Nvol2, const int ex) {
    int ist2 = ist/2;
    int ieo = (ist + leo) % 2;
    int offset = (ieo + 2*ex)* Nvol2;
    return (in % 2) + 2*(ist2 % VLEND2)
           + VLEND*( (in/2) + (Nin/2)*((ist2 + offset)/VLEND2) );
  }

  template<typename REALTYPE>
  inline int idxh(const int in, const int Nin,
                  const int ist2, const int Nvol2, const int ex) {
    return (in % 2) + 2*(ist2 % VLEND2)
           + VLEND*( (in/2) + (Nin/2)*((ist2 + Nvol2*ex)/VLEND2) );
  }

  template<>
  inline int idx<float>(const int in,  const int Nin,   const int ist,
		 const int leo, const int Nvol2, const int ex) {
    int ist2 = ist/2;
    int ieo = (ist + leo) % 2;
    int offset = (ieo + 2*ex)* Nvol2;
    return (in % 2) + 2*(ist2 % VLENS2)
           + VLENS*( (in/2) + (Nin/2)*((ist2 + offset)/VLENS2) );
  }

  template<>
  inline int idxh<float>(const int in, const int Nin,
                  const int ist2, const int Nvol2, const int ex) {
    return (in % 2) + 2*(ist2 % VLENS2)
           + VLENS*( (in/2) + (Nin/2)*((ist2 + Nvol2*ex)/VLENS2) );
  }

}

//! Even-odd site index.

/*!
    This class defines even-odd site index for alternative
    implementation.
    For use of QWS version.
                                      [17 Jan 2017 H.Matsufuru]
*/
template<typename REALTYPE>
class Index_eo_alt<REALTYPE, SIMD2>{

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
  Index_eo_alt(){ init(); }

  int site(const int x, const int y, const int z, const int t) const
  { int ieo = (x + leo(y,z,t)) % 2;
    return (x/2) + Nx2 * (y + Ny * (z + Nz * t)) + ieo * Nvol2; }

  int idx(const int in, const int Nin, const int ist, const int ex) const
  { int ist2 = ist/2;
    int leo = Leo[ist2/Nx2];
    return Index_eo_alt_simd2::idx<REALTYPE>(in, Nin, ist, leo, Nvol2, ex);
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
  { return Index_eo_alt_simd2::idxh<REALTYPE>(in, Nin, ist2, Nvol2, ex); }

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

  template <typename AFIELD>
  void split(AFIELD& v_e, AFIELD& v_o, const AFIELD& v);

  template <typename AFIELD>
  void merge(AFIELD& v, const AFIELD& v_e, const AFIELD& v_o);

};

#endif
