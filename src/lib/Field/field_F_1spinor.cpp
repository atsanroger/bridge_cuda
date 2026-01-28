/*!
        @file    field_F_1spinor.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-12-07 20:54:14 #$

        @version $LastChangedRevision: 2555 $
*/

#include "Field/field_F_1spinor.h"
#include "Field/field_thread-inc.h"
#include "ResourceManager/threadManager.h"

#if defined USE_GROUP_SU3
#include "Imp/field_F_imp_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "Imp/field_F_imp_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "Imp/field_F_imp_SU_N-inc.h"
#endif

//====================================================================
void mult_Field_Gn(Field_F_1spinor& y, const int ex,
                   const Field_G& U, int ex1,
                   const Field_F_1spinor& x, int ex2)
{
  int Nvol = y.nvol();
  assert(U.nvol() == Nvol);
  assert(x.nvol() == Nvol);
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());

  const int Nd   = 1;
  const int Nc   = CommonParameters::Nc();
  const int Nc2  = 2 * CommonParameters::Nc();
  const int Ncd2 = Nc2 * Nd;
  const int Ndf  = Nc2 * Nc;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *v = y.ptr(0, is, ex);
  const double *g = U.ptr(0, is, ex1);
  const double *w = x.ptr(0, is, ex2);

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;
    int iv = Ncd2 * site;
    for (int s = 0; s < Nd; ++s) {
      for (int ic = 0; ic < Nc; ++ic) {
        int ig2 = ic * Nc2 + ig;
        int iv2 = s * Nc2 + iv;
        v[2 * ic + iv2]     = mult_Gn_r(&g[ig2], &w[iv2], Nc);
        v[2 * ic + 1 + iv2] = mult_Gn_i(&g[ig2], &w[iv2], Nc);
      }
    }
  }
}


//====================================================================
void mult_Field_Gd(Field_F_1spinor& y, const int ex,
                   const Field_G& U, int ex1,
                   const Field_F_1spinor& x, int ex2)
{
  int Nvol = y.nvol();
  assert(U.nvol() == Nvol);
  assert(x.nvol() == Nvol);
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());

  const int Nd   = 1;
  const int Nc   = CommonParameters::Nc();
  const int Nc2  = 2 * CommonParameters::Nc();
  const int Ncd2 = Nc2 * Nd;
  const int Ndf  = Nc2 * Nc;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *v = y.ptr(0, is, ex);
  const double *g = U.ptr(0, is, ex1);
  const double *w = x.ptr(0, is, ex2);

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;
    int iv = Ncd2 * site;
    for (int s = 0; s < Nd; ++s) {
      for (int ic = 0; ic < Nc; ++ic) {
        int ig2 = ic * 2 + ig;
        int iv2 = s * Nc2 + iv;
        v[2 * ic + iv2]     = mult_Gd_r(&g[ig2], &w[iv2], Nc);
        v[2 * ic + 1 + iv2] = mult_Gd_i(&g[ig2], &w[iv2], Nc);
      }
    }
  }
}


//============================================================END=====
