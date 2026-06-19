/*!
        @file    afield_4spinor-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef ACCEL_AFIELD_4SPINOR_INC_H
#define ACCEL_AFIELD_4SPINOR_INC_H

#include "lib/Tools/gammaMatrix.h"

//====================================================================
template <typename REALTYPE>
void mult_Gn(AField<REALTYPE>& u, const int exu,
             const AField<REALTYPE>& v, const int exv,
             const AField<REALTYPE>& w, const int exw)
{
  typedef AField<REALTYPE> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t> index;

  real_t* up = u.ptr(index.idx(0, NVCD, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NVCD, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NCD];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;
    int ix = VLEN * NCD  * site;

    load_vec(wt, &wp[ix], NCD);

    for(int ic2 = 0; ic2 < NC; ++ic2){
      int iv1 = iv + VLEN * NC * ic2;
      load_vec(vt, &vp[iv1], NC);
      for (int id = 0; id < ND; ++id){
        mult_Vnn(ut, vt, &wt[NC*id], NC);
        int ix2 = ix + VLEN * (ic2 + NC * id);
        save_vec(&up[ix2], &ut, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void mult_Gd(AField<REALTYPE>& u, const int exu,
             const AField<REALTYPE>& v, const int exv,
             const AField<REALTYPE>& w, const int exw)
{
  typedef AField<REALTYPE> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t> index;

  real_t* up = u.ptr(index.idx(0, NVCD, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NVCD, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NCD];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;
    int ix = VLEN * NCD  * site;

    load_vec(wt, &wp[ix], NCD);

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv1 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&vt[ic1], &vp[iv1], 1);
      }
      for(int id = 0; id < ND; ++id){
        mult_Vdn(ut, vt, &wt[NC*id], NC);
        int ix2 = ix + VLEN * (ic2 + NC * id);
        save_vec(&up[ix2], &ut, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void tensorProd_4spinor(AField<REALTYPE>& u, const int exu,
                        const AField<REALTYPE>& v, const int exv,
                        const AField<REALTYPE>& w, const int exw)
{
  typedef AField<REALTYPE> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NVCD, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NVCD, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt, wt;

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;
    int ix = VLEN * NCD  * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){

        int iv1 = iv + VLEN * (ic1 + NC * ic2);

        clear_vec(&ut, 1);
        for(int id = 0; id < ND; ++id){
          int ix1 = ix + VLEN * (ic1 + NC * id);
          int ix2 = ix + VLEN * (ic2 + NC * id);
          load_vec(&vt, &vp[ix1], 1);
          load_vec(&wt, &wp[ix2], 1);
          multadd_dn_vec(ut, vt, wt);
        }
        save_vec(&up[iv1], &ut, 1);

      }
    }
  }

}

//====================================================================
template <typename REALTYPE>
void mult_iGM(AField<REALTYPE>& v, const int exv,
              const GammaMatrix& gm,
              const AField<REALTYPE>& w, const int exw)
{
  typedef AField<REALTYPE> AFIELD;
  typedef REALTYPE real_t;

  assert(v.check_size(w));

  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  int id[ND];

  Index_lex_alt<real_t> index;

  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NVCD, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NVCD, 0, exw));

  real_t vigm1[VLEN];
  Vsimd_t vigm[ND];

  for (int s = 0; s < ND; ++s) {
    id[s] = gm.index(s);
    if(gm.index_c(s) == 0){  // GM component is pure real
      for(int k = 0; k < VLEN2; ++k){
        vigm1[2*k]   = 0.0;
        vigm1[2*k+1] = gm.value_r(s);
      }
    }else{                 // GM component is pure imaginary
      for(int k = 0; k < VLEN2; ++k){
        vigm1[2*k]   = gm.value_r(s);
        vigm1[2*k+1] = 0.0;
      }
    }
    load_vec(&vigm[s], &vigm1[0], 1);
  }

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t vt, wt;

  for(int site = is; site < ns; ++site){
    int ix = VLEN * NCD  * site;

    for(int ic = 0; ic < NC; ++ic){
      for(int s = 0; s < ND; ++s){
        int iv = ix + VLEN * (ic + NC * s);
        int iw = ix + VLEN * (ic + NC * id[s]);
        load_vec(&wt, &wp[iw], 1);
        mult_nn_vec(vt, vigm[s], wt);
        save_vec(&vp[iv], &vt, 1);
      }
    }

  }

}

//============================================================END=====
#endif
