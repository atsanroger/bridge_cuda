/*!
        @file    vsimd_Wilson_SU3_double-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_WILSON_SU3_DOUBLE_INC_INCLUDED
#define VSIMD_WILSON_SU3_DOUBLE_INC_INCLUDED

#include "lib_alt_SIMD2/inline/define_params_SU3.h"

#define  ID1     0
#define  ID2     3
#define  ID3     6
#define  ID4     9

namespace {
  void check_Nc(int Nc)
  {
    vout.general("Fopr_Wilson_impl: implementation for SU(3).\n");
    if(Nc != NC){
      vout.crucial("Nc %d != macro definition: %d\n", Nc, NC);
      exit(EXIT_FAILURE);
    }
  }


  inline void mult_uv(Vsimd_t& v, Vsimd_t *g, Vsimd_t *w, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] =
          g[0].val[kr] * w[0].val[kr] - g[0].val[ki] * w[0].val[ki]
        + g[1].val[kr] * w[1].val[kr] - g[1].val[ki] * w[1].val[ki]
        + g[2].val[kr] * w[2].val[kr] - g[2].val[ki] * w[2].val[ki];
      v.val[ki] =
          g[0].val[kr] * w[0].val[ki] + g[0].val[ki] * w[0].val[kr]
        + g[1].val[kr] * w[1].val[ki] + g[1].val[ki] * w[1].val[kr]
        + g[2].val[kr] * w[2].val[ki] + g[2].val[ki] * w[2].val[kr];
    }
  }

  inline void mult_udagv(Vsimd_t& v, Vsimd_t *g, Vsimd_t *w, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] =
          g[0].val[kr] * w[0].val[kr] + g[0].val[ki] * w[0].val[ki]
        + g[3].val[kr] * w[1].val[kr] + g[3].val[ki] * w[1].val[ki]
        + g[6].val[kr] * w[2].val[kr] + g[6].val[ki] * w[2].val[ki];
      v.val[ki] =
          g[0].val[kr] * w[0].val[ki] - g[0].val[ki] * w[0].val[kr]
        + g[3].val[kr] * w[1].val[ki] - g[3].val[ki] * w[1].val[kr]
        + g[6].val[kr] * w[2].val[ki] - g[6].val[ki] * w[2].val[kr];
    }
  }

  template <typename REALTYPE>
  inline void mult_uv(REALTYPE *v, REALTYPE *g, REALTYPE *w, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr] =  g[kr+VLEN*0] * w[kr+VLEN*0] - g[ki+VLEN*0] * w[ki+VLEN*0]
             + g[kr+VLEN*1] * w[kr+VLEN*1] - g[ki+VLEN*1] * w[ki+VLEN*1]
             + g[kr+VLEN*2] * w[kr+VLEN*2] - g[ki+VLEN*2] * w[ki+VLEN*2];
      v[ki] =  g[kr+VLEN*0] * w[ki+VLEN*0] + g[ki+VLEN*0] * w[kr+VLEN*0]
             + g[kr+VLEN*1] * w[ki+VLEN*1] + g[ki+VLEN*1] * w[kr+VLEN*1]
             + g[kr+VLEN*2] * w[ki+VLEN*2] + g[ki+VLEN*2] * w[kr+VLEN*2];
    }
  }

  template <typename REALTYPE>
  inline void mult_udagv(REALTYPE *v, REALTYPE *g, REALTYPE *w, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr] =  g[kr+VLEN*0] * w[kr+VLEN*0] + g[ki+VLEN*0] * w[ki+VLEN*0]
             + g[kr+VLEN*3] * w[kr+VLEN*1] + g[ki+VLEN*3] * w[ki+VLEN*1]
             + g[kr+VLEN*6] * w[kr+VLEN*2] + g[ki+VLEN*6] * w[ki+VLEN*2];
      v[ki] =  g[kr+VLEN*0] * w[ki+VLEN*0] - g[ki+VLEN*0] * w[kr+VLEN*0]
             + g[kr+VLEN*3] * w[ki+VLEN*1] - g[ki+VLEN*3] * w[kr+VLEN*1]
             + g[kr+VLEN*6] * w[ki+VLEN*2] - g[ki+VLEN*6] * w[kr+VLEN*2];
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_xp(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k + 1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID4)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID4)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID3)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
    inline void load_set_sp2_xp(Vsimd_t *vt1, Vsimd_t *vt2,
                                REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k + 1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID4)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID4)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID3)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_xp1(REALTYPE *vt, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
      vt[icr]       = w[icr + 2*ID1] - w[ici + 2*ID4];
      vt[ici]       = w[ici + 2*ID1] + w[icr + 2*ID4];
      vt[icr + NVC] = w[icr + 2*ID2] - w[ici + 2*ID3];
      vt[ici + NVC] = w[ici + 2*ID2] + w[icr + 2*ID3];
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_xm(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID4)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID4)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID3)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_xm(Vsimd_t *vt1, Vsimd_t *vt2,
                              REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID4)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID4)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID3)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_yp(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID4)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID4)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID3)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_yp(Vsimd_t *vt1, Vsimd_t *vt2,
                              REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID4)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID4)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID3)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_ym(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID4)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID4)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID3)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_ym(Vsimd_t *vt1, Vsimd_t *vt2,
                              REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID4)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID4)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID3)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID3)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_zp(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID3)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID3)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID4)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_zp(Vsimd_t *vt1, Vsimd_t *vt2,
                              REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] - w[ki + VLEN*(ic+ID3)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] + w[kr + VLEN*(ic+ID3)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] + w[ki + VLEN*(ic+ID4)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] - w[kr + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_zm(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = w[kr + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID3)];
       vt1[ki + VLEN*ic] = w[ki + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID3)];
       vt2[kr + VLEN*ic] = w[kr + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID4)];
       vt2[ki + VLEN*ic] = w[ki + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_zm(Vsimd_t *vt1, Vsimd_t *vt2,
                              REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = w[kr + VLEN*(ic+ID1)] + w[ki + VLEN*(ic+ID3)];
       vt1[ic].val[ki] = w[ki + VLEN*(ic+ID1)] - w[kr + VLEN*(ic+ID3)];
       vt2[ic].val[kr] = w[kr + VLEN*(ic+ID2)] - w[ki + VLEN*(ic+ID4)];
       vt2[ic].val[ki] = w[ki + VLEN*(ic+ID2)] + w[kr + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_tp_dirac(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = 2.0 * w[kr + VLEN*(ic+ID3)];
       vt1[ki + VLEN*ic] = 2.0 * w[ki + VLEN*(ic+ID3)];
       vt2[kr + VLEN*ic] = 2.0 * w[kr + VLEN*(ic+ID4)];
       vt2[ki + VLEN*ic] = 2.0 * w[ki + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_tp_dirac(Vsimd_t *vt1, Vsimd_t *vt2,
                                    REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = 2.0 * w[kr + VLEN*(ic+ID3)];
       vt1[ic].val[ki] = 2.0 * w[ki + VLEN*(ic+ID3)];
       vt2[ic].val[kr] = 2.0 * w[kr + VLEN*(ic+ID4)];
       vt2[ic].val[ki] = 2.0 * w[ki + VLEN*(ic+ID4)];
     }
    }
  }

  template <typename REALTYPE>
  inline void set_sp2_tm_dirac(REALTYPE *vt1, REALTYPE *vt2,
                         REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[kr + VLEN*ic] = 2.0 * w[kr + VLEN*(ic+ID1)];
       vt1[ki + VLEN*ic] = 2.0 * w[ki + VLEN*(ic+ID1)];
       vt2[kr + VLEN*ic] = 2.0 * w[kr + VLEN*(ic+ID2)];
       vt2[ki + VLEN*ic] = 2.0 * w[ki + VLEN*(ic+ID2)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_set_sp2_tm_dirac(Vsimd_t *vt1, Vsimd_t *vt2,
                                    REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN2; ++k){
       int kr = 2*k;
       int ki = 2*k+1;
       vt1[ic].val[kr] = 2.0 * w[kr + VLEN*(ic+ID1)];
       vt1[ic].val[ki] = 2.0 * w[ki + VLEN*(ic+ID1)];
       vt2[ic].val[kr] = 2.0 * w[kr + VLEN*(ic+ID2)];
       vt2[ic].val[ki] = 2.0 * w[ki + VLEN*(ic+ID2)];
     }
    }
  }

  inline void set_sp4_xp(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] +=  wt2.val[ki];
      v[ID3].val[ki] += -wt2.val[kr];
      v[ID4].val[kr] +=  wt1.val[ki];
      v[ID4].val[ki] += -wt1.val[kr];
   }
 }

  inline void set_sp4_xm(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] += -wt2.val[ki];
      v[ID3].val[ki] +=  wt2.val[kr];
      v[ID4].val[kr] += -wt1.val[ki];
      v[ID4].val[ki] +=  wt1.val[kr];
    }
  }

  inline void set_sp4_yp(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] += -wt2.val[kr];
      v[ID3].val[ki] += -wt2.val[ki];
      v[ID4].val[kr] +=  wt1.val[kr];
      v[ID4].val[ki] +=  wt1.val[ki];
    }
  }

  inline void set_sp4_ym(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] +=  wt2.val[kr];
      v[ID3].val[ki] +=  wt2.val[ki];
      v[ID4].val[kr] += -wt1.val[kr];
      v[ID4].val[ki] += -wt1.val[ki];
    }
  }

  inline void set_sp4_zp(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] +=  wt1.val[ki];
      v[ID3].val[ki] += -wt1.val[kr];
      v[ID4].val[kr] += -wt2.val[ki];
      v[ID4].val[ki] +=  wt2.val[kr];
    }
  }

  inline void set_sp4_zm(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] +=  wt1.val[kr];
      v[ID1].val[ki] +=  wt1.val[ki];
      v[ID2].val[kr] +=  wt2.val[kr];
      v[ID2].val[ki] +=  wt2.val[ki];
      v[ID3].val[kr] += -wt1.val[ki];
      v[ID3].val[ki] +=  wt1.val[kr];
      v[ID4].val[kr] +=  wt2.val[ki];
      v[ID4].val[ki] += -wt2.val[kr];
    }
  }

  inline void set_sp4_tp_dirac(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID3].val[kr] += wt1.val[kr];
      v[ID3].val[ki] += wt1.val[ki];
      v[ID4].val[kr] += wt2.val[kr];
      v[ID4].val[ki] += wt2.val[ki];
    }
  }

  inline void set_sp4_tm_dirac(Vsimd_t *v,
                         Vsimd_t& wt1, Vsimd_t& wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[ID1].val[kr] += wt1.val[kr];
      v[ID1].val[ki] += wt1.val[ki];
      v[ID2].val[kr] += wt2.val[kr];
      v[ID2].val[ki] += wt2.val[ki];
    }
  }


  template<typename REALTYPE>
  inline void set_aPp5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = a * (w[ID1+ic].val[k] + w[ID3+ic].val[k]);
       v[ID2+ic].val[k] = a * (w[ID2+ic].val[k] + w[ID4+ic].val[k]);
       v[ID3+ic].val[k] = a * (w[ID3+ic].val[k] + w[ID1+ic].val[k]);
       v[ID4+ic].val[k] = a * (w[ID4+ic].val[k] + w[ID2+ic].val[k]);
     }
    }
  }

  template<typename REALTYPE>
  inline void set_aPm5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = a * (w[ID1+ic].val[k] - w[ID3+ic].val[k]);
       v[ID2+ic].val[k] = a * (w[ID2+ic].val[k] - w[ID4+ic].val[k]);
       v[ID3+ic].val[k] = a * (w[ID3+ic].val[k] - w[ID1+ic].val[k]);
       v[ID4+ic].val[k] = a * (w[ID4+ic].val[k] - w[ID2+ic].val[k]);
     }
    }
  }

  template<typename REALTYPE>
  inline void add_aPp5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] += a * (w[ID1+ic].val[k] + w[ID3+ic].val[k]);
       v[ID2+ic].val[k] += a * (w[ID2+ic].val[k] + w[ID4+ic].val[k]);
       v[ID3+ic].val[k] += a * (w[ID3+ic].val[k] + w[ID1+ic].val[k]);
       v[ID4+ic].val[k] += a * (w[ID4+ic].val[k] + w[ID2+ic].val[k]);
     }
    }
  }

  template<typename REALTYPE>
  inline void add_aPm5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] += a * (w[ID1+ic].val[k] - w[ID3+ic].val[k]);
       v[ID2+ic].val[k] += a * (w[ID2+ic].val[k] - w[ID4+ic].val[k]);
       v[ID3+ic].val[k] += a * (w[ID3+ic].val[k] - w[ID1+ic].val[k]);
       v[ID4+ic].val[k] += a * (w[ID4+ic].val[k] - w[ID2+ic].val[k]);
     }
    }
  }

  template<typename REALTYPE>
  inline void set_aPp5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = 2.0 * a * w[ID1+ic].val[k];
       v[ID2+ic].val[k] = 2.0 * a * w[ID2+ic].val[k];
       v[ID3+ic].val[k] = 0.0;
       v[ID4+ic].val[k] = 0.0;
     }
    }
  }

  template<typename REALTYPE>
  inline void set_aPm5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = 0.0;
       v[ID2+ic].val[k] = 0.0;
       v[ID3+ic].val[k] = 2.0 * a * w[ID3+ic].val[k];
       v[ID4+ic].val[k] = 2.0 * a * w[ID4+ic].val[k];
     }
    }
  }

  template<typename REALTYPE>
  inline void add_aPp5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] += 2.0 * a * w[ID1+ic].val[k];
       v[ID2+ic].val[k] += 2.0 * a * w[ID2+ic].val[k];
     }
    }
  }

  template<typename REALTYPE>
  inline void add_aPm5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID3+ic].val[k] += 2.0 * a * w[ID3+ic].val[k];
       v[ID4+ic].val[k] += 2.0 * a * w[ID4+ic].val[k];
     }
    }
  }

  inline void mult_gm5_dirac_vec(Vsimd_t *v, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = w[ID3+ic].val[k];
       v[ID2+ic].val[k] = w[ID4+ic].val[k];
       v[ID3+ic].val[k] = w[ID1+ic].val[k];
       v[ID4+ic].val[k] = w[ID2+ic].val[k];
     }
    }
  }

  inline void mult_gm5_chiral_vec(Vsimd_t *v, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] =   w[ID1+ic].val[k];
       v[ID2+ic].val[k] =   w[ID2+ic].val[k];
       v[ID3+ic].val[k] = - w[ID3+ic].val[k];
       v[ID4+ic].val[k] = - w[ID4+ic].val[k];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_dirac_vec(Vsimd_t *v, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = w[k + VLEN*(ID3 + ic)];
       v[ID2+ic].val[k] = w[k + VLEN*(ID4 + ic)];
       v[ID3+ic].val[k] = w[k + VLEN*(ID1 + ic)];
       v[ID4+ic].val[k] = w[k + VLEN*(ID2 + ic)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_chiral_vec(Vsimd_t *v, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] =   w[k + VLEN*(ID1 + ic)];
       v[ID2+ic].val[k] =   w[k + VLEN*(ID2 + ic)];
       v[ID3+ic].val[k] = - w[k + VLEN*(ID3 + ic)];
       v[ID4+ic].val[k] = - w[k + VLEN*(ID4 + ic)];
     }
    }
  }

  template <typename REALTYPE>
    inline void load_mult_gm5_dirac_vec(Vsimd_t *v,
                                        REALTYPE a, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] = a * w[k + VLEN*(ID3 + ic)];
       v[ID2+ic].val[k] = a * w[k + VLEN*(ID4 + ic)];
       v[ID3+ic].val[k] = a * w[k + VLEN*(ID1 + ic)];
       v[ID4+ic].val[k] = a * w[k + VLEN*(ID2 + ic)];
     }
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_chiral_vec(Vsimd_t *v,
                                       REALTYPE a, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
     for(int k = 0; k < VLEN; ++k){
       v[ID1+ic].val[k] =   a * w[k + VLEN*(ID1 + ic)];
       v[ID2+ic].val[k] =   a * w[k + VLEN*(ID2 + ic)];
       v[ID3+ic].val[k] = - a * w[k + VLEN*(ID3 + ic)];
       v[ID4+ic].val[k] = - a * w[k + VLEN*(ID4 + ic)];
     }
    }
  }

} // end of nameless namespace

#endif
//============================================================END=====
