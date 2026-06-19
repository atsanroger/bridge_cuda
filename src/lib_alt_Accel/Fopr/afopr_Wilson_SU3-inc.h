/*!
        @file    afopr_Wilson_SU3-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $
        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$
        @version $LastChangedRevision: 1928 $
*/

#ifndef ACCEL_AFOPR_WILSON_SU3_INC_INCLUDED
#define ACCEL_AFOPR_WILSON_SU3_INC_INCLUDED

#define  NC      3
#define  NVC     6
#define  NDF    18
#define  NDF2    9
#define  ND      4
#define  ND2     2
#define  NCD    12
#define  NVCD   24

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

  template <typename REALTYPE>
  inline void load_vec(REALTYPE *vt, REALTYPE *v, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        vt[k + VLEN * in] = v[k + VLEN * in];
      }
    }
  }

  template <typename REALTYPE>
  inline void load_vec1(REALTYPE *vt, REALTYPE *v, int k, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      vt[2*in  ] = v[2*k   + VLEN * in];
      vt[2*in+1] = v[2*k+1 + VLEN * in];
    }
  }

  template <typename REALTYPE>
  inline void save_vec(REALTYPE *v, REALTYPE *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[k + VLEN * in] = vt[k + VLEN * in];
      }
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
  inline void shift_vec0_bw(REALTYPE *v, REALTYPE *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[k + VLEN*in] = w[k+2 + VLEN*in];
      }
      v[VLEN-2 + VLEN*in] = 0.0;
      v[VLEN-1 + VLEN*in] = 0.0;
    }
  }

  template <typename REALTYPE>
  inline void shift_vec0_fw(REALTYPE *v, REALTYPE *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[k + VLEN*in] = w[k-2 + VLEN*in];
      }
      v[0 + VLEN*in] = 0.0;
      v[1 + VLEN*in] = 0.0;
    }
  }

  template <typename REALTYPE>
  inline void shift_vec1_bw(REALTYPE *v, REALTYPE *w, REALTYPE *buf, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[k + VLEN*in] = w[k+2 + VLEN*in];
      }
      v[VLEN-2 + VLEN*in] = buf[2*in];
      v[VLEN-1 + VLEN*in] = buf[2*in+1];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec1_fw(REALTYPE *v, REALTYPE *w, REALTYPE *buf, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[k + VLEN*in] = w[k-2 + VLEN*in];
      }
      v[0 + VLEN*in] = buf[2*in];
      v[1 + VLEN*in] = buf[2*in+1];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec2_bw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[k + VLEN*in] = w[k+2 + VLEN*in];
      }
      v[VLEN-2 + VLEN*in] = y[0 + VLEN*in];
      v[VLEN-1 + VLEN*in] = y[1 + VLEN*in];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec2_fw(REALTYPE *v, REALTYPE *w, REALTYPE *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[k + VLEN*in] = w[k-2 + VLEN*in];
      }
      v[0 + VLEN*in] = y[VLEN-2 + VLEN*in];
      v[1 + VLEN*in] = y[VLEN-1 + VLEN*in];
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
  inline void set_sp4_xp(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] +=  wt2[ki];
      v[ki + VLEN * ID3] += -wt2[kr];
      v[kr + VLEN * ID4] +=  wt1[ki];
      v[ki + VLEN * ID4] += -wt1[kr];
   }
 }

  template <typename REALTYPE>
  inline void set_sp4_xm(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] += -wt2[ki];
      v[ki + VLEN * ID3] +=  wt2[kr];
      v[kr + VLEN * ID4] += -wt1[ki];
      v[ki + VLEN * ID4] +=  wt1[kr];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_xm1(REALTYPE *v, int k,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    int kr = 2*k;
    int ki = 2*k+1;
    v[kr +VLEN * ID1] +=  wt1[0];
    v[ki +VLEN * ID1] +=  wt1[1];
    v[kr +VLEN * ID2] +=  wt2[0];
    v[ki +VLEN * ID2] +=  wt2[1];
    v[kr +VLEN * ID3] += -wt2[1];
    v[ki +VLEN * ID3] +=  wt2[0];
    v[kr +VLEN * ID4] += -wt1[1];
    v[ki +VLEN * ID4] +=  wt1[0];
  }

  template <typename REALTYPE>
  inline void set_sp4_yp(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] += -wt2[kr];
      v[ki + VLEN * ID3] += -wt2[ki];
      v[kr + VLEN * ID4] +=  wt1[kr];
      v[ki + VLEN * ID4] +=  wt1[ki];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_ym(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] +=  wt2[kr];
      v[ki + VLEN * ID3] +=  wt2[ki];
      v[kr + VLEN * ID4] += -wt1[kr];
      v[ki + VLEN * ID4] += -wt1[ki];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_zp(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] +=  wt1[ki];
      v[ki + VLEN * ID3] += -wt1[kr];
      v[kr + VLEN * ID4] += -wt2[ki];
      v[ki + VLEN * ID4] +=  wt2[kr];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_zm(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] +=  wt1[kr];
      v[ki + VLEN * ID1] +=  wt1[ki];
      v[kr + VLEN * ID2] +=  wt2[kr];
      v[ki + VLEN * ID2] +=  wt2[ki];
      v[kr + VLEN * ID3] += -wt1[ki];
      v[ki + VLEN * ID3] +=  wt1[kr];
      v[kr + VLEN * ID4] +=  wt2[ki];
      v[ki + VLEN * ID4] += -wt2[kr];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_tp_dirac(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID3] += wt1[kr];
      v[ki + VLEN * ID3] += wt1[ki];
      v[kr + VLEN * ID4] += wt2[kr];
      v[ki + VLEN * ID4] += wt2[ki];
    }
  }

  template <typename REALTYPE>
  inline void set_sp4_tm_dirac(REALTYPE *v,
                         REALTYPE *wt1, REALTYPE *wt2, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v[kr + VLEN * ID1] += wt1[kr];
      v[ki + VLEN * ID1] += wt1[ki];
      v[kr + VLEN * ID2] += wt2[kr];
      v[ki + VLEN * ID2] += wt2[ki];
    }
  }

} // end of nameless namespace

#endif
//============================================================END=====
