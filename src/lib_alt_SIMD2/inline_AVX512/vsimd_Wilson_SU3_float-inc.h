/*!
        @file    vsimd_Wilson_SU3_float-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_WILSON_SU3_FLOAT_INC_INCLUDED
#define VSIMD_WILSON_SU3_FLOAT_INC_INCLUDED

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

  inline void mult_uv(__m512& wt, __m512* ut,  __m512* vt, int Nc)
  {
    __m512 ur, ui, vtp, yt;

    ur  = _mm512_permute_ps(ut[0], 0b10100000);
    ui  = _mm512_permute_ps(ut[0], 0b11110101);
    vtp = _mm512_permute_ps(vt[0], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmaddsub_ps(ur, vt[0], yt);

    ur  = _mm512_permute_ps(ut[1], 0b10100000);
    ui  = _mm512_permute_ps(ut[1], 0b11110101);
    vtp = _mm512_permute_ps(vt[1], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmaddsub_ps(ur, vt[1], yt);
    wt  = _mm512_add_ps(wt, yt);

    ur  = _mm512_permute_ps(ut[2], 0b10100000);
    ui  = _mm512_permute_ps(ut[2], 0b11110101);
    vtp = _mm512_permute_ps(vt[2], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmaddsub_ps(ur, vt[2], yt);
    wt  = _mm512_add_ps(wt, yt);
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

  inline void mult_udagv(__m512& wt, __m512 *ut, __m512 *vt, int Nc)
  {
    __m512 ur, ui, vtp, yt;

    ur  = _mm512_permute_ps(ut[0], 0b10100000);
    ui  = _mm512_permute_ps(ut[0], 0b11110101);
    vtp = _mm512_permute_ps(vt[0], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmsubadd_ps(ur, vt[0], yt);

    ur  = _mm512_permute_ps(ut[3], 0b10100000);
    ui  = _mm512_permute_ps(ut[3], 0b11110101);
    vtp = _mm512_permute_ps(vt[1], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt[1], yt);
    wt  = _mm512_add_ps(wt, yt);

    ur  = _mm512_permute_ps(ut[6], 0b10100000);
    ui  = _mm512_permute_ps(ut[6], 0b11110101);
    vtp = _mm512_permute_ps(vt[2], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt[2], yt);
    wt  = _mm512_add_ps(wt, yt);
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

  inline void load_set_sp2_xp(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC * id)]);
      }
      __m512 a1  = _mm512_set1_ps(1.0);
      __m512 vtp = _mm512_permute_ps(vt[3], 0b10110001);
      vt1[ic] = _mm512_fmaddsub_ps(a1, vt[0], vtp);
      vtp = _mm512_permute_ps(vt[2], 0b10110001);
      vt2[ic] = _mm512_fmaddsub_ps(a1, vt[1], vtp);
    }
  }

  inline void load_set_sp2_xm(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC * id)]);
      }
      __m512 a1  = _mm512_set1_ps(1.0);
      __m512 vtp = _mm512_permute_ps(vt[3], 0b10110001);
      vt1[ic] = _mm512_fmsubadd_ps(a1, vt[0], vtp);
      vtp = _mm512_permute_ps(vt[2], 0b10110001);
      vt2[ic] = _mm512_fmsubadd_ps(a1, vt[1], vtp);
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

  inline void load_set_sp2_yp(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC * id)]);
      }
      vt1[ic] = _mm512_add_ps(vt[0], vt[3]);
      vt2[ic] = _mm512_sub_ps(vt[1], vt[2]);
    }
  }

  inline void load_set_sp2_ym(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC * id)]);
      }
      vt1[ic] = _mm512_sub_ps(vt[0], vt[3]);
      vt2[ic] = _mm512_add_ps(vt[1], vt[2]);
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

  inline void load_set_sp2_zp(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC*id)]);
      }
      __m512 a1  = _mm512_set1_ps(1.0);
      __m512 vtp = _mm512_permute_ps(vt[2], 0b10110001);
      vt1[ic] = _mm512_fmaddsub_ps(a1, vt[0], vtp);
      vtp = _mm512_permute_ps(vt[3], 0b10110001);
      vt2[ic] = _mm512_fmsubadd_ps(a1, vt[1], vtp);
    }
  }

  inline void load_set_sp2_zm(__m512* vt1, __m512* vt2, float* v, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND];
      for(int id = 0; id < ND; ++id){
	vt[id] = _mm512_load_ps(&v[VLEN*(ic + NC*id)]);
      }
      __m512 a1  = _mm512_set1_ps(1.0);
      __m512 vtp = _mm512_permute_ps(vt[2], 0b10110001);
      vt1[ic] = _mm512_fmsubadd_ps(a1, vt[0], vtp);
      vtp = _mm512_permute_ps(vt[3], 0b10110001);
      vt2[ic] = _mm512_fmaddsub_ps(a1, vt[1], vtp);
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

  inline void load_set_sp2_tp_dirac(__m512* vt1, __m512* vt2,
                                    float* v, int Nc)
  {
    __m512 a2 = _mm512_set1_ps(2.0);
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND2];
      vt[0] = _mm512_load_ps(&v[VLEN*(ic + ID3)]);
      vt[1] = _mm512_load_ps(&v[VLEN*(ic + ID4)]);
      vt1[ic] = _mm512_mul_ps(a2, vt[0]);
      vt2[ic] = _mm512_mul_ps(a2, vt[1]);
    }
  }

  inline void load_set_sp2_tm_dirac(__m512* vt1, __m512* vt2,
                                    float* v, int Nc)
  {
    __m512 a2 = _mm512_set1_ps(2.0);
    for(int ic = 0; ic < NC; ++ic){
      __m512 vt[ND2];
      vt[0] = _mm512_load_ps(&v[VLEN*(ic + ID1)]);
      vt[1] = _mm512_load_ps(&v[VLEN*(ic + ID2)]);
      vt1[ic] = _mm512_mul_ps(a2, vt[0]);
      vt2[ic] = _mm512_mul_ps(a2, vt[1]);
    }
  }

  inline void set_sp4_xp(__m512* v2, __m512& wt1, __m512& wt2, int Nc)
  {
    __m512 a1 = _mm512_set1_ps(1.0);
    __m512 wtp;

    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    wtp = _mm512_permute_ps(wt2, 0b10110001);
    v2[ID3] = _mm512_fmsubadd_ps(a1, v2[ID3], wtp);
    wtp = _mm512_permute_ps(wt1, 0b10110001);
    v2[ID4] = _mm512_fmsubadd_ps(a1, v2[ID4], wtp);
  }

  inline void set_sp4_xm(__m512* v2, __m512& wt1, __m512& wt2, int Nc)
  {
    __m512 a1 = _mm512_set1_ps(1.0);
    __m512 wtp;
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    wtp = _mm512_permute_ps(wt2, 0b10110001);
    v2[ID3] = _mm512_fmaddsub_ps(a1, v2[ID3], wtp);
    wtp = _mm512_permute_ps(wt1, 0b10110001);
    v2[ID4] = _mm512_fmaddsub_ps(a1, v2[ID4], wtp);
  }

  inline void set_sp4_yp(__m512* v2, __m512& wt1, __m512& wt2, int Nc)
  {
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    v2[ID3] = _mm512_sub_ps(v2[ID3], wt2);
    v2[ID4] = _mm512_add_ps(v2[ID4], wt1);
  }

  inline void set_sp4_ym(__m512* v2,  __m512& wt1, __m512& wt2, int Nc)
  {
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    v2[ID3] = _mm512_add_ps(v2[ID3], wt2);
    v2[ID4] = _mm512_sub_ps(v2[ID4], wt1);
  }

  inline void set_sp4_zp(__m512* v2, __m512& wt1, __m512& wt2, int Nc)
  {
    __m512 a1 = _mm512_set1_ps(1.0);
    __m512 wtp;
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    wtp = _mm512_permute_ps(wt1, 0b10110001);
    v2[ID3] = _mm512_fmsubadd_ps(a1, v2[ID3], wtp);
    wtp = _mm512_permute_ps(wt2, 0b10110001);
    v2[ID4] = _mm512_fmaddsub_ps(a1, v2[ID4], wtp);
  }

  inline void set_sp4_zm(__m512* v2,  __m512& wt1, __m512& wt2, int Nc)
  {
    __m512 a1 = _mm512_set1_ps(1.0);
    __m512 wtp;
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
    wtp = _mm512_permute_ps(wt1, 0b10110001);
    v2[ID3] = _mm512_fmaddsub_ps(a1, v2[ID3], wtp);
    wtp = _mm512_permute_ps(wt2, 0b10110001);
    v2[ID4] = _mm512_fmsubadd_ps(a1, v2[ID4], wtp);
  }

  inline void set_sp4_tp_dirac(__m512* v2, __m512& wt1, __m512& wt2, int Nc)
  {
    v2[ID3] = _mm512_add_ps(v2[ID3], wt1);
    v2[ID4] = _mm512_add_ps(v2[ID4], wt2);
  }

  inline void set_sp4_tm_dirac(__m512* v2,  __m512& wt1, __m512& wt2, int Nc)
  {
    v2[ID1] = _mm512_add_ps(v2[ID1], wt1);
    v2[ID2] = _mm512_add_ps(v2[ID2], wt2);
  }

  template<typename REALTYPE>
  inline void set_aPp5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av = _mm512_set1_ps(a);
    __m512 vt;
    for(int ic = 0; ic < NC; ++ic){
      vt = _mm512_add_ps(w[ID1+ic], w[ID3+ic]);
      v[ID1+ic] = _mm512_mul_ps(av, vt);
      v[ID3+ic] = _mm512_mul_ps(av, vt);
      vt = _mm512_add_ps(w[ID2+ic], w[ID4+ic]);
      v[ID2+ic] = _mm512_mul_ps(av, vt);
      v[ID4+ic] = _mm512_mul_ps(av, vt);
    }
  }

  template<typename REALTYPE>
  inline void set_aPm5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av  = _mm512_set1_ps(a);
    __m512 avm = _mm512_set1_ps(-a);
    __m512 vt;
    for(int ic = 0; ic < NC; ++ic){
      vt = _mm512_sub_ps(w[ID1+ic], w[ID3+ic]);
      v[ID1+ic] = _mm512_mul_ps(av,  vt);
      v[ID3+ic] = _mm512_mul_ps(avm, vt);
      vt = _mm512_sub_ps(w[ID2+ic], w[ID4+ic]);
      v[ID2+ic] = _mm512_mul_ps(av,  vt);
      v[ID4+ic] = _mm512_mul_ps(avm, vt);
    }
  }

  template<typename REALTYPE>
  inline void add_aPp5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av = _mm512_set1_ps(a);
    __m512 vt;
    for(int ic = 0; ic < NC; ++ic){
      vt = _mm512_add_ps(w[ID1+ic], w[ID3+ic]);
      v[ID1+ic] = _mm512_fmadd_ps(av, vt, v[ID1+ic]);
      v[ID3+ic] = _mm512_fmadd_ps(av, vt, v[ID3+ic]);
      vt = _mm512_add_ps(w[ID2+ic], w[ID4+ic]);
      v[ID2+ic] = _mm512_fmadd_ps(av, vt, v[ID2+ic]);
      v[ID4+ic] = _mm512_fmadd_ps(av, vt, v[ID4+ic]);
    }
  }

  template<typename REALTYPE>
  inline void add_aPm5_dirac_vec(Vsimd_t *v,
                                 REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av  = _mm512_set1_ps(a);
    __m512 avm = _mm512_set1_ps(-a);
    __m512 vt;
    for(int ic = 0; ic < NC; ++ic){
      vt = _mm512_sub_ps(w[ID1+ic], w[ID3+ic]);
      v[ID1+ic] = _mm512_fmadd_ps(av,  vt, v[ID1+ic]);
      v[ID3+ic] = _mm512_fmadd_ps(avm, vt, v[ID3+ic]);
      vt = _mm512_sub_ps(w[ID2+ic], w[ID4+ic]);
      v[ID2+ic] = _mm512_fmadd_ps(av,  vt, v[ID2+ic]);
      v[ID4+ic] = _mm512_fmadd_ps(avm, vt, v[ID4+ic]);
    }
  }

  template<typename REALTYPE>
  inline void set_aPp5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av2  = _mm512_set1_ps(2.0*a);
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = _mm512_mul_ps(av2, w[ID1+ic]);
      v[ID2+ic] = _mm512_mul_ps(av2, w[ID2+ic]);
      v[ID3+ic] = _mm512_setzero_ps();
      v[ID4+ic] = _mm512_setzero_ps();
    }
  }

  template<typename REALTYPE>
  inline void set_aPm5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av2  = _mm512_set1_ps(2.0*a);
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = _mm512_setzero_ps();
      v[ID2+ic] = _mm512_setzero_ps();
      v[ID3+ic] = _mm512_mul_ps(av2, w[ID1+ic]);
      v[ID4+ic] = _mm512_mul_ps(av2, w[ID2+ic]);
    }
  }

  template<typename REALTYPE>
  inline void add_aPp5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av2  = _mm512_set1_ps(2.0*a);
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = _mm512_fmadd_ps(av2, w[ID1+ic], v[ID1+ic]);
      v[ID2+ic] = _mm512_fmadd_ps(av2, w[ID2+ic], v[ID2+ic]);
    }
  }

  template<typename REALTYPE>
  inline void add_aPm5_chiral_vec(Vsimd_t *v,
                                  REALTYPE a, Vsimd_t *w, int Nc)
  {
    __m512 av2  = _mm512_set1_ps(2.0*a);
    for(int ic = 0; ic < NC; ++ic){
      v[ID3+ic] = _mm512_fmadd_ps(av2, w[ID3+ic], v[ID3+ic]);
      v[ID4+ic] = _mm512_fmadd_ps(av2, w[ID4+ic], v[ID4+ic]);
    }
  }

  inline void mult_gm5_dirac_vec(Vsimd_t *v, Vsimd_t *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = w[ID3+ic];
      v[ID2+ic] = w[ID4+ic];
      v[ID3+ic] = w[ID1+ic];
      v[ID4+ic] = w[ID2+ic];
    }
  }

  inline void mult_gm5_chiral_vec(Vsimd_t *v, Vsimd_t *w, int Nc)
  {
    __m512 am  = _mm512_set1_ps(-1.0);
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = w[ID1+ic];
      v[ID2+ic] = w[ID2+ic];
      v[ID3+ic] = _mm512_mul_ps(am, w[ID3+ic]);
      v[ID4+ic] = _mm512_mul_ps(am, w[ID4+ic]);
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_dirac_vec(Vsimd_t *v, REALTYPE *w, int Nc)
  {
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = _mm512_load_ps(&w[VLEN*(ID3+ic)]);
      v[ID2+ic] = _mm512_load_ps(&w[VLEN*(ID4+ic)]);
      v[ID3+ic] = _mm512_load_ps(&w[VLEN*(ID1+ic)]);
      v[ID4+ic] = _mm512_load_ps(&w[VLEN*(ID2+ic)]);
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_chiral_vec(Vsimd_t *v, REALTYPE *w, int Nc)
  {
    __m512 am  = _mm512_set1_ps(-1.0);
    __m512 wt;
    for(int ic = 0; ic < NC; ++ic){
      v[ID1+ic] = _mm512_load_ps(&w[VLEN*(ID1+ic)]);
      v[ID2+ic] = _mm512_load_ps(&w[VLEN*(ID2+ic)]);
      wt = _mm512_load_ps(&w[VLEN*(ID3+ic)]);
      v[ID3+ic] = _mm512_mul_ps(am, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID4+ic)]);
      v[ID4+ic] = _mm512_mul_ps(am, wt);
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_dirac_vec(Vsimd_t *v,
                                      REALTYPE a, REALTYPE *w, int Nc)
  {
    __m512 av = _mm512_set1_ps(a);
    __m512 wt;
    for(int ic = 0; ic < NC; ++ic){
      wt = _mm512_load_ps(&w[VLEN*(ID3+ic)]);
      v[ID1+ic] = _mm512_mul_ps(av, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID4+ic)]);
      v[ID2+ic] = _mm512_mul_ps(av, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID1+ic)]);
      v[ID3+ic] = _mm512_mul_ps(av, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID2+ic)]);
      v[ID4+ic] = _mm512_mul_ps(av, wt);
    }
  }

  template <typename REALTYPE>
  inline void load_mult_gm5_chiral_vec(Vsimd_t *v,
                                       REALTYPE a, REALTYPE *w, int Nc)
  {
    __m512 av = _mm512_set1_ps(a);
    __m512 am = _mm512_set1_ps(-a);
    __m512 wt;
    for(int ic = 0; ic < NC; ++ic){
      wt = _mm512_load_ps(&w[VLEN*(ID1+ic)]);
      v[ID1+ic] = _mm512_mul_ps(av, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID2+ic)]);
      v[ID2+ic] = _mm512_mul_ps(av, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID3+ic)]);
      v[ID3+ic] = _mm512_mul_ps(am, wt);
      wt = _mm512_load_ps(&w[VLEN*(ID4+ic)]);
      v[ID4+ic] = _mm512_mul_ps(am, wt);
    }
  }

} // end of nameless namespace

//============================================================END=====

#endif
