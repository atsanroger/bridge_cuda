/*!
        @file    vsimd_common_SU3_double-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_COMMON_SU3_DOUBLE_INC_INCLUDED
#define VSIMD_COMMON_SU3_DOUBLE_INC_INCLUDED

namespace {

  inline void mult_Vnn(__m512d& wt, __m512d* ut,  __m512d* vt, int Nc)
  {
    __m512d ur, ui, vtp, yt;

    ur  = _mm512_permute_pd(ut[0], 0b00000000);
    ui  = _mm512_permute_pd(ut[0], 0b11111111);
    vtp = _mm512_permute_pd(vt[0], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    wt  = _mm512_fmaddsub_pd(ur, vt[0], yt);

    ur  = _mm512_permute_pd(ut[1], 0b00000000);
    ui  = _mm512_permute_pd(ut[1], 0b11111111);
    vtp = _mm512_permute_pd(vt[1], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmaddsub_pd(ur, vt[1], yt);
    wt  = _mm512_add_pd(wt, yt);

    ur  = _mm512_permute_pd(ut[2], 0b00000000);
    ui  = _mm512_permute_pd(ut[2], 0b11111111);
    vtp = _mm512_permute_pd(vt[2], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmaddsub_pd(ur, vt[2], yt);
    wt  = _mm512_add_pd(wt, yt);
  }

  inline void mult_Vdn(__m512d& wt, __m512d* ut,  __m512d* vt, int Nc)
  {
    __m512d ur, ui, vtp, yt;

    ur  = _mm512_permute_pd(ut[0], 0b00000000);
    ui  = _mm512_permute_pd(ut[0], 0b11111111);
    vtp = _mm512_permute_pd(vt[0], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    wt  = _mm512_fmsubadd_pd(ur, vt[0], yt);

    ur  = _mm512_permute_pd(ut[1], 0b00000000);
    ui  = _mm512_permute_pd(ut[1], 0b11111111);
    vtp = _mm512_permute_pd(vt[1], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmsubadd_pd(ur, vt[1], yt);
    wt  = _mm512_add_pd(wt, yt);

    ur  = _mm512_permute_pd(ut[2], 0b00000000);
    ui  = _mm512_permute_pd(ut[2], 0b11111111);
    vtp = _mm512_permute_pd(vt[2], 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmsubadd_pd(ur, vt[2], yt);
    wt  = _mm512_add_pd(wt, yt);
  }

} // end of nameless namespace

//============================================================END=====
#endif
