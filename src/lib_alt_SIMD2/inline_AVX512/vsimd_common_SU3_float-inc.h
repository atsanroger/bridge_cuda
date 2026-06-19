/*!
        @file    vsimd_common_SU3_float-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_COMMON_SU3_FLOAT_INC_INCLUDED
#define VSIMD_COMMON_SU3_FLOAT_INC_INCLUDED

namespace {

  inline void mult_Vnn(__m512& wt, __m512* ut,  __m512* vt, int Nc)
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

  inline void mult_Vdn(__m512& wt, __m512 *ut, __m512 *vt, int Nc)
  {
    __m512 ur, ui, vtp, yt;

    ur  = _mm512_permute_ps(ut[0], 0b10100000);
    ui  = _mm512_permute_ps(ut[0], 0b11110101);
    vtp = _mm512_permute_ps(vt[0], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmsubadd_ps(ur, vt[0], yt);

    ur  = _mm512_permute_ps(ut[1], 0b10100000);
    ui  = _mm512_permute_ps(ut[1], 0b11110101);
    vtp = _mm512_permute_ps(vt[1], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt[1], yt);
    wt  = _mm512_add_ps(wt, yt);

    ur  = _mm512_permute_ps(ut[2], 0b10100000);
    ui  = _mm512_permute_ps(ut[2], 0b11110101);
    vtp = _mm512_permute_ps(vt[2], 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt[2], yt);
    wt  = _mm512_add_ps(wt, yt);
  }


} // end of nameless namespace

//============================================================END=====

#endif
