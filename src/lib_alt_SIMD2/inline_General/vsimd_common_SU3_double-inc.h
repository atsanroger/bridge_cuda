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

  inline void mult_Vnn(Vsimd_t& v, Vsimd_t *g, Vsimd_t *w, int Nc)
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

  inline void mult_Vdn(Vsimd_t& v, Vsimd_t *g, Vsimd_t *w, int Nc)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] =
          g[0].val[kr] * w[0].val[kr] + g[0].val[ki] * w[0].val[ki]
        + g[1].val[kr] * w[1].val[kr] + g[1].val[ki] * w[1].val[ki]
        + g[2].val[kr] * w[2].val[kr] + g[2].val[ki] * w[2].val[ki];
      v.val[ki] =
          g[0].val[kr] * w[0].val[ki] - g[0].val[ki] * w[0].val[kr]
        + g[1].val[kr] * w[1].val[ki] - g[1].val[ki] * w[1].val[kr]
        + g[2].val[kr] * w[2].val[ki] - g[2].val[ki] * w[2].val[kr];
    }
  }

} // end of nameless namespace

#endif
//============================================================END=====
