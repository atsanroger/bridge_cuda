/*!
        @file    vsimd_common_double-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_COMMON_DOUBLE_INC_INCLUDED
#define VSIMD_COMMON_DOUBLE_INC_INCLUDED

namespace {

  inline void clear_vec(Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        vt[in].val[k] = 0.0;
      }
    }
  }

  template <typename REALTYPE>
  inline void clear_vec(REALTYPE *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        vt[k + VLEN * in] = 0.0;
      }
    }
  }

  inline void copy_vec(Vsimd_t *v, Vsimd_t *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[in].val[k] = w[in].val[k];
      }
    }
  }

  inline void add_vec(Vsimd_t *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[in].val[k] += vt[in].val[k];
      }
    }
  }

  template <typename REALTYPE>
  inline void add_vec(REALTYPE *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[k + VLEN * in] += vt[in].val[k];
      }
    }
  }

  inline void sub_vec(Vsimd_t *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[in].val[k] -= vt[in].val[k];
      }
    }
  }

  template <typename REALTYPE>
  inline void sub_vec(REALTYPE *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[k + VLEN * in] -= vt[in].val[k];
      }
    }
  }

  inline void multadd_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN; ++k){
      v.val[k] += w1.val[k] * w2.val[k];
    }
  }

  inline void mult_nn_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] = w1.val[kr] * w2.val[kr] - w1.val[ki] * w2.val[ki];
      v.val[ki] = w1.val[kr] * w2.val[ki] + w1.val[ki] * w2.val[kr];
    }
  }

  inline void mult_nd_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] =   w1.val[kr] * w2.val[kr] + w1.val[ki] * w2.val[ki];
      v.val[ki] = - w1.val[kr] * w2.val[ki] + w1.val[ki] * w2.val[kr];
    }
  }

  inline void mult_dn_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] = w1.val[kr] * w2.val[kr] + w1.val[ki] * w2.val[ki];
      v.val[ki] = w1.val[kr] * w2.val[ki] - w1.val[ki] * w2.val[kr];
    }
  }

  inline void multadd_nn_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] += w1.val[kr] * w2.val[kr] - w1.val[ki] * w2.val[ki];
      v.val[ki] += w1.val[kr] * w2.val[ki] + w1.val[ki] * w2.val[kr];
    }
  }

  inline void multadd_nd_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] +=   w1.val[kr] * w2.val[kr] + w1.val[ki] * w2.val[ki];
      v.val[ki] += - w1.val[kr] * w2.val[ki] + w1.val[ki] * w2.val[kr];
    }
  }

  inline void multadd_dn_vec(Vsimd_t& v, Vsimd_t& w1, Vsimd_t& w2)
  {
    for(int k = 0; k < VLEN2; ++k){
      int kr = 2*k;
      int ki = 2*k+1;
      v.val[kr] += w1.val[kr] * w2.val[kr] + w1.val[ki] * w2.val[ki];
      v.val[ki] += w1.val[kr] * w2.val[ki] - w1.val[ki] * w2.val[kr];
    }
  }

  template <typename REALTYPE>
  inline void axpy_vec(Vsimd_t *y, REALTYPE a, Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].val[k] += a * x[in].val[k];
      }
    }
  }

  template <typename REALTYPE>
  inline void axpy_vec(Vsimd_t *y, REALTYPE ar, REALTYPE ai,
                                           Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        int kr = 2*k;
        int ki = 2*k + 1;
        y[in].val[kr] += ar * x[in].val[kr] - ai * x[in].val[ki];
        y[in].val[ki] += ar * x[in].val[ki] + ai * x[in].val[kr];
      }
    }
  }

  template <typename REALTYPE>
  inline void aypx_vec(REALTYPE a, Vsimd_t *y, Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].val[k] = a * y[in].val[k] + x[in].val[k];
      }
    }
  }

  template <typename REALTYPE>
  inline void aypx_vec(REALTYPE ar, REALTYPE ai, Vsimd_t *y,
                                           Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        int kr = 2 * k;
        int ki = 2 * k + 1;
        REALTYPE yr = y[in].val[kr];
        REALTYPE yi = y[in].val[ki];
        y[in].val[kr] = ar * yr - ai * yi + x[in].val[kr];
        y[in].val[ki] = ar * yi + ai * yr + x[in].val[ki];
      }
    }
  }

  template <typename REALTYPE>
  inline void scal_vec(Vsimd_t *y, REALTYPE a, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].val[k] = a * y[in].val[k];
      }
    }
  }

  template <typename REALTYPE>
    inline void scal_vec(Vsimd_t *y, REALTYPE ar, REALTYPE ai, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        int kr = 2 * k;
        int ki = 2 * k + 1;
        REALTYPE yr = y[in].val[kr];
        REALTYPE yi = y[in].val[ki];
        y[in].val[kr] = ar * yr - ai * yi;
        y[in].val[ki] = ar * yi + ai * yr;
      }
    }
  }

  template <typename REALTYPE>
  inline void set_vec(Vsimd_t *y, REALTYPE a, Vsimd_t *x, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        y[in].val[k] = a * x[in].val[k];
      }
    }
  }

  inline void unit_vec(Vsimd_t *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        y[in].val[2*k]   = 1.0;
        y[in].val[2*k+1] = 0.0;
      }
    }
  }

  inline void conjg_vec(Vsimd_t *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        y[in].val[2*k+1] = - y[in].val[2*k+1];
      }
    }
  }

  inline void real_vec(Vsimd_t *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        y[in].val[2*k+1] = 0.0;
      }
    }
  }

  inline void imag_vec(Vsimd_t *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        y[in].val[2*k] = 0.0;
      }
    }
  }

  template <typename REALTYPE>
  inline void xI_vec(Vsimd_t *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN2; ++k){
        REALTYPE tmp = y[in].val[2*k];
        y[in].val[2*k]   = - y[in].val[2*k+1];
        y[in].val[2*k+1] = tmp;
      }
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
  inline void load_vec(Vsimd_t *vt, REALTYPE *v, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        vt[in].val[k] = v[k + VLEN * in];
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
  inline void save_vec(REALTYPE *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        v[k + VLEN * in] = vt[in].val[k];
      }
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
  inline void shift_vec0_bw(Vsimd_t *v, REALTYPE *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[in].val[k] = w[k+2 + VLEN*in];
      }
      v[in].val[VLEN-2] = 0.0;
      v[in].val[VLEN-1] = 0.0;
    }
  }

  template <typename REALTYPE>
  inline void shift_vec0_fw(Vsimd_t *v, REALTYPE *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[in].val[k] = w[k-2 + VLEN*in];
      }
      v[in].val[0] = 0.0;
      v[in].val[1] = 0.0;
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
  inline void shift_vec1_bw(Vsimd_t *v, REALTYPE *w, REALTYPE *buf, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[in].val[k] = w[k+2 + VLEN*in];
      }
      v[in].val[VLEN-2] = buf[2*in];
      v[in].val[VLEN-1] = buf[2*in+1];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec1_fw(Vsimd_t *v, REALTYPE *w, REALTYPE *buf, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[in].val[k] = w[k-2 + VLEN*in];
      }
      v[in].val[0] = buf[2*in];
      v[in].val[1] = buf[2*in+1];
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
  inline void shift_vec2_bw(Vsimd_t *v, REALTYPE *w, REALTYPE *y,
                            REALTYPE bc, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[in].val[k] = w[k+2 + VLEN*in];
      }
      v[in].val[VLEN-2] = bc * y[0 + VLEN*in];
      v[in].val[VLEN-1] = bc * y[1 + VLEN*in];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec2_bw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN-2; ++k){
        v[in].val[k] = w[k+2 + VLEN*in];
      }
      v[in].val[VLEN-2] = y[0 + VLEN*in];
      v[in].val[VLEN-1] = y[1 + VLEN*in];
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
  inline void shift_vec2_fw(Vsimd_t *v, REALTYPE *w, REALTYPE *y,
                            REALTYPE bc, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[in].val[k] = w[k-2 + VLEN*in];
      }
      v[in].val[0] = bc * y[VLEN-2 + VLEN*in];
      v[in].val[1] = bc * y[VLEN-1 + VLEN*in];
    }
  }

  template <typename REALTYPE>
  inline void shift_vec2_fw(Vsimd_t *v, REALTYPE *w, REALTYPE *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 2; k < VLEN; ++k){
        v[in].val[k] = w[k-2 + VLEN*in];
      }
      v[in].val[0] = y[VLEN-2 + VLEN*in];
      v[in].val[1] = y[VLEN-1 + VLEN*in];
    }
  }



} // end of nameless namespace

#endif
//============================================================END=====
