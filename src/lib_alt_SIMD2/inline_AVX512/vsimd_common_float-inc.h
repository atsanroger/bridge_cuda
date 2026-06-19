/*!
        @file    vsimd_common_float-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef VSIMD_COMMON_FLOAT_INC_INCLUDED
#define VSIMD_COMMON_FLOAT_INC_INCLUDED

// This file is Intel AVX-512 version.

namespace {

  template <typename REALTYPE>
  inline void clear_vec(REALTYPE *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      for(int k = 0; k < VLEN; ++k){
        vt[k + VLEN * in] = 0.0;
      }
    }
  }

  inline void clear_vec(__m512 *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      vt[in] = _mm512_setzero_ps();
    }
  }

  inline void copy_vec(__m512 *v, __m512 *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = w[in];
    }
  }

  template <typename REALTYPE>
  inline void add_vec(REALTYPE *v, __m512 *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      __m512 wt;
      wt = _mm512_load_ps(&v[VLEN*in]);
      wt = _mm512_add_ps(wt, vt[in]);
      _mm512_store_ps(&v[VLEN*in], wt);
    }
  }

  inline void add_vec(__m512 *v, __m512 *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = _mm512_add_ps(v[in], vt[in]);
    }
  }

  template <typename REALTYPE>
  inline void sub_vec(REALTYPE *v, __m512 *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      __m512 wt;
      wt = _mm512_load_ps(&v[VLEN*in]);
      wt = _mm512_sub_ps(wt, vt[in]);
      _mm512_store_ps(&v[VLEN*in], wt);
    }
  }

  inline void sub_vec(__m512 *v, __m512 *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = _mm512_sub_ps(v[in], vt[in]);
    }
  }

  inline void multadd_vec(__m512& wt, __m512& ut, __m512& vt)
  {
    __m512 yt;
    yt = _mm512_mul_ps(ut, vt);
    wt = _mm512_add_ps(wt, yt);
  }

  inline void mult_nn_vec(__m512& wt, __m512& ut,  __m512& vt)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmaddsub_ps(ur, vt, yt);
  }

  inline void mult_nd_vec(__m512& wt, __m512& vt,  __m512& ut)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmsubadd_ps(ur, vt, yt);
  }

  inline void mult_dn_vec(__m512& wt, __m512& ut,  __m512& vt)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    wt  = _mm512_fmsubadd_ps(ur, vt, yt);
  }

  inline void multadd_nn_vec(__m512& wt, __m512& ut,  __m512& vt)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmaddsub_ps(ur, vt, yt);
    wt  = _mm512_add_ps(wt, yt);
  }

  inline void multadd_nd_vec(__m512& wt, __m512& vt,  __m512& ut)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt, yt);
    wt  = _mm512_add_ps(wt, yt);
  }

  inline void multadd_dn_vec(__m512& wt, __m512& ut,  __m512& vt)
  {
    __m512 ur, ui, vtp, yt;
    ur  = _mm512_permute_ps(ut, 0b10100000);
    ui  = _mm512_permute_ps(ut, 0b11110101);
    vtp = _mm512_permute_ps(vt, 0b10110001);
    yt  = _mm512_mul_ps(ui, vtp);
    yt  = _mm512_fmsubadd_ps(ur, vt, yt);
    wt  = _mm512_add_ps(wt, yt);
  }

  inline void axpy_vec(__m512 *y, float a, __m512 *x, int Nin)
  {
    __m512 av = _mm512_set1_ps(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmadd_ps(av, x[in], y[in]);
    }
  }

  inline void axpy_vec(__m512 *y, float ar, float ai, __m512 *x, int Nin)
  {
    __m512 avr = _mm512_set1_ps(ar);
    __m512 avi = _mm512_set1_ps(ai);
    for(int in = 0; in < Nin; ++in){
      __m512 xp  = _mm512_permute_ps(x[in], 0b10110001);
      __m512 yt  = _mm512_mul_ps(avi, xp);
      yt = _mm512_fmaddsub_ps(avr, x[in], yt);
      y[in] = _mm512_add_ps(y[in], yt);
    }
  }

  inline void aypx_vec(float a, __m512 *y, __m512 *x, int Nin)
  {
    __m512 av = _mm512_set1_ps(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmadd_ps(av, y[in], x[in]);
    }
  }

  inline void aypx_vec(float ar, float ai, __m512 *y, __m512 *x, int Nin)
  {
    __m512 avr = _mm512_set1_ps(ar);
    __m512 avi = _mm512_set1_ps(ai);
    for(int in = 0; in < Nin; ++in){
      __m512 yp  = _mm512_permute_ps(y[in], 0b10110001);
      __m512 yt  = _mm512_mul_ps(avi, yp);
      yt = _mm512_fmaddsub_ps(avr, y[in], yt);
      y[in] = _mm512_add_ps(yt, x[in]);
    }
  }

  inline void scal_vec(__m512 *y, float a, int Nin)
  {
    __m512 av = _mm512_set1_ps(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_ps(av, y[in]);
    }
  }

  inline void scal_vec(__m512 *y, float ar, float ai, int Nin)
  {
    __m512 avr = _mm512_set1_ps(ar);
    __m512 avi = _mm512_set1_ps(ai);
    for(int in = 0; in < Nin; ++in){
      __m512 yp  = _mm512_permute_ps(y[in], 0b10110001);
      __m512 yt  = _mm512_mul_ps(avi, yp);
      y[in] = _mm512_fmaddsub_ps(avr, y[in], yt);
    }
  }

  inline void set_vec(__m512 *y, float a, __m512 *x, int Nin)
  {
    __m512 av = _mm512_set1_ps(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_ps(av, x[in]);
    }
  }

  inline void unit_vec(__m512 *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_set4_ps(0.0, 1.0, 0.0, 1.0);
    }
  }

  inline void conjg_vec(__m512 *y, int Nin)
  {
    __m512 av = _mm512_set1_ps(0.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmsubadd_ps(av, av, y[in]);
    }
  }

  inline void real_vec(__m512 *y, int Nin)
  {
    __m512 av = _mm512_set4_ps(0.0, 1.0, 0.0, 1.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_ps(av, y[in]);
    }
  }

  inline void imag_vec(__m512 *y, int Nin)
  {
    __m512 av = _mm512_set4_ps(1.0, 0.0, 1.0, 0.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_ps(av, y[in]);
    }
  }

  template<typename REALTYPE>
    inline void xI_vec(__m512 *y, int Nin)
    {
      __m512 av = _mm512_set4_ps(1.0, -1.0, 1.0, -1.0);
      for(int in = 0; in < Nin; ++in){
	__m512 yp = _mm512_permute_ps(y[in], 0b10110001);
	y[in] = _mm512_mul_ps(av, yp);
      }
    }

  inline void load_vec(__m512 *vt, float *v, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      vt[in] = _mm512_load_ps(&v[VLEN*in]);
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
    inline void save_vec(REALTYPE *v, __m512 *vt, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	_mm512_store_ps(&v[VLEN*in], vt[in]);
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
    inline void shift_vec0_bw(__m512 *v, REALTYPE *w, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in],
			       w[14 + VLEN*in], w[15 + VLEN*in],
			       0.0,             0.0            );
      }
    }

  template <typename REALTYPE>
    inline void shift_vec0_fw(__m512 *v, REALTYPE *w, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(0.0,           0.0,
				 w[0  + VLEN*in], w[1  + VLEN*in],
				 w[2  + VLEN*in], w[3  + VLEN*in],
				 w[4  + VLEN*in], w[5  + VLEN*in],
				 w[6  + VLEN*in], w[7  + VLEN*in],
				 w[8  + VLEN*in], w[9  + VLEN*in],
				 w[10 + VLEN*in], w[11 + VLEN*in],
				 w[12 + VLEN*in], w[13 + VLEN*in]);
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
    inline void shift_vec1_bw(__m512 *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in],
			       w[14 + VLEN*in], w[15 + VLEN*in],
			       y[2*in],         y[2*in+1]);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec1_fw(__m512 *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(y[2*in],         y[2*in+1],
			       w[0  + VLEN*in], w[1  + VLEN*in],
			       w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in]);
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
    inline void shift_vec2_bw(__m512 *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in],
			       w[14 + VLEN*in], w[15 + VLEN*in],
			       y[0  + VLEN*in], y[1  + VLEN*in]);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_bw(__m512 *v, REALTYPE *w, REALTYPE *y,
			      REALTYPE bc, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
        REALTYPE y0 = y[0  + VLEN*in] * bc;
        REALTYPE y1 = y[1  + VLEN*in] * bc;
	v[in] = _mm512_setr_ps(w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in],
			       w[14 + VLEN*in], w[15 + VLEN*in],
			       y0,              y1);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_fw(__m512 *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_ps(y[14 + VLEN*in], y[15 + VLEN*in],
			       w[0  + VLEN*in], w[1  + VLEN*in],
			       w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in]);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_fw(__m512 *v, REALTYPE *w, REALTYPE *y,
			      REALTYPE bc, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
        REALTYPE y0 = y[14 + VLEN*in] * bc;
        REALTYPE y1 = y[15 + VLEN*in] * bc;
	v[in] = _mm512_setr_ps(y0,              y1,
			       w[0  + VLEN*in], w[1  + VLEN*in],
			       w[2  + VLEN*in], w[3  + VLEN*in],
			       w[4  + VLEN*in], w[5  + VLEN*in],
			       w[6  + VLEN*in], w[7  + VLEN*in],
			       w[8  + VLEN*in], w[9  + VLEN*in],
			       w[10 + VLEN*in], w[11 + VLEN*in],
			       w[12 + VLEN*in], w[13 + VLEN*in]);
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


} // end of nameless namespace

//============================================================END=====

#endif
