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

  inline void clear_vec(__m512d *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      vt[in] = _mm512_setzero_pd();
    }
  }

  inline void copy_vec(Vsimd_t *v, Vsimd_t *w, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = w[in];
    }
  }

  template <typename REALTYPE>
    inline void add_vec(REALTYPE *v, __m512d *vt, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	__m512d wt;
	wt = _mm512_load_pd(&v[VLEN*in]);
	wt = _mm512_add_pd(wt, vt[in]);
	_mm512_store_pd(&v[VLEN*in], wt);
      }
    }

  template <typename REALTYPE>
    inline void sub_vec(REALTYPE *v, __m512d *vt, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	__m512d wt;
	wt = _mm512_load_pd(&v[VLEN*in]);
	wt = _mm512_sub_pd(wt, vt[in]);
	_mm512_store_pd(&v[VLEN*in], wt);
      }
    }

  inline void add_vec(Vsimd_t *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = _mm512_add_pd(v[in], vt[in]);
    }
  }

  inline void sub_vec(Vsimd_t *v, Vsimd_t *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      v[in] = _mm512_sub_pd(v[in], vt[in]);
    }
  }

  inline void multadd_vec(__m512d& wt, __m512d& ut,  __m512d& vt)
  {
    __m512d yt;
    yt  = _mm512_mul_pd(ut, vt);
    wt  = _mm512_add_pd(wt, yt);
  }

  inline void mult_nn_vec(__m512d& wt, __m512d& ut,  __m512d& vt)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    wt  = _mm512_fmaddsub_pd(ur, vt, yt);
  }

  inline void mult_nd_vec(__m512d& wt, __m512d& vt,  __m512d& ut)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    wt  = _mm512_fmsubadd_pd(ur, vt, yt);
  }

  inline void mult_dn_vec(__m512d& wt, __m512d& ut,  __m512d& vt)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    wt  = _mm512_fmsubadd_pd(ur, vt, yt);
  }

  inline void multadd_nn_vec(__m512d& wt, __m512d& ut,  __m512d& vt)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmaddsub_pd(ur, vt, yt);
    wt  = _mm512_add_pd(wt, yt);
  }


  inline void multadd_nd_vec(__m512d& wt, __m512d& vt,  __m512d& ut)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmsubadd_pd(ur, vt, yt);
    wt  = _mm512_add_pd(wt, yt);
  }

  inline void multadd_dn_vec(__m512d& wt, __m512d& ut,  __m512d& vt)
  {
    __m512d ur, ui, vtp, yt;
    ur  = _mm512_permute_pd(ut, 0b00000000);
    ui  = _mm512_permute_pd(ut, 0b11111111);
    vtp = _mm512_permute_pd(vt, 0b01010101);
    yt  = _mm512_mul_pd(ui, vtp);
    yt  = _mm512_fmsubadd_pd(ur, vt, yt);
    wt  = _mm512_add_pd(wt, yt);
  }

  inline void axpy_vec(__m512d *y, double a, __m512d *x, int Nin)
  {
    __m512d av = _mm512_set1_pd(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmadd_pd(av, x[in], y[in]);
    }
  }

  inline void axpy_vec(__m512d *y, double ar, double ai, __m512d *x, int Nin)
  {
    __m512d avr = _mm512_set1_pd(ar);
    __m512d avi = _mm512_set1_pd(ai);
    for(int in = 0; in < Nin; ++in){
      __m512d xp  = _mm512_permute_pd(x[in], 0b01010101);
      __m512d yt  = _mm512_mul_pd(avi, xp);
      yt = _mm512_fmaddsub_pd(avr, x[in], yt);
      y[in] = _mm512_add_pd(y[in], yt);
    }
  }

  /*
  inline void aypx_vec(__m512d *y, double a, __m512d *x, int Nin)
  {
    __m512d av = _mm512_set1_pd(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmadd_pd(av, y[in], x[in]);
    }
  }
  */

  inline void aypx_vec(double a, __m512d *y, __m512d *x, int Nin)
  {
    __m512d av = _mm512_set1_pd(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmadd_pd(av, y[in], x[in]);
    }
  }

  inline void aypx_vec(double ar, double ai, __m512d *y, __m512d *x, int Nin)
  {
    __m512d avr = _mm512_set1_pd(ar);
    __m512d avi = _mm512_set1_pd(ai);

    for(int in = 0; in < Nin; ++in){
      __m512d yp  = _mm512_permute_pd(y[in], 0b01010101);
      __m512d yt  = _mm512_mul_pd(avi, yp);
      yt = _mm512_fmaddsub_pd(avr, y[in], yt);
      y[in] = _mm512_add_pd(yt, x[in]);
    }
  }

  inline void scal_vec(__m512d *y, double a, int Nin)
  {
    __m512d av = _mm512_set1_pd(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_pd(av, y[in]);
    }
  }

  inline void scal_vec(__m512d *y, double ar, double ai, int Nin)
  {
    __m512d avr = _mm512_set1_pd(ar);
    __m512d avi = _mm512_set1_pd(ai);

    for(int in = 0; in < Nin; ++in){
      __m512d yp  = _mm512_permute_pd(y[in], 0b01010101);
      __m512d yt  = _mm512_mul_pd(avi, yp);
      y[in] = _mm512_fmaddsub_pd(avr, y[in], yt);
    }
  }

  inline void set_vec(__m512d *y, double a, __m512d *x, int Nin)
  {
    __m512d av = _mm512_set1_pd(a);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_pd(av, x[in]);
    }
  }

  inline void unit_vec(__m512d *y, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_set4_pd(0.0, 1.0, 0.0, 1.0);
    }
  }

  inline void conjg_vec(__m512d *y, int Nin)
  {
    __m512d av = _mm512_set1_pd(0.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_fmsubadd_pd(av, av, y[in]);
    }
  }

  inline void real_vec(__m512d *y, int Nin)
  {
    __m512d av = _mm512_set4_pd(0.0, 1.0, 0.0, 1.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_pd(av, y[in]);
    }
  }

  inline void imag_vec(__m512d *y, int Nin)
  {
    __m512d av = _mm512_set4_pd(1.0, 0.0, 1.0, 0.0);
    for(int in = 0; in < Nin; ++in){
      y[in] = _mm512_mul_pd(av, y[in]);
    }
  }

  template<typename REALTYPE>
  inline void xI_vec(__m512d *y, int Nin)
  {
    __m512d av = _mm512_set4_pd(1.0, -1.0, 1.0, -1.0);
    for(int in = 0; in < Nin; ++in){
      __m512d yp = _mm512_permute_pd(y[in], 0b01010101);
      y[in] = _mm512_mul_pd(av, yp);
    }
  }

  template <typename REALTYPE>
  inline void load_vec(__m512d *vt, REALTYPE *v, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      vt[in] = _mm512_load_pd(&v[VLEN*in]);
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
  inline void save_vec(REALTYPE *v, __m512d *vt, int Nin)
  {
    for(int in = 0; in < Nin; ++in){
      _mm512_store_pd(&v[VLEN*in], vt[in]);
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
    inline void shift_vec0_bw(__m512d *v, REALTYPE *w, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in],
			       w[6 + VLEN*in], w[7 + VLEN*in],
			       0.0,  0.0);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec0_fw(__m512d *v, REALTYPE *w, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(0.0, 0.0,
			       w[0 + VLEN*in], w[1 + VLEN*in],
			       w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in]);

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
    inline void shift_vec1_bw(__m512d *v, REALTYPE *w, REALTYPE *buf, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in],
			       w[6 + VLEN*in], w[7 + VLEN*in],
			       buf[2*in],      buf[2*in+1]);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec1_fw(__m512d *v, REALTYPE *w, REALTYPE *buf, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(buf[2*in], buf[2*in+1],
			       w[0 + VLEN*in], w[1 + VLEN*in],
			       w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in]);

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
    inline void shift_vec2_bw(__m512d *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in],
			       w[6 + VLEN*in], w[7 + VLEN*in],
			       y[0 + VLEN*in], y[1 + VLEN*in]);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_bw(__m512d *v, REALTYPE *w, REALTYPE *y,
                              REALTYPE bc, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
        REALTYPE y0 = y[0 + VLEN*in] * bc;
        REALTYPE y1 = y[1 + VLEN*in] * bc;
	v[in] = _mm512_setr_pd(w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in],
			       w[6 + VLEN*in], w[7 + VLEN*in],
			       y0,             y1);
      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_fw(__m512d *v, REALTYPE *w, REALTYPE *y,
                              REALTYPE bc, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
        REALTYPE y0 = y[VLEN-2 + VLEN*in] * bc;
        REALTYPE y1 = y[VLEN-1 + VLEN*in] * bc;
	v[in] = _mm512_setr_pd(y0,             y1,
			       w[0 + VLEN*in], w[1 + VLEN*in],
			       w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in]);

      }
    }

  template <typename REALTYPE>
    inline void shift_vec2_fw(__m512d *v, REALTYPE *w, REALTYPE *y, int Nin)
    {
      for(int in = 0; in < Nin; ++in){
	v[in] = _mm512_setr_pd(y[VLEN-2+VLEN*in], y[VLEN-1 + VLEN*in],
			       w[0 + VLEN*in], w[1 + VLEN*in],
			       w[2 + VLEN*in], w[3 + VLEN*in],
			       w[4 + VLEN*in], w[5 + VLEN*in]);

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
