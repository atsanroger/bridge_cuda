/*!
      @file    afield_fp16_acle-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-01-13 22:19:26 #$
      @version $LastChangedRevision: 2569 $
*/

#include <cassert>

#include "lib_alt_QXS/inline/afield_th-inc.h"


//====================================================================
namespace {

  inline void convert_to32(svbool_t& pg32,
                           svfloat32_t& vt1, svfloat32_t& vt2,
                           svfloat16_t wt)
  {
     svfloat16_t wr = svrev_f16(wt);
     vt1 = svcvt_f32_f16_z(pg32, wt);
     vt2 = svcvt_f32_f16_z(pg32, wr);
  }
}

//====================================================================
template<>
//real_t AField<real_t, QXS>::norm2() const
void AField<real_t, QXS>::norm2_double(double& rr) const
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  real_t *vp = const_cast<AField<real_t, QXS> *>(this)->ptr(0);

#pragma omp barrier

  int Nouter = (Nstv > Nin) ? Nstv : Nin;
  int Ninner = (Nstv * Nin) / Nouter;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg   = set_predicate();
  svbool_t pg32 = svptrue_b32();

  svfloat32_t yt;
  yt = svdup_f32_m(yt, pg32, 0.0);  // setting vector to 0.0

  for (int ex = 0; ex < Nex; ++ex) {
    svfloat32_t tmp_yt;
    tmp_yt = svdup_f32_m(tmp_yt, pg32, 0.0);  // setting vector to 0.0

    for (int i = is; i < ns; ++i) {
      real_t   *v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      svfloat32_t xt1, xt2;
      xt1 = svdup_f32_m(xt1, pg32, 0.0);  // setting vector to 0.0
      xt2 = svdup_f32_m(xt2, pg32, 0.0);  // setting vector to 0.0
      for (int in = 0; in < Ninner; ++in) {
        svfloat16_t vt;
        load_vec(pg, vt, &v[VLEN * in]);
        svfloat32_t vt1, vt2;
        convert_to32(pg32, vt1, vt2, vt);
        xt1 = svmla_m(pg32, xt1, vt1, vt1);
        xt2 = svmla_m(pg32, xt2, vt2, vt2);
      }
      tmp_yt = svadd_m(pg32, tmp_yt, xt1);
      tmp_yt = svadd_m(pg32, tmp_yt, xt2);
    }
    yt = svadd_m(pg32, yt, tmp_yt);
  }

  float at = svaddv(pg32, yt);   // reduction of elements

#pragma omp barrier

  double ad = double(at);
  ThreadManager::reduce_sum_global(ad, ith, nth);

  rr = ad;
}

//====================================================================
template<>
real_t AField<real_t, QXS>::norm2() const
{
  double rr;
  norm2_double(rr);

  return real_t(rr);
}

//====================================================================
template<>
void AField<real_t, QXS>::dot_double(double& prod,
				     const AField<real_t, QXS>& w) const
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  real_t *RESTRICT vp = const_cast<AField<real_t, QXS> *>(this)->ptr(0);
  real_t *RESTRICT wp = const_cast<AField<real_t, QXS> *>(&w)->ptr(0);
  assert(vp != wp);

#pragma omp barrier

  svbool_t pg = set_predicate();
  svbool_t pg32 = svptrue_b32();

  svfloat32_t yt;
  yt = svdup_f32_m(yt, pg32, 0.0);  // setting vector to 0.0


  for (int ex = 0; ex < Nex; ++ex) {
    svfloat32_t tmp_yt;
    tmp_yt = svdup_f32_m(tmp_yt, pg32, 0.0);  // setting vector to 0.0

    for (int site = is; site < ns; ++site) {
      real_t* v = &vp[VLEN * Nin * (site + Nstv * ex)];
      real_t* w = &wp[VLEN * Nin * (site + Nstv * ex)];
      svfloat32_t xt1, xt2;
      xt1 = svdup_f32_m(xt1, pg32, 0.0);  // setting vector to 0.0
      xt2 = svdup_f32_m(xt2, pg32, 0.0);  // setting vector to 0.0

      for (int in = 0; in < Nin; ++in) {
        svfloat16_t vt, wt;
        load_vec(pg, wt, &w[VLEN * in]);
        load_vec(pg, vt, &v[VLEN * in]);
        svfloat32_t wt1, wt2, vt1, vt2;
        convert_to32(pg32, wt1, wt2, wt);
        convert_to32(pg32, vt1, vt2, vt);
        xt1 = svmla_m(pg32, xt1, vt1, wt1);
        xt2 = svmla_m(pg32, xt2, vt2, wt2);
      }
      tmp_yt = svadd_m(pg32, tmp_yt, xt1);
      tmp_yt = svadd_m(pg32, tmp_yt, xt2);
    }

    yt = svadd_m(pg32, yt, tmp_yt);

  }

  float at = svaddv(pg32, yt);   // reduction of elements

#pragma omp barrier

  double ad = double(at);
  ThreadManager::reduce_sum_global(ad, ith, nth);

  prod = ad;
}

//====================================================================
template<>
real_t AField<real_t, QXS>::dot(const AField<real_t, QXS>& w) const
{
  double prod;
  dot_double(prod, w);

  return real_t(prod);
  
}

//====================================================================
template<>
void AField<real_t, QXS>::dotc(real_t& ar, real_t& ai,
                                 const AField<real_t, QXS>& w) const
{
  assert(check_size(w));

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;
  int Nin2 = Nin / 2;

  real_t *RESTRICT vp = const_cast<AField<real_t, QXS> *>(this)->ptr(0);
  real_t *RESTRICT wp = const_cast<AField<real_t, QXS> *>(&w)->ptr(0);
  assert(vp != wp);

#pragma omp barrier

  int Nouter  = (Nstv > Nin2) ? Nstv : Nin2;
  int Ninner2 = (Nstv * Nin2) / Nouter;
  int Ninner  = 2 * Ninner2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg = set_predicate();
  svbool_t pg32 = svptrue_b32();

  svfloat32_t ytr, yti;
  ytr = svdup_f32_m(ytr, pg32, 0.0);  // setting vector to 0.0
  yti = svdup_f32_m(yti, pg32, 0.0);  // setting vector to 0.0

  for (int ex = 0; ex < Nex; ++ex) {
    svfloat32_t tmp_ytr, tmp_yti;
    tmp_ytr = svdup_f32_m(tmp_ytr, pg32, 0.0);
    tmp_yti = svdup_f32_m(tmp_yti, pg32, 0.0);

    for (int i = is; i < ns; ++i) {
      real_t* v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      real_t* w = &wp[VLEN * Ninner * (i + Nouter * ex)];

      svfloat32_t xtr, xti;
      xtr = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0
      xti = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0

      for (int in = 0; in < Ninner2; ++in) {
        svfloat16_t vtr, vti, wtr, wti;
        int inr = 2 * in;
        int ini = 2 * in + 1;
        load_vec(pg, wtr, &w[VLEN * inr]);
        load_vec(pg, wti, &w[VLEN * ini]);
        load_vec(pg, vtr, &v[VLEN * inr]);
        load_vec(pg, vti, &v[VLEN * ini]);
        svfloat32_t vtr1, vtr2, vti1, vti2;
        svfloat32_t wtr1, wtr2, wti1, wti2;
        convert_to32(pg32, vtr1, vtr2, vtr);
        convert_to32(pg32, vti1, vti2, vti);
        convert_to32(pg32, wtr1, wtr2, wtr);
        convert_to32(pg32, wti1, wti2, wti);

        xtr = svmla_m(pg32, xtr, vtr1, wtr1);
        xtr = svmla_m(pg32, xtr, vti1, wti1);
        xtr = svmla_m(pg32, xtr, vtr2, wtr2);
        xtr = svmla_m(pg32, xtr, vti2, wti2);

        xti = svmla_m(pg32, xti, vtr1, wti1);
        xti = svmls_m(pg32, xti, vti1, wtr1);
        xti = svmla_m(pg32, xti, vtr2, wti2);
        xti = svmls_m(pg32, xti, vti2, wtr2);

      }
      tmp_ytr = svadd_m(pg32, tmp_ytr, xtr);
      tmp_yti = svadd_m(pg32, tmp_yti, xti);
    }

    ytr = svadd_m(pg32, ytr, tmp_ytr);
    yti = svadd_m(pg32, yti, tmp_yti);
  }

  float atr = svaddv(pg32, ytr);   // reduction of elements
  float ati = svaddv(pg32, yti);   // reduction of elements

#pragma omp barrier

  double sum[2] = { double(atr), double(ati) };
  ThreadManager::reduce_sum_global(sum, 2, ith, nth);
  ar = real_t(sum[0]);
  ai = real_t(sum[1]);

}

//====================================================================
template<>
void AField<real_t, QXS>::dotc_double(std::complex<double>& a,
                                 const AField<real_t, QXS>& w) const
{
  assert(check_size(w));

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;
  int Nin2 = Nin / 2;

  real_t *RESTRICT vp = const_cast<AField<real_t, QXS> *>(this)->ptr(0);
  real_t *RESTRICT wp = const_cast<AField<real_t, QXS> *>(&w)->ptr(0);
  assert(vp != wp);

#pragma omp barrier

  int Nouter  = (Nstv > Nin2) ? Nstv : Nin2;
  int Ninner2 = (Nstv * Nin2) / Nouter;
  int Ninner  = 2 * Ninner2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg = set_predicate();
  svbool_t pg32 = svptrue_b32();

  svfloat32_t ytr, yti;
  ytr = svdup_f32_m(ytr, pg32, 0.0);  // setting vector to 0.0
  yti = svdup_f32_m(yti, pg32, 0.0);  // setting vector to 0.0

  for (int ex = 0; ex < Nex; ++ex) {
    svfloat32_t tmp_ytr, tmp_yti;
    tmp_ytr = svdup_f32_m(tmp_ytr, pg32, 0.0);
    tmp_yti = svdup_f32_m(tmp_yti, pg32, 0.0);

    for (int i = is; i < ns; ++i) {
      real_t* v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      real_t* w = &wp[VLEN * Ninner * (i + Nouter * ex)];

      svfloat32_t xtr, xti;
      xtr = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0
      xti = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0

      for (int in = 0; in < Ninner2; ++in) {
        svfloat16_t vtr, vti, wtr, wti;
        int inr = 2 * in;
        int ini = 2 * in + 1;
        load_vec(pg, wtr, &w[VLEN * inr]);
        load_vec(pg, wti, &w[VLEN * ini]);
        load_vec(pg, vtr, &v[VLEN * inr]);
        load_vec(pg, vti, &v[VLEN * ini]);
        svfloat32_t vtr1, vtr2, vti1, vti2;
        svfloat32_t wtr1, wtr2, wti1, wti2;
        convert_to32(pg32, vtr1, vtr2, vtr);
        convert_to32(pg32, vti1, vti2, vti);
        convert_to32(pg32, wtr1, wtr2, wtr);
        convert_to32(pg32, wti1, wti2, wti);

        xtr = svmla_m(pg32, xtr, vtr1, wtr1);
        xtr = svmla_m(pg32, xtr, vti1, wti1);
        xtr = svmla_m(pg32, xtr, vtr2, wtr2);
        xtr = svmla_m(pg32, xtr, vti2, wti2);

        xti = svmla_m(pg32, xti, vtr1, wti1);
        xti = svmls_m(pg32, xti, vti1, wtr1);
        xti = svmla_m(pg32, xti, vtr2, wti2);
        xti = svmls_m(pg32, xti, vti2, wtr2);

      }
      tmp_ytr = svadd_m(pg32, tmp_ytr, xtr);
      tmp_yti = svadd_m(pg32, tmp_yti, xti);
    }

    ytr = svadd_m(pg32, ytr, tmp_ytr);
    yti = svadd_m(pg32, yti, tmp_yti);
  }

  float atr = svaddv(pg32, ytr);   // reduction of elements
  float ati = svaddv(pg32, yti);   // reduction of elements

#pragma omp barrier

  double sum[2] = { double(atr), double(ati) };
  ThreadManager::reduce_sum_global(sum, 2, ith, nth);
  //ar = real_t(sum[0]);
  //ai = real_t(sum[1]);

  a = cmplx(sum[0], sum[1]);
  
}

//====================================================================
/*
template<>
typename AField<real_t, QXS>::complex_t
AField<real_t, QXS>::dotc_and_norm2(real_t& norm2, real_t& w_norm2,
                                    const AField<real_t, QXS>& w) const
{
  if (m_element_type == Element_type::REAL) {
    vout.crucial("%s: dotc_and_norm2 for real field is unsupported\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  assert(check_size(w));

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;
  int Nin2 = Nin / 2;

  real_t *RESTRICT vp = const_cast<AField<real_t, QXS> *>(this)->ptr(0);
  real_t *RESTRICT wp = const_cast<AField<real_t, QXS> *>(&w)->ptr(0);
  assert(vp != wp);

#pragma omp barrier

  int Nouter  = (Nstv > Nin2) ? Nstv : Nin2;
  int Ninner2 = (Nstv * Nin2) / Nouter;
  int Ninner  = 2 * Ninner2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg = set_predicate();
  svbool_t pg32 = svptrue_b32();

  svfloat32_t ytr, yti, w2, v2;
  ytr = svdup_f32_m(ytr, pg32, 0.0);  // setting vector to 0.0
  yti = svdup_f32_m(yti, pg32, 0.0);  // setting vector to 0.0
  w2  = svdup_f32_m(w2,  pg32, 0.0);  // setting vector to 0.0
  v2  = svdup_f32_m(v2,  pg32, 0.0);  // setting vector to 0.0

  for (int ex = 0; ex < Nex; ++ex) {

    svfloat32_t tmp_ytr, tmp_yti, tmp_vt2, tmp_wt2;
    tmp_ytr = svdup_f32_m(tmp_ytr, pg32, 0.0);
    tmp_yti = svdup_f32_m(tmp_yti, pg32, 0.0);
    tmp_vt2 = svdup_f32_m(tmp_vt2, pg32, 0.0);
    tmp_wt2 = svdup_f32_m(tmp_wt2, pg32, 0.0);

    for (int iout = is; iout < ns; ++iout) {
      real_t   *v = &vp[VLEN * Ninner * (iout + Nouter * ex)];
      real_t   *w = &wp[VLEN * Ninner * (iout + Nouter * ex)];

      svfloat32_t xtr, xti, vt2, wt2;
      xtr = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0
      xti = svdup_f32_m(xtr, pg32, 0.0);  // setting vector to 0.0
      vt2 = svdup_f32_m(vt2, pg32, 0.0);  // setting vector to 0.0
      wt2 = svdup_f32_m(wt2, pg32, 0.0);  // setting vector to 0.0

      svfloat16_t vtr, vti, wtr, wti;
      for (int in = 0; in < Ninner2; ++in) {
        int inr = 2 * in;
        int ini = 2 * in + 1;
        load_vec(pg, wtr, &w[VLEN * inr]);
        load_vec(pg, wti, &w[VLEN * ini]);
        load_vec(pg, vtr, &v[VLEN * inr]);
        load_vec(pg, vti, &v[VLEN * ini]);

        svfloat32_t vtr1, vtr2, vti1, vti2;
        svfloat32_t wtr1, wtr2, wti1, wti2;
        convert_to32(pg32, vtr1, vtr2, vtr);
        convert_to32(pg32, vti1, vti2, vti);
        convert_to32(pg32, wtr1, wtr2, wtr);
        convert_to32(pg32, wti1, wti2, wti);

        xtr = svmla_m(pg32, xtr, vtr1, wtr1);
        xtr = svmla_m(pg32, xtr, vti1, wti1);
        xtr = svmla_m(pg32, xtr, vtr2, wtr2);
        xtr = svmla_m(pg32, xtr, vti2, wti2);

        xti = svmla_m(pg32, xti, vtr1, wti1);
        xti = svmls_m(pg32, xti, vti1, wtr1);
        xti = svmla_m(pg32, xti, vtr2, wti2);
        xti = svmls_m(pg32, xti, vti2, wtr2);

        vt2 = svmla_m(pg32, vt2, vtr1, vtr1);
        vt2 = svmla_m(pg32, vt2, vti1, vti1);
        vt2 = svmla_m(pg32, vt2, vtr2, vtr2);
        vt2 = svmla_m(pg32, vt2, vti2, vti2);

        wt2 = svmla_m(pg32, wt2, wtr1, wtr1);
        wt2 = svmla_m(pg32, wt2, wti1, wti1);
        wt2 = svmla_m(pg32, wt2, wtr2, wtr2);
        wt2 = svmla_m(pg32, wt2, wti2, wti2);
      }
      tmp_ytr = svadd_m(pg32, tmp_ytr, xtr);
      tmp_yti = svadd_m(pg32, tmp_yti, xti);
      tmp_vt2 = svadd_m(pg32, tmp_vt2, vt2);
      tmp_wt2 = svadd_m(pg32, tmp_wt2, wt2);
    }
    ytr = svadd_m(pg32, ytr, tmp_ytr);
    yti = svadd_m(pg32, yti, tmp_yti);
    v2  = svadd_m(pg32, v2,  tmp_vt2);
    w2  = svadd_m(pg32, w2,  tmp_wt2);
  }

  float atr = svaddv(pg32, ytr);   // reduction of elements
  float ati = svaddv(pg32, yti);   // reduction of elements
  float wwt = svaddv(pg32, w2);    // reduction of elements
  float vvt = svaddv(pg32, v2);    // reduction of elements

  double sum[4] = { double(atr), double(ati),
                    double(wwt), double(vvt) };

#pragma omp barrier

  ThreadManager::reduce_sum_global(sum, 4, ith, nth);

  complex_t result = { real_t(sum[0]), real_t(sum[1]) };

  w_norm2 = real_t(sum[2]);
  norm2   = real_t(sum[3]);

  return result;
}
*/
//============================================================END=====
