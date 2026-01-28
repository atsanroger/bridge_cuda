/*!
      @file    bridgeQXS_Staggered.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_STAGGERED_INCLUDED
#define BRIDGEQXS_STAGGERED_INCLUDED

namespace BridgeQXS {
  // real_t = double

  void mult_staggered_clear(double *v, int *Nsize, int Nin);

  // y = b * y + a *x
  void mult_staggered_axpby(double b, double *y, double a, double *x,
                            int *Nsize, int Nin);

  void mult_staggered_phase(double *v, double *ph, int *Nsize, int Nin);

  // v = U * w
  void mult_staggered_mult_Gn(double *v, double *u, double *w, int *Nsize);

  // v = U^d * w
  void mult_staggered_mult_Gd(double *v, double *u, double *w, int *Nsize);

  void mult_staggered_bulk(double *v2, double *up, double *v1,
                           double qm, int jd,
                           int *Nsize, int *do_comm);

  void mult_staggered_1(double *buf_xp, double *buf_xm,
                        double *buf_yp, double *buf_ym,
                        double *buf_zp, double *buf_zm,
                        double *buf_tp, double *buf_tm,
                        double *up, double *v1,
                        int *Nsize, int *do_comm);

  void mult_staggered_2(double *v2, double *up, double *v1,
                        double *buf_xp, double *buf_xm,
                        double *buf_yp, double *buf_ym,
                        double *buf_zp, double *buf_zm,
                        double *buf_tp, double *buf_tm,
                        double qm, int jd,
                        int *Nsize, int *do_comm);

  void mult_staggered_eo_bulk(double *v2, double *up,
                              double *v1, double *xp,
                              double qm, int jd,
                              int *Nsize, int *do_comm,
                              int *Leo, int ieo, int iflag);

  void mult_staggered_eo_1(double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double *up, double *v1,
                           int *Nsize, int *do_comm,
                           int *Leo, int ieo);

  void mult_staggered_eo_2(double *v2, double *up, double *v1,
                           double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double qm, int *Nsize, int *do_comm,
                           int *Leo, int ieo, int iflag);

  // real_t = float

  void mult_staggered_clear(float *v, int *Nsize, int Nin);

  void mult_staggered_axpby(float b, float *v, float a, float *w,
                            int *Nsize, int Nin);

  void mult_staggered_phase(float *v, float *ph, int *Nsize, int Nin);

  // v = U * w
  void mult_staggered_mult_Gn(float *v, float *u, float *w, int *Nsize);

  // v = U^d * w
  void mult_staggered_mult_Gd(float *v, float *u, float *w, int *Nsize);

  void mult_staggered_bulk(float *v2, float *up, float *v1,
                           float qm, int jd,
                           int *Nsize, int *do_comm);

  void mult_staggered_1(float *buf_xp, float *buf_xm,
                        float *buf_yp, float *buf_ym,
                        float *buf_zp, float *buf_zm,
                        float *buf_tp, float *buf_tm,
                        float *up, float *v1,
                        int *Nsize, int *do_comm);

  void mult_staggered_2(float *v2, float *up, float *v1,
                        float *buf_xp, float *buf_xm,
                        float *buf_yp, float *buf_ym,
                        float *buf_zp, float *buf_zm,
                        float *buf_tp, float *buf_tm,
                        float qm, int jd,
                        int *Nsize, int *do_comm);

  void mult_staggered_eo_bulk(float *v2, float *up,
                              float *v1, float *xp,
                              float qm, int jd,
                              int *Nsize, int *do_comm,
                              int *Leo, int ieo, int iflag);

  void mult_staggered_eo_1(float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float *up, float *v1,
                           int *Nsize, int *do_comm,
                           int *Leo, int ieo);

  void mult_staggered_eo_2(float *v2, float *up, float *v1,
                           float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float qm, int *Nsize, int *do_comm,
                           int *Leo, int ieo, int iflag);
}
#endif
