/*!
      @file    bridgeQXS_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

#ifndef BRIDGEQXS_WILSON_INCLUDED
#define BRIDGEQXS_WILSON_INCLUDED

namespace BridgeQXS {
  // real_t = double

  void mult_wilson_bulk_dirac(double *v2, double *up, double *v1,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double *up, double *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(double *v2, double *up, double *v1,
                           double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(double *v2, double *up, double *v1,
                                 double *xp,
                                 double kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double *up, double *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(double *v2, double *up, double *v1,
                              double *xp,
                              double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(double *v2, double *v1, int *Nsize);

  // real_t = float

  void mult_wilson_bulk_dirac(float *v2, float *up, float *v1,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float *up, float *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(float *v2, float *up, float *v1,
                           float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(float *v2, float *up, float *v1,
                                 float *xp,
                                 float kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float *up, float *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(float *v2, float *up, float *v1,
                              float *xp,
                              float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(float *v2, float *v1, int *Nsize);

#ifdef USE_FP16
  // real_t = half

  void mult_wilson_bulk_dirac(half *v2, half *up, half *v1,
                              half kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(half *buf_xp, half *buf_xm,
                           half *buf_yp, half *buf_ym,
                           half *buf_zp, half *buf_zm,
                           half *buf_tp, half *buf_tm,
                           half *up, half *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(half *v2, half *up, half *v1,
                           half *buf_xp, half *buf_xm,
                           half *buf_yp, half *buf_ym,
                           half *buf_zp, half *buf_zm,
                           half *buf_tp, half *buf_tm,
                           half kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(half *v2, half *up, half *v1,
                                 half *xp,
                                 half kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(half *buf_xp, half *buf_xm,
                              half *buf_yp, half *buf_ym,
                              half *buf_zp, half *buf_zm,
                              half *buf_tp, half *buf_tm,
                              half *up, half *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(half *v2, half *up, half *v1,
                              half *xp,
                              half *buf_xp, half *buf_xm,
                              half *buf_yp, half *buf_ym,
                              half *buf_zp, half *buf_zm,
                              half *buf_tp, half *buf_tm,
                              half kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(half *v2, half *v1, int *Nsize);
#endif

}
#endif
