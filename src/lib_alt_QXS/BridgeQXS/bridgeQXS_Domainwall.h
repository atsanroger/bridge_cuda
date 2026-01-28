/*!
      @file    bridgeQXS_Domainwall.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

#ifndef BRIDGEQXS_DOMAINWALL_INCLUDED
#define BRIDGEQXS_DOMAINWALL_INCLUDED

namespace BridgeQXS {
  // real_t = double

  void mult_domainwall_5din_5dir_dirac(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void accum_mult_domainwall_5din_5dir_dirac(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(double *vp, double *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_clear(double *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    double *vp, double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    double *vp, double *up, double *wp,
    double *buf2_xp, double *buf2_xm,
    double *buf2_yp, double *buf2_ym,
    double *buf2_zp, double *buf2_zm,
    double *buf2_tp, double *buf2_tm,
    double mq, double M0, int Ns, int *bc,
    int *Nsize, int *do_comm);


  void mult_domainwall_5din_eo_5dir_dirac(double *yp, double *wp,
                                          double mq, double M0, int Ns, int *bc,
                                          double *b, double *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(double *vp, double *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(double *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    double *vp, double *yp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    double *vp, double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    double *buf1_xp, double *buf1_xm,
    double *buf1_yp, double *buf1_ym,
    double *buf1_zp, double *buf1_zm,
    double *buf1_tp, double *buf1_tm,
    double *up, double *wp,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    double *vp, double *up, double *wp,
    double *buf2_xp, double *buf2_xm,
    double *buf2_yp, double *buf2_ym,
    double *buf2_zp, double *buf2_zm,
    double *buf2_tp, double *buf2_tm,
    double mq, double M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    double *vp, double *up, double *wp,
    double *yp,
    double mq, double M0, int Ns, int *bc,
    double *b, double *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  // real_t = float

  void mult_domainwall_5din_5dir_dirac(
    float *vp, float *yp, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void accum_mult_domainwall_5din_5dir_dirac(
    float *vp, float *yp, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    float *vp, float *yp, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(float *vp, float *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_clear(float *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    float *vp, float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    float *vp, float *up, float *wp,
    float *buf2_xp, float *buf2_xm,
    float *buf2_yp, float *buf2_ym,
    float *buf2_zp, float *buf2_zm,
    float *buf2_tp, float *buf2_tm,
    float mq, float M0, int Ns, int *bc,
    int *Nsize, int *do_comm);


  void mult_domainwall_5din_eo_5dir_dirac(float *yp, float *wp,
                                          float mq, float M0, int Ns, int *bc,
                                          float *b, float *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(float *vp, float *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(float *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    float *vp, float *yp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    float *vp, float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    float *buf1_xp, float *buf1_xm,
    float *buf1_yp, float *buf1_ym,
    float *buf1_zp, float *buf1_zm,
    float *buf1_tp, float *buf1_tm,
    float *up, float *wp,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    float *vp, float *up, float *wp,
    float *buf2_xp, float *buf2_xm,
    float *buf2_yp, float *buf2_ym,
    float *buf2_zp, float *buf2_zm,
    float *buf2_tp, float *buf2_tm,
    float mq, float M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    float *vp, float *up, float *wp,
    float *yp,
    float mq, float M0, int Ns, int *bc,
    float *b, float *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  // real_t
  void mult_domainwall_5din_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);


  void mult_domainwall_5din_dd_hopb_dirac(
    real_t *vp, real_t *up, real_t *wp,
    real_t mq, real_t M0, int Ns, int *bc,
    real_t *b, real_t *c,
    int *Nsize, int *block_size,
    int ieo);


  void mult_domainwall_5din_dd_5dir_dirac(
    real_t *vp, real_t *yp, real_t *wp,
    real_t mq, real_t M0, int Ns, int *bc,
    real_t *b, real_t *c,
    int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_5dirdag_dirac(
    real_t *vp, real_t *yp, real_t *wp,
    real_t mq, real_t M0, int Ns, int *bc,
    real_t *b, real_t *c,
    int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_mult_gm5_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    int *block_size, int ieo);

  // real_t = half

#ifdef USE_FP16
  void mult_domainwall_5din_5dir_dirac(
    half *vp, half *yp, half *wp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Nsize, int *do_comm);

  void accum_mult_domainwall_5din_5dir_dirac(
    half *vp, half *yp, half *wp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_5dirdag_dirac(
    half *vp, half *yp, half *wp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_mult_gm5_dirac(half *vp, half *wp,
                                           int Ns, int *Nsize);

  void mult_domainwall_5din_mult_gm5R_dirac(half *vp, half *wp,
                                            int Ns, int *Nsize);

  void mult_domainwall_5din_mult_R(half *vp, half *wp,
                                   int Ns, int *Nsize);

  void mult_domainwall_5din_clear(half *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
    half *vp, half *up, half *wp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop1_dirac(
    half *buf1_xp, half *buf1_xm,
    half *buf1_yp, half *buf1_ym,
    half *buf1_zp, half *buf1_zm,
    half *buf1_tp, half *buf1_tm,
    half *up, half *wp,
    half mq, half M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
    half *vp, half *up, half *wp,
    half *buf2_xp, half *buf2_xm,
    half *buf2_yp, half *buf2_ym,
    half *buf2_zp, half *buf2_zm,
    half *buf2_tp, half *buf2_tm,
    half mq, half M0, int Ns, int *bc,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_5dir_dirac(half *yp, half *wp,
                                          half mq, half M0, int Ns, int *bc,
                                          half *b, half *c,
                                          int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_mult_gm5_dirac(half *vp, half *wp,
                                              int Ns, int *Nsize);

  void mult_domainwall_5din_eo_clear(half *vp, int Ns, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
    half *vp, half *yp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Nsize, int *do_comm);

  void mult_domainwall_5din_eo_hopb_dirac(
    half *vp, half *up, half *wp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac(
    half *buf1_xp, half *buf1_xm,
    half *buf1_yp, half *buf1_ym,
    half *buf1_zp, half *buf1_zm,
    half *buf1_tp, half *buf1_tm,
    half *up, half *wp,
    half mq, half M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop1_dirac_res(
    half *buf1_xp, half *buf1_xm,
    half *buf1_yp, half *buf1_ym,
    half *buf1_zp, half *buf1_zm,
    half *buf1_tp, half *buf1_tm,
    half *up, half *wp,
    half mq, half M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_hop2_dirac(
    half *vp, half *up, half *wp,
    half *buf2_xp, half *buf2_xm,
    half *buf2_yp, half *buf2_ym,
    half *buf2_zp, half *buf2_zm,
    half *buf2_tp, half *buf2_tm,
    half mq, half M0, int Ns, int *bc,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);

  void mult_domainwall_5din_eo_bulk_dirac(
    half *vp, half *up, half *wp,
    half *yp,
    half mq, half M0, int Ns, int *bc,
    half *b, half *c,
    int *Leo, int *Nsize, int *do_comm,
    const int ieo);
#endif

  // real_t
  void mult_domainwall_5din_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_L_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_U_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Ldag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *e, real_t *dpinv, real_t *dm);

  void mult_domainwall_5din_eo_Udag_inv_dirac(
    real_t *vp, real_t *wp,
    int Ns, int *Nsize,
    real_t *f, real_t *dpinv, real_t *dm);

}
#endif
