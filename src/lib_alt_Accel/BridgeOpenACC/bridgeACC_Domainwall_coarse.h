/*!
      @file    mult_Domainwall_coarse_openacc2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2590 $
*/

#ifndef BRIDGEACC_DOMAINWALL_COARSE_INCLUDED
#define BRIDGEACC_DOMAINWALL_COARSE_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_domainwall_coarse_mult_gm5(double *restrict v2,
				   double *restrict v1,
				   int Ncol, int *Nsize);

  void mult_domainwall_coarse_bulk(
                           double *restrict v2, double *restrict u,
                           double *restrict v1,
                           int Ncol, int *Nsize, int *bc);

  void mult_domainwall_coarse_prec(
                           double *restrict v2,
                           double *restrict ct, double *restrict v1,
                           int Ncol, int *Nsize);

  void mult_domainwall_coarse_1(
                   double *restrict buf_xp, double *restrict buf_xm,
                   double *restrict buf_yp, double *restrict buf_ym,
                   double *restrict buf_zp, double *restrict buf_zm,
                   double *restrict buf_tp, double *restrict buf_tm,
                   double *restrict v1,
                   int Ncol, int *Nsize, int *do_comm);

  void mult_domainwall_coarse_2(
                   double *restrict v2, double *restrict u,
                   double *restrict buf_xp, double *restrict buf_xm,
                   double *restrict buf_yp, double *restrict buf_ym,
                   double *restrict buf_zp, double *restrict buf_zm,
                   double *restrict buf_tp, double *restrict buf_tm,
                   int Ncol, int *Nsize, int *do_comm);

  // real_t = float

  void mult_domainwall_coarse_mult_gm5(float *restrict v2,
				   float *restrict v1,
				   int Ncol, int *Nsize);

  void mult_domainwall_coarse_bulk(
                           float *restrict v2, float *restrict u,
                           float *restrict v1,
                           int Ncol, int *Nsize, int *bc);

  void mult_domainwall_coarse_prec(
                           float *restrict v2,
                           float *restrict ct, float *restrict v1,
                           int Ncol, int *Nsize);

  void mult_domainwall_coarse_1(
                   float *restrict buf_xp, float *restrict buf_xm,
                   float *restrict buf_yp, float *restrict buf_ym,
                   float *restrict buf_zp, float *restrict buf_zm,
                   float *restrict buf_tp, float *restrict buf_tm,
                   float *restrict v1,
                   int Ncol, int *Nsize, int *do_comm);

  void mult_domainwall_coarse_2(
                   float *restrict v2, float *restrict u,
                   float *restrict buf_xp, float *restrict buf_xm,
                   float *restrict buf_yp, float *restrict buf_ym,
                   float *restrict buf_zp, float *restrict buf_zm,
                   float *restrict buf_tp, float *restrict buf_tm,
                   int Ncol, int *Nsize, int *do_comm);

}
#endif

//============================================================END=====
