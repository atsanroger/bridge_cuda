/*!
      @file    mult_Clover_coarse.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

#ifndef BRIDGEACC_CLOVER_COARSE_INCLUDED
#define BRIDGEACC_CLOVER_COARSE_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_clover_coarse_mult_gm5(double* v2,
				   double* v1,
				   int Ncol, int *Nsize);

  void mult_clover_coarse_bulk(
                           double* v2, double* u,
                           double* ct, double* v1,
                           int Ncol, int *Nsize, int *bc);

  void mult_clover_coarse_1(
                   double* buf_xp, double* buf_xm,
                   double* buf_yp, double* buf_ym,
                   double* buf_zp, double* buf_zm,
                   double* buf_tp, double* buf_tm,
                   double* u, double* v1,
                   int Ncol, int *Nsize, int *do_comm);

  void mult_clover_coarse_2(
                   double* v2, double* u,
                   double* buf_xp, double* buf_xm,
                   double* buf_yp, double* buf_ym,
                   double* buf_zp, double* buf_zm,
                   double* buf_tp, double* buf_tm,
                   int Ncol, int *Nsize, int *do_comm);

  // real_t = float

  void mult_clover_coarse_mult_gm5(float* v2,
				   float* v1,
				   int Ncol, int *Nsize);

  void mult_clover_coarse_bulk(
                           float* v2, float* u,
                           float* ct, float* v1,
                           int Ncol, int *Nsize, int *bc);

  void mult_clover_coarse_1(
                   float* buf_xp, float* buf_xm,
                   float* buf_yp, float* buf_ym,
                   float* buf_zp, float* buf_zm,
                   float* buf_tp, float* buf_tm,
                   float* u, float* v1,
                   int Ncol, int *Nsize, int *do_comm);

  void mult_clover_coarse_2(
                   float* v2, float* u,
                   float* buf_xp, float* buf_xm,
                   float* buf_yp, float* buf_ym,
                   float* buf_zp, float* buf_zm,
                   float* buf_tp, float* buf_tm,
                   int Ncol, int *Nsize, int *do_comm);

}
#endif

//============================================================END=====
