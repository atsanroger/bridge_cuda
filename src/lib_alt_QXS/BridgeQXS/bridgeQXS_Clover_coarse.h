/*!
      @file    bridgeQXS_Clover.h
      @brief
      @author  Isaku Kanamori (kanamori)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_CLOVER_COARSE_INCLUDED
#define BRIDGEQXS_CLOVER_COARSE_INCLUDED
#include <vector>

namespace BridgeQXS {
  void mult_coarse_1(double *buf1_xp, double *buf1_xm,
                     double *buf1_yp, double *buf1_ym,
                     double *buf1_zp, double *buf1_zm,
                     double *buf1_tp, double *buf1_tm,
                     double *u0, double *v1, const int *Nsize,
                     int ncol, const int *do_comm);


  void mult_coarse_b(double *v2,
                     double *u0, double *c0,
                     double *v1,
                     const int *Nsize, int ncol,
                     const int *do_comm, double *work);


  void mult_coarse_2(double *v2, double *u0, double *v1,
                     double *buf2_xp, double *buf2_xm,
                     double *buf2_yp, double *buf2_ym,
                     double *buf2_zp, double *buf2_zm,
                     double *buf2_tp, double *buf2_tm,
                     const int *Nsize, int ncol, const int *do_comm,
                     double *work,
                     std::vector<int>& list);


  void mult_coarse_1(float *buf1_xp, float *buf1_xm,
                     float *buf1_yp, float *buf1_ym,
                     float *buf1_zp, float *buf1_zm,
                     float *buf1_tp, float *buf1_tm,
                     float *u0, float *v1, const int *Nsize,
                     int ncol, const int *do_comm);


  void mult_coarse_b(float *v2,
                     float *u0, float *c0,
                     float *v1,
                     const int *Nsize, int ncol,
                     const int *do_comm, float *work);


  void mult_coarse_2(float *v2, float *u0, float *v1,
                     float *buf2_xp, float *buf2_xm,
                     float *buf2_yp, float *buf2_ym,
                     float *buf2_zp, float *buf2_zm,
                     float *buf2_tp, float *buf2_tm,
                     const int *Nsize, int ncol, const int *do_comm,
                     float *work,
                     std::vector<int>& list);
}
#endif
