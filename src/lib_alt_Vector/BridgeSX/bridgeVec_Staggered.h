/*!
      @file    afopr_Staggered_vector.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_STAGGERED_INCLUDED
#define BRIDGEVEC_STAGGERED_INCLUDED

namespace BridgeVec {


// double = double

void mult_staggered_D(double *RESTRICT v2, double *RESTRICT u,
                      double *RESTRICT v1,
                      int *Nsize, int *bc, double mq, int jdag);

void mult_staggered_1(
                  double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                  double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                  double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                  double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                  double *RESTRICT u, double *RESTRICT v1,
                  int *Nsize, int *bc, int *do_comm);
 
void mult_staggered_2(double *RESTRICT v2, double *RESTRICT u,
                      double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                      double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                      double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                      double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                      int *Nsize, int *bc,  int *do_comm, double mq, int jdag);


void mult_staggered_phase_dev(double *RESTRICT u, double *RESTRICT ph,
                              int *Nsize, int Nc);

void mult_staggered_gm5(double *RESTRICT v2, double *RESTRICT v1,
                       int *Nsize, int Nc);

void mult_staggered_aypx(double a, double *RESTRICT v2, double *RESTRICT v1,
                         int *Nsize, int Nc);

void mult_staggered_axpy(double *RESTRICT v2, double a, double *RESTRICT v1,
                         int *Nsize, int Nc);

void mult_staggered_scal(double *RESTRICT v, double a,
                         int *Nsize, int Nc);

void mult_staggered_clear(double *RESTRICT v2, int *Nsize, int Nc);

void mult_staggered_xp1(double *RESTRICT buf, double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xp2(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xpb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xm1(double *RESTRICT buf, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xm2(double *RESTRICT v2, double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xmb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_yp1(double *RESTRICT buf, double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_yp2(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ypb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ym1(double *RESTRICT buf, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ym2(double *RESTRICT v2, double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

 void mult_staggered_ymb(double *RESTRICT v2, double *RESTRICT u,
                         double *RESTRICT v1,
                         int *Nsize, int *bc, int Nc);

void mult_staggered_zp1(double *RESTRICT buf, double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zp2(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zpb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zm1(double *RESTRICT buf, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zm2(double *RESTRICT v2, double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zmb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tp1(double *RESTRICT buf, double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tp2(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tpb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tm1(double *RESTRICT buf, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tm2(double *RESTRICT v2, double *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tmb(double *RESTRICT v2, double *RESTRICT u,
                        double *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

// real_t = float

void mult_staggered_D(float *RESTRICT v2, float *RESTRICT u,
                      float *RESTRICT v1,
                      int *Nsize, int *bc, float mq, int jdag);

void mult_staggered_1(
                  float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                  float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                  float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                  float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                  float *RESTRICT u, float *RESTRICT v1,
                  int *Nsize, int *bc, int *do_comm);
 
void mult_staggered_2(float *RESTRICT v2, float *RESTRICT u,
                      float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                      float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                      float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                      float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                      int *Nsize, int *bc,  int *do_comm, float mq, int jdag);


void mult_staggered_phase_dev(float *RESTRICT u, float *RESTRICT ph,
                              int *Nsize, int Nc);

void mult_staggered_gm5(float *RESTRICT v2, float *RESTRICT v1,
                       int *Nsize, int Nc);

void mult_staggered_aypx(float a, float *RESTRICT v2, float *RESTRICT v1,
                         int *Nsize, int Nc);

void mult_staggered_axpy(float *RESTRICT v2, float a, float *RESTRICT v1,
                         int *Nsize, int Nc);

void mult_staggered_scal(float *RESTRICT v, float a,
                         int *Nsize, int Nc);

void mult_staggered_clear(float *RESTRICT v2, int *Nsize, int Nc);

void mult_staggered_xp1(float *RESTRICT buf, float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xp2(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xpb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xm1(float *RESTRICT buf, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xm2(float *RESTRICT v2, float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_xmb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_yp1(float *RESTRICT buf, float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_yp2(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ypb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ym1(float *RESTRICT buf, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_ym2(float *RESTRICT v2, float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

 void mult_staggered_ymb(float *RESTRICT v2, float *RESTRICT u,
                         float *RESTRICT v1,
                         int *Nsize, int *bc, int Nc);

void mult_staggered_zp1(float *RESTRICT buf, float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zp2(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

 void mult_staggered_zpb(float *RESTRICT v2, float *RESTRICT u,
                         float *RESTRICT v1,
                         int *Nsize, int *bc, int Nc);

void mult_staggered_zm1(float *RESTRICT buf, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zm2(float *RESTRICT v2, float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_zmb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tp1(float *RESTRICT buf, float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tp2(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tpb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tm1(float *RESTRICT buf, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tm2(float *RESTRICT v2, float *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc);

void mult_staggered_tmb(float *RESTRICT v2, float *RESTRICT u,
                        float *RESTRICT v1,
                        int *Nsize, int *bc, int Nc);

}

#endif
//============================================================END=====

