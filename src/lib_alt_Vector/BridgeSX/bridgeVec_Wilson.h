/*!
      @file    bridgeVec_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_WILSON_INCLUDED
#define BRIDGEVEC_WILSON_INCLUDED

namespace BridgeVec {

void mult_wilson_H_dirac(double *RESTRICT v2, double *RESTRICT u,
                         double *RESTRICT v1,
                         int *Nsize, int *bc, double kappa);

void mult_wilson_D_dirac(double *RESTRICT v2, double *RESTRICT u,
                         double *RESTRICT v1,
                         int *Nsize, int *bc, double kappa);

void mult_wilson_H_chiral(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, double kappa);

void mult_wilson_D_chiral(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, double kappa);

void mult_wilson_1_dirac(
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double *RESTRICT u, double *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_1_chiral(
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double *RESTRICT u, double *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_dirac(double *RESTRICT v2, double *RESTRICT u, 
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_chiral(double *RESTRICT v2, double *RESTRICT u,
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_dirac(double *RESTRICT v2, double *RESTRICT u, 
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_chiral(double *RESTRICT v2, double *RESTRICT u, 
                    double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                    double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                    double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                    double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_gm5_dirac(double *RESTRICT v2, double *RESTRICT v1,
                           int *Nsize, int Nc);

void mult_wilson_gm5_aypx_dirac(double a, double *RESTRICT v2,
                                double *RESTRICT v1, int *Nsize, int Nc);

void mult_wilson_gm5_chiral(double *RESTRICT v2, double *RESTRICT v1,
                            int *Nsize, int Nc);

void mult_wilson_gm5_aypx_chiral(double a, double *RESTRICT v2,
                                 double *RESTRICT v1, int *Nsize, int Nc);

void mult_wilson_aypx(double a, double *RESTRICT v2, double *RESTRICT v1,
                      int *Nsize, int Nc);

void mult_wilson_axpy(double *RESTRICT v2, double a, double *RESTRICT v1,
                      int *Nsize, int Nc);

void mult_wilson_scal(double *RESTRICT v2, double a,
                      int *Nsize, int Nc);

void mult_wilson_clear(double *RESTRICT v2, int *Nsize, int Nc);

void mult_wilson_xp1(double *RESTRICT buf, double *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xp2(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT buf,  int *Nsize, int *bc, int Nc);

void mult_wilson_xpb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm1(double *RESTRICT buf, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm2(double *RESTRICT v2, double *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xmb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_yp1(double *RESTRICT buf, double *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_yp2(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT buf, int *Nsize, int *bc, int Nc);

void mult_wilson_ypb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym1(double *RESTRICT buf, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym2(double *RESTRICT v2, double *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_ymb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zp1(double *RESTRICT buf, double *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zp2(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT buf, int *Nsize, int *bc, int Nc);

void mult_wilson_zpb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm1(double *RESTRICT buf, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm2(double *RESTRICT v2, double *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zmb(double *RESTRICT v2, double *RESTRICT u,
                     double *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_dirac(double *RESTRICT buf, double *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT buf,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_dirac(double *RESTRICT buf, double *RESTRICT u,
                           double *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_dirac(double *RESTRICT v2, double *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_chiral(double *RESTRICT buf, double *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_chiral(double *RESTRICT buf, double *RESTRICT u,
                            double *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_chiral(double *RESTRICT v2, double *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

//====================================================================

void mult_wilson_H_dirac(float *RESTRICT v2, float *RESTRICT u,
                         float *RESTRICT v1,
                         int *Nsize, int *bc, float kappa);

void mult_wilson_D_dirac(float *RESTRICT v2, float *RESTRICT u,
                         float *RESTRICT v1,
                         int *Nsize, int *bc, float kappa);

void mult_wilson_H_chiral(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, float kappa);

void mult_wilson_D_chiral(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, float kappa);

void mult_wilson_1_dirac(
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float *RESTRICT u, float *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_1_chiral(
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float *RESTRICT u, float *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_dirac(float *RESTRICT v2, float *RESTRICT u, 
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_chiral(float *RESTRICT v2, float *RESTRICT u,
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_dirac(float *RESTRICT v2, float *RESTRICT u, 
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_chiral(float *RESTRICT v2, float *RESTRICT u, 
                    float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                    float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                    float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                    float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_gm5_dirac(float *RESTRICT v2, float *RESTRICT v1,
                           int *Nsize, int Nc);

void mult_wilson_gm5_aypx_dirac(float a, float *RESTRICT v2,
                                float *RESTRICT v1, int *Nsize, int Nc);

void mult_wilson_gm5_chiral(float *RESTRICT v2, float *RESTRICT v1,
                            int *Nsize, int Nc);

void mult_wilson_gm5_aypx_chiral(float a, float *RESTRICT v2,
                                 float *RESTRICT v1, int *Nsize, int Nc);

void mult_wilson_aypx(float a, float *RESTRICT v2, float *RESTRICT v1,
                      int *Nsize, int Nc);

void mult_wilson_axpy(float *RESTRICT v2, float a, float *RESTRICT v1,
                      int *Nsize, int Nc);

void mult_wilson_scal(float *RESTRICT v2, float a,
                      int *Nsize, int Nc);

void mult_wilson_clear(float *RESTRICT v2, int *Nsize, int Nc);

void mult_wilson_xp1(float *RESTRICT buf, float *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xp2(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT buf,  int *Nsize, int *bc, int Nc);

void mult_wilson_xpb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm1(float *RESTRICT buf, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm2(float *RESTRICT v2, float *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xmb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_yp1(float *RESTRICT buf, float *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_yp2(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT buf, int *Nsize, int *bc, int Nc);

void mult_wilson_ypb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym1(float *RESTRICT buf, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym2(float *RESTRICT v2, float *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_ymb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zp1(float *RESTRICT buf, float *RESTRICT v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zp2(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT buf, int *Nsize, int *bc, int Nc);

void mult_wilson_zpb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm1(float *RESTRICT buf, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm2(float *RESTRICT v2, float *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zmb(float *RESTRICT v2, float *RESTRICT u,
                     float *RESTRICT v1, int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_dirac(float *RESTRICT buf, float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT buf,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_dirac(float *RESTRICT buf, float *RESTRICT u,
                           float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_dirac(float *RESTRICT v2, float *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_chiral(float *RESTRICT buf, float *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_chiral(float *RESTRICT buf, float *RESTRICT u,
                            float *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_chiral(float *RESTRICT v2, float *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT v1,
                            int *Nsize, int *bc, int Nc);

}

#endif
//============================================================END=====
