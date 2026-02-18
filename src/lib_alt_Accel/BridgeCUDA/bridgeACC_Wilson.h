/*!
      @file    bridgeACC_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2581 $
*/

#ifndef BRIDGEACC_WILSON_INCLUDED
#define BRIDGEACC_WILSON_INCLUDED

namespace BridgeACC {

// real_t = double

void mult_wilson_H_dirac(double* v2, double* u,
                         double* v1,
                         int *Nsize, int *bc, double kappa);

void mult_wilson_D_dirac(double* v2, double* u,
                         double* v1,
                         int *Nsize, int *bc, double kappa);

void mult_wilson_H_chiral(double* v2, double* u,
                          double* v1,
                          int *Nsize, int *bc, double kappa);

void mult_wilson_D_chiral(double* v2, double* u,
                          double* v1,
                          int *Nsize, int *bc, double kappa);

// QDW Wrappers
void mult_wilson_qdw_D_dirac(double* v2, double* u,
                             double* v1,
                             int *Nsize, int *bc, double kappa);

void mult_wilson_qdw_D_chiral(double* v2, double* u,
                              double* v1,
                              int *Nsize, int *bc, double kappa);


void mult_wilson_1_dirac(
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double* u, double* v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_1_chiral(
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double* u, double* v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_dirac(double* v2, double* u, 
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_chiral(double* v2, double* u,
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_dirac(double* v2, double* u, 
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_chiral(double* v2, double* u, 
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_Meo_dirac(double* v2, double* u,
                           double* v1, double* x1, 
                           const int ieo, const int jeo, const int iflag,
                           int *Nsize, int *bc, double kappa);

void mult_wilson_Meo_chiral(double* v2, double* u,
                            double* v1, double* x1, 
                            const int ieo, const int jeo, const int iflag,
                            int *Nsize, int *bc, double kappa);

void mult_wilson_eo_1_dirac(
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double* u, double* v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_1_chiral(
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double* u, double* v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_2_dirac(double* v2, double* u, 
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_2_chiral(double* v2, double* u,
                    double* buf_xp, double* buf_xm,
                    double* buf_yp, double* buf_ym,
                    double* buf_zp, double* buf_zm,
                    double* buf_tp, double* buf_tm,
                    double kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc);


void mult_wilson_gm5_dirac(double* v2, double* v1,
                           int *Nsize, int Nc);

void mult_wilson_gm5_aypx_dirac(double a, double* v2,
                                double* v1, int *Nsize, int Nc);

void mult_wilson_gm5_chiral(double* v2, double* v1,
                            int *Nsize, int Nc);

void mult_wilson_gm5_aypx_chiral(double a, double* v2,
                                 double* v1, int *Nsize, int Nc);

void mult_wilson_aypx(double a, double* v2, double* v1,
                      int *Nsize, int Nc);

void mult_wilson_axpy(double* v2, double a, double* v1,
                      int *Nsize, int Nc);

void mult_wilson_scal(double* v2, double a,
                      int *Nsize, int Nc);

void mult_wilson_clear(double* v2, int *Nsize, int Nc);

void mult_wilson_xp1(double* buf, double* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xp2(double* v2, double* u,
                     double* buf,  int *Nsize, int *bc, int Nc);

void mult_wilson_xpb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm1(double* buf, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm2(double* v2, double* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xmb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_yp1(double* buf, double* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_yp2(double* v2, double* u,
                     double* buf, int *Nsize, int *bc, int Nc);

void mult_wilson_ypb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym1(double* buf, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym2(double* v2, double* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_ymb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zp1(double* buf, double* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zp2(double* v2, double* u,
                     double* buf, int *Nsize, int *bc, int Nc);

void mult_wilson_zpb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm1(double* buf, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm2(double* v2, double* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zmb(double* v2, double* u,
                     double* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_dirac(double* buf, double* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_dirac(double* v2, double* u,
                           double* buf,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_dirac(double* v2, double* u,
                           double* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_dirac(double* buf, double* u,
                           double* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_dirac(double* v2, double* buf, 
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_dirac(double* v2, double* u,
                           double* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_chiral(double* buf, double* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_chiral(double* v2, double* u,
                            double* buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_chiral(double* v2, double* u,
                            double* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_chiral(double* buf, double* u,
                            double* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_chiral(double* v2, double* buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_chiral(double* v2, double* u,
                            double* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_xp1_eo(double* buf, double* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xp2_eo(double* v2, double* u,
                        double* buf,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xpb_eo(double* v2, double* u,
                        double* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xm1_eo(double* buf, double* u,
                        double* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xm2_eo(double* v2, double* buf, 
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xmb_eo(double* v2, double* u,
                        double* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void set_block_config(double* u, int* Nsize,
                      int* block_size);

void mult_wilson_dd_dirac(double* v2, double* u, double* v1, double kappa,
                          int *bc, int *Nsize, int *block_size, int ieo);


// real_t = float

void mult_wilson_H_dirac(float* v2, float* u,
                         float* v1,
                         int *Nsize, int *bc, float kappa);

void mult_wilson_D_dirac(float* v2, float* u,
                         float* v1,
                         int *Nsize, int *bc, float kappa);

void mult_wilson_H_chiral(float* v2, float* u,
                          float* v1,
                          int *Nsize, int *bc, float kappa);

void mult_wilson_D_chiral(float* v2, float* u,
                          float* v1,
                          int *Nsize, int *bc, float kappa);

// QDW Wrappers (float stubs)
void mult_wilson_qdw_D_dirac(float* v2, float* u,
                             float* v1,
                             int *Nsize, int *bc, float kappa);

void mult_wilson_qdw_D_chiral(float* v2, float* u,
                              float* v1,
                              int *Nsize, int *bc, float kappa);

void mult_wilson_1_dirac(
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float* u, float* v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_1_chiral(
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float* u, float* v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_dirac(float* v2, float* u, 
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_2_chiral(float* v2, float* u,
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_dirac(float* v2, float* u, 
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_H2_chiral(float* v2, float* u, 
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_Meo_dirac(float* v2, float* u,
                           float* v1, float* x1, 
                           const int ieo, const int jeo, const int iflag,
                           int *Nsize, int *bc, float kappa);

void mult_wilson_Meo_chiral(float* v2, float* u,
                            float* v1, float* x1, 
                            const int ieo, const int jeo, const int iflag,
                            int *Nsize, int *bc, float kappa);

void mult_wilson_eo_1_dirac(
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float* u, float* v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_1_chiral(
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float* u, float* v1,
                    const int ieo, const int jeo,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_2_dirac(float* v2, float* u, 
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_eo_2_chiral(float* v2, float* u,
                    float* buf_xp, float* buf_xm,
                    float* buf_yp, float* buf_ym,
                    float* buf_zp, float* buf_zm,
                    float* buf_tp, float* buf_tm,
                    float kappa, const int ieo, const int jeo,
                    const int iflag,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_gm5_dirac(float* v2, float* v1,
                           int *Nsize, int Nc);

void mult_wilson_gm5_aypx_dirac(float a, float* v2,
                                float* v1, int *Nsize, int Nc);

void mult_wilson_gm5_chiral(float* v2, float* v1,
                            int *Nsize, int Nc);

void mult_wilson_gm5_aypx_chiral(float a, float* v2,
                                 float* v1, int *Nsize, int Nc);

void mult_wilson_aypx(float a, float* v2, float* v1,
                      int *Nsize, int Nc);

void mult_wilson_axpy(float* v2, float a, float* v1,
                      int *Nsize, int Nc);

void mult_wilson_scal(float* v2, float a,
                      int *Nsize, int Nc);

void mult_wilson_clear(float* v2, int *Nsize, int Nc);

void mult_wilson_xp1(float* buf, float* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xp2(float* v2, float* u,
                     float* buf,  int *Nsize, int *bc, int Nc);

void mult_wilson_xpb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm1(float* buf, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_xm2(float* v2, float* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_xmb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_yp1(float* buf, float* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_yp2(float* v2, float* u,
                     float* buf, int *Nsize, int *bc, int Nc);

void mult_wilson_ypb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym1(float* buf, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_ym2(float* v2, float* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_ymb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zp1(float* buf, float* v1,
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zp2(float* v2, float* u,
                     float* buf, int *Nsize, int *bc, int Nc);

void mult_wilson_zpb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm1(float* buf, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_zm2(float* v2, float* buf, 
                     int *Nsize, int *bc, int Nc);

void mult_wilson_zmb(float* v2, float* u,
                     float* v1, int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_dirac(float* buf, float* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_dirac(float* v2, float* u,
                           float* buf,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_dirac(float* v2, float* u,
                           float* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_dirac(float* buf, float* u,
                           float* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_dirac(float* v2, float* buf, 
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_dirac(float* v2, float* u,
                           float* v1,
                           int *Nsize, int *bc, int Nc);

void mult_wilson_tp1_chiral(float* buf, float* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tp2_chiral(float* v2, float* u,
                            float* buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tpb_chiral(float* v2, float* u,
                            float* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm1_chiral(float* buf, float* u,
                            float* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tm2_chiral(float* v2, float* buf, 
                            int *Nsize, int *bc, int Nc);

void mult_wilson_tmb_chiral(float* v2, float* u,
                            float* v1,
                            int *Nsize, int *bc, int Nc);

void mult_wilson_xp1_eo(float* buf, float* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xp2_eo(float* v2, float* u,
                        float* buf,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xpb_eo(float* v2, float* u,
                        float* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xm1_eo(float* buf, float* u,
                        float* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xm2_eo(float* v2, float* buf, 
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_xmb_eo(float* v2, float* u,
                        float* v1,
                        int *Nsize, int *bc, int ieo, int Nc);

void mult_wilson_dd_dirac(float* v2, float* u, float* v1, float kappa,
                          int *bc, int *Nsize, int *block_size, int ieo);

void set_block_config(float* u, int* Nsize,
                      int* block_size);


}

#endif
