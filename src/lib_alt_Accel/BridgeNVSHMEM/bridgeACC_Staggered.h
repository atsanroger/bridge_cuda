/*!
      @file    bridgeACC_Staggered.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2534 $
*/

#ifndef BRIDGEACC_STAGGERED_INCLUDED
#define BRIDGEACC_STAGGERED_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_staggered_D(double *restrict v2, double *restrict u_up,
                        double *restrict u_dn, double *restrict v1,
                        double mq, int *bc, int *Nsize, int jdag);

  void mult_staggered_1(double *restrict buf_xp, double *restrict buf_xm,
                        double *restrict buf_yp, double *restrict buf_ym,
                        double *restrict buf_zp, double *restrict buf_zm,
                        double *restrict buf_tp, double *restrict buf_tm,
                        double *restrict u_dn, double *restrict v1,
                        int *bc, int *Nsize, int *do_comm);

  void mult_staggered_2(double *restrict v2, double *restrict u_up,
                        double *restrict buf_xp, double *restrict buf_xm,
                        double *restrict buf_yp, double *restrict buf_ym,
                        double *restrict buf_zp, double *restrict buf_zm,
                        double *restrict buf_tp, double *restrict buf_tm,
                        double mq, int *bc,
                        int *Nsize, int *do_comm, int jdag);

  void mult_staggered_Meo(double *restrict v2,
			  double *restrict u_up, double *restrict u_dn,
			  double *restrict v1, double *restrict x1, 
			  double mq2, int *bc, int *Nsize,
			  const int ieo, const int jeo, const int iflag);

  void mult_staggered_1eo(double *restrict buf_xp, double *restrict buf_xm,
                          double *restrict buf_yp, double *restrict buf_ym,
                          double *restrict buf_zp, double *restrict buf_zm,
                          double *restrict buf_tp, double *restrict buf_tm,
                          double *restrict u, double *restrict v1,
                          int *bc, int *Nsize, int *do_comm,
                          const int ieo, const int jeo);

  void mult_staggered_2eo(double *restrict v2, double *restrict u, 
                          double *restrict buf_xp, double *restrict buf_xm,
                          double *restrict buf_yp, double *restrict buf_ym,
                          double *restrict buf_zp, double *restrict buf_zm,
                          double *restrict buf_tp, double *restrict buf_tm,
                          double mq, int *bc,  int *Nsize, int *do_comm,
                          const int ieo, const int jeo, const int iflag);

  void mult_staggered_phase_dev(double *restrict u, double *restrict ph,
                                int *Nsize, int Nc);

  void mult_staggered_gm5(double *restrict v2, double *restrict v1,
                          int *Nsize, int Nc);

  void mult_staggered_aypx(double a, double *restrict v2, double *restrict v1,
                           int *Nsize, int Nc);

  void mult_staggered_axpy(double *restrict v2, double a, double *restrict v1,
                           int *Nsize, int Nc);

  void mult_staggered_scal(double *restrict v, double a,
                           int *Nsize, int Nc);

  void mult_staggered_clear(double *restrict v2, int *Nsize, int Nc);

  void mult_staggered_xp1(double *restrict buf, double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp2(double *restrict v2, double *restrict u,
                          double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xpb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm1(double *restrict buf, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm2(double *restrict v2, double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xmb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp1(double *restrict buf, double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp2(double *restrict v2, double *restrict u,
                          double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ypb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym1(double *restrict buf, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym2(double *restrict v2, double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_ymb(double *restrict v2, double *restrict u,
                           double *restrict v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zp1(double *restrict buf, double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zp2(double *restrict v2, double *restrict u,
                          double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zpb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm1(double *restrict buf, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm2(double *restrict v2, double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zmb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp1(double *restrict buf, double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp2(double *restrict v2, double *restrict u,
                          double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tpb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm1(double *restrict buf, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm2(double *restrict v2, double *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tmb(double *restrict v2, double *restrict u,
                          double *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp1_eo(double *restrict buf, double *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xp2_eo(double *restrict v2, double *restrict u,
                             double *restrict buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xpb_eo(double *restrict v2, double *restrict u,
                             double *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm1_eo(double *restrict buf, double *restrict u,
                             double *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm2_eo(double *restrict v2, double *restrict buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xmb_eo(double *restrict v2, double *restrict u,
                             double *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  // real_t = float

  void mult_staggered_D(float *restrict v2, float *restrict u_up,
                        float *restrict u_dn, float *restrict v1,
                        float mq, int *bc, int *Nsize, int jdag);

  void mult_staggered_1(float *restrict buf_xp, float *restrict buf_xm,
                        float *restrict buf_yp, float *restrict buf_ym,
                        float *restrict buf_zp, float *restrict buf_zm,
                        float *restrict buf_tp, float *restrict buf_tm,
                        float *restrict u_dn, float *restrict v1,
                        int *bc, int *Nsize, int *do_comm);

  void mult_staggered_2(float *restrict v2, float *restrict u_up,
                        float *restrict buf_xp, float *restrict buf_xm,
                        float *restrict buf_yp, float *restrict buf_ym,
                        float *restrict buf_zp, float *restrict buf_zm,
                        float *restrict buf_tp, float *restrict buf_tm,
                        float mq, int *bc,
                        int *Nsize, int *do_comm, int jdag);

  void mult_staggered_Meo(float *restrict v2,
			  float *restrict u_up, float *restrict u_dn,
			  float *restrict v1, float *restrict x1, 
			  float mq2, int *bc, int *Nsize,
			  const int ieo, const int jeo, const int iflag);

  void mult_staggered_1eo(float *restrict buf_xp, float *restrict buf_xm,
                          float *restrict buf_yp, float *restrict buf_ym,
                          float *restrict buf_zp, float *restrict buf_zm,
                          float *restrict buf_tp, float *restrict buf_tm,
                          float *restrict u, float *restrict v1,
                          int *bc, int *Nsize, int *do_comm,
                          const int ieo, const int jeo);

  void mult_staggered_2eo(float *restrict v2, float *restrict u, 
                          float *restrict buf_xp, float *restrict buf_xm,
                          float *restrict buf_yp, float *restrict buf_ym,
                          float *restrict buf_zp, float *restrict buf_zm,
                          float *restrict buf_tp, float *restrict buf_tm,
                          float mq, int *bc,  int *Nsize, int *do_comm,
                          const int ieo, const int jeo, const int iflag);

  void mult_staggered_phase_dev(float *restrict u, float *restrict ph,
                                int *Nsize, int Nc);

  void mult_staggered_gm5(float *restrict v2, float *restrict v1,
                         int *Nsize, int Nc);

  void mult_staggered_aypx(float a, float *restrict v2, float *restrict v1,
                           int *Nsize, int Nc);

  void mult_staggered_axpy(float *restrict v2, float a, float *restrict v1,
                           int *Nsize, int Nc);

  void mult_staggered_scal(float *restrict v, float a,
                           int *Nsize, int Nc);

  void mult_staggered_clear(float *restrict v2, int *Nsize, int Nc);

  void mult_staggered_xp1(float *restrict buf, float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp2(float *restrict v2, float *restrict u,
                          float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xpb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm1(float *restrict buf, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm2(float *restrict v2, float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xmb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp1(float *restrict buf, float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp2(float *restrict v2, float *restrict u,
                          float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ypb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym1(float *restrict buf, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym2(float *restrict v2, float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_ymb(float *restrict v2, float *restrict u,
                           float *restrict v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zp1(float *restrict buf, float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zp2(float *restrict v2, float *restrict u,
                          float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_zpb(float *restrict v2, float *restrict u,
                           float *restrict v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zm1(float *restrict buf, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm2(float *restrict v2, float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zmb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp1(float *restrict buf, float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp2(float *restrict v2, float *restrict u,
                          float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tpb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm1(float *restrict buf, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm2(float *restrict v2, float *restrict buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tmb(float *restrict v2, float *restrict u,
                          float *restrict v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp1_eo(float *restrict buf, float *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xp2_eo(float *restrict v2, float *restrict u,
                             float *restrict buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xpb_eo(float *restrict v2, float *restrict u,
                             float *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm1_eo(float *restrict buf, float *restrict u,
                             float *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm2_eo(float *restrict v2, float *restrict buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xmb_eo(float *restrict v2, float *restrict u,
                             float *restrict v1,
                             int *Nsize, int *bc, int ieo, int Nc);

}

#endif
