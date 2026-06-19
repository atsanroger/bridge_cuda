/*!
      @file    bridgeACC_Domainwall.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2620 $
*/

#ifndef BRIDGEACC_DOMAINWALL_INCLUDED
#define BRIDGEACC_DOMAINWALL_INCLUDED

#include "bridge_defs.h"

namespace BridgeACC {
  // real_t = double

  void mult_domainwall_5din_5dir_dirac(
        double *RESTRICT vp, double *RESTRICT yp, double *RESTRICT wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  void mult_domainwall_5din_5dirdag_dirac(
        double *RESTRICT vp, double *RESTRICT yp, double *RESTRICT wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  void mult_domainwall_5din_mult_gm5_dirac(
        double *RESTRICT vp, double *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_mult_gm5R_dirac(
        double *RESTRICT vp, double *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_mult_R(
        double *RESTRICT vp, double *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_clear(double *RESTRICT vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
         double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
         int Ns, int *bc, int *Nsize, int *do_comm, int flag);

  void mult_domainwall_5din_hop1_dirac(
        double *RESTRICT buf1_xp, double *RESTRICT buf1_xm,
        double *RESTRICT buf1_yp, double *RESTRICT buf1_ym,
        double *RESTRICT buf1_zp, double *RESTRICT buf1_zm,
        double *RESTRICT buf1_tp, double *RESTRICT buf1_tm,
        double *RESTRICT up, double *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
        double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
        double *RESTRICT buf2_xp, double *RESTRICT buf2_xm,
        double *RESTRICT buf2_yp, double *RESTRICT buf2_ym,
        double *RESTRICT buf2_zp, double *RESTRICT buf2_zm,
        double *RESTRICT buf2_tp, double *RESTRICT buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm);

  void mult_domainwall_5din_ee_5dir_dirac(
        double *RESTRICT vp, double *RESTRICT wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  void mult_domainwall_5din_eo_5dir_dirac(
        double *RESTRICT yp, double *RESTRICT wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  void mult_domainwall_5din_ee_5dirdag_dirac(
        double *RESTRICT vp, double *RESTRICT wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
        double *RESTRICT vp, double *RESTRICT yp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

  //void mult_domainwall_5din_eo_hopb_dirac(
  //      double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
  //      int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo);

  void mult_domainwall_5din_eo_hopb_dirac_4d(
        double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
	int jgm5);

  void mult_domainwall_5din_eo_hopb_dirac_5d(
        double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
	int jgm5);

  void mult_domainwall_5din_eo_hop1_dirac(
        double *RESTRICT buf1_xp, double *RESTRICT buf1_xm,
        double *RESTRICT buf1_yp, double *RESTRICT buf1_ym,
        double *RESTRICT buf1_zp, double *RESTRICT buf1_zm,
        double *RESTRICT buf1_tp, double *RESTRICT buf1_tm,
        double *RESTRICT up, double *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
        int jgm5);

  void mult_domainwall_5din_eo_hop2_dirac(
        double *RESTRICT vp, double *RESTRICT up, double *RESTRICT wp,
        double *RESTRICT buf2_xp, double *RESTRICT buf2_xm,
        double *RESTRICT buf2_yp, double *RESTRICT buf2_ym,
        double *RESTRICT buf2_zp, double *RESTRICT buf2_zm,
        double *RESTRICT buf2_tp, double *RESTRICT buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo);

  void mult_domainwall_5din_L_inv_dirac(
         double *RESTRICT vp, double *RESTRICT wp,
         int Ns, int *Nsize, double *e, double *dpinv, double *dm);

  void mult_domainwall_5din_U_inv_dirac(
         double *RESTRICT vp, double *RESTRICT wp,
         int Ns, int *Nsize, double *f, double *dpinv, double *dm);

  void mult_domainwall_5din_Ldag_inv_dirac(
         double *RESTRICT vp, double *RESTRICT wp,
         int Ns, int *Nsize, double *e, double *dpinv, double *dm);

  void mult_domainwall_5din_Udag_inv_dirac(
         double *RESTRICT vp, double *RESTRICT wp,
         int Ns, int *Nsize, double *f, double *dpinv, double *dm);

  void mult_domainwall_5din_dd_5dir_dirac(
         double *RESTRICT vp, double *RESTRICT yp, double *RESTRICT wp,
         double mq, double M0, int Ns, double *b, double *c,
         int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_5dirdag_dirac(
         double *RESTRICT vp, double *RESTRICT yp, double *RESTRICT wp,
         double mq, double M0, int Ns, double *b, double *c,
         int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_hopb_dirac(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *block_size, int ieo, int flag);

  void mult_domainwall_5din_xpb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_xmb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_ypb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_ymb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_zpb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_zmb(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_tpb_dirac(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_tmb_dirac(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_xp1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xm1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_yp1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ym1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zp1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zm1(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tp1_dirac(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tm1_dirac(
                            double *RESTRICT buf, double *RESTRICT up,
                            double *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xp2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xm2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_yp2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ym2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zp2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zm2(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tp2_dirac(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tm2_dirac(
                            double *RESTRICT vp, double *RESTRICT up,
                            double *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ee_inv_dirac_4d(
                        double *RESTRICT vp, double *RESTRICT wp,
			int jd,
                        int Ns, double *RESTRICT mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_inv_dirac_5d(
                        double *RESTRICT vp, double *RESTRICT wp,
			int jd,
                        int Ns, double *RESTRICT mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_LUinv_dirac(
                        double *RESTRICT vp, double *RESTRICT wp,
                        int Ns, int *Nsize, double *e, double *f,
                        double *dpinv, double *dm);

  void mult_domainwall_5din_ee_LUdaginv_dirac(
                        double *RESTRICT vp, double *RESTRICT wp,
                        int Ns, int *Nsize, double *e, double *f,
                        double *dpinv, double *dm);

  // real_t = float

  void mult_domainwall_5din_5dir_dirac(
         float *RESTRICT vp, float *RESTRICT yp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  void mult_domainwall_5din_5dirdag_dirac(
         float *RESTRICT vp, float *RESTRICT yp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  void mult_domainwall_5din_mult_gm5_dirac(
         float *RESTRICT vp, float *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_mult_gm5R_dirac(
         float *RESTRICT vp, float *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_mult_R(
         float *RESTRICT vp, float *RESTRICT wp, int Ns, int *Nsize);

  void mult_domainwall_5din_clear(float *RESTRICT vp, int Ns, int *Nsize);

  void mult_domainwall_5din_hopb_dirac(
         float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
         int Ns, int *bc, int *Nsize, int *do_comm, int flag);

  void mult_domainwall_5din_hop1_dirac(
         float *RESTRICT buf1_xp, float *RESTRICT buf1_xm,
         float *RESTRICT buf1_yp, float *RESTRICT buf1_ym,
         float *RESTRICT buf1_zp, float *RESTRICT buf1_zm,
         float *RESTRICT buf1_tp, float *RESTRICT buf1_tm,
         float *RESTRICT up, float *RESTRICT wp,
         int Ns, int *bc, int *Nsize, int *do_comm);

  void mult_domainwall_5din_hop2_dirac(
         float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
         float *RESTRICT buf2_xp, float *RESTRICT buf2_xm,
         float *RESTRICT buf2_yp, float *RESTRICT buf2_ym,
         float *RESTRICT buf2_zp, float *RESTRICT buf2_zm,
         float *RESTRICT buf2_tp, float *RESTRICT buf2_tm,
         int Ns, int *bc, int *Nsize, int *do_comm);

  void mult_domainwall_5din_ee_5dir_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  void mult_domainwall_5din_eo_5dir_dirac(
         float *RESTRICT yp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  void mult_domainwall_5din_ee_5dirdag_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  void mult_domainwall_5din_eo_5dirdag_dirac(
         float *RESTRICT vp, float *RESTRICT yp,
         float mq, float M0, int Ns, float *b, float *c, int *Nsize);

  //void mult_domainwall_5din_eo_hopb_dirac(
  //      float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
  //      int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo);

  void mult_domainwall_5din_eo_hopb_dirac_4d(
        float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
        int jgm5);

  void mult_domainwall_5din_eo_hopb_dirac_5d(
        float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
        int jgm5);

  void mult_domainwall_5din_eo_hop1_dirac(
         float *RESTRICT buf1_xp, float *RESTRICT buf1_xm,
         float *RESTRICT buf1_yp, float *RESTRICT buf1_ym,
         float *RESTRICT buf1_zp, float *RESTRICT buf1_zm,
         float *RESTRICT buf1_tp, float *RESTRICT buf1_tm,
         float *RESTRICT up, float *RESTRICT wp,
         int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
         int jgm5);

  void mult_domainwall_5din_eo_hop2_dirac(
         float *RESTRICT vp, float *RESTRICT up, float *RESTRICT wp,
         float *RESTRICT buf2_xp, float *RESTRICT buf2_xm,
         float *RESTRICT buf2_yp, float *RESTRICT buf2_ym,
         float *RESTRICT buf2_zp, float *RESTRICT buf2_zm,
         float *RESTRICT buf2_tp, float *RESTRICT buf2_tm,
         int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo);

  void mult_domainwall_5din_L_inv_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         int Ns, int *Nsize, float *e, float *dpinv, float *dm);

  void mult_domainwall_5din_U_inv_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         int Ns, int *Nsize, float *f, float *dpinv, float *dm);

  void mult_domainwall_5din_Ldag_inv_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         int Ns, int *Nsize, float *e, float *dpinv, float *dm);

  void mult_domainwall_5din_Udag_inv_dirac(
         float *RESTRICT vp, float *RESTRICT wp,
         int Ns, int *Nsize, float *f, float *dpinv, float *dm);

  void mult_domainwall_5din_dd_5dir_dirac(
         float *RESTRICT vp, float *RESTRICT yp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c,
         int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_5dirdag_dirac(
         float *RESTRICT vp, float *RESTRICT yp, float *RESTRICT wp,
         float mq, float M0, int Ns, float *b, float *c,
         int *Nsize, int *block_size, int ieo);

  void mult_domainwall_5din_dd_hopb_dirac(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *block_size, int ieo, int flag);

  void mult_domainwall_5din_xpb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_xmb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_ypb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_ymb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_zpb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_zmb(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_tpb_dirac(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_tmb_dirac(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

  void mult_domainwall_5din_xp1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xm1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_yp1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ym1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zp1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zm1(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tp1_dirac(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tm1_dirac(
                            float *RESTRICT buf, float *RESTRICT up,
                            float *RESTRICT wp,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xp2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_xm2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_yp2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ym2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zp2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_zm2(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tp2_dirac(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_tm2_dirac(
                            float *RESTRICT vp, float *RESTRICT up,
                            float *RESTRICT buf,
                            int Ns, int *bc, int *Nsize);

  void mult_domainwall_5din_ee_inv_dirac_4d(
                        float *RESTRICT vp, float *RESTRICT wp,
			int jd,
                        int Ns, float *RESTRICT mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_inv_dirac_5d(
                        float *RESTRICT vp, float *RESTRICT wp,
			int jd,
                        int Ns, float *RESTRICT mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_LUinv_dirac(
                        float *RESTRICT vp, float *RESTRICT wp,
                        int Ns, int *Nsize, float *e, float *f,
                        float *dpinv, float *dm);

  void mult_domainwall_5din_ee_LUdaginv_dirac(
                        float *RESTRICT vp, float *RESTRICT wp,
                        int Ns, int *Nsize, float *e, float *f,
                        float *dpinv, float *dm);

}
#endif
