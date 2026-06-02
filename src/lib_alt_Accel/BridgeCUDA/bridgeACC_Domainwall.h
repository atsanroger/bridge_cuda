/*!
      @file    bridgeACC_Domainwall_5din.h
      @brief
      @author  Wei Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-01-01 13:51:53 #$
      @version $LastChangedRevision: 2621 $
*/

#ifndef BRIDGEACC_DOMAINWALL_5DIN_ACC_INCLUDED
#define BRIDGEACC_DOMAINWALL_5DIN_ACC_INCLUDED

namespace BridgeACC {

// real_t = double

void initDomainwallConstantMemory(
    double* e, double* f, double* dpinv,
    double* dm, double* b, double* c, int Ns);
    
void initDomainwallConstantMemory(
    float* e, float* f, float* dpinv,
    float* dm, float* b, float* c, int Ns);

// Push float-float pairs from double coefficients into const_*_float (hi)
// and const_*_ff_lo (lo). Used by extended_precision=true on a float base
// so kernels can read coefficients at ~48-bit precision and do FP32-only
// dw_add/dw_mul on them.
void initDomainwallConstantMemoryFF(
    double* e, double* f, double* dpinv,
    double* dm, double* b, double* c, int Ns);

// Triple-word (QTW) variant: pushes (hi, mid, lo) FP32 decomposition of each
// double coefficient. const_*_float gets the hi word; const_*_tw_mid the mid;
// const_*_tw_lo the lo. Together they encode each coefficient at ~72-bit
// precision using only FP32 storage / fma arithmetic.
void initDomainwallConstantMemoryTW(
    double* e, double* f, double* dpinv,
    double* dm, double* b, double* c, int Ns);

void mult_domainwall_5din_5dir_dirac(
        double *vp, double *yp, double *wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

void mult_domainwall_5din_5dirdag_dirac(
        double *vp, double *yp, double *wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

void mult_domainwall_5din_mult_gm5_dirac(
        double *vp, double *wp, int Ns, int *Nsize);

void mult_domainwall_5din_mult_gm5R_dirac(
        double *vp, double *wp, int Ns, int *Nsize );

void mult_domainwall_5din_mult_R(
        double *vp, double *wp, int Ns, int *Nsize);

void mult_domainwall_5din_clear(
        double *vp, int Ns, int *Nsize);

void mult_domainwall_5din_hopb_dirac(
        double *vp, double *up, double *wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int flag);

void mult_domainwall_5din_hop1_dirac(
        double * buf1_xp, double * buf1_xm,
        double * buf1_yp, double * buf1_ym,
        double * buf1_zp, double * buf1_zm,
        double * buf1_tp, double * buf1_tm,
        double * up, double * wp,
        int Ns, int *bc, int *Nsize, int *do_comm);

void mult_domainwall_5din_hop2_dirac(
        double * vp, double * up, double * wp,
        double * buf2_xp, double * buf2_xm,
        double * buf2_yp, double * buf2_ym,
        double * buf2_zp, double * buf2_zm,
        double * buf2_tp, double * buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm);

void mult_domainwall_5din_ee_5dir_dirac(
    double *vp, double *wp, double mq, double M0, int Ns, 
    double *b, double *c, double alpha, int *Nsize);

void mult_domainwall_5din_eo_5dir_dirac(
    double *yp, double *wp, double mq, double M0, int Ns, 
    double *b, double *c, double alpha, int *Nsize);

void mult_domainwall_5din_ee_5dirdag_dirac(
    double *vp, double *yp, double mq, double M0, int Ns, 
    double *b, double *c, double alpha, int *Nsize);

void mult_domainwall_5din_eo_5dirdag_dirac(
    double *vp, double *yp, double mq, double M0, int Ns, 
    double *b, double *c, double alpha, int *Nsize);

void mult_domainwall_5din_eo_hopb_dirac_4d(
    double *vp, double *up, double *wp, int Ns, int *bc, 
    int *Nsize, int *do_comm, int ieo, int jeo,
    int jgm5);

void mult_domainwall_5din_eo_hopb_dirac_5d(
    double *vp, double *up, double *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo,
    int jgm5, bool recon);

void mult_domainwall_5din_eo_hop1_dirac(
        double * buf1_xp, double * buf1_xm,
        double * buf1_yp, double * buf1_ym,
        double * buf1_zp, double * buf1_zm,
        double * buf1_tp, double * buf1_tm,
        double * up     , double * wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5);

void mult_domainwall_5din_eo_hop2_dirac(
        double * vp     , double * up, double * wp,
        double * buf2_xp, double * buf2_xm,
        double * buf2_yp, double * buf2_ym,
        double * buf2_zp, double * buf2_zm,
        double * buf2_tp, double * buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo);

void mult_domainwall_5din_L_inv_dirac(
        double *vp, double *wp, int Ns, int *Nsize, double *e, double *dpinv, double *dm);

void mult_domainwall_5din_U_inv_dirac(
        double *vp, double *wp, int Ns, int *Nsize, double *f, double *dpinv, double *dm);

void mult_domainwall_5din_Ldag_inv_dirac(
        double *vp, double *wp, int Ns, int *Nsize, double *e, double *dpinv, double *dm);

void mult_domainwall_5din_Udag_inv_dirac(
        double *vp, double *wp, int Ns, int *Nsize, double *f, double *dpinv, double *dm);

void mult_domainwall_5din_dd_5dir_dirac(
           double * vp, double * yp, double * wp,
           double mq, double M0, int Ns, double *b, double *c,
           int *Nsize, int *block_size, int ieo);

void mult_domainwall_5din_dd_5dirdag_dirac(
           double * vp, double * yp, double * wp,
           double mq, double M0, int Ns, double *b, double *c,
           int *Nsize, int *block_size, int ieo);

void mult_domainwall_5din_dd_hopb_dirac(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *block_size, int ieo, int flag);

void mult_domainwall_5din_xpb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_xmb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_ypb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_ymb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_zpb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_zmb(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_tpb_dirac(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_tmb_dirac(
                            double * vp, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_xp1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xm1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_yp1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_ym1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zp1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zm1(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tp1_dirac(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tm1_dirac(
                            double * buf, double * up,
                            double * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xp2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xm2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_yp2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_ym2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zp2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zm2(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tp2_dirac(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tm2_dirac(
                            double * vp, double * up,
                            double * buf,
                            int Ns, int *bc, int *Nsize);

// QDW Domain Wall kernels (double precision only)
void mult_domainwall_5din_5dir_dirac_qdw(
        double *vp, double *yp, double *wp,
        double mq, double M0, int Ns, double *b, double *c, int *Nsize);

void mult_domainwall_5din_hopb_qdw_dirac(
        double *vp, double *up, double *wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int flag);

// QDW EO Domain Wall kernels (double overloads)
void mult_domainwall_5din_ee_5dir_dirac_qdw(
        double *vp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);

void mult_domainwall_5din_eo_5dir_dirac_qdw(
        double *yp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);

void mult_domainwall_5din_ee_5dirdag_dirac_qdw(
        double *vp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);

void mult_domainwall_5din_eo_5dirdag_dirac_qdw(
        double *vp, double *yp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);

void mult_domainwall_5din_eo_hopb_qdw_dirac_5d(
        double *vp, double *up, double *up_lo, double *wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5, bool recon);

// QDW EO Domain Wall kernels — float overloads (vector is float-float;
// scalar args stay double so extended_precision can promote coefficients).
void mult_domainwall_5din_ee_5dir_dirac_qdw(
        float *vp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dir_dirac_qdw(
        float *yp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_ee_5dirdag_dirac_qdw(
        float *vp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dirdag_dirac_qdw(
        float *vp, float *yp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_hopb_qdw_dirac_5d(
        float *vp, float *up, float *up_lo, float *wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5, bool recon);

// QDW LU inverse (alpha promotable to double via extended_precision)
void mult_domainwall_5din_ee_LUinv_dirac_qdw(
        double* vp, double* wp, int Ns, int *Nsize, double alpha, bool ext);
void mult_domainwall_5din_ee_LUdaginv_dirac_qdw(
        double* vp, double* wp, int Ns, int *Nsize, double alpha, bool ext);

void mult_domainwall_5din_ee_LUinv_dirac_qdw(
        float* vp, float* wp, int Ns, int *Nsize, double alpha, bool ext);
void mult_domainwall_5din_ee_LUdaginv_dirac_qdw(
        float* vp, float* wp, int Ns, int *Nsize, double alpha, bool ext);

// QTW (triple-word, 6 reals/cplx) EO Domain Wall kernels.
// `ext` is reserved for future TW gauge-link toggling; today the FP-gauge
// path is used (up_mid==up_lo==nullptr in hopb).
void mult_domainwall_5din_ee_5dir_dirac_qtw(
        double *vp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dir_dirac_qtw(
        double *yp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_ee_5dirdag_dirac_qtw(
        double *vp, double *wp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dirdag_dirac_qtw(
        double *vp, double *yp, double mq, double M0, int Ns,
        double *b, double *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_hopb_qtw_dirac_5d(
        double *vp, double *up, double *up_mid, double *up_lo, double *wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5, bool recon);
void mult_domainwall_5din_ee_LUinv_dirac_qtw(
        double *vp, double *wp, int Ns, int *Nsize, double alpha, bool ext);
void mult_domainwall_5din_ee_LUdaginv_dirac_qtw(
        double *vp, double *wp, int Ns, int *Nsize, double alpha, bool ext);

void mult_domainwall_5din_ee_5dir_dirac_qtw(
        float *vp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dir_dirac_qtw(
        float *yp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_ee_5dirdag_dirac_qtw(
        float *vp, float *wp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_5dirdag_dirac_qtw(
        float *vp, float *yp, double mq, double M0, int Ns,
        float *b, float *c, double alpha, int *Nsize, bool ext);
void mult_domainwall_5din_eo_hopb_qtw_dirac_5d(
        float *vp, float *up, float *up_mid, float *up_lo, float *wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5, bool recon);
void mult_domainwall_5din_ee_LUinv_dirac_qtw(
        float *vp, float *wp, int Ns, int *Nsize, double alpha, bool ext);
void mult_domainwall_5din_ee_LUdaginv_dirac_qtw(
        float *vp, float *wp, int Ns, int *Nsize, double alpha, bool ext);

// anchor1
  void set_block_config(double* u, int* Nsize,
                        int* block_size);

  void mult_domainwall_5din_ee_inv_dirac_4d(
                        double * vp, double * wp, int jd,
                        int Ns, double * mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_inv_dirac_5d(
                        double * vp, double * wp, int jd,
                        int Ns, double * mat_inv, int *Nsize);

  void mult_domainwall_5din_ee_LUdaginv_dirac(
                        double* vp, double* wp,
                        int Ns, int *Nsize,
                        double* e, double* f,
                        double* dpinv, double* dm, double alpha);

  void mult_domainwall_5din_ee_LUinv_dirac(
                        double* vp, double* wp,
                        int Ns, int *Nsize,
                        double* e, double* f,
                        double* dpinv, double* dm, double alpha);

//-----------------float

// real_t = float

void mult_domainwall_5din_5dir_dirac(
        float *vp, float *yp, float *wp,
        float mq, float M0, int Ns, float *b, float *c, int *Nsize);

void mult_domainwall_5din_5dirdag_dirac(
        float *vp, float *yp, float *wp,
        float mq, float M0, int Ns, float *b, float *c, int *Nsize);

void mult_domainwall_5din_mult_gm5_dirac(
        float *vp, float *wp, int Ns, int *Nsize);

void mult_domainwall_5din_mult_gm5R_dirac(
        float *vp, float *wp, int Ns, int *Nsize );

void mult_domainwall_5din_mult_R(
        float *vp, float *wp, int Ns, int *Nsize);

void mult_domainwall_5din_clear(
        float *vp, int Ns, int *Nsize);

void mult_domainwall_5din_hopb_dirac(
        float *vp, float *up, float *wp,
        int Ns, int *bc, int *Nsize, int *do_comm, int flag);

void mult_domainwall_5din_hop1_dirac(
         float * buf1_xp, float * buf1_xm,
         float * buf1_yp, float * buf1_ym,
         float * buf1_zp, float * buf1_zm,
         float * buf1_tp, float * buf1_tm,
         float * up, float * wp,
         int Ns, int *bc, int *Nsize, int *do_comm);

void mult_domainwall_5din_hop2_dirac(
         float * vp, float * up, float * wp,
         float * buf2_xp, float * buf2_xm,
         float * buf2_yp, float * buf2_ym,
         float * buf2_zp, float * buf2_zm,
         float * buf2_tp, float * buf2_tm,
         int Ns, int *bc, int *Nsize, int *do_comm);

void mult_domainwall_5din_ee_5dir_dirac(
    float *vp, float *wp, float mq, float M0, int Ns, 
    float *b, float *c, float alpha, int *Nsize);

void mult_domainwall_5din_eo_5dir_dirac(
    float *yp, float *wp, float mq, float M0, int Ns, 
    float *b, float *c, float alpha, int *Nsize);

void mult_domainwall_5din_ee_5dirdag_dirac(
    float *vp, float *yp, float mq, float M0, int Ns, 
    float *b, float *c, float alpha, int *Nsize);

void mult_domainwall_5din_eo_5dirdag_dirac(
    float *vp, float *yp, float mq, float M0, int Ns, 
    float *b, float *c, float alpha, int *Nsize);

void mult_domainwall_5din_eo_hopb_dirac_4d(
    float *vp, float *up, float *wp, int Ns, int *bc, 
    int *Nsize, int *do_comm, int ieo, int jeo,
    int jgm5);

void mult_domainwall_5din_eo_hopb_dirac_5d(
    float *vp, float *up, float *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo,
    int jgm5, bool recon);

void mult_domainwall_5din_eo_hop1_dirac(
        float * buf1_xp, float * buf1_xm,
        float * buf1_yp, float * buf1_ym,
        float * buf1_zp, float * buf1_zm,
        float * buf1_tp, float * buf1_tm,
        float * up     , float * wp,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo, int jgm5);

void mult_domainwall_5din_eo_hop2_dirac(
        float * vp     , float * up, float * wp,
        float * buf2_xp, float * buf2_xm,
        float * buf2_yp, float * buf2_ym,
        float * buf2_zp, float * buf2_zm,
        float * buf2_tp, float * buf2_tm,
        int Ns, int *bc, int *Nsize, int *do_comm,
        int ieo, int jeo);

void mult_domainwall_5din_L_inv_dirac(
        float *vp, float *wp, int Ns, int *Nsize, float *e, float *dpinv, float *dm);

void mult_domainwall_5din_U_inv_dirac(
        float *vp, float *wp, int Ns, int *Nsize, float *f, float *dpinv, float *dm);

void mult_domainwall_5din_Ldag_inv_dirac(
        float *vp, float *wp, int Ns, int *Nsize, float *e, float *dpinv, float *dm);

void mult_domainwall_5din_Udag_inv_dirac(
        float *vp, float *wp, int Ns, int *Nsize, float *f, float *dpinv, float *dm);

void mult_domainwall_5din_dd_5dir_dirac(
                float * vp, float * yp, float * wp,
                float mq, float M0, int Ns, float *b, float *c,
                int *Nsize, int *block_size, int ieo);

void mult_domainwall_5din_dd_5dirdag_dirac(
                float * vp, float * yp, float * wp,
                float mq, float M0, int Ns, float *b, float *c,
                int *Nsize, int *block_size, int ieo);

void mult_domainwall_5din_dd_hopb_dirac(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *block_size, int ieo, int flag);

void mult_domainwall_5din_xpb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_xmb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_ypb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_ymb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_zpb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_zmb(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_tpb_dirac(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_tmb_dirac(
                            float * vp, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize,
                            int *do_comm, int flag);

void mult_domainwall_5din_xp1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xm1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_yp1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_ym1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zp1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zm1(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tp1_dirac(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tm1_dirac(
                            float * buf, float * up,
                            float * wp,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xp2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_xm2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_yp2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_ym2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zp2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_zm2(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tp2_dirac(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

void mult_domainwall_5din_tm2_dirac(
                            float * vp, float * up,
                            float * buf,
                            int Ns, int *bc, int *Nsize);

// anchor2
void set_block_config(float* u, int* Nsize,
                        int* block_size);
//

void mult_domainwall_5din_ee_inv_dirac_4d(
                        float * vp, float * wp, int jd,
                        int Ns, float * mat_inv, int *Nsize);

void mult_domainwall_5din_ee_inv_dirac_5d(
                        float * vp, float * wp, int jd,
                        int Ns, float * mat_inv, int *Nsize);


void mult_domainwall_5din_ee_LUdaginv_dirac(
                        float* vp, float* wp,
                        int Ns, int *Nsize,
                        float* e, float* f,
                        float* dpinv, float* dm, float alpha);

void mult_domainwall_5din_ee_LUinv_dirac(
                        float* vp, float* wp,
                        int Ns, int *Nsize,
                        float* e, float* f,
                        float* dpinv, float* dm, float alpha);

}
#endif
