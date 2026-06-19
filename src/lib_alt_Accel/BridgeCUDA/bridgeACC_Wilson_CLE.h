/*!
      @file    bridgeACC_Wilson_CLE.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2534 $
*/

#ifndef BRIDGEACC_WILSON_CLE_INCLUDED
#define BRIDGEACC_WILSON_CLE_INCLUDED

namespace BridgeACC {

// real_t = double

void mult_wilson_cle_D_dirac(
                    double *restrict v2, double *restrict u_up,
                    double *restrict u_dn, double *restrict v1,
                    int *Nsize, int *bc, double kappa);

void mult_wilson_cle_1_dirac(
                    double *restrict buf_xp, double *restrict buf_xm,
                    double *restrict buf_yp, double *restrict buf_ym,
                    double *restrict buf_zp, double *restrict buf_zm,
                    double *restrict buf_tp, double *restrict buf_tm,
                    double *restrict u, double *restrict v1,
                    int *Nsize, int *bc, int *do_comm, int Nc);

void mult_wilson_cle_2_dirac(
                    double *restrict v2, double *restrict u, 
                    double *restrict buf_xp, double *restrict buf_xm,
                    double *restrict buf_yp, double *restrict buf_ym,
                    double *restrict buf_zp, double *restrict buf_zm,
                    double *restrict buf_tp, double *restrict buf_tm,
                    double kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc);

}

#endif
