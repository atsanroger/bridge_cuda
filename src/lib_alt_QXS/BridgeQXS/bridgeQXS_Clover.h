/*!
      @file    bridgeQXS_Clover.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_CLOVER_INCLUDED
#define BRIDGEQXS_CLOVER_INCLUDED

namespace BridgeQXS {
  void mult_clover_bulk_dirac(double *v2, double *up, double *ct,
                              double *v1, double kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_clover_bulk_dirac_chrot(
    double *v2, double *up, double *ct,
    double *v1, double kappa, int *bc,
    int *Nsize, int *do_comm);

  void mult_clover_dd_dirac(double *v2, double *up, double *ct,
                            double *v1, double kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void mult_clover_dd_dirac_chrot(
    double *v2, double *up, double *ct,
    double *v1, double kappa, int *bc,
    int *Nsize, int *block_size, int ieo);

  void mult_clover_cswinv_dirac_chrot(
    double *v2, double *ct, double *v1,
    int *Nsize);


  void mult_clover_bulk_dirac(float *v2, float *up, float *ct,
                              float *v1, float kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_clover_bulk_dirac_chrot(
    float *v2, float *up, float *ct,
    float *v1, float kappa, int *bc,
    int *Nsize, int *do_comm);

  void mult_clover_dd_dirac(float *v2, float *up, float *ct,
                            float *v1, float kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void mult_clover_dd_dirac_chrot(
    float *v2, float *up, float *ct,
    float *v1, float kappa, int *bc,
    int *Nsize, int *block_size, int ieo);

  void mult_clover_cswinv_dirac_chrot(
    float *v2, float *ct, float *v1,
    int *Nsize);
}
#endif
