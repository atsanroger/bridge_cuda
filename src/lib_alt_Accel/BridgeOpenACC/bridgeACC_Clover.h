/*!
      @file    mult_Clover_openacc2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

#ifndef BRIDGEACC_CLOVER_INCLUDED
#define BRIDGEACC_CLOVER_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_clover_H_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT ct, double *RESTRICT v1,
                           int *Nsize, int *bc, double kappa);

  void mult_clover_D_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT ct, double *RESTRICT v1,
                           int *Nsize, int *bc, double kappa);

  void mult_clover_H_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT ct, double *RESTRICT v1,
                            int *Nsize, int *bc, double kappa);

  void mult_clover_D_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT ct, double *RESTRICT v1,
                            int *Nsize, int *bc, double kappa);

  void mult_csw_dirac(double *RESTRICT v2, double *RESTRICT u,
                      double *RESTRICT v1, int *Nsize, int iflag, int Nc);

  void mult_csw_chiral(double *RESTRICT v2, double *RESTRICT u,
                       double *RESTRICT v1, int *Nsize, int iflag, int Nc);


  void mult_clover_dd_dirac(double *RESTRICT v2, double *RESTRICT up,
                            double *RESTRICT ct, double *RESTRICT v1,
                            double kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  // void set_block_config(double *RESTRICT u, int* Nsize,
  //                      std::vector<int>& block_size);

  // real_t = float

  void mult_clover_H_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT ct, float *RESTRICT v1,
                           int *Nsize, int *bc, float kappa);

  void mult_clover_D_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT ct, float *RESTRICT v1,
                           int *Nsize, int *bc, float kappa);

  void mult_clover_H_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT ct, float *RESTRICT v1,
                            int *Nsize, int *bc, float kappa);

  void mult_clover_D_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT ct, float *RESTRICT v1,
                            int *Nsize, int *bc, float kappa);

  void mult_csw_dirac(float *RESTRICT v2, float *RESTRICT u,
                      float *RESTRICT v1, int *Nsize, int iflag, int Nc);

  void mult_csw_chiral(float *RESTRICT v2, float *RESTRICT u,
                       float *RESTRICT v1, int *Nsize, int iflag, int Nc);

  void mult_clover_dd_dirac(float *RESTRICT v2, float *RESTRICT up,
                            float *RESTRICT ct, float *RESTRICT v1,
                            float kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  //void set_block_config(float *RESTRICT u, int* Nsize,
  //                      std::vector<int>& block_size);

}
#endif

//============================================================END=====
