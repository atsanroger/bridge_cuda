/*!
      @file    mult_Clover_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2606 $
*/

#ifndef BRIDGEACC_CLOVER_INCLUDED
#define BRIDGEACC_CLOVER_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_clover_H_dirac(double* v2, double* u,
                           double* ct, double* v1,
                           int *Nsize, int *bc, double kappa);

  void mult_clover_D_dirac(double* v2, double* u,
                           double* ct, double* v1,
                           int *Nsize, int *bc, double kappa);

  void mult_clover_H_chiral(double* v2, double* u,
                            double* ct, double* v1,
                            int *Nsize, int *bc, double kappa);

  void mult_clover_D_chiral(double* v2, double* u,
                            double* ct, double* v1,
                            int *Nsize, int *bc, double kappa);

  void mult_csw_dirac(double* v2, double* u,
                      double* v1, int *Nsize, int iflag, int Nc);

  void mult_csw_chiral(double* v2, double* u,
                       double* v1, int *Nsize, int iflag, int Nc);


  void mult_clover_dd_dirac(double* v2, double* up,
                            double* ct, double* v1,
                            double kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void set_block_config(double* u, int* Nsize,
                        int* block_size);
  //                        std::vector<int>& block_size);

  // real_t = float

  void mult_clover_H_dirac(float* v2, float* u,
                           float* ct, float* v1,
                           int *Nsize, int *bc, float kappa);

  void mult_clover_D_dirac(float* v2, float* u,
                           float* ct, float* v1,
                           int *Nsize, int *bc, float kappa);

  void mult_clover_H_chiral(float* v2, float* u,
                            float* ct, float* v1,
                            int *Nsize, int *bc, float kappa);

  void mult_clover_D_chiral(float* v2, float* u,
                            float* ct, float* v1,
                            int *Nsize, int *bc, float kappa);

  void mult_csw_dirac(float* v2, float* u,
                      float* v1, int *Nsize, int iflag, int Nc);

  void mult_csw_chiral(float* v2, float* u,
                       float* v1, int *Nsize, int iflag, int Nc);

  void mult_clover_dd_dirac(float* v2, float* up,
                            float* ct, float* v1,
                            float kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void set_block_config(float* u, int* Nsize,
                        int* block_size);
  //                      std::vector<int>& block_size);

}
#endif

//============================================================END=====
