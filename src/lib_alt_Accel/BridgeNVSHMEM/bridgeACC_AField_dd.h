/*!
      @file    bridgeACC_AField_dd.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-03-21 21:09:15 #$
      @version $LastChangedRevision: 2587 $
*/

#ifndef BRIDGEACC_AFIELD_DD_INCLUDED
#define BRIDGEACC_AFIELD_DD_INCLUDED

namespace BridgeACC {

// real_t = double

  void block_axpy_eo(double* vp, double* out,
		     double* wp, double fac,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_scal_eo(double* vp, double* out,
                     int Nin, int Nex,
                     int *Nsize, int *block_size, int ieo);


  void block_dotc_eo(double* out,
		     double* vp, double* wp,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_norm2_eo(double* out, double* vp,
		      int Nin, int Nex,
		      int *Nsize, int *block_size, int ieo);

  void reduce_block(double* out, double* out_red,
                    int* Nsize, int* block_size, int ieo);

// real_t = float

  void block_axpy_eo(float* vp, float* out,
		     float* wp, float fac,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_scal_eo(float* vp, float* out,
                     int Nin, int Nex,
                     int *Nsize, int *block_size, int ieo);


  void block_dotc_eo(float* out,
		     float* vp, float* wp,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_norm2_eo(float* out, float* vp,
		      int Nin, int Nex,
		      int *Nsize, int *block_size, int ieo);

  void reduce_block(float* out, float* out_red,
                    int* Nsize, int* block_size, int ieo);

//====================================================================

}

#endif
