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

  void block_axpy_eo(double *restrict vp, double *restrict out,
		     double *restrict wp, double fac,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_scal_eo(double *restrict vp, double *restrict out,
                     int Nin, int Nex,
                     int *Nsize, int *block_size, int ieo);


  void block_dotc_eo(double *restrict out,
		     double *restrict vp, double *restrict wp,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_norm2_eo(double *restrict out, double *restrict vp,
		      int Nin, int Nex,
		      int *Nsize, int *block_size, int ieo);

  void reduce_block(double *restrict out, double *restrict out_red,
		    int *Nsize, int* block_size, int ieo);

// real_t = float

  void block_axpy_eo(float *restrict vp, float *restrict out,
		     float *restrict wp, float fac,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_scal_eo(float *restrict vp, float *restrict out,
                     int Nin, int Nex,
                     int *Nsize, int *block_size, int ieo);


  void block_dotc_eo(float *restrict out,
		     float *restrict vp, float *restrict wp,
		     int Nin, int Nex,
		     int *Nsize, int *block_size, int ieo);

  void block_norm2_eo(float *restrict out, float *restrict vp,
		      int Nin, int Nex,
		      int *Nsize, int *block_size, int ieo);

  void reduce_block(float *restrict out, float *restrict out_red,
		    int *Nsize, int* block_size, int ieo);

//====================================================================

}

#endif
