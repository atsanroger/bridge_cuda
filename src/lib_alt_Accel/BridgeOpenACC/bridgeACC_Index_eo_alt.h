/*!
      @file    index_eo_alt_openacc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef BRIDGEACC_INDEX_EO_ALT_INCLUDED
#define BRIDGEACC_INDEX_EO_ALT_INCLUDED

namespace BridgeACC {

// real_t = double
void split(double* restrict ve, double* restrict vo,
           double* restrict w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge(double* restrict v,
           double* restrict we, double* restrict wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

// real_t = float
void split(float* restrict ve, float* restrict vo,
           float* restrict w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge(float* restrict v,
           float* restrict we, float* restrict wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

}

//============================================================END=====
#endif
