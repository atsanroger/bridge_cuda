/*!
      @file    index_eo_alt_openacc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-12 13:12:27 #$
      @version $LastChangedRevision: 2539 $
*/

#ifndef BRIDGEACC_INDEX_EO_ALT_INCLUDED
#define BRIDGEACC_INDEX_EO_ALT_INCLUDED

namespace BridgeACC {

// real_t = double
void split(double* ve, double* vo,
           double* w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge(double* v,
           double* we, double* wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

// real_t = float
void split(float* ve, float* vo,
           float* w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge(float* v,
           float* we, float* wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad);

}

//============================================================END=====
#endif
