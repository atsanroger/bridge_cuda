/*!
      @file    bridgeVec_aindex_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_AINDEX_EO_INCLUDED
#define BRIDGEVEC_AINDEX_EO_INCLUDED

namespace BridgeVec {

// real_t = double
void split_dev(double* RESTRICT ve, double* RESTRICT vo,
               double* RESTRICT w, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge_dev(double* RESTRICT v,
               double* RESTRICT we, double* RESTRICT wo, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad);

// real_t = float
void split_dev(float* RESTRICT ve, float* RESTRICT vo,
               float* RESTRICT w, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad);

void merge_dev(float* RESTRICT v,
               float* RESTRICT we, float* RESTRICT wo, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad);

}

#endif

//============================================================END=====
