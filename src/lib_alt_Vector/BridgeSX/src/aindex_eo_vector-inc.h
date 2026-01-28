/*!
      @file    aindex_eo_vector-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

//====================================================================
void split_dev(real_t* RESTRICT ve, real_t* RESTRICT vo,
               real_t* RESTRICT w, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst2 = Nst/2;

  int nv2 = nin * nvol2_pad;
  int nv  = nin * nvol_pad;

  for(int ist = 0; ist < Nst; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int ieo = (ix + iy + iz + it + ieo_origin) % 2;
    int ist2 = ist/2;
    if(ieo == 0){
      for(int in = 0; in < nin; ++in){
        ve[IDXV(nin, Nst2, in, ist2, 0)] = w[IDXV(nin, Nst, in, ist, 0)];
      }
    }else{
      for(int in = 0; in < nin; ++in){
        vo[IDXV(nin, Nst2, in, ist2, 0)] = w[IDXV(nin, Nst, in, ist, 0)];
      }
    }
  }

}

//====================================================================
void merge_dev(real_t* RESTRICT v,
               real_t* RESTRICT we, real_t* RESTRICT wo, int ieo_origin,
               int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst2 = Nst/2;

  int nv2 = nin * nvol2_pad;
  int nv  = nin * nvol_pad;

  for(int ist = 0; ist < Nst; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int ieo = (ix + iy + iz + it + ieo_origin) % 2;
    int ist2 = ist/2;
    if(ieo == 0){
      for(int in = 0; in < nin; ++in){
        v[IDXV(nin, Nst, in, ist, 0)] = we[IDXV(nin, Nst2, in, ist2, 0)];
      }
    }else{
      for(int in = 0; in < nin; ++in){
        v[IDXV(nin, Nst, in, ist, 0)] = wo[IDXV(nin, Nst2, in, ist2, 0)];
      }
    }
  }

}

//============================================================END=====
