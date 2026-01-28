/*!
      @file    index_eo_alt_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-29 15:46:01 #$
      @version $LastChangedRevision: 2542 $
*/

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

//====================================================================
void split(real_t* restrict ve, real_t* restrict vo,
           real_t* restrict w, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nvol  = Nx * Ny * Nz * Nt;
  int Nvol2 = Nvol/2;

  int nv2 = nin * nvol2_pad;
  int nv  = nin * nvol_pad;

#pragma acc data present(ve[0:nv2], vo[0:nv2], w[0:nv]), \
                 copyin(ieo_origin, Nx, Ny, Nz, Nt, Nvol, Nvol2)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop
    for(int ist2 = nvol2_pad - NWP; ist2 < nvol2_pad; ++ist2){
      for(int in = 0; in < nin; ++in){
        ve[IDX2(nin, in, ist2)] = 0.0;
        vo[IDX2(nin, in, ist2)] = 0.0;
      }
    }
  }

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop
    for(int ist = 0; ist < Nvol; ++ist){
      int ix = ist % Nx;
      int iy = (ist/Nx) % Ny;
      int iz = (ist/(Nx*Ny)) % Nz;
      int it = ist/(Nx*Ny*Nz);
      int ieo = (ix + iy + iz + it + ieo_origin) % 2;
      int ist2 = ist/2;
      if(ieo == 0){
        for(int in = 0; in < nin; ++in){
          ve[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
        }
      }else{
        for(int in = 0; in < nin; ++in){
          vo[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
        }
      }
    }
  }

 } // acc data

}

//====================================================================
void merge(real_t* restrict v,
           real_t* restrict we, real_t* restrict wo, int ieo_origin,
           int nin, int* Nsize, int nvol2_pad, int nvol_pad)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nvol  = Nx * Ny * Nz * Nt;
  int Nvol2 = Nvol/2;

  int nv2 = nin * nvol2_pad;
  int nv  = nin * nvol_pad;

#pragma acc data present(we[0:nv2], wo[0:nv2], v[0:nv]), \
                 copyin(ieo_origin, Nx, Ny, Nz, Nt, Nvol, Nvol2)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < nvol_pad; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int ieo = (ix + iy + iz + it + ieo_origin) % 2;
    int ist2 = ist/2;
    if(ist < Nvol){
      if(ieo == 0){
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = we[IDX2(nin, in, ist2)];
        }
      }else{
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = wo[IDX2(nin, in, ist2)];
        }
      }

    }else{
      for(int in = 0; in < nin; ++in){
        v[IDX2(nin, in, ist)] = 0.0;
      }
    }
  }
 }

}

//============================================================END=====
