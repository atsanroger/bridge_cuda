/*!
      @file    bridgeACC_Wilson_double.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"
#include "bridgeACC.h"

#include "bridgeACC_AField.h"


typedef double real_t;


namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "src/mult_Wilson_cuda-inc.h"
#include "src/mult_Wilson_cuda2-inc.h"
#include "src/mult_Wilson_dd_cuda-inc.h"
#include "src/mult_Wilson_cuda_qdw-inc.h"

//====================================================================
void mult_wilson_qdw_D_dirac(double* v2, double* u, double* v1,
                             int* Nsize, int* bc, double kappa)
{
  // Cast to specific types for QDW kernel
  // v1, v2 are double* but hold double4 data (QDW).
  // u is standard double* (Gauge).
  // The kernel expects double4* for v1, v2.
  
  double4* v2_dev = (double4*)dev_ptr(v2);
  double4* v1_dev = (double4*)dev_ptr(v1);
  real_t* u_dev  = (real_t*)dev_ptr(u);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_qdw_D_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                              Nx, Ny, Nz, Nt,
                                              bc_x, bc_y, bc_z, bc_t,
                                              kappa);

  CHECK(cudaDeviceSynchronize());
}

void mult_wilson_qdw_D_chiral(double* v2, double* u, double* v1,
                              int* Nsize, int* bc, double kappa)
{
  double4* v2_dev = (double4*)dev_ptr(v2);
  double4* v1_dev = (double4*)dev_ptr(v1);
  real_t* u_dev  = (real_t*)dev_ptr(u);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_qdw_D_chiral_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                               Nx, Ny, Nz, Nt,
                                               bc_x, bc_y, bc_z, bc_t,
                                               kappa);

  CHECK(cudaDeviceSynchronize());
}

void mult_wilson_qdw_gm5_dirac(double* v2, double* v1, int* Nsize, int Nc)
{
  double4* v2_dev = (double4*)dev_ptr(v2);
  double4* v1_dev = (double4*)dev_ptr(v1);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_qdw_gm5_dirac_kernel<<<nbl,nth>>>(v2_dev, v1_dev, size);

  CHECK(cudaDeviceSynchronize());
}

void mult_wilson_qdw_gm5_chiral(double* v2, double* v1, int* Nsize, int Nc)
{
  double4* v2_dev = (double4*)dev_ptr(v2);
  double4* v1_dev = (double4*)dev_ptr(v1);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_qdw_gm5_chiral_kernel<<<nbl,nth>>>(v2_dev, v1_dev, size);

  CHECK(cudaDeviceSynchronize());
}



}

//============================================================END=====
