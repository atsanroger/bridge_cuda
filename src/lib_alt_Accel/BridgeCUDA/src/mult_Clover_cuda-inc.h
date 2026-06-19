/*!
      @file    mult_Clover_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2538 $
*/

// This code explicitly assumes SU(3) gauge group.

#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)
#define MULT_GDXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v0+u1*v1 + u2*v2+u3*v3 + u4*v4+u5*v5)
#define MULT_GDXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v1-u1*v0 + u2*v3-u3*v2 + u4*v5-u5*v4)

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)


//====================================================================
__global__ void mult_clover_D_dirac_kernel(
                              real_t* v2, real_t* u,
                              real_t* ct, real_t* v1,
                              int Nx, int Ny, int Nz, int Nt,
                              int bc_x, int bc_y, int bc_z, int bc_t,
                              real_t kappa)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  real_t* u_up = u;
  real_t* u_dn = u;

  real_t bc2;

  real_t u_0, u_1, u_2, u_3, u_4, u_5;
  real_t u_6, u_7, u_8, u_9, u10, u11;
  real_t u12, u13, u14, u15, u16, u17;
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
  real_t wt1r, wt1i, wt2r, wt2i;

  real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
  real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
  real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
  real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

  int ix  = site % Nx;
  int izt = site/(Nx * Ny);
  int iy  = (site/Nx) % Ny;
  int it  = site/(Nx * Ny * Nz);
  int iz  = (site/(Nx * Ny)) % Nz;
  int ixy = site % (Nx * Ny);


#include "inc/mult_Wilson_xyz_cuda-inc.h"

#include "inc/mult_Wilson_t_dirac_cuda-inc.h"

#include "inc/mult_Clover_csw_dirac_cuda-inc.h"

  // write back to global memory
  v2[IDX2_SP_R(0,0,site)] = v2_01;
  v2[IDX2_SP_I(0,0,site)] = v2_11;
  v2[IDX2_SP_R(1,0,site)] = v2_21;
  v2[IDX2_SP_I(1,0,site)] = v2_31;
  v2[IDX2_SP_R(2,0,site)] = v2_41;
  v2[IDX2_SP_I(2,0,site)] = v2_51;

  v2[IDX2_SP_R(0,1,site)] = v2_02;
  v2[IDX2_SP_I(0,1,site)] = v2_12;
  v2[IDX2_SP_R(1,1,site)] = v2_22;
  v2[IDX2_SP_I(1,1,site)] = v2_32;
  v2[IDX2_SP_R(2,1,site)] = v2_42;
  v2[IDX2_SP_I(2,1,site)] = v2_52;

  v2[IDX2_SP_R(0,2,site)] = v2_03;
  v2[IDX2_SP_I(0,2,site)] = v2_13;
  v2[IDX2_SP_R(1,2,site)] = v2_23;
  v2[IDX2_SP_I(1,2,site)] = v2_33;
  v2[IDX2_SP_R(2,2,site)] = v2_43;
  v2[IDX2_SP_I(2,2,site)] = v2_53;

  v2[IDX2_SP_R(0,3,site)] = v2_04;
  v2[IDX2_SP_I(0,3,site)] = v2_14;
  v2[IDX2_SP_R(1,3,site)] = v2_24;
  v2[IDX2_SP_I(1,3,site)] = v2_34;
  v2[IDX2_SP_R(2,3,site)] = v2_44;
  v2[IDX2_SP_I(2,3,site)] = v2_54;

}

//====================================================================
void mult_clover_D_dirac(real_t* v2, real_t* u, real_t* ct, real_t* v1,
                         int* Nsize, int* bc, real_t kappa)
{

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
  real_t* u_dev  = (real_t*)dev_ptr(u);
  real_t* ct_dev  = (real_t*)dev_ptr(ct);

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  mult_clover_D_dirac_kernel<<<nbl,nth>>>(
                                     v2_dev, u_dev, ct_dev, v1_dev,
                                     Nx, Ny, Nz, Nt,
                                     bc_x, bc_y, bc_z, bc_t, kappa);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
void mult_clover_D_chiral(real_t* v2, real_t* u, real_t* ct, real_t* v1,
                          int* Nsize, int* bc, real_t kappa)
{
  printf("mult_clover_D_chiral not implemented.\n");
  exit(1);
}

//============================================================END=====
