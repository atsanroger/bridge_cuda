/*!
      @file    mult_Wilson_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)
#define MULT_GDXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v0+u1*v1 + u2*v2+u3*v3 + u4*v4+u5*v5)
#define MULT_GDXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v1-u1*v0 + u2*v3-u3*v2 + u4*v5-u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)


//====================================================================
__global__
void mult_wilson_gm5_dirac_kernel(real_t* v2, real_t* v1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  for(int ic = 0; ic < NC; ++ic){
    v2[IDX2_SP_R(ic, 0, ist)] = v1[IDX2_SP_R(ic, 2, ist)];
    v2[IDX2_SP_I(ic, 0, ist)] = v1[IDX2_SP_I(ic, 2, ist)];
    v2[IDX2_SP_R(ic, 1, ist)] = v1[IDX2_SP_R(ic, 3, ist)];
    v2[IDX2_SP_I(ic, 1, ist)] = v1[IDX2_SP_I(ic, 3, ist)];
    v2[IDX2_SP_R(ic, 2, ist)] = v1[IDX2_SP_R(ic, 0, ist)];
    v2[IDX2_SP_I(ic, 2, ist)] = v1[IDX2_SP_I(ic, 0, ist)];
    v2[IDX2_SP_R(ic, 3, ist)] = v1[IDX2_SP_R(ic, 1, ist)];
    v2[IDX2_SP_I(ic, 3, ist)] = v1[IDX2_SP_I(ic, 1, ist)];
  }

}

//====================================================================
void mult_wilson_gm5_dirac(real_t* v2, real_t* v1, int *Nsize, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_gm5_dirac_kernel<<<nbl,nth>>>(v2_dev, v1_dev);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_gm5_chiral_kernel(real_t* v2, real_t* v1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  for(int ic = 0; ic < NC; ++ic){
    v2[IDX2_SP_R(ic, 0, ist)] =  v1[IDX2_SP_R(ic, 0, ist)];
    v2[IDX2_SP_I(ic, 0, ist)] =  v1[IDX2_SP_I(ic, 0, ist)];
    v2[IDX2_SP_R(ic, 1, ist)] =  v1[IDX2_SP_R(ic, 1, ist)];
    v2[IDX2_SP_I(ic, 1, ist)] =  v1[IDX2_SP_I(ic, 1, ist)];
    v2[IDX2_SP_R(ic, 2, ist)] = -v1[IDX2_SP_R(ic, 2, ist)];
    v2[IDX2_SP_I(ic, 2, ist)] = -v1[IDX2_SP_I(ic, 2, ist)];
    v2[IDX2_SP_R(ic, 3, ist)] = -v1[IDX2_SP_R(ic, 3, ist)];
    v2[IDX2_SP_I(ic, 3, ist)] = -v1[IDX2_SP_I(ic, 3, ist)];
  }

}

//====================================================================
void mult_wilson_gm5_chiral(real_t* v2, real_t* v1, int *Nsize, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int size = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = size/nth;

  mult_wilson_gm5_chiral_kernel<<<nbl,nth>>>(v2_dev, v1_dev);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_D_dirac_kernel(
                             real_t* v2, real_t* u, real_t* v1,
                             int Nx, int Ny, int Nz, int Nt,
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             real_t kappa)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad  = CEIL_NWP(Nst);

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

  // aypx and write back to global memory
  v2[IDX2_SP_R(0,0,site)] = v1[IDX2_SP_R(0,0,site)] - kappa * v2_01;
  v2[IDX2_SP_I(0,0,site)] = v1[IDX2_SP_I(0,0,site)] - kappa * v2_11;
  v2[IDX2_SP_R(1,0,site)] = v1[IDX2_SP_R(1,0,site)] - kappa * v2_21;
  v2[IDX2_SP_I(1,0,site)] = v1[IDX2_SP_I(1,0,site)] - kappa * v2_31;
  v2[IDX2_SP_R(2,0,site)] = v1[IDX2_SP_R(2,0,site)] - kappa * v2_41;
  v2[IDX2_SP_I(2,0,site)] = v1[IDX2_SP_I(2,0,site)] - kappa * v2_51;

  v2[IDX2_SP_R(0,1,site)] = v1[IDX2_SP_R(0,1,site)] - kappa * v2_02;
  v2[IDX2_SP_I(0,1,site)] = v1[IDX2_SP_I(0,1,site)] - kappa * v2_12;
  v2[IDX2_SP_R(1,1,site)] = v1[IDX2_SP_R(1,1,site)] - kappa * v2_22;
  v2[IDX2_SP_I(1,1,site)] = v1[IDX2_SP_I(1,1,site)] - kappa * v2_32;
  v2[IDX2_SP_R(2,1,site)] = v1[IDX2_SP_R(2,1,site)] - kappa * v2_42;
  v2[IDX2_SP_I(2,1,site)] = v1[IDX2_SP_I(2,1,site)] - kappa * v2_52;

  v2[IDX2_SP_R(0,2,site)] = v1[IDX2_SP_R(0,2,site)] - kappa * v2_03;
  v2[IDX2_SP_I(0,2,site)] = v1[IDX2_SP_I(0,2,site)] - kappa * v2_13;
  v2[IDX2_SP_R(1,2,site)] = v1[IDX2_SP_R(1,2,site)] - kappa * v2_23;
  v2[IDX2_SP_I(1,2,site)] = v1[IDX2_SP_I(1,2,site)] - kappa * v2_33;
  v2[IDX2_SP_R(2,2,site)] = v1[IDX2_SP_R(2,2,site)] - kappa * v2_43;
  v2[IDX2_SP_I(2,2,site)] = v1[IDX2_SP_I(2,2,site)] - kappa * v2_53;

  v2[IDX2_SP_R(0,3,site)] = v1[IDX2_SP_R(0,3,site)] - kappa * v2_04;
  v2[IDX2_SP_I(0,3,site)] = v1[IDX2_SP_I(0,3,site)] - kappa * v2_14;
  v2[IDX2_SP_R(1,3,site)] = v1[IDX2_SP_R(1,3,site)] - kappa * v2_24;
  v2[IDX2_SP_I(1,3,site)] = v1[IDX2_SP_I(1,3,site)] - kappa * v2_34;
  v2[IDX2_SP_R(2,3,site)] = v1[IDX2_SP_R(2,3,site)] - kappa * v2_44;
  v2[IDX2_SP_I(2,3,site)] = v1[IDX2_SP_I(2,3,site)] - kappa * v2_54;

}

//====================================================================
void mult_wilson_D_dirac(real_t* v2, real_t* u, real_t* v1,
                         int* Nsize, int* bc, real_t kappa)
{

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
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

  mult_wilson_D_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                          Nx, Ny, Nz, Nt,
					  bc_x, bc_y, bc_z, bc_t,
                                          kappa);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_D_chiral_kernel(
                             real_t* v2, real_t* u, real_t* v1,
                             int Nx, int Ny, int Nz, int Nt,
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             real_t kappa)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad  = CEIL_NWP(Nst);

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

#include "inc/mult_Wilson_t_chiral_cuda-inc.h"

  // aypx and write back to global memory
  v2[IDX2_SP_R(0,0,site)] = v1[IDX2_SP_R(0,0,site)] - kappa * v2_01;
  v2[IDX2_SP_I(0,0,site)] = v1[IDX2_SP_I(0,0,site)] - kappa * v2_11;
  v2[IDX2_SP_R(1,0,site)] = v1[IDX2_SP_R(1,0,site)] - kappa * v2_21;
  v2[IDX2_SP_I(1,0,site)] = v1[IDX2_SP_I(1,0,site)] - kappa * v2_31;
  v2[IDX2_SP_R(2,0,site)] = v1[IDX2_SP_R(2,0,site)] - kappa * v2_41;
  v2[IDX2_SP_I(2,0,site)] = v1[IDX2_SP_I(2,0,site)] - kappa * v2_51;

  v2[IDX2_SP_R(0,1,site)] = v1[IDX2_SP_R(0,1,site)] - kappa * v2_02;
  v2[IDX2_SP_I(0,1,site)] = v1[IDX2_SP_I(0,1,site)] - kappa * v2_12;
  v2[IDX2_SP_R(1,1,site)] = v1[IDX2_SP_R(1,1,site)] - kappa * v2_22;
  v2[IDX2_SP_I(1,1,site)] = v1[IDX2_SP_I(1,1,site)] - kappa * v2_32;
  v2[IDX2_SP_R(2,1,site)] = v1[IDX2_SP_R(2,1,site)] - kappa * v2_42;
  v2[IDX2_SP_I(2,1,site)] = v1[IDX2_SP_I(2,1,site)] - kappa * v2_52;

  v2[IDX2_SP_R(0,2,site)] = v1[IDX2_SP_R(0,2,site)] - kappa * v2_03;
  v2[IDX2_SP_I(0,2,site)] = v1[IDX2_SP_I(0,2,site)] - kappa * v2_13;
  v2[IDX2_SP_R(1,2,site)] = v1[IDX2_SP_R(1,2,site)] - kappa * v2_23;
  v2[IDX2_SP_I(1,2,site)] = v1[IDX2_SP_I(1,2,site)] - kappa * v2_33;
  v2[IDX2_SP_R(2,2,site)] = v1[IDX2_SP_R(2,2,site)] - kappa * v2_43;
  v2[IDX2_SP_I(2,2,site)] = v1[IDX2_SP_I(2,2,site)] - kappa * v2_53;

  v2[IDX2_SP_R(0,3,site)] = v1[IDX2_SP_R(0,3,site)] - kappa * v2_04;
  v2[IDX2_SP_I(0,3,site)] = v1[IDX2_SP_I(0,3,site)] - kappa * v2_14;
  v2[IDX2_SP_R(1,3,site)] = v1[IDX2_SP_R(1,3,site)] - kappa * v2_24;
  v2[IDX2_SP_I(1,3,site)] = v1[IDX2_SP_I(1,3,site)] - kappa * v2_34;
  v2[IDX2_SP_R(2,3,site)] = v1[IDX2_SP_R(2,3,site)] - kappa * v2_44;
  v2[IDX2_SP_I(2,3,site)] = v1[IDX2_SP_I(2,3,site)] - kappa * v2_54;

}

//====================================================================
void mult_wilson_D_chiral(real_t* v2, real_t* u, real_t* v1,
                          int* Nsize, int* bc, real_t kappa)
{

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
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

  mult_wilson_D_chiral_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                           Nx, Ny, Nz, Nt,
                                           bc_x, bc_y, bc_z, bc_t,
                                           kappa);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_1x_kernel(real_t* buf_xp, real_t* buf_xm,
                           real_t* u, real_t* v1,
                           int Nx, int Ny, int Nz, int Nt, 
                           int bc_x)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

 {
  int ix  = 0;
  int ist = ix + Nx * iyzt;
  real_t bc2 = bc_x;
  for(int ic = 0; ic < NC; ++ic){
    real_t vt1[2], vt2[2];
    vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
    vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
    vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
    vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
    buf_xp[IDXBF_R(ic, 0, iyzt)] = bc2 * vt1[0];
    buf_xp[IDXBF_I(ic, 0, iyzt)] = bc2 * vt1[1];
    buf_xp[IDXBF_R(ic, 1, iyzt)] = bc2 * vt2[0];
    buf_xp[IDXBF_I(ic, 1, iyzt)] = bc2 * vt2[1];
  }
 }

 {
  int ix  = Nx-1;
  int ist = ix + Nx * iyzt;

  real_t vt1[NVC], vt2[NVC], ut[NDF];
  for(int ic = 0; ic < NC; ++ic){
    int icr = 2*ic;
    int ici = 2*ic + 1;
    vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
    vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
    vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
    vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
  }

  real_t wt1[2], wt2[2];

  for(int ic = 0; ic < NC; ++ic){
     for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    buf_xm[IDXBF_R(ic, 0, iyzt)] = wt1[0];
    buf_xm[IDXBF_I(ic, 0, iyzt)] = wt1[1];
    buf_xm[IDXBF_R(ic, 1, iyzt)] = wt2[0];
    buf_xm[IDXBF_I(ic, 1, iyzt)] = wt2[1];
  }
 }

}

//====================================================================
__global__
void mult_wilson_1y_kernel(real_t* buf_yp, real_t* buf_ym,
                           real_t* u, real_t* v1,
                           int Nx, int Ny, int Nz, int Nt, 
                           int bc_y)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;

 {
  int iy   = 0;
  int ist  = ix + Nx * (iy + Ny * izt);
  int ixzt = ix + Nx * izt;
  real_t bc2 = bc_y;

  for(int ic = 0; ic < NC; ++ic){
    real_t vt1[2], vt2[2];
    vt1[0] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
    vt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
    vt2[0] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
    vt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
    buf_yp[IDXBF_R(ic, 0, ixzt)] = bc2 * vt1[0];
    buf_yp[IDXBF_I(ic, 0, ixzt)] = bc2 * vt1[1];
    buf_yp[IDXBF_R(ic, 1, ixzt)] = bc2 * vt2[0];
    buf_yp[IDXBF_I(ic, 1, ixzt)] = bc2 * vt2[1];
  }
 }

 {
  int iy   = Ny-1;
  int ist  = ix + Nx * (iy + Ny*izt);
  int istu = ist + Nst;
  int ixzt = ix + Nx * izt;

  real_t vt1[NVC], vt2[NVC];
  for(int ic = 0; ic < NC; ++ic){
    int icr = 2*ic;
    int ici = 2*ic + 1;
    vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
    vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
    vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
    vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
  }

  for(int ic = 0; ic < NC; ++ic){
    real_t ut[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
    }

    real_t wt1[2], wt2[2];
    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    buf_ym[IDXBF_R(ic, 0, ixzt)] = wt1[0];
    buf_ym[IDXBF_I(ic, 0, ixzt)] = wt1[1];
    buf_ym[IDXBF_R(ic, 1, ixzt)] = wt2[0];
    buf_ym[IDXBF_I(ic, 1, ixzt)] = wt2[1];
  }
 }

}

//====================================================================
__global__
void mult_wilson_1z_kernel(real_t* buf_zp, real_t* buf_zm,
                           real_t* u, real_t* v1,
                           int Nx, int Ny, int Nz, int Nt, 
                           int bc_z)
{
  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst = Nx * Ny * Nz * Nt;
  int Nxy = Nx * Ny;

  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;

 {
  int iz   = 0;
  int ist  = ixy + Nxy * (iz + Nz * it);
  int ixyt = ixy + Nxy * it;
  real_t bc2 = bc_z;

  for(int ic = 0; ic < NC; ++ic){
    real_t wt1[2], wt2[2];
    wt1[0] = v1[IDX2_SP_R(ic, 0, ist)] - v1[IDX2_SP_I(ic, 2, ist)];
    wt1[1] = v1[IDX2_SP_I(ic, 0, ist)] + v1[IDX2_SP_R(ic, 2, ist)];
    wt2[0] = v1[IDX2_SP_R(ic, 1, ist)] + v1[IDX2_SP_I(ic, 3, ist)];
    wt2[1] = v1[IDX2_SP_I(ic, 1, ist)] - v1[IDX2_SP_R(ic, 3, ist)];
    buf_zp[IDXBF_R(ic, 0, ixyt)] = bc2 * wt1[0];
    buf_zp[IDXBF_I(ic, 0, ixyt)] = bc2 * wt1[1];
    buf_zp[IDXBF_R(ic, 1, ixyt)] = bc2 * wt2[0];
    buf_zp[IDXBF_I(ic, 1, ixyt)] = bc2 * wt2[1];
  }
 }

 {
  int iz = Nz-1;
  int ist  = ixy + Nxy * (iz + Nz * it);
  int istu = ist + Nst * 2;
  int ixyt = ixy + Nxy * it;

  real_t vt1[NVC], vt2[NVC];
  for(int ic = 0; ic < NC; ++ic){
    int icr = 2*ic;
    int ici = 2*ic + 1;
    vt1[icr] = v1[IDX2_SP_R(ic, 0, ist)] + v1[IDX2_SP_I(ic, 2, ist)];
    vt1[ici] = v1[IDX2_SP_I(ic, 0, ist)] - v1[IDX2_SP_R(ic, 2, ist)];
    vt2[icr] = v1[IDX2_SP_R(ic, 1, ist)] - v1[IDX2_SP_I(ic, 3, ist)];
    vt2[ici] = v1[IDX2_SP_I(ic, 1, ist)] + v1[IDX2_SP_R(ic, 3, ist)];
  }

  for(int ic = 0; ic < NC; ++ic){
    real_t ut[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
    }

    real_t wt1[2], wt2[2];
    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    buf_zm[IDXBF_R(ic, 0, ixyt)] = wt1[0];
    buf_zm[IDXBF_I(ic, 0, ixyt)] = wt1[1];
    buf_zm[IDXBF_R(ic, 1, ixyt)] = wt2[0];
    buf_zm[IDXBF_I(ic, 1, ixyt)] = wt2[1];
  }

 }

}

//====================================================================
__global__
void mult_wilson_1t_dirac_kernel(real_t* buf_tp, real_t* buf_tm,
                                 real_t* u, real_t* v1,
                                 int Nx, int Ny, int Nz, int Nt, 
                                 int bc_t)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int Nst  = Nxyz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

 {
  int it = 0;
  int ist  = ixyz + Nxyz * it;
  real_t bc2 = 2.0 * bc_t;
  for(int ic = 0; ic < NC; ++ic){
    buf_tp[IDXBF_R(ic, 0, ixyz)] = bc2 * v1[IDX2_SP_R(ic, 2, ist)];
    buf_tp[IDXBF_I(ic, 0, ixyz)] = bc2 * v1[IDX2_SP_I(ic, 2, ist)];
    buf_tp[IDXBF_R(ic, 1, ixyz)] = bc2 * v1[IDX2_SP_R(ic, 3, ist)];
    buf_tp[IDXBF_I(ic, 1, ixyz)] = bc2 * v1[IDX2_SP_I(ic, 3, ist)];
  }
 }

 {
  int it   = Nt-1;
  int ist  = ixyz + Nxyz * it;
  int istu = ist + Nst_pad * 3;

  real_t vt1[NVC], vt2[NVC];
  for(int ic = 0; ic < NC; ++ic){
    int icr = 2*ic;
    int ici = 2*ic + 1;
    vt1[icr] = 2.0 * v1[IDX2_SP_R(ic, 0, ist)];
    vt1[ici] = 2.0 * v1[IDX2_SP_I(ic, 0, ist)];
    vt2[icr] = 2.0 * v1[IDX2_SP_R(ic, 1, ist)];
    vt2[ici] = 2.0 * v1[IDX2_SP_I(ic, 1, ist)];
  }

  for(int ic = 0; ic < NC; ++ic){
    real_t ut[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
    }

    real_t wt1[2], wt2[2];
    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    buf_tm[IDXBF_R(ic, 0, ixyz)] = wt1[0];
    buf_tm[IDXBF_I(ic, 0, ixyz)] = wt1[1];
    buf_tm[IDXBF_R(ic, 1, ixyz)] = wt2[0];
    buf_tm[IDXBF_I(ic, 1, ixyz)] = wt2[1];
  }

  /*
  real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
  real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
  real_t u_0, u_1, u_2, u_3, u_4, u_5;
  real_t u_6, u_7, u_8, u_9, u10, u11;
  real_t u12, u13, u14, u15, u16, u17;
  real_t wt1r, wt1i, wt2r, wt2i;

  vt1_0 = 2.0 * v1[IDX2_SP_R(0,0,ist)];
  vt1_1 = 2.0 * v1[IDX2_SP_I(0,0,ist)];
  vt1_2 = 2.0 * v1[IDX2_SP_R(1,0,ist)];
  vt1_3 = 2.0 * v1[IDX2_SP_I(1,0,ist)];
  vt1_4 = 2.0 * v1[IDX2_SP_R(2,0,ist)];
  vt1_5 = 2.0 * v1[IDX2_SP_I(2,0,ist)];

  vt2_0 = 2.0 * v1[IDX2_SP_R(0,1,ist)];
  vt2_1 = 2.0 * v1[IDX2_SP_I(0,1,ist)];
  vt2_2 = 2.0 * v1[IDX2_SP_R(1,1,ist)];
  vt2_3 = 2.0 * v1[IDX2_SP_I(1,1,ist)];
  vt2_4 = 2.0 * v1[IDX2_SP_R(2,1,ist)];
  vt2_5 = 2.0 * v1[IDX2_SP_I(2,1,ist)];

  u_0 = u[IDX2_G_R(0,0,istu)];
  u_1 = u[IDX2_G_I(0,0,istu)];
  u_2 = u[IDX2_G_R(0,1,istu)];
  u_3 = u[IDX2_G_I(0,1,istu)];
  u_4 = u[IDX2_G_R(0,2,istu)];
  u_5 = u[IDX2_G_I(0,2,istu)];

  wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

  buf_tm[IDXBF_R(0, 0, ixyz)] = wt1r;
  buf_tm[IDXBF_I(0, 0, ixyz)] = wt1i;
  buf_tm[IDXBF_R(0, 1, ixyz)] = wt2r;
  buf_tm[IDXBF_I(0, 1, ixyz)] = wt2i;

  u_6 = u[IDX2_G_R(1,0,istu)];
  u_7 = u[IDX2_G_I(1,0,istu)];
  u_8 = u[IDX2_G_R(1,1,istu)];
  u_9 = u[IDX2_G_I(1,1,istu)];
  u10 = u[IDX2_G_R(1,2,istu)];
  u11 = u[IDX2_G_I(1,2,istu)];

  wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

  buf_tm[IDXBF_R(1, 0, ixyz)] = wt1r;
  buf_tm[IDXBF_I(1, 0, ixyz)] = wt1i;
  buf_tm[IDXBF_R(1, 1, ixyz)] = wt2r;
  buf_tm[IDXBF_I(1, 1, ixyz)] = wt2i;

#ifdef SU3_3RD_ROW_RECONST
  u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
  u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
  u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
  u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
  u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
  u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
  u12 = u[IDX2_G_R(2,0,istu)];
  u13 = u[IDX2_G_I(2,0,istu)];
  u14 = u[IDX2_G_R(2,1,istu)];
  u15 = u[IDX2_G_I(2,1,istu)];
  u16 = u[IDX2_G_R(2,2,istu)];
  u17 = u[IDX2_G_I(2,2,istu)];
#endif

  wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                   vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
  wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
  wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                   vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

  buf_tm[IDXBF_R(2, 0, ixyz)] = wt1r;
  buf_tm[IDXBF_I(2, 0, ixyz)] = wt1i;
  buf_tm[IDXBF_R(2, 1, ixyz)] = wt2r;
  buf_tm[IDXBF_I(2, 1, ixyz)] = wt2i;
  */
 }

}

//====================================================================
__global__
void mult_wilson_1t_chiral_kernel(real_t* buf_tp, real_t* buf_tm,
                                 real_t* u, real_t* v1,
                                 int Nx, int Ny, int Nz, int Nt, 
                                 int bc_t)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int Nst  = Nxyz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

 {
  int it  = 0;
  int ist = ixyz + Nxyz * it;
  real_t bc2 = bc_t;
  for(int ivc = 0; ivc < NVC; ++ivc){
    buf_tp[IDXBF(ivc, 0, ixyz)]
      = bc2 * (v1[IDX2_SP(ivc, 0, ist)] + v1[IDX2_SP(ivc, 2, ist)]);
    buf_tp[IDXBF(ivc, 1, ixyz)]
      = bc2 * (v1[IDX2_SP(ivc, 1, ist)] + v1[IDX2_SP(ivc, 3, ist)]);
  }
 }

 {
  int it   = Nt-1;
  int ist  = ixyz + Nxyz * it;
  int istu = ist + Nst_pad * 3;

  real_t vt1[NVC], vt2[NVC];
  for(int ivc = 0; ivc < NVC; ++ivc){
    vt1[ivc] = v1[IDX2_SP(ivc, 0, ist)] - v1[IDX2_SP(ivc, 2, ist)];
    vt2[ivc] = v1[IDX2_SP(ivc, 1, ist)] - v1[IDX2_SP(ivc, 3, ist)];
  }

  for(int ic = 0; ic < NC; ++ic){
    real_t ut[NVC];
    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, istu)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, istu)];
    }

    real_t wt1[2], wt2[2];
    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    buf_tm[IDXBF_R(ic, 0, ixyz)] = wt1[0];
    buf_tm[IDXBF_I(ic, 0, ixyz)] = wt1[1];
    buf_tm[IDXBF_R(ic, 1, ixyz)] = wt2[0];
    buf_tm[IDXBF_I(ic, 1, ixyz)] = wt2[1];
  }

 }

}

//====================================================================
void mult_wilson_1_dirac(real_t* buf_xp, real_t* buf_xm,
                         real_t* buf_yp, real_t* buf_ym,
                         real_t* buf_zp, real_t* buf_zm,
                         real_t* buf_tp, real_t* buf_tm,
                         real_t* u, real_t* v1,
                         int *Nsize, int *bc, int *do_comm, int Nc)
{

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];

  real_t* u_dev  = (real_t*)dev_ptr(u);
  real_t* v1_dev = (real_t*)dev_ptr(v1);

  // real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
  // real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
  // real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
  // real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
  // real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
  // real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
  // real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
  // real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  if(do_comm[0] > 0){
    int size = Ny * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1x_kernel<<<nbl,nth>>>(buf_xp, buf_xm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_x);
  }

  if(do_comm[1] > 0){
    int size = Nx * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1y_kernel<<<nbl,nth>>>(buf_yp, buf_ym,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_y);
  }

  if(do_comm[2] > 0){
    int size = Nx * Ny * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1z_kernel<<<nbl,nth>>>(buf_zp, buf_zm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_z);
  }

  if(do_comm[3] > 0){
    int size = Nx * Ny * Nz;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1t_dirac_kernel<<<nbl,nth>>>(
                                       buf_tp, buf_tm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_t);
  }

  CHECK(cudaDeviceSynchronize());

  // if(do_comm[0] > 0){
  //   int size_bx = NVC * 2 * Ny * Nz * Nt;
  //   update_host(buf_xp, 0, size_bx);
  //   update_host(buf_xm, 0, size_bx);
  // }

  // if(do_comm[1] > 0){
  //   int size_by = NVC * 2 * Nx * Nz * Nt;
  //   update_host(buf_yp, 0, size_by);
  //   update_host(buf_ym, 0, size_by);
  // }

  // if(do_comm[2] > 0){
  //   int size_bz = NVC * 2 * Nx * Ny * Nt;
  //   update_host(buf_zp, 0, size_bz);
  //   update_host(buf_zm, 0, size_bz);
  // }

  // if(do_comm[3] > 0){
  //   int size_bt = NVC * 2 * Nx * Ny * Nz;
  //   update_host(buf_tp, 0, size_bt);
  //   update_host(buf_tm, 0, size_bt);
  // }


}

//====================================================================
void mult_wilson_1_chiral(real_t* buf_xp, real_t* buf_xm,
                          real_t* buf_yp, real_t* buf_ym,
                          real_t* buf_zp, real_t* buf_zm,
                          real_t* buf_tp, real_t* buf_tm,
                          real_t* u, real_t* v1,
                          int *Nsize, int *bc, int *do_comm, int Nc)
{

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  real_t* u_dev  = (real_t*)dev_ptr(u);
  real_t* v1_dev = (real_t*)dev_ptr(v1);

  if(do_comm[0] > 0){
    // real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
    // real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
    int size = Ny * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1x_kernel<<<nbl,nth>>>(buf_xp, buf_xm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_x);
  }

  if(do_comm[1] > 0){
    // real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
    // real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
    int size = Nx * Nz * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1y_kernel<<<nbl,nth>>>(buf_yp, buf_ym,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_y);
  }

  if(do_comm[2] > 0){
    // real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
    // real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
    int size = Nx * Ny * Nt;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1z_kernel<<<nbl,nth>>>(buf_zp, buf_zm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_z);
  }

  if(do_comm[3] > 0){
    // real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
    // real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);
    int size = Nx * Ny * Nz;
    int nth = VECTOR_LENGTH;
    int nbl = size/nth;
    mult_wilson_1t_chiral_kernel<<<nbl,nth>>>(
                                       buf_tp, buf_tm,
                                       u_dev, v1_dev,
                                       Nx, Ny, Nz, Nt, bc_t);
  }

  CHECK(cudaDeviceSynchronize());

  // if(do_comm[0] > 0){
  //   int size_bx = NVC * 2 * Ny * Nz * Nt;
  //   update_host(buf_xp, 0, size_bx);
  //   update_host(buf_xm, 0, size_bx);
  // }

  // if(do_comm[1] > 0){
  //   int size_by = NVC * 2 * Nx * Nz * Nt;
  //   update_host(buf_yp, 0, size_by);
  //   update_host(buf_ym, 0, size_by);
  // }

  // if(do_comm[2] > 0){
  //   int size_bz = NVC * 2 * Nx * Ny * Nt;
  //   update_host(buf_zp, 0, size_bz);
  //   update_host(buf_zm, 0, size_bz);
  // }

  // if(do_comm[3] > 0){
  //   int size_bt = NVC * 2 * Nx * Ny * Nz;
  //   update_host(buf_tp, 0, size_bt);
  //   update_host(buf_tm, 0, size_bt);
  // }

}

//====================================================================
__global__
void mult_wilson_2_dirac_kernel(
                             real_t* v2, real_t* u, 
                             real_t* buf_xp, real_t* buf_xm,
                             real_t* buf_yp, real_t* buf_ym,
                             real_t* buf_zp, real_t* buf_zm,
                             real_t* buf_tp, real_t* buf_tm,
                             real_t kappa,
                             int Nx, int Ny, int Nz, int Nt, 
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             int do_comm_x, int do_comm_y,
                             int do_comm_z, int do_comm_t)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx * Ny)) % Nz;
  int it = ist/(Nx * Ny * Nz);

  real_t v2L[NVC * ND];
  for(int id = 0; id < ND; ++id){
    for(int ivc = 0; ivc < NVC; ++ivc){
      v2L[ivc + NVC * id] = 0.0;
    }
  }

  int ncomm = 0;

#include "inc/mult_Wilson_2xyz_cuda-inc.h"

  // idir = 3
  if(do_comm_t > 0){

    if(it == Nt-1){
      int ixyz = ix + Nx * (iy + Ny * iz);
      int istu = ist + Nst_pad * 3;

      /*
      real_t vt1[NVC], vt2[NVC];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_tp[IDXBF(ivc, 0, ixyz)];
        vt2[ivc] = buf_tp[IDXBF(ivc, 1, ixyz)];
      }

      for(int ic = 0; ic < NC; ++ic){

        real_t ut[NVC];
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
          ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
        }

        real_t wt1[2], wt2[2];
        wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        v2L[ic*2   + ID3] += wt1[0];
        v2L[ic*2+1 + ID3] += wt1[1];
        v2L[ic*2   + ID4] += wt2[0];
        v2L[ic*2+1 + ID4] += wt2[1];
      }
      */
      real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
      real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1_0 = buf_tp[IDXBF(0, 0, ixyz)];
        vt1_1 = buf_tp[IDXBF(1, 0, ixyz)];
        vt1_2 = buf_tp[IDXBF(2, 0, ixyz)];
        vt1_3 = buf_tp[IDXBF(3, 0, ixyz)];
        vt1_4 = buf_tp[IDXBF(4, 0, ixyz)];
        vt1_5 = buf_tp[IDXBF(5, 0, ixyz)];
        vt2_0 = buf_tp[IDXBF(0, 1, ixyz)];
        vt2_1 = buf_tp[IDXBF(1, 1, ixyz)];
        vt2_2 = buf_tp[IDXBF(2, 1, ixyz)];
        vt2_3 = buf_tp[IDXBF(3, 1, ixyz)];
        vt2_4 = buf_tp[IDXBF(4, 1, ixyz)];
        vt2_5 = buf_tp[IDXBF(5, 1, ixyz)];
      }

      for(int ic = 0; ic < NC; ++ic){

        real_t ut_0, ut_1, ut_2, ut_3, ut_4, ut_5;
        ut_0 = u[IDX2_G_R(0, ic, istu)];
        ut_1 = u[IDX2_G_I(0, ic, istu)];
        ut_2 = u[IDX2_G_R(1, ic, istu)];
        ut_3 = u[IDX2_G_I(1, ic, istu)];
        ut_4 = u[IDX2_G_R(2, ic, istu)];
        ut_5 = u[IDX2_G_I(2, ic, istu)];

        real_t wt1_r, wt1_i, wt2_r, wt2_i;
        wt1_r = MULT_UV_R(ut_0, ut_1, ut_2, ut_3, ut_4, ut_5,
                         vt1_0,vt1_1,vt1_2,vt1_3,vt1_4,vt1_5);
        wt1_i = MULT_UV_I(ut_0, ut_1, ut_2, ut_3, ut_4, ut_5,
                         vt1_0,vt1_1,vt1_2,vt1_3,vt1_4,vt1_5);
        wt2_r = MULT_UV_R(ut_0, ut_1, ut_2, ut_3, ut_4, ut_5,
                         vt2_0,vt2_1,vt2_2,vt2_3,vt2_4,vt2_5);
        wt2_i = MULT_UV_I(ut_0, ut_1, ut_2, ut_3, ut_4, ut_5,
                         vt2_0,vt2_1,vt2_2,vt2_3,vt2_4,vt2_5);
        v2L[ic*2   + ID3] += wt1_r;
        v2L[ic*2+1 + ID3] += wt1_i;
        v2L[ic*2   + ID4] += wt2_r;
        v2L[ic*2+1 + ID4] += wt2_i;
      }
      ++ncomm;
    }

    if(it == 0){
      int ixyz = ix + Nx * (iy + Ny * iz);
      real_t bc2 = bc_t;

      for(int ic = 0; ic < NC; ++ic){
        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = bc2 * buf_tm[IDXBF_R(ic, 0, ixyz)];
        wt1i = bc2 * buf_tm[IDXBF_I(ic, 0, ixyz)];
        wt2r = bc2 * buf_tm[IDXBF_R(ic, 1, ixyz)];
        wt2i = bc2 * buf_tm[IDXBF_I(ic, 1, ixyz)];
        v2L[ic*2   + ID1] += wt1r;
        v2L[ic*2+1 + ID1] += wt1i;
        v2L[ic*2   + ID2] += wt2r;
        v2L[ic*2+1 + ID2] += wt2i;
      }
      ++ncomm;
    }

  }

  //axpy
  if(ncomm > 0){
    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        v2[IDX2_SP_R(ic, id, ist)] += -kappa * v2L[ic*2   + NVC*id];
        v2[IDX2_SP_I(ic, id, ist)] += -kappa * v2L[ic*2+1 + NVC*id];
      }
    }
  }

}

//====================================================================
void mult_wilson_2_dirac(real_t* v2, real_t* u, 
                         real_t* buf_xp, real_t* buf_xm,
                         real_t* buf_yp, real_t* buf_ym,
                         real_t* buf_zp, real_t* buf_zm,
                         real_t* buf_tp, real_t* buf_tm,
                         real_t kappa,
                         int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  // if(do_comm[0] > 0){
  //   int size_bx = NVC * 2 * Ny * Nz * Nt;
  //   update_device(buf_xp, 0, size_bx);
  //   update_device(buf_xm, 0, size_bx);
  // }

  // if(do_comm[1] > 0){
  //   int size_by = NVC * 2 * Nx * Nz * Nt;
  //   update_device(buf_yp, 0, size_by);
  //   update_device(buf_ym, 0, size_by);
  // }

  // if(do_comm[2] > 0){
  //   int size_bz = NVC * 2 * Nx * Ny * Nt;
  //   update_device(buf_zp, 0, size_bz);
  //   update_device(buf_zm, 0, size_bz);
  // }

  // if(do_comm[3] > 0){
  //   int size_bt = NVC * 2 * Nx * Ny * Nz;
  //   update_device(buf_tp, 0, size_bt);
  //   update_device(buf_tm, 0, size_bt);
  // }

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* u_dev  = (real_t*)dev_ptr(u);

  // real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
  // real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
  // real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
  // real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
  // real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
  // real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
  // real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
  // real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  mult_wilson_2_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev,
					  buf_xp, buf_xm,
					  buf_yp, buf_ym,
					  buf_zp, buf_zm,
					  buf_tp, buf_tm,
                                          kappa,
                                          Nx, Ny, Nz, Nt,
					  bc_x, bc_y, bc_z, bc_t,
                                          do_comm_x, do_comm_y,
                                          do_comm_z, do_comm_t);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_2_chiral_kernel(
                             real_t* v2, real_t* u, 
                             real_t* buf_xp, real_t* buf_xm,
                             real_t* buf_yp, real_t* buf_ym,
                             real_t* buf_zp, real_t* buf_zm,
                             real_t* buf_tp, real_t* buf_tm,
                             real_t kappa,
                             int Nx, int Ny, int Nz, int Nt, 
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             int do_comm_x, int do_comm_y,
                             int do_comm_z, int do_comm_t)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz * Nt;
  int Nst  = Nxyz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int ix = ist % Nx;
  int iy = (ist/Nx) % Ny;
  int iz = (ist/(Nx * Ny)) % Nz;
  int it = ist/Nxyz;

  real_t v2L[NVC * ND];
  for(int id = 0; id < ND; ++id){
    for(int ivc = 0; ivc < NVC; ++ivc){
      v2L[ivc + NVC * id] = 0.0;
    }
  }

  int ncomm = 0;

#include "inc/mult_Wilson_2xyz_cuda-inc.h"

  // idir = 3
  if(do_comm_t > 0){

    if(it == Nt-1){
      int ixyz = ix + Nx * (iy + Ny * iz);
      int istu = ist + Nst_pad * 3;

      real_t vt1[NVC], vt2[NVC];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_tp[IDXBF(ivc, 0, ixyz)];
        vt2[ivc] = buf_tp[IDXBF(ivc, 1, ixyz)];
      }

      for(int ic = 0; ic < NC; ++ic){

        real_t ut[NVC];
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
          ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
        }

        real_t wt1[2], wt2[2];
        wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                          vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2L[ic*2   + ID1] += wt1[0];
        v2L[ic*2+1 + ID1] += wt1[1];
        v2L[ic*2   + ID2] += wt2[0];
        v2L[ic*2+1 + ID2] += wt2[1];
        v2L[ic*2   + ID3] += wt1[0];
        v2L[ic*2+1 + ID3] += wt1[1];
        v2L[ic*2   + ID4] += wt2[0];
        v2L[ic*2+1 + ID4] += wt2[1];
      }
      ++ncomm;
    }

    if(it == 0){
      int ixyz = ix + Nx * (iy + Ny * iz);
      real_t bc2 = bc_t;

      for(int ic = 0; ic < NC; ++ic){
        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = bc2 * buf_tm[IDXBF_R(ic, 0, ixyz)];
        wt1i = bc2 * buf_tm[IDXBF_I(ic, 0, ixyz)];
        wt2r = bc2 * buf_tm[IDXBF_R(ic, 1, ixyz)];
        wt2i = bc2 * buf_tm[IDXBF_I(ic, 1, ixyz)];

	v2L[ic*2   + ID1] += wt1r;
	v2L[ic*2+1 + ID1] += wt1i;
	v2L[ic*2   + ID2] += wt2r;
	v2L[ic*2+1 + ID2] += wt2i;
	v2L[ic*2   + ID3] -= wt1r;
	v2L[ic*2+1 + ID3] -= wt1i;
	v2L[ic*2   + ID4] -= wt2r;
	v2L[ic*2+1 + ID4] -= wt2i;

      }
      ++ncomm;
    }

  }

  //axpy
  if(ncomm > 0){
    for(int id = 0; id < ND; ++id){
      for(int ivc = 0; ivc < NVC; ++ivc){
        v2[IDX2_SP(ivc, id, ist)] += -kappa * v2L[ivc + NVC * id];
      }
    }
  }

}

//====================================================================
void mult_wilson_2_chiral(real_t* v2, real_t* u, 
                          real_t* buf_xp, real_t* buf_xm,
                          real_t* buf_yp, real_t* buf_ym,
                          real_t* buf_zp, real_t* buf_zm,
                          real_t* buf_tp, real_t* buf_tm,
                          real_t kappa,
                          int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  // if(do_comm[0] > 0){
  //   int size_bx = NVC * 2 * Ny * Nz * Nt;
  //   update_device(buf_xp, 0, size_bx);
  //   update_device(buf_xm, 0, size_bx);
  // }

  // if(do_comm[1] > 0){
  //   int size_by = NVC * 2 * Nx * Nz * Nt;
  //   update_device(buf_yp, 0, size_by);
  //   update_device(buf_ym, 0, size_by);
  // }

  // if(do_comm[2] > 0){
  //   int size_bz = NVC * 2 * Nx * Ny * Nt;
  //   update_device(buf_zp, 0, size_bz);
  //   update_device(buf_zm, 0, size_bz);
  // }

  // if(do_comm[3] > 0){
  //   int size_bt = NVC * 2 * Nx * Ny * Nz;
  //   update_device(buf_tp, 0, size_bt);
  //   update_device(buf_tm, 0, size_bt);
  // }

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* u_dev  = (real_t*)dev_ptr(u);

  // real_t* buf_xp_dev = (real_t*)dev_ptr(buf_xp);
  // real_t* buf_xm_dev = (real_t*)dev_ptr(buf_xm);
  // real_t* buf_yp_dev = (real_t*)dev_ptr(buf_yp);
  // real_t* buf_ym_dev = (real_t*)dev_ptr(buf_ym);
  // real_t* buf_zp_dev = (real_t*)dev_ptr(buf_zp);
  // real_t* buf_zm_dev = (real_t*)dev_ptr(buf_zm);
  // real_t* buf_tp_dev = (real_t*)dev_ptr(buf_tp);
  // real_t* buf_tm_dev = (real_t*)dev_ptr(buf_tm);

  int bc_x = bc[0];
  int bc_y = bc[1];
  int bc_z = bc[2];
  int bc_t = bc[3];

  int do_comm_x = do_comm[0];
  int do_comm_y = do_comm[1];
  int do_comm_z = do_comm[2];
  int do_comm_t = do_comm[3];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;

  mult_wilson_2_chiral_kernel<<<nbl,nth>>>(v2_dev, u_dev,
                                           buf_xp, buf_xm,
                                           buf_yp, buf_ym,
                                           buf_zp, buf_zm,
                                           buf_tp, buf_tm,
                                           kappa,
                                           Nx, Ny, Nz, Nt,
                                           bc_x, bc_y, bc_z, bc_t,
                                           do_comm_x, do_comm_y,
                                           do_comm_z, do_comm_t);

  CHECK(cudaDeviceSynchronize());

}

//============================================================END=====
