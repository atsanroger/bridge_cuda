/*!
      @file    mult_Wilson_cuda_qdw-inc.h
      @brief   QDW Wilson fermion kernels
*/

// Definition of QDW macros for projection if not in AField
#define PROJ_X_P(res, a, b) \
  dw_add(a.x, a.z, -b.y, -b.w, res.x, res.z); \
  dw_add(a.y, a.w, b.x, b.z, res.y, res.w);

#define PROJ_X_M(res, a, b) \
  dw_add(a.x, a.z, b.y, b.w, res.x, res.z); \
  dw_add(a.y, a.w, -b.x, -b.z, res.y, res.w);

#define PROJ_Y_P(res, a, b) \
  dw_add(a.x, a.z, b.x, b.z, res.x, res.z); \
  dw_add(a.y, a.w, b.y, b.w, res.y, res.w);

#define PROJ_Y_M(res, a, b) \
  dw_add(a.x, a.z, -b.x, -b.z, res.x, res.z); \
  dw_add(a.y, a.w, -b.y, -b.w, res.y, res.w);

#define PROJ_Z_P(res, a, b) \
  dw_add(a.x, a.z, -b.y, -b.w, res.x, res.z); \
  dw_add(a.y, a.w, b.x, b.z, res.y, res.w); // Same as X? No, Z gamma is different.
// Gamma matrices:
// X: 0 0 0 -i
//    0 0 -i 0
//    0 i 0 0
//    i 0 0 0
// Z: 0 0 -i 0
//    0 0 0 i
//    i 0 0 0
//    0 -i 0 0
// 1-gZ:
// s0 += i s2
// s1 -= i s3
// s2 += i s0
// s3 -= i s1
// PROJ_Z_P (1-gZ) inputs:
// For s0: (s0, s2). res = s0 + i s2.
// For s2: (s2, s0). res = s2 + i s0.
// So same macro PROJ_X_P (a + i b) works if inputs are swapped correctly.

#define PROJ_T_P(res, a) \
  dw_scal(2.0, QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(res), QDW_R_LO(res)); \
  dw_scal(2.0, QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(res), QDW_I_LO(res));

#define MULT_MI(res, a) \
  res.x = a.y; res.y = -a.x; \
  res.z = a.w; res.w = -a.z;

#define MULT_PI(res, a) \
  res.x = -a.y; res.y = a.x; \
  res.z = -a.w; res.w = a.z;

#define QDW_SCAL(res, s, a) \
  dw_scal((double)s, QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(res), QDW_R_LO(res)); \
  dw_scal((double)s, QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(res), QDW_I_LO(res));


//====================================================================
__global__
void mult_wilson_qdw_D_dirac_kernel(
                             double4* v2, real_t* u, double4* v1,
                             int Nx, int Ny, int Nz, int Nt,
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             real_t kappa)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad  = CEIL_NWP(Nst);

  real_t* u_up = u;
  real_t* u_dn = u;

  // QDW Accumulators
  // Color 0
  double4 v2_c0_s0 = {0,0,0,0}, v2_c0_s1 = {0,0,0,0}, v2_c0_s2 = {0,0,0,0}, v2_c0_s3 = {0,0,0,0};
  // Color 1
  double4 v2_c1_s0 = {0,0,0,0}, v2_c1_s1 = {0,0,0,0}, v2_c1_s2 = {0,0,0,0}, v2_c1_s3 = {0,0,0,0};
  // Color 2
  double4 v2_c2_s0 = {0,0,0,0}, v2_c2_s1 = {0,0,0,0}, v2_c2_s2 = {0,0,0,0}, v2_c2_s3 = {0,0,0,0};

  int ix  = site % Nx;
  int ly  = (site/Nx) % Ny;
  int lz  = (site/(Nx * Ny)) % Nz;
  int it  = site/(Nx * Ny * Nz);
  int ixy = site % (Nx * Ny); // Used in T kernel? Check.
  // Standard var names recalculation to match include files
  int izt = site/(Nx * Ny);
  int iy = ly;
  int iz = lz;

  #include "inc/mult_Wilson_xyz_cuda_qdw-inc.h"
  #include "inc/mult_Wilson_t_dirac_cuda_qdw-inc.h"

  // Final update: v2 = v1 - kappa * temp
  // v1, v2 are double4. kappa is double.
  // v2[site] = QDW_SUB(v1[site], QDW_SCAL(temp, kappa))
  
  // Color 0
  double4 v1_val;
  double4 tmp_k;
  
  // c0 s0
  v1_val = v1[IDX2_QDW(0, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s0);
  QDW_SUB(v2[IDX2_QDW(0, 0, site)], v1_val, tmp_k);
  
  // c0 s1
  v1_val = v1[IDX2_QDW(0, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s1);
  QDW_SUB(v2[IDX2_QDW(0, 1, site)], v1_val, tmp_k);
  
  // c0 s2
  v1_val = v1[IDX2_QDW(0, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s2);
  QDW_SUB(v2[IDX2_QDW(0, 2, site)], v1_val, tmp_k);
  
  // c0 s3
  v1_val = v1[IDX2_QDW(0, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s3);
  QDW_SUB(v2[IDX2_QDW(0, 3, site)], v1_val, tmp_k);

  // ... Repeat for Color 1, 2 ...
  // c1 s0
  v1_val = v1[IDX2_QDW(1, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s0);
  QDW_SUB(v2[IDX2_QDW(1, 0, site)], v1_val, tmp_k);
  // c1 s1
  v1_val = v1[IDX2_QDW(1, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s1);
  QDW_SUB(v2[IDX2_QDW(1, 1, site)], v1_val, tmp_k);
  // c1 s2
  v1_val = v1[IDX2_QDW(1, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s2);
  QDW_SUB(v2[IDX2_QDW(1, 2, site)], v1_val, tmp_k);
  // c1 s3
  v1_val = v1[IDX2_QDW(1, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s3);
  QDW_SUB(v2[IDX2_QDW(1, 3, site)], v1_val, tmp_k);

  // c2 s0
  v1_val = v1[IDX2_QDW(2, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s0);
  QDW_SUB(v2[IDX2_QDW(2, 0, site)], v1_val, tmp_k);
  // c2 s1
  v1_val = v1[IDX2_QDW(2, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s1);
  QDW_SUB(v2[IDX2_QDW(2, 1, site)], v1_val, tmp_k);
  // c2 s2
  v1_val = v1[IDX2_QDW(2, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s2);
  QDW_SUB(v2[IDX2_QDW(2, 2, site)], v1_val, tmp_k);
  // c2 s3
  v1_val = v1[IDX2_QDW(2, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s3);
  QDW_SUB(v2[IDX2_QDW(2, 3, site)], v1_val, tmp_k);

}

//====================================================================
__global__
void mult_wilson_qdw_D_chiral_kernel(
                             double4* v2, real_t* u, double4* v1,
                             int Nx, int Ny, int Nz, int Nt,
                             int bc_x, int bc_y, int bc_z, int bc_t,
                             real_t kappa)
{
  int site = blockIdx.x * blockDim.x + threadIdx.x;

  int Nst  = Nx * Ny * Nz * Nt;
  int Nst_pad  = CEIL_NWP(Nst);

  real_t* u_up = u;
  real_t* u_dn = u;

  // QDW Accumulators
  // Color 0
  double4 v2_c0_s0 = {0,0,0,0}, v2_c0_s1 = {0,0,0,0}, v2_c0_s2 = {0,0,0,0}, v2_c0_s3 = {0,0,0,0};
  // Color 1
  double4 v2_c1_s0 = {0,0,0,0}, v2_c1_s1 = {0,0,0,0}, v2_c1_s2 = {0,0,0,0}, v2_c1_s3 = {0,0,0,0};
  // Color 2
  double4 v2_c2_s0 = {0,0,0,0}, v2_c2_s1 = {0,0,0,0}, v2_c2_s2 = {0,0,0,0}, v2_c2_s3 = {0,0,0,0};

  int ix  = site % Nx;
  int ly  = (site/Nx) % Ny;
  int lz  = (site/(Nx * Ny)) % Nz;
  int it  = site/(Nx * Ny * Nz);
  int ixy = site % (Nx * Ny); 
  // Standard var names recalculation to match include files
  int izt = site/(Nx * Ny);
  int iy = ly;
  int iz = lz;

  #include "inc/mult_Wilson_xyz_cuda_qdw-inc.h"

  #include "inc/mult_Wilson_t_chiral_cuda_qdw-inc.h"

  // Final update: v2 = v1 - kappa * temp
  // Same as Dirac
  
  // Color 0
  double4 v1_val;
  double4 tmp_k;
  
  // c0 s0
  v1_val = v1[IDX2_QDW(0, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s0);
  QDW_SUB(v2[IDX2_QDW(0, 0, site)], v1_val, tmp_k);
  
  // c0 s1
  v1_val = v1[IDX2_QDW(0, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s1);
  QDW_SUB(v2[IDX2_QDW(0, 1, site)], v1_val, tmp_k);
  
  // c0 s2
  v1_val = v1[IDX2_QDW(0, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s2);
  QDW_SUB(v2[IDX2_QDW(0, 2, site)], v1_val, tmp_k);
  
  // c0 s3
  v1_val = v1[IDX2_QDW(0, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c0_s3);
  QDW_SUB(v2[IDX2_QDW(0, 3, site)], v1_val, tmp_k);

  // ... Repeat for Color 1, 2 ...
  // c1 s0
  v1_val = v1[IDX2_QDW(1, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s0);
  QDW_SUB(v2[IDX2_QDW(1, 0, site)], v1_val, tmp_k);
  // c1 s1
  v1_val = v1[IDX2_QDW(1, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s1);
  QDW_SUB(v2[IDX2_QDW(1, 1, site)], v1_val, tmp_k);
  // c1 s2
  v1_val = v1[IDX2_QDW(1, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s2);
  QDW_SUB(v2[IDX2_QDW(1, 2, site)], v1_val, tmp_k);
  // c1 s3
  v1_val = v1[IDX2_QDW(1, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c1_s3);
  QDW_SUB(v2[IDX2_QDW(1, 3, site)], v1_val, tmp_k);

  // c2 s0
  v1_val = v1[IDX2_QDW(2, 0, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s0);
  QDW_SUB(v2[IDX2_QDW(2, 0, site)], v1_val, tmp_k);
  // c2 s1
  v1_val = v1[IDX2_QDW(2, 1, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s1);
  QDW_SUB(v2[IDX2_QDW(2, 1, site)], v1_val, tmp_k);
  // c2 s2
  v1_val = v1[IDX2_QDW(2, 2, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s2);
  QDW_SUB(v2[IDX2_QDW(2, 2, site)], v1_val, tmp_k);
  // c2 s3
  v1_val = v1[IDX2_QDW(2, 3, site)];
  QDW_SCAL(tmp_k, kappa, v2_c2_s3);
  QDW_SUB(v2[IDX2_QDW(2, 3, site)], v1_val, tmp_k);

}

//====================================================================
__global__ void mult_wilson_qdw_gm5_dirac_kernel(double4* v2, double4* v1, int size) {
  int site = blockIdx.x * blockDim.x + threadIdx.x;
  if(site >= size) return;
  
  for(int ic=0; ic<3; ic++){
    v2[IDX2_QDW(ic, 0, site)] = v1[IDX2_QDW(ic, 0, site)];
    v2[IDX2_QDW(ic, 1, site)] = v1[IDX2_QDW(ic, 1, site)];
    
    double4 s2 = v1[IDX2_QDW(ic, 2, site)];
    double4 s3 = v1[IDX2_QDW(ic, 3, site)];
    s2.x = -s2.x; s2.y = -s2.y; s2.z = -s2.z; s2.w = -s2.w;
    s3.x = -s3.x; s3.y = -s3.y; s3.z = -s3.z; s3.w = -s3.w;
    v2[IDX2_QDW(ic, 2, site)] = s2;
    v2[IDX2_QDW(ic, 3, site)] = s3;
  }
}

//====================================================================
__global__ void mult_wilson_qdw_gm5_chiral_kernel(double4* v2, double4* v1, int size) {
  int site = blockIdx.x * blockDim.x + threadIdx.x;
  if(site >= size) return;

  for(int ic=0; ic<3; ic++){
    double4 s0 = v1[IDX2_QDW(ic, 0, site)];
    double4 s1 = v1[IDX2_QDW(ic, 1, site)];
    double4 s2 = v1[IDX2_QDW(ic, 2, site)];
    double4 s3 = v1[IDX2_QDW(ic, 3, site)];
    
    v2[IDX2_QDW(ic, 0, site)] = s2;
    v2[IDX2_QDW(ic, 1, site)] = s3;
    v2[IDX2_QDW(ic, 2, site)] = s0;
    v2[IDX2_QDW(ic, 3, site)] = s1;
  }
}
