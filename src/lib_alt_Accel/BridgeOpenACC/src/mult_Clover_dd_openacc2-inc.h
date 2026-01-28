/*!
      @file    mult_Clover_dd_openacc2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
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
void mult_clover_dd_dirac(real_t *RESTRICT v2, real_t *RESTRICT u,
                          real_t *RESTRICT ct, real_t *RESTRICT v1,
			  real_t kappa, int *bc,
			  int *Nsize, int *block_size, int ieo)
{

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;
  int nvst = NVC * ND * Nst;
  int ngst = NDF * Nst * 4;
  int ncst = NDF * Nst * 8;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  real_t *RESTRICT u_up = u;
  real_t *RESTRICT u_dn = u;

#pragma acc data present(v1[0:nvst], u_up[0:ngst], u_dn[0:ngst], \
                         ct[0:ncst], v2[0:nvst]) \
                 copyin(bc[0:4], kappa, Nx, Ny, Nz, Nt, Nst, \
                        Bx, By, Bz, Bt, Bsize, NBx, NBy, NBz, Nblock)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

    int Nst2 = Nst/2;
#pragma acc loop
    for(int site2 = 0; site2 < Nst2; ++site2){
      int bsite = site2 % Bsize;
      int block = (site2/Bsize) * 2;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);
      int jeo = (ieo + ibx + iby + ibz + ibt) % 2;
      if(jeo == 1){
        if(NBx > 1){
          ++ibx;
	}else if(NBy > 1){
          ++iby;
	}else if(NBz > 1){
          ++ibz;
	}else{
          ++ibt;
	}
      }

      /*
    for(int site2 = 0; site2 < Nst; ++site2){
      int bsite = site2 % Bsize;
      int block = site2/Bsize;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);
      int jeo = (ieo + ibx + iby + ibz + ibt) % 2;
      if(jeo == 1) continue;
      */

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite/Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt/By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt/Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

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

      {
#include "inc/mult_Wilson_openacc2_xyz-inc.h"

#include "inc/mult_Wilson_openacc2_t_dirac-inc.h"

#include "inc/mult_Clover_openacc2_csw_dirac-inc.h"
      }

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

  }

}

//====================================================================
void set_block_config(real_t *RESTRICT u, int* Nsize,
                      std::vector<int>& block_size)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);
  int Ngst = NDF * Nst_pad * 4;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

#pragma acc data present(u[0:Ngst]) \
                 copyin(Nx, Ny, Nz, Nt, Nst_pad, Bx, By, Bz, Bt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int site = 0; site < Nst_pad; ++site){

    int ix = site % Nx;
    int iy = (site/Nx) % Ny;
    int iz = (site/(Nx * Ny)) % Nz;
    int it = site/(Nx * Ny * Nz);

    int mu = 0;
    if((ix + 1) % Bx == 0){
      int istg = site + Nst_pad * mu;
      for(int in = 0; in < NDF; ++in){
        u[IDX2(NDF, in, istg)] = 0.0;
      }
    }

    mu = 1;
    if((iy + 1) % By == 0){
      int istg = site + Nst_pad * mu;
      for(int in = 0; in < NDF; ++in){
        u[IDX2(NDF, in, istg)] = 0.0;
      }
    }

    mu = 2;
    if((iz + 1) % Bz == 0){
      int istg = site + Nst_pad * mu;
      for(int in = 0; in < NDF; ++in){
        u[IDX2(NDF, in, istg)] = 0.0;
      }
    }

    mu = 3;
    if((it + 1) % Bt == 0){
      int istg = site + Nst_pad * mu;
      for(int in = 0; in < NDF; ++in){
        u[IDX2(NDF, in, istg)] = 0.0;
      }
    }

  }
 
 }

}

//============================================================END=====
