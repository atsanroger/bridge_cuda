/*!
      @file    afopr_Wilson_simd2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "mult_common_th_simd2-inc.h"
#include "mult_Wilson_simd2_parts-inc.h"

namespace QYS{

//====================================================================
void mult_wilson_bulk_dirac(real_t *v2, real_t *u, real_t *v1,
                     real_t kappa, int *bc, int *Nsize, int *do_comm)
{
  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  int Nxy = Nxv * Ny;
  int Nxyz = Nxv * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  // bulk
  for(int site = is; site < ns; ++site){
    int ix   = site % Nxv;
    int iyzt = site/Nxv;
    int iy  = iyzt % Ny;
    int izt = site/Nxy;
    int iz  = izt % Nz;
    int it  = izt/Nz;
    int ixy = ix + Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    clear_vec(v2v, NCD);

    if(ix < Nxv-1 || do_comm[0] == 0){
      real_t *ut = &u[NDF * Nst * 0];
      int nei = ix+1 + Nxv * iyzt;
      if(ix == Nxv-1) nei = 0 + Nxv * iyzt;
      shift_vec2_bw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NC*ND);
      mult_wilson_xpb(v2v, &ut[VLEN*NDF2*site], zL, NC);
    }

    if(ix > 0 ||  do_comm[0] == 0){
      real_t *ut = &u[NDF * Nst * 0];
      int nei = ix-1 + Nxv * iyzt;
      if(ix == 0) nei = Nxv-1 + Nxv * iyzt;
      shift_vec2_fw(zL, &v1[VLEN*NCD*site], &v1[VLEN*NCD*nei], NCD);
      shift_vec2_fw(uL, &ut[VLEN*NDF2*site], &ut[VLEN*NDF2*nei], NDF2);
      mult_wilson_xmb(v2v, uL, zL, NC);
    }

    if(iy < Ny-1 ||  do_comm[1] == 0){
      real_t *ut = &u[NDF * Nst * 1];
      int iy2 = (iy + 1) % Ny;
      int nei = ix + Nxv * (iy2 + Ny * izt);
      mult_wilson_ypb(v2v, &ut[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }

    if(iy > 0 ||  do_comm[1] == 0){
      real_t *ut = &u[NDF * Nst * 1];
      int iy2 = (iy - 1 + Ny) % Ny;
      int nei = ix + Nxv * (iy2 + Ny * izt);
      mult_wilson_ymb(v2v, &ut[VLEN*NDF2*nei], &v1[VLEN*NCD*nei], NC);
    }

    if(iz < Nz-1 ||  do_comm[2] == 0){
      real_t *ut = &u[NDF * Nst * 2];
      int iz2 = (iz + 1) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz*it);
      mult_wilson_zpb(v2v, &ut[VLEN*NDF2*site], &v1[VLEN*NCD*nei], NC);
    }

    if(iz > 0 ||  do_comm[2] == 0){
      real_t *ut = &u[NDF * Nst * 2];
      int iz2 = (iz - 1 + Nz) % Nz;
      int nei = ixy + Nxy * (iz2 + Nz*it);
      mult_wilson_zmb(v2v, &ut[VLEN*NDF2*nei],
                      &v1[VLEN*NCD*nei], NC);
    }

    if(it < Nt-1 ||  do_comm[3] == 0){
      real_t *ut = &u[NDF * Nst * 3];
      int it2 = (it + 1) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tpb_dirac(v2v, &ut[VLEN*NDF2*site],
                            &v1[VLEN*NCD*nei], NC);
    }

    if(it > 0 ||  do_comm[3] == 0){
      real_t *ut = &u[NDF * Nst * 3];
      int it2 = (it - 1 + Nt) % Nt;
      int nei = ixyz + Nxyz * it2;
      mult_wilson_tmb_dirac(v2v, &ut[VLEN*NDF2*nei],
                            &v1[VLEN*NCD*nei], NC);
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

}

//====================================================================
void mult_wilson_1_dirac(
                    real_t *buf_xp, real_t *buf_xm,
                    real_t *buf_yp, real_t *buf_ym,
                    real_t *buf_zp, real_t *buf_zm,
                    real_t *buf_tp, real_t *buf_tm,
                    real_t *u, real_t *v1,
                    int *bc, int *Nsize, int *do_comm)
{
  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  int Nxy = Nxv * Ny;
  int Nxyz = Nxv * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int site = is; site < ns; ++site){
    int ix   = site % Nxv;
    int iyzt = site/Nxv;
    int iy  = iyzt % Ny;
    int izt = site/Nxy;
    int iz  = izt % Nz;
    int it  = izt/Nz;
    int ixy = ix + Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    if(do_comm[0] == 1){
      if(ix == 0){
        mult_wilson_xp1(&buf_xp[NVC*ND2*iyzt], &v1[VLEN * NCD*site], NC);
      }
      if(ix == Nxv-1 ){
        real_t *ut = &u[NDF * Nst * 0];
        mult_wilson_xm1(&buf_xm[NVC*ND2*iyzt], &ut[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == 0){
        int ibf = VLEN * NC * ND2 * (ix + Nxv * izt);
        mult_wilson_yp1(&buf_yp[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iy == Ny-1){
        real_t *ut = &u[NDF * Nst * 1];
        int ibf = VLEN * NC * ND2 * (ix + Nxv * izt);
        mult_wilson_ym1(&buf_ym[ibf], &ut[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == 0){
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp1(&buf_zp[ibf], &v1[VLEN*NCD*site], NC);
      }
      if(iz == Nz-1){
        real_t *ut = &u[NDF * Nst * 2];
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zm1(&buf_zm[ibf], &ut[VLEN*NDF2*site],
                        &v1[VLEN*NCD*site], NC);
      }
    }

    if(do_comm[3] == 1){
      if(it == 0){
        mult_wilson_tp1_dirac(&buf_tp[VLEN*NC*ND2*ixyz],
                              &v1[VLEN*NCD*site], NC);
      }
      if(it == Nt-1){
        real_t *ut = &u[NDF * Nst * 3];
        mult_wilson_tm1_dirac(&buf_tm[VLEN*NC*ND2*ixyz], 
                       &ut[VLEN*NDF2*site], &v1[VLEN*NCD*site], NC);
      }
    }

  }

}

//====================================================================
void mult_wilson_2_dirac(real_t *v2, real_t *u, real_t *v1,
                         real_t *buf_xp, real_t *buf_xm,
                         real_t *buf_yp, real_t *buf_ym,
                         real_t *buf_zp, real_t *buf_zm,
                         real_t *buf_tp, real_t *buf_tm,
                         real_t kappa, int *bc, int *Nsize, int *do_comm)
{
  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  int Nxy = Nxv * Ny;
  int Nxyz = Nxv * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  // boundary
  for(int site = is; site < ns; ++site){
    int ix   = site % Nxv;
    int iyzt = site/Nxv;
    int iy  = iyzt % Ny;
    int izt = site/Nxy;
    int iz  = izt % Nz;
    int it  = izt/Nz;
    int ixy = ix + Nxv * iy;
    int ixyz = ixy + Nxy * iz;

    load_vec(v2v, &v2[VLEN*NCD*site], NCD);

    if(do_comm[0] == 1){
      if(ix == Nxv-1){
        real_t *ut = &u[NDF * Nst * 0];
        shift_vec0_bw(zL, &v1[VLEN*NCD*site], NCD);
        mult_wilson_xpb(v2v, &ut[VLEN*NDF2*site], zL, NC);
        mult_wilson_xp2(v2v, &ut[VLEN*NDF2*site],
                             &buf_xp[NVC*ND2*iyzt], NC);
      }
      if(ix == 0){
        real_t *ut = &u[NDF * Nst * 0];
        shift_vec0_fw(zL, &v1[VLEN*NCD*site], NCD);
        shift_vec0_fw(uL, &ut[VLEN*NDF2*site], NDF2);
        mult_wilson_xmb(v2v, uL, zL, NC);
        mult_wilson_xm2(v2v, &buf_xm[NVC*ND2*iyzt], NC);
      }
    }

    if(do_comm[1] == 1){
      if(iy == Ny-1){
        real_t *ut = &u[NDF * Nst * 1];
        int ibf = VLEN * NC * ND2 * (ix + Nxv * izt);
        mult_wilson_yp2(v2v, &ut[VLEN*NDF2*site], &buf_yp[ibf], NC);
      }
      if(iy == 0){
        int ibf = VLEN * NC * ND2 * (ix + Nxv * izt);
        mult_wilson_ym2(v2v, &buf_ym[ibf], NC);
      }
    }

    if(do_comm[2] == 1){
      if(iz == Nz-1){
        real_t *ut = &u[NDF * Nst * 2];
        int ibf = VLEN * NC * ND2 * (ixy + Nxy * it);
        mult_wilson_zp2(v2v, &ut[VLEN*NDF2*site], &buf_zp[ibf], NC);
      }
      if(iz == 0){
        int ibf = VLEN * NC *ND2 * (ixy + Nxy * it);
        mult_wilson_zm2(v2v, &buf_zm[ibf], NC);
      }
    }

    if(do_comm[3] == 1){
      if(it == Nt-1){
        real_t *ut = &u[NDF * Nst * 3];
        mult_wilson_tp2_dirac(v2v, &ut[VLEN*NDF2*site],
                              &buf_tp[VLEN*NC*ND2*ixyz], NC);
      }
      if(it == 0){
        mult_wilson_tm2_dirac(v2v, &buf_tm[VLEN*NC*ND2*ixyz], NC);
      }
    }

    Vsimd_t v1v[NCD];
    load_vec(v1v, &v1[VLEN*NCD*site], NCD);
    aypx_vec(-kappa, v2v, v1v, NCD);
    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

}

//============================================================END=====

} // end of namespace QYS
