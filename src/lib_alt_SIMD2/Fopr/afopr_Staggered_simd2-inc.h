/*!
      @file    afopr_Staggered_simd2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

// This code explicitly assumes SU(3) gauge group.
#include "lib_alt_SIMD2/inline/define_params_SU3.h"

namespace QYS{

//====================================================================
void mult_staggered_parity(real_t *v, real_t *ph, int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    Vsimd_t pht, vt, wt;
    load_vec(&pht, &ph[VLEN * ist], 1);
    for(int ic = 0; ic < NC; ++ic){
      load_vec(&vt, &v[VLEN*(ic + NC*ist)], 1);
      mult_nn_vec(wt, pht, vt);
      save_vec(&v[VLEN*(ic + NC*ist)], &wt, 1);
    }
  }

}

//====================================================================
void mult_staggered_stg_phase(real_t *v, real_t *ph, int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    Vsimd_t pht, vt, wt;
    load_vec(&pht, &ph[VLEN * ist], 1);
    for(int idf = 0; idf < NDF2; ++idf){
      load_vec(&vt, &v[VLEN*(idf + NDF2*ist)], 1);
      mult_nn_vec(wt, pht, vt);
      save_vec(&v[VLEN*(idf + NDF2*ist)], &wt, 1);
    }
  }

}

//====================================================================
void mult_staggered_aypx(real_t a, real_t *v2, real_t *v1,
                         int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    Vsimd_t vt2[NC], vt1[NC];
    load_vec(vt1, &v1[VLEN * NC * ist], NC);
    load_vec(vt2, &v2[VLEN * NC * ist], NC);
    aypx_vec(a, vt2, vt1, NC);
    save_vec(&v2[VLEN * NC * ist], vt2, NC);
  }

}

//====================================================================
void mult_staggered_axpy(real_t *v2, real_t a, real_t *v1,
                         int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    Vsimd_t vt2[NC], vt1[NC];
    load_vec(vt1, &v1[VLEN * NC * ist], NC);
    load_vec(vt2, &v2[VLEN * NC * ist], NC);
    axpy_vec(vt2, a, vt1, NC);
    save_vec(&v2[VLEN * NC * ist], vt2, NC);
  }

}

//====================================================================
void mult_staggered_scal(real_t *v, real_t a, int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    Vsimd_t vt[NC];
    load_vec(vt, &v[VLEN * NC * ist], NC);
    scal_vec(vt, a, NC);
    save_vec(&v[VLEN * NC * ist], vt, NC);
  }

}

//====================================================================
void mult_staggered_clear(real_t *v, int *Nsize, int Nc)
{
  int Nstv = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  Vsimd_t vt[NC];
  clear_vec(vt, NC);
  for(int ist = is; ist < ns; ++ist){
    save_vec(&v[VLEN * NC * ist], vt, NC);
  }

}

//====================================================================
void mult_staggered_xp1(real_t *buf, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nyzt);

  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = 0;
    int ist = ix + Nxv * iyzt;

    real_t vt[NVC];
    load_vec1(vt, &v1[VLEN*NC*ist], 0, NC);

    for(int ivc = 0; ivc < NVC; ++ivc){
      buf[ivc + NVC*iyzt] = bc2 * vt[ivc];
    }

  }

}

//====================================================================
void mult_staggered_xp2(real_t *v2, real_t *u, real_t *buf, 
                        int *bc, int *Nsize, int* do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  real_t v0[VLEN*NC];
  for(int ic = 0; ic < NC; ++ic){
    for(int k = 0; k < VLEN; ++k){
      v0[k + VLEN*ic] = 0.0;
    }
  }

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nxv;
    int iyzt = ist/Nxv;

    if(ix == Nxv-1){
      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      //shift_vec1_bw(vt, &v1[VLEN*NC*ist], &buf[NVC*iyzt], NC);
      shift_vec1_bw(vt, v0, &buf[NVC*iyzt], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],  NDF2);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);
    }

  }

}

//====================================================================
void mult_staggered_xpb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int* do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nxv;
    int iyzt = ist/Nxv;

    Vsimd_t v2v[NC], vt[NC], ut[NDF2];

    if(ix < Nxv-1){
      int nei = ix+1 + Nxv * iyzt;
      shift_vec2_bw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], NC);
    }else{
      if(do_comm[idir] == 0){
        int nei = 0 + Nxv * iyzt;
        shift_vec2_bw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], bc2, NC);
      }else{
        shift_vec0_bw(vt, &v1[VLEN*NC*ist], NC);
      }
    }

    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],  NDF2);

    for(int ic = 0; ic < NC; ++ic){
      int ic2 = NC * ic;
      mult_uv(v2v[ic], &ut[ic2], vt, NC);
    }

    add_vec(&v2[VLEN*NC*ist], v2v, NC);

  }

}

//====================================================================
void mult_staggered_xm1(real_t *buf, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nyzt);

  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = Nxv-1;
    int ist = ix + Nxv * iyzt;

    Vsimd_t vt[NC], ut[NDF2], wt;
    real_t  wt2[VLEN];

    load_vec(vt, &v1[VLEN*NC*ist], NC);
    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],  NDF2);

    for(int ic = 0; ic < NC; ++ic){
      mult_udagv(wt, &ut[ic], vt, NC);
      save_vec(wt2, &wt, 1);
      buf[2*ic   + NVC*iyzt] = wt2[VLEN-2];
      buf[2*ic+1 + NVC*iyzt] = wt2[VLEN-1];
    }

  }

}

//====================================================================
void mult_staggered_xm2(real_t *v2, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  real_t bc2 = real_t(bc[idir]);

  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nxv;
    int iyzt = ist/Nxv;

    if(ix == 0){
      for(int ic = 0; ic < NC; ++ic){
        v2[0 + VLEN*(ic + NC*ist)] += -bc2 * buf[2*ic   + NVC*iyzt];
        v2[1 + VLEN*(ic + NC*ist)] += -bc2 * buf[2*ic+1 + NVC*iyzt];
      }
    }

  }

}

//====================================================================
void mult_staggered_xmb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 0;

  int Nxv  = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nstv = Nxv * Nyzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nxv;
    int iyzt = ist/Nxv;

    Vsimd_t v2v[NC];
    Vsimd_t vt[NC], ut[NDF2];

    if(ix > 0){
      int nei = ix-1 + Nxv * iyzt;
      shift_vec2_fw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], NC);
      shift_vec2_fw(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],
                        &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(v2v[ic], &ut[ic], vt, NC);
      }

    }else{  // ix == 0
      int nei = Nxv-1 + Nxv * iyzt;
      if(do_comm[idir] == 0){
        shift_vec2_fw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], bc2, NC);
      }else{
        shift_vec0_fw(vt, &v1[VLEN*NC*ist], NC);
      }
      shift_vec2_fw(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],
                        &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(v2v[ic], &ut[ic], vt, NC);
      }

    }

    sub_vec(&v2[VLEN*NC*ist], v2v, NC);

  }

}

//====================================================================
void mult_staggered_yp1(real_t *buf, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;
  int Nxzt = Nxv * Nzt;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxzt);

  real_t bc2 = real_t(bc[idir]);

  for(int ixzt = is; ixzt < ns; ++ixzt){
    int ix  = ixzt % Nxv;
    int izt = ixzt/Nxv;
    int iy  = 0;
    int ist = ix + Nxv * (iy + Ny * izt);

    for(int ic = 0; ic < NC; ++ic){
      for(int k = 0; k < VLEN; ++k){
        buf[k + VLEN*(ic + NC*ixzt)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }

  }

}

//====================================================================
void mult_staggered_yp2(real_t *v2, real_t *u, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nxv;
    int iy  = (ist/Nxv) % Ny;
    int izt = ist/(Nxv*Ny);
    int ixzt = ix + Nxv * izt;

    if(iy == Ny-1){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &buf[VLEN*NC*ixzt], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_ypb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nxv;
    int iy  = (ist/Nxv) % Ny;
    int izt = ist/(Nxv*Ny);
    int iy2 = (iy+1) % Ny;
    int nei = ix + Nxv * (iy2 + Ny * izt);

    if(iy < Ny-1 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(iy == Ny-1) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_ym1(real_t *buf, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;
  int Nxzt = Nxv * Nzt;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxzt);

  for(int ixzt = is; ixzt < ns; ++ixzt){
    int ix  = ixzt % Nxv;
    int izt = ixzt/Nxv;
    int iy  = Ny-1;
    int ist = ix + Nxv * (iy + Ny * izt);

    Vsimd_t v2v[NC], vt[NC], ut[NDF2];
    load_vec(vt, &v1[VLEN*NC*ist], NC);
    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

    for(int ic = 0; ic < NC; ++ic){
      mult_udagv(v2v[ic], &ut[ic], vt, NC);
    }

    save_vec(&buf[VLEN*NC*ixzt], v2v, NC);

  }

}

//====================================================================
void mult_staggered_ym2(real_t *v2, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  real_t bc2 = real_t(bc[idir]);

  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nxv;
    int iy  = (ist/Nxv) % Ny;
    int izt = ist/(Nxv*Ny);
    int ixzt = ix + Nxv * izt;

    if(iy == 0){
      Vsimd_t vt[NC];
      load_vec(vt, &buf[VLEN*NC*ixzt], NC);
      scal_vec(vt, bc2, NC);
      sub_vec(&v2[VLEN*NC*ist], vt, NC);
    }

  }

}

//====================================================================
void mult_staggered_ymb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 1;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nstv = Nxv * Ny * Nzt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nxv;
    int iy  = (ist/Nxv) % Ny;
    int izt = ist/(Nxv*Ny);
    int iy2 = (iy-1+Ny) % Ny;
    int nei = ix + Nxv * (iy2 + Ny * izt);

    if(iy > 0 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(iy == 0) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(v2v[ic], &ut[ic], vt, NC);
      }

      sub_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_zp1(real_t *buf, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyt = Nxv * Ny * Nt;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxyt);

  real_t bc2 = real_t(bc[idir]);

  int Nxy = Nxv * Ny;
  for(int ixzt = is; ixzt < ns; ++ixzt){
    int ixy = ixzt % Nxy;
    int it  = ixzt/Nxy;
    int iz   = 0;
    int ist  = ixy + Nxy * (iz + Nz * it);

    for(int ic = 0; ic < NC; ++ic){
      for(int k = 0; k < VLEN; ++k){
        buf[k + VLEN*(ic + NC*ixzt)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }

  }

}

//====================================================================
void mult_staggered_zp2(real_t *v2, real_t *u, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  int Nxy = Nxv * Ny;
  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int ixyt = ixy + Nxy * it;

    if(iz == Nz-1){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &buf[VLEN*NC*ixyt], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_zpb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  int Nxy = Nxv * Ny;
  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz+1) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);

    if(iz < Nz-1 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(iz == Nz-1) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_zm1(real_t *buf, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyt = Nxv * Ny * Nt;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxyt);

  int Nxy = Nxv * Ny;
  for(int ixzt = is; ixzt < ns; ++ixzt){
    int ixy = ixzt % Nxy;
    int it  = ixzt/Nxy;
    int iz   = Nz-1;
    int ist  = ixy + Nxy * (iz + Nz * it);

    Vsimd_t v2v[NC], vt[NC], ut[NDF2];
    load_vec(vt, &v1[VLEN*NC*ist], NC);
    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

    for(int ic = 0; ic < NC; ++ic){
      mult_udagv(v2v[ic], &ut[ic], vt, NC);
    }

    save_vec(&buf[VLEN*NC*ixzt], v2v, NC);

  }

}

//====================================================================
void mult_staggered_zm2(real_t *v2, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  real_t bc2 = real_t(bc[idir]);

  int Nxy = Nxv * Ny;
  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int ixyt = ixy + Nxy * it;

    if(iz == 0){
      Vsimd_t vt[NC];
      load_vec(vt, &buf[VLEN*NC*ixyt], NC);
      scal_vec(vt, bc2, NC);
      sub_vec(&v2[VLEN*NC*ist], vt, NC);
    }

  }

}

//==================================================================== 
void mult_staggered_zmb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 2;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  int Nxy = Nxv * Ny;
  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz-1+Nz) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);

    if(iz > 0 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(iz == 0) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(v2v[ic], &ut[ic], vt, NC);
      }

      sub_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_tp1(real_t *buf, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyz = Nxv * Ny * Nz;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxyz);

  int Nxy = Nxv * Ny;
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = 0;
    int ist  = ixyz + Nxyz * it;

    for(int ic = 0; ic < NC; ++ic){
      for(int k = 0; k < VLEN; ++k){
        buf[k + VLEN*(ic + NC*ixyz)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }

  }

}

//====================================================================
void mult_staggered_tp2(real_t *v2, real_t *u, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyz = Nxv * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;

    if(it == Nt-1){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &buf[VLEN*NC*ixyz], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_tpb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  int Nxyz = Nxv * Ny * Nz;
  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;

    if(it < Nt-1 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(it == Nt-1) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(v2v[ic], &ut[ic2], vt, NC);
      }

      add_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//====================================================================
void mult_staggered_tm1(real_t *buf, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyz = Nxv * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxyz);

  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    Vsimd_t v2v[NC], vt[NC], ut[NDF2];
    load_vec(vt, &v1[VLEN*NC*ist], NC);
    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

    for(int ic = 0; ic < NC; ++ic){
      mult_udagv(v2v[ic], &ut[ic], vt, NC);
    }

    save_vec(&buf[VLEN*NC*ixyz], v2v, NC);

  }

}

//====================================================================
void mult_staggered_tm2(real_t *v2, real_t *buf, 
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyz = Nxv * Ny * Nz;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;

    if(it == 0){
      Vsimd_t vt[NC];
      load_vec(vt, &buf[VLEN*NC*ixyz], NC);
      scal_vec(vt, bc2, NC);
      sub_vec(&v2[VLEN*NC*ist], vt, NC);
    }

  }

}

//==================================================================== 
void mult_staggered_tmb(real_t *v2, real_t *u, real_t *v1,
                        int *bc, int *Nsize, int *do_comm)
{
  int idir = 3;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyz = Nxv * Ny * Nz;

  real_t bc2 = real_t(bc[idir]);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nstv);

  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;

    if(it > 0 || do_comm[idir] == 0){

      Vsimd_t v2v[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(it == 0) scal_vec(vt, bc2, NC);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(v2v[ic], &ut[ic], vt, NC);
      }

      sub_vec(&v2[VLEN*NC*ist], v2v, NC);

    }

  }

}

//============================================================END=====

} // end of namespace QYS
