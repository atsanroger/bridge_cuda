/*!
      @file    afopr_Staggered_simd2R-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

// This code explicitly assumes SU(3) gauge group.
#include "lib_alt_SIMD2/inline/define_params_SU3.h"

namespace QYS{

//====================================================================
void mult_staggered_1(real_t *buf_xp, real_t *buf_xm,
                      real_t *buf_yp, real_t *buf_ym,
                      real_t *buf_zp, real_t *buf_zm,
                      real_t *buf_tp, real_t *buf_tm,
                      real_t *u, real_t *v1,
                      int *bc, int *Nsize, int *do_comm)
{

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
      buf_xp[ivc + NVC*iyzt] = bc2 * vt[ivc];
    }
  }

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
      buf_xm[2*ic   + NVC*iyzt] = wt2[VLEN-2];
      buf_xm[2*ic+1 + NVC*iyzt] = wt2[VLEN-1];
    }
  }

 }

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
        buf_yp[k + VLEN*(ic + NC*ixzt)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }
  }

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
    save_vec(&buf_ym[VLEN*NC*ixzt], v2v, NC);
  }
 }

 {
  int idir = 2;

  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;
  int Nxyt = Nxv * Ny * Nt;
  int Nxy  = Nxv * Ny;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, Nxyt);

  real_t bc2 = real_t(bc[idir]);

  for(int ixzt = is; ixzt < ns; ++ixzt){
    int ixy = ixzt % Nxy;
    int it  = ixzt/Nxy;
    int iz   = 0;
    int ist  = ixy + Nxy * (iz + Nz * it);

    for(int ic = 0; ic < NC; ++ic){
      for(int k = 0; k < VLEN; ++k){
        buf_zp[k + VLEN*(ic + NC*ixzt)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }
  }

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

    save_vec(&buf_zm[VLEN*NC*ixzt], v2v, NC);
  }

 }

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
        buf_tp[k + VLEN*(ic + NC*ixyz)]
                          = bc2 * v1[k + VLEN*(ic + NC*ist)];
      }
    }
  }

  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    Vsimd_t v2v[NC], vt[NC], ut[NDF2];
    load_vec(vt, &v1[VLEN*NC*ist], NC);
    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

    for(int ic = 0; ic < NC; ++ic){
      mult_udagv(v2v[ic], &ut[ic], vt, NC);
    }
    save_vec(&buf_tm[VLEN*NC*ixyz], v2v, NC);
  }
 }

}

//====================================================================
void mult_staggered_2(real_t *v2, real_t *u,
                      real_t *buf_xp, real_t *buf_xm,
                      real_t *buf_yp, real_t *buf_ym,
                      real_t *buf_zp, real_t *buf_zm,
                      real_t *buf_tp, real_t *buf_tm,
                      int *bc, int *Nsize, int* do_comm,
                      real_t mq, int jd)
{
  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

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
    int iy   = iyzt % Ny;
    int izt  = iyzt/Ny;
    int iz   = izt % Nz;
    int it   = izt/Nz;
    int ixzt = ix + Nxv * izt;
    int ixyz = ix + Nxv * (iy + Ny * iz);
    int ixyt = ix + Nxv * (iy + Ny * it);

    Vsimd_t v2v[NC];
    clear_vec(v2v, NC);

    int idir = 0;
    if(do_comm[idir] == 1){
      real_t bc2 = real_t(bc[idir]);

      if(ix == Nxv-1){
        Vsimd_t vt2[NC], vt[NC], ut[NDF2];
        shift_vec1_bw(vt, v0, &buf_xp[NVC*iyzt], NC);
        load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],  NDF2);

        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          mult_uv(vt2[ic], &ut[ic2], vt, NC);
        }
        add_vec(v2v, vt2, NC);
      }

      if(ix == 0){
        Vsimd_t vt[NC];
        shift_vec1_fw(vt, v0, &buf_xm[NVC*iyzt], NC);
        axpy_vec(v2v, -bc2, vt, NC);
      }

    }

    idir = 1;
    if(do_comm[idir] == 1){
      real_t bc2 = real_t(bc[idir]);

      if(iy == Ny-1){
        Vsimd_t vt2[NC], vt[NC], ut[NDF2];
        load_vec(vt, &buf_yp[VLEN*NC*ixzt], NC);
        load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);
        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          mult_uv(vt2[ic], &ut[ic2], vt, NC);
        }
        add_vec(v2v, vt2, NC);
      }

      if(iy == 0){
        Vsimd_t vt[NC];
        load_vec(vt, &buf_ym[VLEN*NC*ixzt], NC);
        axpy_vec(v2v, -bc2, vt, NC);
      }

    }

    idir = 2;
    if(do_comm[idir] == 1){
      real_t bc2 = real_t(bc[idir]);

      if(iz == Nz-1){
        Vsimd_t vt2[NC], vt[NC], ut[NDF2];
        load_vec(vt, &buf_zp[VLEN*NC*ixyt], NC);
        load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);
        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          mult_uv(vt2[ic], &ut[ic2], vt, NC);
        }
        add_vec(v2v, vt2, NC);
      }

      if(iz == 0){
        Vsimd_t vt[NC];
        load_vec(vt, &buf_zm[VLEN*NC*ixyt], NC);
        axpy_vec(v2v, -bc2, vt, NC);
      }

    }

    idir = 3;
    if(do_comm[idir] == 1){
      real_t bc2 = real_t(bc[idir]);

      if(it == Nt-1){
        Vsimd_t vt2[NC], vt[NC], ut[NDF2];
        load_vec(vt, &buf_tp[VLEN*NC*ixyz], NC);
        load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);
        for(int ic = 0; ic < NC; ++ic){
          int ic2 = NC * ic;
          mult_uv(vt2[ic], &ut[ic2], vt, NC);
        }
        add_vec(v2v, vt2, NC);
      }

      if(it == 0){
        Vsimd_t vt[NC];
        load_vec(vt, &buf_tm[VLEN*NC*ixyz], NC);
        axpy_vec(v2v, -bc2, vt, NC);
      }

    }

    Vsimd_t v[NC];
    load_vec(v, &v2[VLEN*NC*ist], NC);
    // real_t fac = real_t(jd) * 0.5/mq;  // hopping normalization
    real_t fac = real_t(jd) * 0.5;      // mass normalization
    axpy_vec(v, fac, v2v, NC);

    save_vec(&v2[VLEN*NC*ist], v, NC);

  }

}

//====================================================================
void mult_staggered_bulk(real_t *v2, real_t *u, real_t *v1,
                         int *bc, int *Nsize, int* do_comm,
                         real_t mq, int jd)
{
  int Nxv  = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nstv = Nxv * Ny * Nz * Nt;
  int Nst  = Nstv * VLEN2;

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
    int iy   = iyzt % Ny;
    int izt  = iyzt/Ny;
    int iz   = izt % Nz;
    int it   = izt/Nz;

    Vsimd_t v2v[NC];
    clear_vec(v2v, NC);

    Vsimd_t vt2[NC], vt[NC], ut[NDF2];

    int idir = 0;

    if(ix < Nxv-1){
      int nei = ix+1 + Nxv * iyzt;
      shift_vec2_bw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], NC);
    }else{
      if(do_comm[idir] == 0){
        int nei = 0 + Nxv * iyzt;
        real_t bc2 = real_t(bc[idir]);
        shift_vec2_bw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], bc2, NC);
      }else{
        shift_vec0_bw(vt, &v1[VLEN*NC*ist], NC);
      }
    }

    load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],  NDF2);

    for(int ic = 0; ic < NC; ++ic){
      int ic2 = NC * ic;
      mult_uv(vt2[ic], &ut[ic2], vt, NC);
    }

    add_vec(v2v, vt2, NC);

    if(ix > 0){
      int nei = ix-1 + Nxv * iyzt;
      real_t bc2 = real_t(bc[idir]);
      shift_vec2_fw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], NC);
      shift_vec2_fw(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],
                        &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(vt2[ic], &ut[ic], vt, NC);
      }
    }else{  // ix == 0
      int nei = Nxv-1 + Nxv * iyzt;
      real_t bc2 = real_t(bc[idir]);
      if(do_comm[idir] == 0){
        shift_vec2_fw(vt, &v1[VLEN*NC*ist], &v1[VLEN*NC*nei], bc2, NC);
      }else{
        shift_vec0_fw(vt, &v1[VLEN*NC*ist], NC);
      }
      shift_vec2_fw(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir],
                        &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(vt2[ic], &ut[ic], vt, NC);
      }
    }

    sub_vec(v2v, vt2, NC);

    idir = 1;

    if(iy < Ny-1 || do_comm[idir] == 0){
      int iy2 = (iy+1) % Ny;
      int nei = ix + Nxv * (iy2 + Ny * izt);

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(iy == Ny-1){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(vt2[ic], &ut[ic2], vt, NC);
      }

      add_vec(v2v, vt2, NC);

    }

    if(iy > 0 || do_comm[idir] == 0){
      int iy2 = (iy-1+Ny) % Ny;
      int nei = ix + Nxv * (iy2 + Ny * izt);

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(iy == 0){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(vt2[ic], &ut[ic], vt, NC);
      }

      sub_vec(v2v, vt2, NC);
    }

    idir = 2;

    if(iz < Nz-1 || do_comm[idir] == 0){
      int iz2 = (iz+1) % Nz;
      int nei = ix + Nxv * (iy + Ny * (iz2 + Nz * it));
      real_t bc2 = real_t(bc[idir]);

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(iz == Nz-1){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(vt2[ic], &ut[ic2], vt, NC);
      }

      add_vec(v2v, vt2, NC);
    }

    if(iz > 0 || do_comm[idir] == 0){
      int iz2 = (iz-1+Nz) % Nz;
      int nei = ix + Nxv * (iy + Ny * (iz2 + Nz * it));

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(iz == 0){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(vt2[ic], &ut[ic], vt, NC);
      }

      sub_vec(v2v, vt2, NC);

    }

    idir = 3;

    if(it < Nt-1 || do_comm[idir] == 0){
      int it2 = (it+1) % Nt;
      int nei = ix + Nxv * (iy + Ny * (iz + Nz * it2));

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*ist + NDF*Nst*idir], NDF2);

      if(it == Nt-1){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        int ic2 = NC * ic;
        mult_uv(vt2[ic], &ut[ic2], vt, NC);
      }

      add_vec(v2v, vt2, NC);

    }

    if(it > 0 || do_comm[idir] == 0){
      int it2 = (it-1+Nt) % Nt;
      int nei = ix + Nxv * (iy + Ny * (iz + Nz * it2));
      real_t bc2 = real_t(bc[idir]);

      Vsimd_t vt2[NC], vt[NC], ut[NDF2];
      load_vec(vt, &v1[VLEN*NC*nei], NC);
      load_vec(ut, &u[VLEN*NDF2*nei + NDF*Nst*idir], NDF2);

      if(it == 0){
        real_t bc2 = real_t(bc[idir]);
        scal_vec(vt, bc2, NC);
      }

      for(int ic = 0; ic < NC; ++ic){
        mult_udagv(vt2[ic], &ut[ic], vt, NC);
      }

      sub_vec(v2v, vt2, NC);

    }

    Vsimd_t w[NC];
    load_vec(w, &v1[VLEN*NC*ist], NC);

    // real_t fac = real_t(jd) * 0.5/mq;
    // aypx_vec(fac, v2v, w, NC);

    real_t fac = real_t(jd) * 0.5;
    scal_vec(v2v, fac, NC);
    axpy_vec(v2v, mq, w, NC);

    save_vec(&v2[VLEN*NC*ist], v2v, NC);

  }

}

//====================================================================

} // end of namespace QYS
