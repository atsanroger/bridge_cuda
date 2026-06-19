/*!
      @file    afopr_Staggered_vector-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

//====================================================================
void mult_staggered_gm5(real_t *RESTRICT v2, real_t *RESTRICT v1,
                           int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  // this is incorrect!
  for(int ist = 0; ist < Nst; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] = v1[IDXV_1SP_R(Nst, ic, ist)];
      v2[IDXV_1SP_I(Nst, ic, ist)] = v1[IDXV_1SP_I(Nst, ic, ist)];
    }
  }

}

//====================================================================
void mult_staggered_phase_dev(real_t *RESTRICT u, real_t *RESTRICT ph,
                              int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  for(int mu = 0; mu < 4; ++mu){
    for(int ist = 0; ist < Nst; ++ist){
      real_t ph1 = ph[IDXV(1, Nst, 0, ist, mu)];
#pragma _NEC unroll(NDF)
      for(int idf = 0; idf < NDF; ++idf){
        u[IDXV(NDF, Nst, idf, ist, mu)] *= ph1;
      }

    }
  }

}

//====================================================================
void mult_staggered_aypx(real_t a, real_t *RESTRICT v2, real_t *RESTRICT v1,
                         int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  for(int ist = 0; ist < Nst; ++ist){
#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] = a * v2[IDXV_1SP(Nst, ivc, ist)]
                                      + v1[IDXV_1SP(Nst, ivc, ist)];
    }
  }

}

//====================================================================
void mult_staggered_axpy(real_t *RESTRICT v2, real_t a, real_t *RESTRICT v1,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  for(int ist = 0; ist < Nst; ++ist){
#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] += a * v1[IDXV_1SP(Nst, ivc, ist)];
    }
  }

}

//====================================================================
void mult_staggered_scal(real_t *RESTRICT v, real_t a,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  for(int ist = 0; ist < Nst; ++ist){
#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      v[IDXV_1SP(Nst, ivc, ist)] *= a;
    }
  }

}

//====================================================================
void mult_staggered_clear(real_t *RESTRICT v2, int *Nsize, int Nc)
{
  int Nst = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  for(int ist = 0; ist < Nst; ++ist){
#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] = 0.0;
    }
  }

}

//====================================================================
void mult_staggered_xp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst = Nx * Nyzt;

  real_t bc2 = bc[0];

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      buf[IDXV_1SP(Nyzt, ivc, iyzt)] = bc2 * v1[IDXV_1SP(Nst, ivc, ist)];
    }

  }

}

//====================================================================
void mult_staggered_xp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                     real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC], ut[NDF];

#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = buf[IDXV_1SP(Nyzt, ivc, iyzt)];
    }

#pragma _NEC unroll(NC)
    for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
 void mult_staggered_xpb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ((ix+1) % Nx) + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == Nx-1) bc2 = bc[0];

    real_t vt[NVC], ut[NVC];

#pragma _NEC unroll(NVC)
    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = bc2 * v1[IDXV_1SP(Nst, ivc, nei)];
    }

#pragma _NEC unroll(NC)
    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
void mult_staggered_xm1(real_t *RESTRICT buf, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      buf[IDXV_1SP_R(Nyzt, ic, iyzt)]
               = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      buf[IDXV_1SP_I(Nyzt, ic, iyzt)]
               = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }

}

//====================================================================
void mult_staggered_xm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  real_t bc2 = bc[0];

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] += -bc2 * buf[IDXV_1SP(Nyzt, ivc, iyzt)];
    }

  }

}

//====================================================================
void mult_staggered_xmb(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ix-1 + Nx * iyzt;
    if(ix == 0) nei = Nx-1 + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == 0) bc2 = bc[0];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] += -bc2 * wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] += -bc2 * wti;

    }

  }

}

//====================================================================
void mult_staggered_yp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nxzt = Nx * Nzt;
  int Nst  = Nxzt * Ny;

  real_t bc2 = bc[1];

  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    int iy   = 0;
    int ix   = ixzt % Nx;
    int izt  = ixzt/Nx;
    int ist  = ix + Nx * (iy + Ny * izt);

    for(int ivc = 0; ivc < NVC; ++ivc){
      buf[IDXV_1SP(Nxzt, ivc, ixzt)] = bc2 * v1[IDXV_1SP(Nst, ivc, ist)];
    }

  }

}

//====================================================================
void mult_staggered_yp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nxzt = Nx * Nzt;
  int Nst  = Nxzt * Ny;

  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    int iy   = Ny-1;
    int ix = ixzt % Nx;
    int izt = ixzt/Nx;
    int ist  = ix + Nx * (iy + Ny * izt);

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = buf[IDXV_1SP(Nxzt, ivc, ixzt)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
void mult_staggered_ypb(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy+1) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);
    real_t bc2 = 1.0;
    if(iy == Ny-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = bc2 * v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
void mult_staggered_ym1(real_t *RESTRICT buf, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nxzt = Nx * Nzt;
  int Nst  = Nxzt * Ny;

  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    int iy  = Ny-1;
    int ix  = ixzt % Nx;
    int izt = ixzt/Nx;
    int ist = ix + Nx * (iy + Ny*izt);

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      buf[IDXV_1SP_R(Nxzt, ic, ixzt)]
               = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      buf[IDXV_1SP_I(Nxzt, ic, ixzt)]
               = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }

}

//====================================================================
void mult_staggered_ym2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nxzt = Nx * Nzt;
  int Nst  = Nxzt * Ny;

  real_t bc2 = bc[idir];

  for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
    int iy = 0;
    int ix = ixzt % Nx;
    int izt = ixzt/Nx;
    int ist  = ix + Nx*(iy + Ny*izt);

     for(int ivc = 0; ivc < NVC; ++ivc){
       v2[IDXV_1SP(Nst, ivc, ist)] += -bc2 * buf[IDXV_1SP(Nxzt, ivc, ixzt)];
     }

  }

}

//====================================================================
void mult_staggered_ymb(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst  = Nx * Ny * Nzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy-1+Ny) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);

    real_t bc2 = 1.0;
    if(iy == 0) bc2 = bc[1];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  -bc2 * wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  -bc2 * wti;

    }

  }

}

//====================================================================
void mult_staggered_zp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nxyt = Nxy * Nt;
  int Nst  = Nxyt * Nz;

  real_t bc2 = bc[idir];

  for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    int iz  = 0;
    int ixy = ixyt % Nxy;
    int it  = ixyt/Nxy;
    int ist = ixy + Nxy * (iz + Nz * it);

    for(int ivc = 0; ivc < NVC; ++ivc){
      buf[IDXV_1SP(Nxyt, ivc, ixyt)]
                           = bc2 * v1[IDXV_1SP(Nst, ivc, ist)];
    }

  }

}

//====================================================================
void mult_staggered_zp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT buf, 
                        int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nxyt = Nxy * Nt;
  int Nst  = Nxyt * Nz;

  for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    int iz  = Nz-1;
    int ixy = ixyt % Nxy;
    int it  = ixyt/Nxy;
    int ist = ixy + Nxy * (iz + Nz * it);

     real_t vt[NVC], ut[NDF];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt[ivc] = buf[IDXV_1SP(Nxyt, ivc, ixyt)];
     }

     for(int ic = 0; ic < NC; ++ic){

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
         ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
       }

       real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
       real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

       v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
       v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

     }

  }

}

//====================================================================
void mult_staggered_zpb(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz+1) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == Nz-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = bc2 * v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
void mult_staggered_zm1(real_t *RESTRICT buf, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nxyt = Nxy * Nt;
  int Nst  = Nxyt * Nz;

  for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    int iz = Nz-1;
    int ixy = ixyt % Nxy;
    int it  = ixyt/Nxy;
    int ist  = ixy + Nxy * (iz + Nz * it);
 
    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      buf[IDXV_1SP_R(Nxyt, ic, ixyt)]
               = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      buf[IDXV_1SP_I(Nxyt, ic, ixyt)]
               = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }

}

//====================================================================
void mult_staggered_zm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nxyt = Nxy * Nt;
  int Nst  = Nxyt * Nz;

  real_t bc2 = bc[idir];

  for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
    int iz  = 0;
    int ixy = ixyt % Nxy;
    int it  = ixyt/Nxy;
    int ist = ixy + Nxy * (iz + Nz * it);

    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] += -bc2 * buf[IDXV_1SP(Nxyt, ivc, ixyt)];
     }

  }

}

//==================================================================== 
void mult_staggered_zmb(real_t *RESTRICT v2, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz-1+Nz) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == 0) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  -bc2 * wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  -bc2 * wti;

    }

  }

}

//====================================================================
void mult_staggered_tp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  real_t bc2 = bc[idir];

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;
    real_t bc2 = bc[3];

    real_t *wt = buf;
    for(int ic = 0; ic < NC; ++ic){
      wt[IDXV_1SP_R(Nxyz, ic, ixyz)] = bc2 * v1[IDXV_1SP_R(Nst, ic, ist)];
      wt[IDXV_1SP_I(Nxyz, ic, ixyz)] = bc2 * v1[IDXV_1SP_I(Nst, ic, ist)];
    }

  }

  /*
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;

    for(int ivc = 0; ivc < NVC; ++ivc){
      buf[IDXV_1SP(Nxyz, ivc, ixyz)] = bc2 * v1[IDXV_1SP(Nst, ivc, ist)];
    }

  }
  */
}

//====================================================================
void mult_staggered_tp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                           real_t *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = Nt-1;
    int ist = ixyz + Nxyz * it;

    real_t *v1 = buf;
    int Nbf = Nxyz;
    int ibf = ixyz;

    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_1SP_R(Nst, ic, ist)] +=
         u[IDXV_G_R(Nst, 0, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
       - u[IDXV_G_I(Nst, 0, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
       + u[IDXV_G_R(Nst, 1, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
       - u[IDXV_G_I(Nst, 1, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
       + u[IDXV_G_R(Nst, 2, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 2, ibf)]
       - u[IDXV_G_I(Nst, 2, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 2, ibf)];
      v2[IDXV_1SP_I(Nst, ic, ist)] +=
         u[IDXV_G_R(Nst, 0, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 0, ibf)]
       + u[IDXV_G_I(Nst, 0, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 0, ibf)]
       + u[IDXV_G_R(Nst, 1, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 1, ibf)]
       + u[IDXV_G_I(Nst, 1, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 1, ibf)]
       + u[IDXV_G_R(Nst, 2, ic, ist, 0)] * v1[IDXV_1SP_I(Nbf, 2, ibf)]
       + u[IDXV_G_I(Nst, 2, ic, ist, 0)] * v1[IDXV_1SP_R(Nbf, 2, ibf)];
    }

  }

  /*
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = buf[IDXV_1SP(Nxyz, ivc, ixyz)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }
  */

}

//====================================================================
 void mult_staggered_tpb(real_t *RESTRICT v2, real_t *RESTRICT u,
                            real_t *RESTRICT v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == Nt-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = bc2 * v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  wti;

    }

  }

}

//====================================================================
void mult_staggered_tm1(real_t *RESTRICT buf, real_t *RESTRICT u,
                        real_t *RESTRICT v1,
                        int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nxyz * Nt;

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t *wt = buf;
    int ibf = ixyz;
    int Nbf = Nxyz;

    for(int ic = 0; ic < NC; ++ic){
      wt[IDXV_1SP_R(Nbf, ic, ibf)] =
          u[IDXV_G_R(Nst, ic, 0, ist, 0)] * v1[IDXV_1SP_R(Nst, 0, ist)]
        + u[IDXV_G_I(Nst, ic, 0, ist, 0)] * v1[IDXV_1SP_I(Nst, 0, ist)]
        + u[IDXV_G_R(Nst, ic, 1, ist, 0)] * v1[IDXV_1SP_R(Nst, 1, ist)]
        + u[IDXV_G_I(Nst, ic, 1, ist, 0)] * v1[IDXV_1SP_I(Nst, 1, ist)]
        + u[IDXV_G_R(Nst, ic, 2, ist, 0)] * v1[IDXV_1SP_R(Nst, 2, ist)]
        + u[IDXV_G_I(Nst, ic, 2, ist, 0)] * v1[IDXV_1SP_I(Nst, 2, ist)];
      wt[IDXV_1SP_I(Nbf, ic, ibf)] =
	  u[IDXV_G_R(Nst, ic, 0, ist, 0)] * v1[IDXV_1SP_I(Nst, 0, ist)]
        - u[IDXV_G_I(Nst, ic, 0, ist, 0)] * v1[IDXV_1SP_R(Nst, 0, ist)]
        + u[IDXV_G_R(Nst, ic, 1, ist, 0)] * v1[IDXV_1SP_I(Nst, 1, ist)]
        - u[IDXV_G_I(Nst, ic, 1, ist, 0)] * v1[IDXV_1SP_R(Nst, 1, ist)]
        + u[IDXV_G_R(Nst, ic, 2, ist, 0)] * v1[IDXV_1SP_I(Nst, 2, ist)]
        - u[IDXV_G_I(Nst, ic, 2, ist, 0)] * v1[IDXV_1SP_R(Nst, 2, ist)];
    }

  }

  /*
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = Nt-1;
    int ist = ixyz + Nxyz * it;

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, ist)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      buf[IDXV_1SP_R(Nxyz, ic, ixyz)] 
               = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      buf[IDXV_1SP_I(Nxyz, ic, ixyz)] 
               = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }
  */
}

//====================================================================
void mult_staggered_tm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  real_t bc2 = bc[idir];

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;

    for(int ivc = 0; ivc < NVC; ++ ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] -= bc2 * buf[IDXV_1SP(Nxyz, ivc, ixyz)];
    }
  }

  /*
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;

    for(int ivc = 0; ivc < NVC; ++ivc){
      v2[IDXV_1SP(Nst, ivc, ist)] += -bc2 * buf[IDXV_1SP(Nxyz, ivc, ixyz)];
    }

  }
  */
}

//==================================================================== 
void mult_staggered_tmb(real_t *RESTRICT v2, real_t *RESTRICT u,
                           real_t *RESTRICT v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == 0) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = v1[IDXV_1SP(Nst, ivc, nei)];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDXV_1SP_R(Nst, ic, ist)] +=  -bc2 * wtr;
      v2[IDXV_1SP_I(Nst, ic, ist)] +=  -bc2 * wti;

    }

  }

}

//====================================================================

