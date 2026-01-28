/*!
      @file    mult_Staggered_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2542 $
*/


#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

//====================================================================
void mult_staggered_gm5(real_t *restrict v2, real_t *restrict v1,
                           int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size = NVC * Nst_pad;

#pragma acc data present(v2[0:size], v1[0:size]) copyin(Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst_pad; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] = v1[IDX2_1SP_R(ic, ist)];
      v2[IDX2_1SP_I(ic, ist)] = v1[IDX2_1SP_I(ic, ist)];
    }
  }

 }

}

//====================================================================
void mult_staggered_phase_dev(real_t *restrict u, real_t *restrict ph,
                              int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size_u  = NDF * Nst_pad * 4;
  int size_ph = Nst_pad * 4;

#pragma acc data present(u[0:size_u], ph[0:size_ph]) copyin(Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int mu = 0; mu < 4; ++mu){
    for(int ist = 0; ist < Nst_pad; ++ist){
      int istex = ist + Nst_pad * mu;
      real_t ph1 = 0.0;
      if(ist < Nst) ph1 = ph[IDX2(1, 0, istex)];
      for(int idf = 0; idf < NDF; ++idf){
        u[IDX2(NDF, idf, istex)] *= ph1;
      }
    }
  }

 }

}

//====================================================================
void mult_staggered_aypx(real_t a, real_t *restrict v2, real_t *restrict v1,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size = NVC * Nst_pad;

#pragma acc data present(v2[0:size], v1[0:size]) copyin(a, Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] = a * v2[IDX2_1SP_R(ic, ist)]
                                  + v1[IDX2_1SP_R(ic, ist)];
      v2[IDX2_1SP_I(ic, ist)] = a * v2[IDX2_1SP_I(ic, ist)]
                                  + v1[IDX2_1SP_I(ic, ist)];
    }
  }

  }

}

//====================================================================
void mult_staggered_axpy(real_t *restrict v2, real_t a, real_t *restrict v1,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size = NVC * Nst_pad;

#pragma acc data present(v2[0:size], v1[0:size]) copyin(a, Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] += a * v1[IDX2_1SP_R(ic, ist)];
      v2[IDX2_1SP_I(ic, ist)] += a * v1[IDX2_1SP_I(ic, ist)];
    }
  }

 }

}

//====================================================================
void mult_staggered_scal(real_t *restrict v, real_t a,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size = NVC * Nst_pad;

#pragma acc data present(v[0:size]) copyin(a, Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v[IDX2_1SP_R(ic, ist)] *= a;
      v[IDX2_1SP_I(ic, ist)] *= a;
    }
  }

  }

}

//====================================================================
void mult_staggered_clear(real_t *restrict v2, int *Nsize, int Nc)
{
  int Nst = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int Nst_pad = CEIL_NWP(Nst);
  int size = NVC * Nst_pad;

#pragma acc data present(v2[0:size]), copyin(Nst, Nst_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] = 0.0;
      v2[IDX2_1SP_I(ic, ist)] = 0.0;
    }
  }

 }

}

//====================================================================
void mult_staggered_xp1(real_t *restrict buf, real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];

  real_t bc2 = bc[0];

  int size   = NVC * Nx * Nyzt;
  int size_b = NVC * Nyzt;

#pragma acc data present(v1[0:size], buf[0:size_b]), copyin(bc2, Nx, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
    }

    real_t *vt1 = &buf[NVC*iyzt];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = bc2 * vt[icr];
      vt1[ici] = bc2 * vt[ici];
    }

  }

 }

}

//====================================================================
void mult_staggered_xp2(real_t *restrict v2, real_t *restrict u,
                     real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int size   = NVC * Nst;
  int size_b = NVC * Nyzt;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(Nx, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = buf[ivc + NVC*iyzt];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
 void mult_staggered_xpb(real_t *restrict v2, real_t *restrict u,
                      real_t *restrict v1,
                      int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], Nx, Nyzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ((ix+1) % Nx) + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == Nx-1) bc2 = bc[0];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = bc2 * v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = bc2 * v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_xm1(real_t *restrict buf, real_t *restrict u,
                     real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int size   = NVC * Nst;
  int size_b = NVC * Nyzt;
  int size_u = NDF * Nst;

#pragma acc data present(v1[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(Nx, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
    }

    real_t *wt1 = &buf[NVC*iyzt];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }

 }

}

//====================================================================
void mult_staggered_xm2(real_t *restrict v2, real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  real_t bc2 = bc[0];

  int size   = NVC * Nst;
  int size_b = NVC * Nyzt;

#pragma acc data present(v2[0:size], buf[0:size_b]), copyin(bc2, Nx, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] += -bc2 * buf[2*ic   + NVC * iyzt];
      v2[IDX2_1SP_I(ic, ist)] += -bc2 * buf[2*ic+1 + NVC * iyzt];
    }

  }

 }

}

//====================================================================
 void mult_staggered_xmb(real_t *restrict v2, real_t *restrict u,
                      real_t *restrict v1,
                      int *Nsize, int *bc, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], Nx, Nyzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ix-1 + Nx * iyzt;
    if(ix == 0) nei = Nx-1 + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == 0) bc2 = bc[0];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  -bc2 * wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  -bc2 * wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_yp1(real_t *restrict buf, real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];

  real_t bc2 = bc[1];

  int size   = NVC * Nx * Ny * Nzt;
  int size_b = NVC * Nx * Nzt;

#pragma acc data present(v1[0:size], buf[0:size_b]), \
                 copyin(bc2, Nx, Ny, Nzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;

     real_t vt[NVC];

     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
       vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
     }

     real_t *vt1 = &buf[NVC*ixzt];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt1[ivc] = bc2 * vt[ivc];
     }

   }
  }

 }

}

//====================================================================
void mult_staggered_yp2(real_t *restrict v2, real_t *restrict u,
                     real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int size   = NVC * Nst;
  int size_b = NVC * Nx * Nzt;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(Nx, Ny, Nzt, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;

    real_t vt[NVC], ut[NDF];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ivc] = buf[ivc + NVC*ixzt];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

   }
  }

 }

}

//====================================================================
 void mult_staggered_ypb(real_t *restrict v2, real_t *restrict u,
                      real_t *restrict v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], idir, Nx, Ny, Nzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy+1) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);
    real_t bc2 = 1.0;
    if(iy == Ny-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = bc2 * v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = bc2 * v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_ym1(real_t *restrict buf, real_t *restrict u,
                     real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int size   = NVC * Nst;
  int size_b = NVC * Nx * Nzt;
  int size_u = NDF * Nst;

#pragma acc data present(v1[0:size], buf[0:size_b], u[0:size_u]), \
            copyin(Nx, Ny, Nzt, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny*izt);
     int ixzt = ix + Nx * izt;

    real_t vt[NVC], ut[NDF];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
    }

    real_t *wt1 = &buf[NVC*ixzt];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

   }
  }

 }

}

//====================================================================
void mult_staggered_ym2(real_t *restrict v2, real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  real_t bc2 = bc[idir];

  int size   = NVC * Nst;
  int size_b = NVC * Nx * Nzt;

#pragma acc data present(v2[0:size], buf[0:size_b]), \
                 copyin(bc2, Nx, Ny, Nzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy = 0;
     int ist  = ix + Nx*(iy + Ny*izt);
     int ixzt = ix + Nx * izt;

     for(int ic = 0; ic < NC; ++ic){
       v2[IDX2_1SP_R(ic, ist)] += -bc2 * buf[2*ic   + NVC * ixzt];
       v2[IDX2_1SP_I(ic, ist)] += -bc2 * buf[2*ic+1 + NVC * ixzt];
     }

   }
  }

 }

}

//====================================================================
 void mult_staggered_ymb(real_t *restrict v2, real_t *restrict u,
                      real_t *restrict v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst  = Nx * Ny * Nzt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], idir, Nx, Ny, Nzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy-1+Ny) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);

    real_t bc2 = 1.0;
    if(iy == 0) bc2 = bc[1];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  -bc2 * wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  -bc2 * wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_zp1(real_t *restrict buf, real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  real_t bc2 = bc[idir];

  int size   = NVC * Nxy * Nz * Nt;
  int size_b = NVC * Nxy * Nt;

#pragma acc data present(v1[0:size], buf[0:size_b]), \
                 copyin(bc2, Nxy, Nz, Nt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt[NVC];

     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
       vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
     }

     real_t *vt1 = &buf[NVC*ixyt];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt1[ivc] = bc2 * vt[ivc];
     }

   }
  }

 }

}

//====================================================================
void mult_staggered_zp2(real_t *restrict v2, real_t *restrict u,
                     real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int size   = NVC * Nst;
  int size_b = NVC * Nxy * Nt;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], buf[0:size_b], u[0:size_u]), \
            copyin(Nxy, Nz, Nt, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt[NVC], ut[NDF];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt[ivc] = buf[ivc + NVC*ixyt];
     }

     for(int ic = 0; ic < NC; ++ic){

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
         ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
       }

       real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
       real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

       v2[IDX2_1SP_R(ic, ist)] +=  wtr;
       v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

   }
  }

 }

}

//====================================================================
 void mult_staggered_zpb(real_t *restrict v2, real_t *restrict u,
                      real_t *restrict v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], idir, Nxy, Nz, Nt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz+1) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == Nz-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = bc2 * v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = bc2 * v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_zm1(real_t *restrict buf, real_t *restrict u,
                     real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int size   = NVC * Nst;
  int size_b = NVC * Nxy * Nt;
  int size_u = NDF * Nst;

#pragma acc data present(v1[0:size], buf[0:size_b], u[0:size_u]), \
            copyin(Nxy, Nz, Nt, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt[NVC], ut[NVC];

     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
       vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
     }

     real_t *wt1 = &buf[NVC*ixyt];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
         ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
       }

       wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
       wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

   }
  }

 }

}

//====================================================================
void mult_staggered_zm2(real_t *restrict v2, real_t *restrict buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  real_t bc2 = bc[idir];

  int size   = NVC * Nst;
  int size_b = NVC * Nxy * Nt;

#pragma acc data present(v2[0:size], buf[0:size_b]), \
  copyin(bc2, Nxy, Nz, Nt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     for(int ic = 0; ic < NC; ++ic){
       v2[IDX2_1SP_R(ic, ist)] += -bc2 * buf[2*ic   + NVC * ixyt];
       v2[IDX2_1SP_I(ic, ist)] += -bc2 * buf[2*ic+1 + NVC * ixyt];
     }

   }
  }

 }

}

//==================================================================== 
void mult_staggered_zmb(real_t *restrict v2, real_t *restrict u,
                     real_t *restrict v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                         copyin(bc[0:4], idir, Nxy, Nz, Nt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz-1+Nz) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == 0) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  -bc2 * wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  -bc2 * wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_tp1(real_t *restrict buf, real_t *restrict v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];

  real_t bc2 = bc[idir];

  int size   = NVC * Nxyz * Nt;
  int size_b = NVC * Nxyz;

#pragma acc data present(v1[0:size], buf[0:size_b]), copyin(bc2, Nxyz)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;

     real_t vt[NVC];

     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
       vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
     }

     real_t *vt1 = &buf[NVC*ixyz];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt1[ivc] = bc2 * vt[ivc];
     }

  }

 }

}

//====================================================================
void mult_staggered_tp2(real_t *restrict v2, real_t *restrict u,
                           real_t *restrict buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int size   = NVC * Nst;
  int size_b = NVC * Nxyz;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(Nxyz, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

     real_t vt[NVC], ut[NDF];

     for(int ivc = 0; ivc < NVC; ++ivc){
       vt[ivc] = buf[ivc + NVC*ixyz];
     }

     for(int ic = 0; ic < NC; ++ic){

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
         ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
       }

       real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
       real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                              vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

       v2[IDX2_1SP_R(ic, ist)] +=  wtr;
       v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
 void mult_staggered_tpb(real_t *restrict v2, real_t *restrict u,
                            real_t *restrict v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                 copyin(bc[0:4], idir, Nxyz, Nt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == Nt-1) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = bc2 * v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = bc2 * v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  wti;

    }

  }

 }

}

//====================================================================
void mult_staggered_tm1(real_t *restrict buf, real_t *restrict u,
                           real_t *restrict v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int size   = NVC * Nst;
  int size_b = NVC * Nxyz;
  int size_u = NDF * Nst;

#pragma acc data present(v1[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(Nxyz, Nt, Nst, idir)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt[NVC], ut[NDF];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
    }

    real_t *wt1 = &buf[NVC*ixyz];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

    }

  }

 }

}

//====================================================================
void mult_staggered_tm2(real_t *restrict v2, real_t *restrict buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  real_t bc2 = bc[idir];

  int size   = NVC * Nst;
  int size_b = NVC * Nxyz;

#pragma acc data present(v2[0:size], buf[0:size_b]), copyin(bc2, Nxyz)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;

    for(int ic = 0; ic < NC; ++ic){
      v2[IDX2_1SP_R(ic, ist)] += -bc2 * buf[2*ic   + NVC * ixyz];
      v2[IDX2_1SP_I(ic, ist)] += -bc2 * buf[2*ic+1 + NVC * ixyz];
    }

  }

 }

}

//==================================================================== 
void mult_staggered_tmb(real_t *restrict v2, real_t *restrict u,
                           real_t *restrict v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int size   = NVC * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                 copyin(bc[0:4], idir, Nxyz, Nt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == 0) bc2 = bc[idir];

    real_t vt[NVC], ut[NVC];

    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic]   = v1[ IDX2_1SP_R(ic, nei) ];
      vt[2*ic+1] = v1[ IDX2_1SP_I(ic, nei) ];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
      }

      real_t wtr = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      real_t wti = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);

      v2[IDX2_1SP_R(ic, ist)] +=  -bc2 * wtr;
      v2[IDX2_1SP_I(ic, ist)] +=  -bc2 * wti;

    }

  }

 }

}

//====================================================================

