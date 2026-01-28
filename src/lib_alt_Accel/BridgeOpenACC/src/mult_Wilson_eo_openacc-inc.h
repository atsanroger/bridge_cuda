/*!
      @file    mult_Wilson_eo_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/


#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)


//====================================================================
/*
inline void check_setup()
{
  if(VLEN2 < 2){
    vout.crucial("VLEN2 = %d is too small for this implementation\n",
                 VLEN2);
    exit(EXIT_FAILURE);
  }
}
*/

//====================================================================
void mult_wilson_xp1_eo(real_t *restrict buf, real_t *restrict v1,
                        int *Nsize, int *bc, int ieo, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  real_t bc2 = bc[0];

  int size   = NVC * ND * Nx * Nyzt;
  int size_b = NVC * 2 * ((Nyzt+1)/2);

#pragma acc data present(v1[0:size], buf[0:size_b]), \
                 copyin(ieo, bc2, Nx, Ny, Nz, Nt, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int keo = (ieo + iy + iz + it) % 2;

    if(keo == 1){
      int iyzt2 = iyzt/2;

      real_t vt[NVC*ND];
      for(int id = 0; id < ND; ++id){
       for(int ic = 0; ic < NC; ++ic){
         vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
         vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
       }
      }

      real_t *vt1 = &buf[      NVC*2*iyzt2];
      real_t *vt2 = &buf[NVC + NVC*2*iyzt2];

      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID4] );
        vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID4] );
        vt2[icr] = bc2 *( vt[icr + ID2] - vt[ici + ID3] );
        vt2[ici] = bc2 *( vt[ici + ID2] + vt[icr + ID3] );
      }

    }

  }

 }

}

//====================================================================
void mult_wilson_xp2_eo(real_t *restrict v2, real_t *restrict u,
                        real_t *restrict buf, 
                        int *Nsize, int *bc, int ieo, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;
  int Nst  = Nx * Nyzt;

  int size   = NVC * ND * Nst;
  int size_b = NVC * 2 * ((Nyzt+1)/2);
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(ieo, Nx, Ny, Nz, Nt, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int keo = (ieo + iy + iz + it) % 2;

    if(keo == 1){
      int iyzt2 = iyzt/2;

      real_t vt1[NVC], vt2[NVC], ut[NDF];
      real_t wt1[2], wt2[2];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf[ivc +       2*NVC*iyzt2];
        vt2[ivc] = buf[ivc + NVC + 2*NVC*iyzt2];
      }

      for(int ic = 0; ic < NC; ++ic){

        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
          ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
        }

        wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
        v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
        v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
        v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
        v2[IDX2_SP_R(ic, 2, ist)] +=  wt2[1];
        v2[IDX2_SP_I(ic, 2, ist)] += -wt2[0];
        v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[1];
        v2[IDX2_SP_I(ic, 3, ist)] += -wt1[0];

      }

    }

  }

 }

}

//====================================================================
 void mult_wilson_xpb_eo(real_t *restrict v2, real_t *restrict u,
                         real_t *restrict v1,
                         int *Nsize, int *bc, int ieo, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;
  int Nst  = Nx * Nyzt;

  int size   = NVC * ND * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                 copyin(ieo, bc[0:4], Nx, Nyzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  for(int ist = 0; ist < Nst; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int iyzt = ist/Nx;
    int keo  = (ieo + iy + iz + it) % 2;
    int nei  = ((ix + keo) % Nx) + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == Nx-1 && keo == 1) bc2 = bc[0];

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID4] );
      vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID4] );
      vt2[icr] = bc2 *( vt[icr + ID2] - vt[ici + ID3] );
      vt2[ici] = bc2 *( vt[ici + ID2] + vt[icr + ID3] );
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
        ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
      v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
      v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
      v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
      v2[IDX2_SP_R(ic, 2, ist)] +=  wt2[1];
      v2[IDX2_SP_I(ic, 2, ist)] += -wt2[0];
      v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[1];
      v2[IDX2_SP_I(ic, 3, ist)] += -wt1[0];

    }

  }

 }

}

//====================================================================
void mult_wilson_xm1_eo(real_t *restrict buf, real_t *restrict u,
                        real_t *restrict v1,
                        int *Nsize, int *bc, int ieo, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;
  int Nst  = Nx * Nyzt;

  int size   = NVC * ND * Nst;
  int size_b = NVC * 2 * ((Nyzt+1)/2);
  int size_u = NDF * Nst;

#pragma acc data present(v1[0:size], buf[0:size_b], u[0:size_u]), \
                 copyin(ieo, Nx, Ny, Nz, Nt, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int keo = (ieo + iy + iz + it) % 2;

    if(keo == 0){
      int iyzt2 = iyzt/2;

      real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];

      for(int id = 0; id < ND; ++id){
       for(int ic = 0; ic < NC; ++ic){
         vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
         vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
       }
      }

      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = vt[icr + ID1] + vt[ici + ID4];
        vt1[ici] = vt[ici + ID1] - vt[icr + ID4];
        vt2[icr] = vt[icr + ID2] + vt[ici + ID3];
        vt2[ici] = vt[ici + ID2] - vt[icr + ID3];
      }

      real_t *wt1 = &buf[      NVC*2*iyzt2];
      real_t *wt2 = &buf[NVC + NVC*2*iyzt2];

      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;

        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
          ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
        }

        wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                             vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      }

    }

  }

 }

}

//====================================================================
void mult_wilson_xm2_eo(real_t *restrict v2, real_t *restrict buf, 
                        int *Nsize, int *bc, int ieo, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;
  int Nst  = Nx * Nyzt;

  real_t bc2 = bc[0];

  int size   = NVC * ND * Nst;
  int size_b = NVC * 2 * ((Nyzt+1)/2);

#pragma acc data present(v2[0:size], buf[0:size_b]), \
                 copyin(bc2, ieo, Nx, Ny, Nz, Nt, Nyzt)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int iy = iyzt % Ny;
    int iz = (iyzt/Ny) % Nz;
    int it = iyzt/(Ny * Nz);
    int ist = ix + Nx * iyzt;
    int keo = (ieo + iy + iz + it) % 2;

    if(keo == 0){
      int iyzt2 = iyzt/2;

      real_t wt1[2], wt2[2];

      for(int ic = 0; ic < NC; ++ic){
        int icr = 2 * ic;
        int ici = 2 * ic + 1;
        wt1[0] = bc2 * buf[icr       + 2*NVC*iyzt2];
        wt1[1] = bc2 * buf[ici       + 2*NVC*iyzt2];
        wt2[0] = bc2 * buf[icr + NVC + 2*NVC*iyzt2];
        wt2[1] = bc2 * buf[ici + NVC + 2*NVC*iyzt2];

        v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
        v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
        v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
        v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
        v2[IDX2_SP_R(ic, 2, ist)] += -wt2[1];
        v2[IDX2_SP_I(ic, 2, ist)] +=  wt2[0];
        v2[IDX2_SP_R(ic, 3, ist)] += -wt1[1];
        v2[IDX2_SP_I(ic, 3, ist)] +=  wt1[0];
      }

    }

  }

 }

}

//====================================================================
 void mult_wilson_xmb_eo(real_t *restrict v2, real_t *restrict u,
                         real_t *restrict v1,
                         int *Nsize, int *bc, int ieo, int Nc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;
  int Nst  = Nx * Nyzt;

  int size   = NVC * ND * Nst;
  int size_u = NDF * Nst;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]), \
                 copyin(bc[0:4], ieo, Nx, Ny, Nz, Nt, Nyzt, Nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);
    int iyzt = ist/Nx;
    int keo = (ieo + iy + iz + it) % 2;
    int ix2 = (ix - 1 + keo + Nx) % Nx;
    int nei  = ix2 + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == 0 && keo == 0) bc2 = bc[0];

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = vt[icr + ID1] + vt[ici + ID4];
      vt1[ici] = vt[ici + ID1] - vt[icr + ID4];
      vt2[icr] = vt[icr + ID2] + vt[ici + ID3];
      vt2[ici] = vt[ici + ID2] - vt[icr + ID3];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
        ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDX2_SP_R(ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDX2_SP_I(ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDX2_SP_R(ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDX2_SP_I(ic, 1, ist)] +=  bc2 * wt2[1];
      v2[IDX2_SP_R(ic, 2, ist)] += -bc2 * wt2[1];
      v2[IDX2_SP_I(ic, 2, ist)] +=  bc2 * wt2[0];
      v2[IDX2_SP_R(ic, 3, ist)] += -bc2 * wt1[1];
      v2[IDX2_SP_I(ic, 3, ist)] +=  bc2 * wt1[0];

    }

  }

 }

}

//====================================================================

