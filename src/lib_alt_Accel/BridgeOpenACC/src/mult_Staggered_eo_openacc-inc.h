/*!
      @file    mult_Staggered_eo_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)


//====================================================================
void mult_staggered_xp1_eo(real_t *restrict buf, real_t *restrict v1,
                           int *Nsize, int *bc, int ieo, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  real_t bc2 = bc[0];

  int size   = NVC * Nx * Nyzt;
  int size_b = NVC * ((Nyzt+1)/2);

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
      real_t *vt1 = &buf[NVC*iyzt2];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = bc2 * v1[IDX2_1SP(ivc, ist)];
      }
    }

  }

 }

}

//====================================================================
void mult_staggered_xp2_eo(real_t *restrict v2, real_t *restrict u,
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

  int size   = NVC * Nst;
  int size_b = NVC * ((Nyzt+1)/2);
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

      real_t vt[NVC], ut[NVC];

      for(int ivc = 0; ivc < NVC; ++ivc){
	vt[ivc] = buf[ivc + NVC*iyzt2];
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
 void mult_staggered_xpb_eo(real_t *restrict v2, real_t *restrict u,
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

  int size   = NVC * Nst;
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
void mult_staggered_xm1_eo(real_t *restrict buf, real_t *restrict u,
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

  int size   = NVC * Nst;
  int size_b = NVC * ((Nyzt+1)/2);
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

      real_t vt[NVC], ut[NVC];

      for(int ic = 0; ic < NC; ++ic){
	vt[2*ic  ] = v1[ IDX2_1SP_R(ic, ist) ];
	vt[2*ic+1] = v1[ IDX2_1SP_I(ic, ist) ];
      }

      real_t *wt1 = &buf[NVC*iyzt2];

      for(int ic = 0; ic < NC; ++ic){
	for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
          ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
	}
        wt1[2*ic  ] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                                vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
        wt1[2*ic+1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                                vt[0], vt[1], vt[2], vt[3], vt[4], vt[5]);
      }

    }

  }

 }

}

//====================================================================
void mult_staggered_xm2_eo(real_t *restrict v2, real_t *restrict buf, 
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

  int size   = NVC * Nst;
  int size_b = NVC * ((Nyzt+1)/2);

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

      for(int ic = 0; ic < NC; ++ic){
	v2[IDX2_1SP_R(ic, ist)] += -bc2 * buf[2*ic   + NVC*iyzt2];
	v2[IDX2_1SP_I(ic, ist)] += -bc2 * buf[2*ic+1 + NVC*iyzt2];
      }

    }

  }

 }

}

//====================================================================
void mult_staggered_xmb_eo(real_t *restrict v2, real_t *restrict u,
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

  int size   = NVC * Nst;
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

