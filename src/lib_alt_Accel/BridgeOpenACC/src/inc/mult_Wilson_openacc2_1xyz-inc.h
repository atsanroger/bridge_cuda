/*!
      @file    mult_Wilson_openacc2_1xyz-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

 // int idir = 0;
 if(do_comm[0] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nyzt = Ny * Nz * Nt;

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;
    real_t bc2 = bc[0];
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

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nyzt = Ny * Nz * Nt;

#pragma acc loop
  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
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

 } // do_comm[0]

 // idir = 1;
 if(do_comm[1] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nzt = Nz * Nt;

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;
     real_t bc2 = bc[1];

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
  }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nzt = Nz * Nt;
   int Nst = Nx * Ny * Nzt;

#pragma acc loop collapse(2)
  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
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

 }

 } // do_comm[1]

 //  idir = 2;
 if(do_comm[2] > 0){

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxy = Nx * Ny;

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;
     real_t bc2 = bc[2];

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
  }

 }

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
   int Nxy = Nx * Ny;
   int Nst = Nxy * Nz * Nt;

#pragma acc loop collapse(2)
  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
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

 }

 } // do_comm[2]

//============================================================END=====
