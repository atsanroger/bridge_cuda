/*!
      @file    afopr_Wilson_vector2_1xyz-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

 if(do_comm[0] > 0){

  int idir = 0;

  int Nyzt = Ny * Nz * Nt;
  int Nbf  = Nyzt;

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;
    real_t bc2 = bc[0];
    for(int ic = 0; ic < NC; ++ic){
      real_t vt1[2], vt2[2];
      vt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
      vt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
      vt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
      vt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
      buf_xp[IDXV_2SP_R(Nbf, ic, 0, iyzt)] = bc2 * vt1[0];
      buf_xp[IDXV_2SP_I(Nbf, ic, 0, iyzt)] = bc2 * vt1[1];
      buf_xp[IDXV_2SP_R(Nbf, ic, 1, iyzt)] = bc2 * vt2[0];
      buf_xp[IDXV_2SP_I(Nbf, ic, 1, iyzt)] = bc2 * vt2[1];
    }
  }

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt1[NVC], vt2[NVC], ut[NDF];
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
      vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
      vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
      vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
    }

    real_t wt1[2], wt2[2];

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, idir)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, idir)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      buf_xm[IDXV_2SP_R(Nbf, ic, 0, iyzt)] = wt1[0];
      buf_xm[IDXV_2SP_I(Nbf, ic, 0, iyzt)] = wt1[1];
      buf_xm[IDXV_2SP_R(Nbf, ic, 1, iyzt)] = wt2[0];
      buf_xm[IDXV_2SP_I(Nbf, ic, 1, iyzt)] = wt2[1];

    }

  }

 } // do_comm[0]

 // idir = 1;
 if(do_comm[1] > 0){

   int Nzt = Nz * Nt;
   int Nbf = Nx * Nzt;

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;
     real_t bc2 = bc[1];

    for(int ic = 0; ic < NC; ++ic){
      real_t vt1[2], vt2[2];
      vt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
      vt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
      vt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
      vt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
      buf_yp[IDXV_2SP_R(Nbf, ic, 0, ixzt)] = bc2 * vt1[0];
      buf_yp[IDXV_2SP_I(Nbf, ic, 0, ixzt)] = bc2 * vt1[1];
      buf_yp[IDXV_2SP_R(Nbf, ic, 1, ixzt)] = bc2 * vt2[0];
      buf_yp[IDXV_2SP_I(Nbf, ic, 1, ixzt)] = bc2 * vt2[1];
    }

   }
  }

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny*izt);
     int ixzt = ix + Nx * izt;
     int idir = 1;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
       vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
       vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
       vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, idir)];
         ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, idir)];
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

       buf_ym[IDXV_2SP_R(Nbf, ic, 0, ixzt)] = wt1[0];
       buf_ym[IDXV_2SP_I(Nbf, ic, 0, ixzt)] = wt1[1];
       buf_ym[IDXV_2SP_R(Nbf, ic, 1, ixzt)] = wt2[0];
       buf_ym[IDXV_2SP_I(Nbf, ic, 1, ixzt)] = wt2[1];

     }

   }
  }

 } // do_comm[1]

 //  idir = 2;
 if(do_comm[2] > 0){

   int Nxy = Nx * Ny;
   int Nbf = Nxy * Nt;

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;
     real_t bc2 = bc[2];

     for(int ic = 0; ic < NC; ++ic){
       real_t wt1[2], wt2[2];
       wt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
       wt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
       wt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
       wt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
       buf_zp[IDXV_2SP_R(Nbf, ic, 0, ixyt)] = bc2 * wt1[0];
       buf_zp[IDXV_2SP_I(Nbf, ic, 0, ixyt)] = bc2 * wt1[1];
       buf_zp[IDXV_2SP_R(Nbf, ic, 1, ixyt)] = bc2 * wt2[0];
       buf_zp[IDXV_2SP_I(Nbf, ic, 1, ixyt)] = bc2 * wt2[1];
     }

   }
  }

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;
     int idir = 2;

     real_t vt1[NVC], vt2[NVC];
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
       vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
       vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
       vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
     }

     for(int ic = 0; ic < NC; ++ic){

       real_t ut[NVC];
       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, idir)];
         ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, idir)];
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
       buf_zm[IDXV_2SP_R(Nbf, ic, 0, ixyt)] = wt1[0];
       buf_zm[IDXV_2SP_I(Nbf, ic, 0, ixyt)] = wt1[1];
       buf_zm[IDXV_2SP_R(Nbf, ic, 1, ixyt)] = wt2[0];
       buf_zm[IDXV_2SP_I(Nbf, ic, 1, ixyt)] = wt2[1];

    }

   }
  }

 } // do_comm[2]

//============================================================END=====
