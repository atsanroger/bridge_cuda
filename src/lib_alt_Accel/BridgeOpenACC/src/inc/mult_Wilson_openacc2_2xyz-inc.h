/*!
      @file    mult_Wilson_openacc2_2xyz-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

    // idir = 0
    if(do_comm[0] > 0){

      if(ix == Nx-1){
        int iyzt = iy + Ny * (iz + Nz * it);

        real_t vt1[NVC], vt2[NVC], ut[NDF];
        real_t wt1[2], wt2[2];

        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_xp[IDXBF(ivc, 0, iyzt)];
          vt2[ivc] = buf_xp[IDXBF(ivc, 1, iyzt)];
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

          v2L[2*ic   + ID1] +=  wt1[0];
          v2L[2*ic+1 + ID1] +=  wt1[1];
          v2L[2*ic   + ID2] +=  wt2[0];
          v2L[2*ic+1 + ID2] +=  wt2[1];
          v2L[2*ic   + ID3] +=  wt2[1];
          v2L[2*ic+1 + ID3] += -wt2[0];
          v2L[2*ic   + ID4] +=  wt1[1];
          v2L[2*ic+1 + ID4] += -wt1[0];

	}

      }

      if(ix == 0){
        int iyzt = iy + Ny * (iz + Nz * it);
        real_t bc2 = bc[0];
        real_t wt1[2], wt2[2];

        for(int ic = 0; ic < NC; ++ic){

          wt1[0] = bc2 * buf_xm[IDXBF_R(ic, 0, iyzt)];
          wt1[1] = bc2 * buf_xm[IDXBF_I(ic, 0, iyzt)];
          wt2[0] = bc2 * buf_xm[IDXBF_R(ic, 1, iyzt)];
          wt2[1] = bc2 * buf_xm[IDXBF_I(ic, 1, iyzt)];

          v2L[2*ic   + ID1] +=  wt1[0];
          v2L[2*ic+1 + ID1] +=  wt1[1];
          v2L[2*ic   + ID2] +=  wt2[0];
          v2L[2*ic+1 + ID2] +=  wt2[1];
          v2L[2*ic   + ID3] += -wt2[1];
          v2L[2*ic+1 + ID3] +=  wt2[0];
          v2L[2*ic   + ID4] += -wt1[1];
          v2L[2*ic+1 + ID4] +=  wt1[0];
	}

      }

    }

    // idir = 1
    if(do_comm[1] > 0){

      if(iy == Ny-1){
	int ixzt = ix + Nx * (iz + Nz * it);
        int istu = ist + Nst;

        real_t vt1[NVC], vt2[NVC], ut[NDF];
        real_t wt1[2], wt2[2];

        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_yp[IDXBF(ivc, 0, ixzt)];
          vt2[ivc] = buf_yp[IDXBF(ivc, 1, ixzt)];
        }

        for(int ic = 0; ic < NC; ++ic){

          for(int ic2 = 0; ic2 < NC; ++ic2){
            ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
            ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
          }

          wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

          v2L[ic*2   + ID1] +=  wt1[0];
          v2L[ic*2+1 + ID1] +=  wt1[1];
          v2L[ic*2   + ID2] +=  wt2[0];
          v2L[ic*2+1 + ID2] +=  wt2[1];
          v2L[ic*2   + ID3] += -wt2[0];
          v2L[ic*2+1 + ID3] += -wt2[1];
          v2L[ic*2   + ID4] +=  wt1[0];
          v2L[ic*2+1 + ID4] +=  wt1[1];

	}

      }

      if(iy == 0){
	int ixzt = ix + Nx * (iz + Nz * it);
        real_t wt1[2], wt2[2];
        real_t bc2 = bc[1];

        for(int ic = 0; ic < NC; ++ic){
          wt1[0] = bc2 * buf_ym[IDXBF_R(ic, 0, ixzt)];
          wt1[1] = bc2 * buf_ym[IDXBF_I(ic, 0, ixzt)];
          wt2[0] = bc2 * buf_ym[IDXBF_R(ic, 1, ixzt)];
          wt2[1] = bc2 * buf_ym[IDXBF_I(ic, 1, ixzt)];
          v2L[ic*2   + ID1] +=  wt1[0];
          v2L[ic*2+1 + ID1] +=  wt1[1];
          v2L[ic*2   + ID2] +=  wt2[0];
          v2L[ic*2+1 + ID2] +=  wt2[1];
          v2L[ic*2   + ID3] +=  wt2[0];
          v2L[ic*2+1 + ID3] +=  wt2[1];
          v2L[ic*2   + ID4] += -wt1[0];
          v2L[ic*2+1 + ID4] += -wt1[1];

        }

      }

    }

    // idir = 2
    if(do_comm[2] > 0){

      if(iz == Nz-1){
        int ixyt = ix + Nx * (iy + Ny * it);
        int istu = ist + Nst * 2;

        real_t vt1[NVC], vt2[NVC];
        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_zp[IDXBF(ivc, 0, ixyt)];
          vt2[ivc] = buf_zp[IDXBF(ivc, 1, ixyt)];
        }

        for(int ic = 0; ic < NC; ++ic){

          real_t ut[NVC], wt1[2], wt2[2];
          for(int ic2 = 0; ic2 < NC; ++ic2){
            ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, istu)];
            ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, istu)];
          }

          wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

          v2L[ic*2   + ID1] +=  wt1[0];
          v2L[ic*2+1 + ID1] +=  wt1[1];
          v2L[ic*2   + ID2] +=  wt2[0];
          v2L[ic*2+1 + ID2] +=  wt2[1];
          v2L[ic*2   + ID3] +=  wt1[1];
          v2L[ic*2+1 + ID3] += -wt1[0];
          v2L[ic*2   + ID4] += -wt2[1];
          v2L[ic*2+1 + ID4] +=  wt2[0];

	}

      }

      if(iz == 0){
        int ixyt = ix + Nx * (iy + Ny * it);
        real_t bc2 = bc[2];
        for(int ic = 0; ic < NC; ++ic){
          real_t wt1[2], wt2[2];
          wt1[0] = bc2 * buf_zm[IDXBF_R(ic, 0, ixyt)];
          wt1[1] = bc2 * buf_zm[IDXBF_I(ic, 0, ixyt)];
          wt2[0] = bc2 * buf_zm[IDXBF_R(ic, 1, ixyt)];
          wt2[1] = bc2 * buf_zm[IDXBF_I(ic, 1, ixyt)];
          v2L[ic*2   + ID1] +=  wt1[0];
          v2L[ic*2+1 + ID1] +=  wt1[1];
          v2L[ic*2   + ID2] +=  wt2[0];
          v2L[ic*2+1 + ID2] +=  wt2[1];
          v2L[ic*2   + ID3] += -wt1[1];
          v2L[ic*2+1 + ID3] +=  wt1[0];
          v2L[ic*2   + ID4] +=  wt2[1];
          v2L[ic*2+1 + ID4] += -wt2[0];
	}

      }

    }

//============================================================END=====
