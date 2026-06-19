/*!
      @file    mult_CloverTerm_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/


#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)


//====================================================================
void mult_csw_dirac(real_t *restrict v2, real_t *restrict u,
                    real_t *restrict v1, int *Nsize, int iflag, int Nc)
  // iflag = 0: mult, iflag = 1: multadd
{
  int Nst    = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int size   = NVC * ND * Nst;
  int size_u = NDF * Nst * ND * ND2;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]) \
                 copyin(Nst, iflag)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){

    real_t ut[NVC], vt[NVC*ND], wt[NVC*ND];

    if(iflag == 0){
      for(int id = 0; id < ND; ++id){
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc + NVC*id] = 0.0;
        }
      }
    }else{
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          wt[2*ic   + NVC*id] =  v2[IDX2_SP_R(ic, id, ist)];
          wt[2*ic+1 + NVC*id] =  v2[IDX2_SP_I(ic, id, ist)];
        }
      }
    }

    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] =  v1[IDX2_SP_R(ic, id, ist)];
        vt[2*ic+1 + NVC*id] =  v1[IDX2_SP_I(ic, id, ist)];
      }
    }

    for(int jd = 0; jd < ND2; ++jd){
     for(int id = 0; id < ND; ++id){
       int igst = ist + Nst * (id + ND * jd);

       for(int ic2 = 0; ic2 < NC; ++ic2){

         for(int ic1 = 0; ic1 < NC; ++ic1){
           ut[2*ic1  ] = u[IDX2_G_R(ic1, ic2, igst)];
           ut[2*ic1+1] = u[IDX2_G_I(ic1, ic2, igst)];
         }

         int id2 = (id + ND2) % ND;
         int j = NVC * id;
         int k = NVC * id2;
         real_t wt1r, wt1i, wt2r, wt2i;
         wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+k], vt[1+k], vt[2+k], vt[3+k], vt[4+k], vt[5+k]);

         wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+k], vt[1+k], vt[2+k], vt[3+k], vt[4+k], vt[5+k]);

         int jd2 = jd + ND2;
         wt[2*ic2   + NVC*jd]  += wt1r;
         wt[2*ic2+1 + NVC*jd]  += wt1i;
         wt[2*ic2   + NVC*jd2] += wt2r;
         wt[2*ic2+1 + NVC*jd2] += wt2i;

       }
     }
    }

    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        v2[IDX2_SP_R(ic,id,ist)] = wt[2*ic   + NVC*id];
        v2[IDX2_SP_I(ic,id,ist)] = wt[2*ic+1 + NVC*id];
      }
    }

  }

 } // acc parallel

 //#pragma omp barrier

}

//====================================================================
void mult_csw_chiral(real_t *restrict v2, real_t *restrict u,
                     real_t *restrict v1, int *Nsize, int iflag, int Nc)
  // iflag = 0: mult, iflag = 1: multadd
{
  int Nst    = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];
  int size   = NVC * ND * Nst;
  int size_u = NDF * Nst * ND * ND2;

#pragma acc data present(v2[0:size], v1[0:size], u[0:size_u]) \
                 copyin(Nst, iflag)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for(int ist = 0; ist < Nst; ++ist){

    real_t ut[NVC], vt[NVC*ND], wt[NVC*ND];

    if(iflag == 0){
      for(int id = 0; id < ND; ++id){
        for(int ivc = 0; ivc < NVC; ++ivc){
          wt[ivc + NVC*id] = 0.0;
        }
      }
    }else{
      for(int id = 0; id < ND; ++id){
        for(int ic = 0; ic < NC; ++ic){
          wt[2*ic   + NVC*id] =  v2[IDX2_SP_R(ic, id, ist)];
          wt[2*ic+1 + NVC*id] =  v2[IDX2_SP_I(ic, id, ist)];
        }
      }
    }

    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] =  v1[IDX2_SP_R(ic, id, ist)];
        vt[2*ic+1 + NVC*id] =  v1[IDX2_SP_I(ic, id, ist)];
      }
    }

    for(int jd = 0; jd < ND2; ++jd){
     for(int id = 0; id < ND2; ++id){
       int igst = ist + Nst * (id + ND2 * jd);

       for(int ic2 = 0; ic2 < NC; ++ic2){

         for(int ic1 = 0; ic1 < NC; ++ic1){
           ut[2*ic1  ] = u[IDX2_G_R(ic1, ic2, igst)];
           ut[2*ic1+1] = u[IDX2_G_I(ic1, ic2, igst)];
         }

         int j = NVC * id;
         real_t wt1r, wt1i;
         wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         wt[2*ic2   + NVC*jd]  += wt1r;
         wt[2*ic2+1 + NVC*jd]  += wt1i;

       }

     }
    }

    for(int jd = 0; jd < ND2; ++jd){
     for(int id = 0; id < ND2; ++id){
       int igst = ist + Nst * (id + ND2*jd + ND);

       for(int ic2 = 0; ic2 < NC; ++ic2){

         for(int ic1 = 0; ic1 < NC; ++ic1){
           ut[2*ic1  ] = u[IDX2_G_R(ic1, ic2, igst)];
           ut[2*ic1+1] = u[IDX2_G_I(ic1, ic2, igst)];
         }

         int j = NVC * (id + ND2);
         real_t wt1r, wt1i;
         wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                   vt[0+j], vt[1+j], vt[2+j], vt[3+j], vt[4+j], vt[5+j]);

         int jd2 = jd + ND2;
         wt[2*ic2   + NVC*jd2]  += wt1r;
         wt[2*ic2+1 + NVC*jd2]  += wt1i;

       }

     }
    }

    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        v2[IDX2_SP_R(ic,id,ist)] = wt[2*ic   + NVC*id];
        v2[IDX2_SP_I(ic,id,ist)] = wt[2*ic+1 + NVC*id];
      }
    }

  }

 } // acc parallel

 //#pragma omp barrier

}

//============================================================END=====
