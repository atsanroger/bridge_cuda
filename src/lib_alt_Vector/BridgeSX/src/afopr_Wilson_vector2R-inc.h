/*!
      @file    afopr_Wilson_vector2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

// This code explicitly assumes SU(3) gauge group.
//#include "lib_alt_Vector/inline/define_params.h"
//#include "lib_alt_Vector/inline/define_index.h"

#define MULT_GXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_GXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)
#define MULT_GDXr(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v0+u1*v1 + u2*v2+u3*v3 + u4*v4+u5*v5)
#define MULT_GDXi(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)  (u0*v1-u1*v0 + u2*v3-u3*v2 + u4*v5-u5*v4)

#define EXT_IMG_R(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (v1r*w2r - v1i*w2i - v2r*w1r + v2i*w1i)
#define EXT_IMG_I(v1r,v1i,v2r,v2i,w1r,w1i,w2r,w2i)  (- v1r*w2i - v1i*w2r + v2r*w1i + v2i*w1r)


//====================================================================
void mult_wilson_D_dirac(real_t *RESTRICT v2, real_t *RESTRICT u,
                         real_t *RESTRICT v1,
                         int *Nsize, int *bc, real_t kappa)
{

  int nx = Nsize[0];
  int ny = Nsize[1];
  int nz = Nsize[2];
  int nt = Nsize[3];
  int nst  = nx*ny*nz*nt;
  int nvst = NVC * ND * nst;
  int ngst = NDF * nst * 4;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, nst);

  //  int nst2 = nst/VLEN;
  //  for(int ist2 = 0; ist2 < nst2; ++ist2){
  //    for(int site2 = 0; site2 < VLEN; ++site2){
  //      int site = site2 + VLEN * ist2;
  int nst2 = (ns-is)/VLEN;
  for(int ist2 = 0; ist2 < nst2; ++ist2){
    for(int site2 = 0; site2 < VLEN; ++site2){
      int site = is + site2 + VLEN * ist2;
 
      real_t bc2;

      real_t u_0, u_1, u_2, u_3, u_4, u_5;
      real_t u_6, u_7, u_8, u_9, u10, u11;
      real_t u12, u13, u14, u15, u16, u17;
      real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
      real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
      real_t wt1r, wt1i, wt2r, wt2i;

      real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
      real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
      real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
      real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;


      //#include "afopr_Wilson_vector2_xyz-inc.h"

      // mult_xp
      int idir = 0;

      int iyzt = site/nx;
      int ix   = site % nx;
      int nn = (ix+1) % nx;

      int isn = nn + iyzt*nx;
      int isg = ix + iyzt*nx;

      /*
      real_t vt1[NVC], vt2[NVC], ut[NVC];
      real_t v2L[NVC][ND+1];

      bc2 = 1.0;
      if(ix == nx-1) bc2 = bc[0];

      for(int ic = 0; ic < NC; ++ic){
        vt1[2*ic]   =  v1[IDXV_SP_R(nst, ic, 0, isn)]
                     - v1[IDXV_SP_I(nst, ic, 3, isn)];
        vt1[2*ic+1] =  v1[IDXV_SP_I(nst, ic, 0, isn)]
                     + v1[IDXV_SP_R(nst, ic, 3, isn)];

        vt2[2*ic]   =  v1[IDXV_SP_R(nst, ic, 1, isn)]
                     - v1[IDXV_SP_I(nst, ic, 2, isn)];
        vt2[2*ic+1] =  v1[IDXV_SP_I(nst, ic, 1, isn)]
                     + v1[IDXV_SP_R(nst, ic, 2, isn)];
      }

      for(int ic = 0; ic < NC; ++ic){

        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2]   = u[IDXV_G_R(nst, ic2, ic, isg, idir)];
          ut[2*ic2+1] = u[IDXV_G_I(nst, ic2, ic, isg, idir)];
	}

        wt1r = MULT_GXr(ut[0],  ut[1],  ut[2],  ut[3],  ut[4],  ut[5],
                        vt1[0], vt1[1], vt1[2], vt1[3], vt1[4], vt1[5]);
        wt1i = MULT_GXi(ut[0],  ut[1],  ut[2],  ut[3],  ut[4],  ut[5],
                        vt1[0], vt1[1], vt1[2], vt1[3], vt1[4], vt1[5]);
        wt2r = MULT_GXr(ut[0],  ut[1],  ut[2],  ut[3],  ut[4],  ut[5],
                        vt2[0], vt2[1], vt2[2], vt2[3], vt2[4], vt2[5]);
        wt2i = MULT_GXi(ut[0],  ut[1],  ut[2],  ut[3],  ut[4],  ut[5],
                        vt2[0], vt2[1], vt2[2], vt2[3], vt2[4], vt2[5]);

        int icr = 2*ic;
        int ici = 2*ic+1;
        v2L[icr][1] =  bc2 * wt1r;
        v2L[ici][1] =  bc2 * wt1i;
        v2L[icr][2] =  bc2 * wt2r;
        v2L[ici][2] =  bc2 * wt2i;
        v2L[icr][3] =  bc2 * wt2i;
        v2L[ici][3] = -bc2 * wt2r;
        v2L[icr][4] =  bc2 * wt1i;
        v2L[ici][4] = -bc2 * wt1r;

      }

      v2_01 = v2L[0][1];
      v2_11 = v2L[1][1];
      v2_02 = v2L[0][2];
      v2_12 = v2L[1][2];
      v2_03 = v2L[0][3];
      v2_13 = v2L[1][3];
      v2_04 = v2L[0][4];
      v2_14 = v2L[1][4];

      v2_21 = v2L[2][1];
      v2_31 = v2L[3][1];
      v2_22 = v2L[2][2];
      v2_32 = v2L[3][2];
      v2_23 = v2L[2][3];
      v2_33 = v2L[3][3];
      v2_24 = v2L[2][4];
      v2_34 = v2L[3][4];

      v2_41 = v2L[4][1];
      v2_51 = v2L[5][1];
      v2_42 = v2L[4][2];
      v2_52 = v2L[5][2];
      v2_43 = v2L[4][3];
      v2_53 = v2L[5][3];
      v2_44 = v2L[4][4];
      v2_54 = v2L[5][4];
      */

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] - v1[IDXV_SP_I(nst, 0,3,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] + v1[IDXV_SP_R(nst, 0,3,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] - v1[IDXV_SP_I(nst, 1,3,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] + v1[IDXV_SP_R(nst, 1,3,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] - v1[IDXV_SP_I(nst, 2,3,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] + v1[IDXV_SP_R(nst, 2,3,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] - v1[IDXV_SP_I(nst, 0,2,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] + v1[IDXV_SP_R(nst, 0,2,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] - v1[IDXV_SP_I(nst, 1,2,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] + v1[IDXV_SP_R(nst, 1,2,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] - v1[IDXV_SP_I(nst, 2,2,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] + v1[IDXV_SP_R(nst, 2,2,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 2,0,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 2,0,isg, idir)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;
      bc2 = 1.0;
      if(ix == nx-1) bc2 = bc[0];
      v2_01 =  bc2 * wt1r;
      v2_11 =  bc2 * wt1i;
      v2_02 =  bc2 * wt2r;
      v2_12 =  bc2 * wt2i;
      v2_03 =  bc2 * wt2i;
      v2_13 = -bc2 * wt2r;
      v2_04 =  bc2 * wt1i;
      v2_14 = -bc2 * wt1r;

      u_6 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 2,1,isg, idir)];
      u11 = u[IDXV_G_I(nst, 2,1,isg, idir)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 =  bc2 * wt1r;
      v2_31 =  bc2 * wt1i;
      v2_22 =  bc2 * wt2r;
      v2_32 =  bc2 * wt2i;
      v2_23 =  bc2 * wt2i;
      v2_33 = -bc2 * wt2r;
      v2_24 =  bc2 * wt1i;
      v2_34 = -bc2 * wt1r;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_41 =  bc2 * wt1r;
      v2_51 =  bc2 * wt1i;
      v2_42 =  bc2 * wt2r;
      v2_52 =  bc2 * wt2i;
      v2_43 =  bc2 * wt2i;
      v2_53 = -bc2 * wt2r;
      v2_44 =  bc2 * wt1i;
      v2_54 = -bc2 * wt1r;



      // mult_xm
      nn = (ix+nx-1) % nx;

      isn = nn + (iyzt)*nx;
      isg = nn + iyzt*nx;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] + v1[IDXV_SP_I(nst, 0,3,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] - v1[IDXV_SP_R(nst, 0,3,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] + v1[IDXV_SP_I(nst, 1,3,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] - v1[IDXV_SP_R(nst, 1,3,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] + v1[IDXV_SP_I(nst, 2,3,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] - v1[IDXV_SP_R(nst, 2,3,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] + v1[IDXV_SP_I(nst, 0,2,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] - v1[IDXV_SP_R(nst, 0,2,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] + v1[IDXV_SP_I(nst, 1,2,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] - v1[IDXV_SP_R(nst, 1,2,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] + v1[IDXV_SP_I(nst, 2,2,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] - v1[IDXV_SP_R(nst, 2,2,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 0,2,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 0,2,isg, idir)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      // ic = 0;
      bc2 = 1.0;
      if(ix == 0) bc2 = bc[0];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt2i;
      v2_13 += +bc2 * wt2r;
      v2_04 += -bc2 * wt1i;
      v2_14 += +bc2 * wt1r;

      u_6 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 1,2,isg, idir)];
      u11 = u[IDXV_G_I(nst, 1,2,isg, idir)];

      wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
      v2_23 += -bc2 * wt2i;
      v2_33 += +bc2 * wt2r;
      v2_24 += -bc2 * wt1i;
      v2_34 += +bc2 * wt1r;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 += -bc2 * wt2i;
      v2_53 += +bc2 * wt2r;
      v2_44 += -bc2 * wt1i;
      v2_54 += +bc2 * wt1r;


      // mult_yp
      int izt = site/(nx*ny);
      int iy  = (site/nx) % ny;
      ix  = site % nx;
      nn = (iy + 1) % ny;
      idir = 1;

      isn = ix + nn*nx + izt*nx*ny;
      isg = ix + iy*nx + izt*nx*ny;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] + v1[IDXV_SP_R(nst, 0,3,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] + v1[IDXV_SP_I(nst, 0,3,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] + v1[IDXV_SP_R(nst, 1,3,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] + v1[IDXV_SP_I(nst, 1,3,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] + v1[IDXV_SP_R(nst, 2,3,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] + v1[IDXV_SP_I(nst, 2,3,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] - v1[IDXV_SP_R(nst, 0,2,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] - v1[IDXV_SP_I(nst, 0,2,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] - v1[IDXV_SP_R(nst, 1,2,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] - v1[IDXV_SP_I(nst, 1,2,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] - v1[IDXV_SP_R(nst, 2,2,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] - v1[IDXV_SP_I(nst, 2,2,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 2,0,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 2,0,isg, idir)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;
      bc2 = 1.0;
      if(iy == ny-1) bc2 = bc[1];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt2r;
      v2_13 += -bc2 * wt2i;
      v2_04 +=  bc2 * wt1r;
      v2_14 +=  bc2 * wt1i;

      u_6 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 2,1,isg, idir)];
      u11 = u[IDXV_G_I(nst, 2,1,isg, idir)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
      v2_23 += -bc2 * wt2r;
      v2_33 += -bc2 * wt2i;
      v2_24 +=  bc2 * wt1r;
      v2_34 +=  bc2 * wt1i;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 += -bc2 * wt2r;
      v2_53 += -bc2 * wt2i;
      v2_44 +=  bc2 * wt1r;
      v2_54 +=  bc2 * wt1i;


      //mult_ym
      nn = (iy + ny - 1) % ny;

      isn = ix + nn*nx + izt*nx*ny;
      isg = ix + nn*nx + izt*nx*ny;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] - v1[IDXV_SP_R(nst, 0,3,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] - v1[IDXV_SP_I(nst, 0,3,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] - v1[IDXV_SP_R(nst, 1,3,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] - v1[IDXV_SP_I(nst, 1,3,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] - v1[IDXV_SP_R(nst, 2,3,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] - v1[IDXV_SP_I(nst, 2,3,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] + v1[IDXV_SP_R(nst, 0,2,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] + v1[IDXV_SP_I(nst, 0,2,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] + v1[IDXV_SP_R(nst, 1,2,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] + v1[IDXV_SP_I(nst, 1,2,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] + v1[IDXV_SP_R(nst, 2,2,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] + v1[IDXV_SP_I(nst, 2,2,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 0,2,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 0,2,isg, idir)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //    ic = 0;
      bc2 = 1.0;
      if(iy == 0) bc2 = bc[1];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 +=  bc2 * wt2r;
      v2_13 +=  bc2 * wt2i;
      v2_04 += -bc2 * wt1r;
      v2_14 += -bc2 * wt1i;

      u_6 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 1,2,isg, idir)];
      u11 = u[IDXV_G_I(nst, 1,2,isg, idir)];

      wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
      v2_23 +=  bc2 * wt2r;
      v2_33 +=  bc2 * wt2i;
      v2_24 += -bc2 * wt1r;
      v2_34 += -bc2 * wt1i;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 +=  bc2 * wt2r;
      v2_53 +=  bc2 * wt2i;
      v2_44 += -bc2 * wt1r;
      v2_54 += -bc2 * wt1i;


      //mult_zp
      idir = 2;
      int it = site/(nx*ny*nz);
      int iz = (site/(nx*ny)) % nz;
      int ixy = site % (nx*ny);

      nn = (iz + 1) % nz;

      isn = ixy + nn*nx*ny + it*nx*ny*nz;
      isg = ixy + iz*nx*ny + it*nx*ny*nz;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] - v1[IDXV_SP_I(nst, 0,2,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] + v1[IDXV_SP_R(nst, 0,2,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] - v1[IDXV_SP_I(nst, 1,2,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] + v1[IDXV_SP_R(nst, 1,2,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] - v1[IDXV_SP_I(nst, 2,2,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] + v1[IDXV_SP_R(nst, 2,2,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] + v1[IDXV_SP_I(nst, 0,3,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] - v1[IDXV_SP_R(nst, 0,3,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] + v1[IDXV_SP_I(nst, 1,3,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] - v1[IDXV_SP_R(nst, 1,3,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] + v1[IDXV_SP_I(nst, 2,3,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] - v1[IDXV_SP_R(nst, 2,3,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 2,0,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 2,0,isg, idir)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //   ic = 0;
      bc2 = 1.0;
      if(iz == nz-1) bc2 = bc[2];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 +=  bc2 * wt1i;
      v2_13 += -bc2 * wt1r;
      v2_04 += -bc2 * wt2i;
      v2_14 +=  bc2 * wt2r;

      u_6 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 2,1,isg, idir)];
      u11 = u[IDXV_G_I(nst, 2,1,isg, idir)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
      v2_23 +=  bc2 * wt1i;
      v2_33 += -bc2 * wt1r;
      v2_24 += -bc2 * wt2i;
      v2_34 +=  bc2 * wt2r;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //   ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 +=  bc2 * wt1i;
      v2_53 += -bc2 * wt1r;
      v2_44 += -bc2 * wt2i;
      v2_54 +=  bc2 * wt2r;


      //mult_zm
      nn = (iz + nz - 1) % nz;

      isn = ixy + nn*nx*ny + it*nx*ny*nz;
      isg = ixy + nn*nx*ny + it*nx*ny*nz;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] + v1[IDXV_SP_I(nst, 0,2,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] - v1[IDXV_SP_R(nst, 0,2,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] + v1[IDXV_SP_I(nst, 1,2,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] - v1[IDXV_SP_R(nst, 1,2,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] + v1[IDXV_SP_I(nst, 2,2,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] - v1[IDXV_SP_R(nst, 2,2,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] - v1[IDXV_SP_I(nst, 0,3,isn)];
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] + v1[IDXV_SP_R(nst, 0,3,isn)];
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] - v1[IDXV_SP_I(nst, 1,3,isn)];
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] + v1[IDXV_SP_R(nst, 1,3,isn)];
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] - v1[IDXV_SP_I(nst, 2,3,isn)];
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] + v1[IDXV_SP_R(nst, 2,3,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 0,2,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 0,2,isg, idir)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //    ic = 0;
      bc2 = 1.0;
      if(iz == 0) bc2 = bc[2];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt1i;
      v2_13 +=  bc2 * wt1r;
      v2_04 +=  bc2 * wt2i;
      v2_14 += -bc2 * wt2r;

      u_6 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 1,2,isg, idir)];
      u11 = u[IDXV_G_I(nst, 1,2,isg, idir)];

      wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
      v2_23 += -bc2 * wt1i;
      v2_33 +=  bc2 * wt1r;
      v2_24 +=  bc2 * wt2i;
      v2_34 += -bc2 * wt2r;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //   ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 += -bc2 * wt1i;
      v2_53 +=  bc2 * wt1r;
      v2_44 +=  bc2 * wt2i;
      v2_54 += -bc2 * wt2r;


      //#include "afopr_Wilson_vector2_t_dirac-inc.h"

      //mult_tp
      idir = 3;
      nn = (it + 1) % nt;

      int ixyz = site % (nx*ny*nz);

      isn = ixyz + nn*nx*ny*nz;
      isg = ixyz + it*nx*ny*nz;

      vt1_0 = 2.0 * v1[IDXV_SP_R(nst, 0,2,isn)];
      vt1_1 = 2.0 * v1[IDXV_SP_I(nst, 0,2,isn)];
      vt1_2 = 2.0 * v1[IDXV_SP_R(nst, 1,2,isn)];
      vt1_3 = 2.0 * v1[IDXV_SP_I(nst, 1,2,isn)];
      vt1_4 = 2.0 * v1[IDXV_SP_R(nst, 2,2,isn)];
      vt1_5 = 2.0 * v1[IDXV_SP_I(nst, 2,2,isn)];

      vt2_0 = 2.0 * v1[IDXV_SP_R(nst, 0,3,isn)];
      vt2_1 = 2.0 * v1[IDXV_SP_I(nst, 0,3,isn)];
      vt2_2 = 2.0 * v1[IDXV_SP_R(nst, 1,3,isn)];
      vt2_3 = 2.0 * v1[IDXV_SP_I(nst, 1,3,isn)];
      vt2_4 = 2.0 * v1[IDXV_SP_R(nst, 2,3,isn)];
      vt2_5 = 2.0 * v1[IDXV_SP_I(nst, 2,3,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 2,0,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 2,0,isg, idir)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;
      bc2 = 1.0;
      if(it == nt-1) bc2 = bc[3];
      v2_03 +=  bc2 * wt1r;
      v2_13 +=  bc2 * wt1i;
      v2_04 +=  bc2 * wt2r;
      v2_14 +=  bc2 * wt2i;

      u_6 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 2,1,isg, idir)];
      u11 = u[IDXV_G_I(nst, 2,1,isg, idir)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_23 +=  bc2 * wt1r;
      v2_33 +=  bc2 * wt1i;
      v2_24 +=  bc2 * wt2r;
      v2_34 +=  bc2 * wt2i;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_43 +=  bc2 * wt1r;
      v2_53 +=  bc2 * wt1i;
      v2_44 +=  bc2 * wt2r;
      v2_54 +=  bc2 * wt2i;


      //mult_tm
      nn = (it + nt - 1) % nt;

      isn = ixyz + nn*nx*ny*nz;
      isg = ixyz + nn*nx*ny*nz;

      vt1_0 = 2.0 * v1[IDXV_SP_R(nst, 0,0,isn)];
      vt1_1 = 2.0 * v1[IDXV_SP_I(nst, 0,0,isn)];
      vt1_2 = 2.0 * v1[IDXV_SP_R(nst, 1,0,isn)];
      vt1_3 = 2.0 * v1[IDXV_SP_I(nst, 1,0,isn)];
      vt1_4 = 2.0 * v1[IDXV_SP_R(nst, 2,0,isn)];
      vt1_5 = 2.0 * v1[IDXV_SP_I(nst, 2,0,isn)];

      vt2_0 = 2.0 * v1[IDXV_SP_R(nst, 0,1,isn)];
      vt2_1 = 2.0 * v1[IDXV_SP_I(nst, 0,1,isn)];
      vt2_2 = 2.0 * v1[IDXV_SP_R(nst, 1,1,isn)];
      vt2_3 = 2.0 * v1[IDXV_SP_I(nst, 1,1,isn)];
      vt2_4 = 2.0 * v1[IDXV_SP_R(nst, 2,1,isn)];
      vt2_5 = 2.0 * v1[IDXV_SP_I(nst, 2,1,isn)];

      u_0 = u[IDXV_G_R(nst, 0,0,isg, idir)];
      u_1 = u[IDXV_G_I(nst, 0,0,isg, idir)];
      u_2 = u[IDXV_G_R(nst, 0,1,isg, idir)];
      u_3 = u[IDXV_G_I(nst, 0,1,isg, idir)];
      u_4 = u[IDXV_G_R(nst, 0,2,isg, idir)];
      u_5 = u[IDXV_G_I(nst, 0,2,isg, idir)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;
      bc2 = 1.0;
      if(it == 0) bc2 = bc[3];
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;

      u_6 = u[IDXV_G_R(nst, 1,0,isg, idir)];
      u_7 = u[IDXV_G_I(nst, 1,0,isg, idir)];
      u_8 = u[IDXV_G_R(nst, 1,1,isg, idir)];
      u_9 = u[IDXV_G_I(nst, 1,1,isg, idir)];
      u10 = u[IDXV_G_R(nst, 1,2,isg, idir)];
      u11 = u[IDXV_G_I(nst, 1,2,isg, idir)];

      wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;

      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);

      wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //   ic = 2;
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;

      //

      // aypx and write back to global memory
      v2[IDXV_SP_R(nst, 0,0,site)] = v1[IDXV_SP_R(nst, 0,0,site)] - kappa * v2_01;
      v2[IDXV_SP_I(nst, 0,0,site)] = v1[IDXV_SP_I(nst, 0,0,site)] - kappa * v2_11;
      v2[IDXV_SP_R(nst, 1,0,site)] = v1[IDXV_SP_R(nst, 1,0,site)] - kappa * v2_21;
      v2[IDXV_SP_I(nst, 1,0,site)] = v1[IDXV_SP_I(nst, 1,0,site)] - kappa * v2_31;
      v2[IDXV_SP_R(nst, 2,0,site)] = v1[IDXV_SP_R(nst, 2,0,site)] - kappa * v2_41;
      v2[IDXV_SP_I(nst, 2,0,site)] = v1[IDXV_SP_I(nst, 2,0,site)] - kappa * v2_51;

      v2[IDXV_SP_R(nst, 0,1,site)] = v1[IDXV_SP_R(nst, 0,1,site)] - kappa * v2_02;
      v2[IDXV_SP_I(nst, 0,1,site)] = v1[IDXV_SP_I(nst, 0,1,site)] - kappa * v2_12;
      v2[IDXV_SP_R(nst, 1,1,site)] = v1[IDXV_SP_R(nst, 1,1,site)] - kappa * v2_22;
      v2[IDXV_SP_I(nst, 1,1,site)] = v1[IDXV_SP_I(nst, 1,1,site)] - kappa * v2_32;
      v2[IDXV_SP_R(nst, 2,1,site)] = v1[IDXV_SP_R(nst, 2,1,site)] - kappa * v2_42;
      v2[IDXV_SP_I(nst, 2,1,site)] = v1[IDXV_SP_I(nst, 2,1,site)] - kappa * v2_52;

      v2[IDXV_SP_R(nst, 0,2,site)] = v1[IDXV_SP_R(nst, 0,2,site)] - kappa * v2_03;
      v2[IDXV_SP_I(nst, 0,2,site)] = v1[IDXV_SP_I(nst, 0,2,site)] - kappa * v2_13;
      v2[IDXV_SP_R(nst, 1,2,site)] = v1[IDXV_SP_R(nst, 1,2,site)] - kappa * v2_23;
      v2[IDXV_SP_I(nst, 1,2,site)] = v1[IDXV_SP_I(nst, 1,2,site)] - kappa * v2_33;
      v2[IDXV_SP_R(nst, 2,2,site)] = v1[IDXV_SP_R(nst, 2,2,site)] - kappa * v2_43;
      v2[IDXV_SP_I(nst, 2,2,site)] = v1[IDXV_SP_I(nst, 2,2,site)] - kappa * v2_53;

      v2[IDXV_SP_R(nst, 0,3,site)] = v1[IDXV_SP_R(nst, 0,3,site)] - kappa * v2_04;
      v2[IDXV_SP_I(nst, 0,3,site)] = v1[IDXV_SP_I(nst, 0,3,site)] - kappa * v2_14;
      v2[IDXV_SP_R(nst, 1,3,site)] = v1[IDXV_SP_R(nst, 1,3,site)] - kappa * v2_24;
      v2[IDXV_SP_I(nst, 1,3,site)] = v1[IDXV_SP_I(nst, 1,3,site)] - kappa * v2_34;
      v2[IDXV_SP_R(nst, 2,3,site)] = v1[IDXV_SP_R(nst, 2,3,site)] - kappa * v2_44;
      v2[IDXV_SP_I(nst, 2,3,site)] = v1[IDXV_SP_I(nst, 2,3,site)] - kappa * v2_54;

    }
  }

}

//====================================================================
void mult_wilson_D_chiral(real_t *RESTRICT v2, real_t *RESTRICT u,
                          real_t *RESTRICT v1,
                          int *Nsize, int *bc, real_t kappa)
{

  int nx = Nsize[0];
  int ny = Nsize[1];
  int nz = Nsize[2];
  int nt = Nsize[3];
  int nst  = nx*ny*nz*nt;
  int nvst = NVC * ND * nst;
  int ngst = NDF * nst * 4;


  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, nst);

  int nst2 = (ns-is)/VLEN;
  for(int ist2 = 0; ist2 < nst2; ++ist2){
    for(int site2 = 0; site2 < VLEN; ++site2){
      int site = is + site2 + VLEN * ist2;

      //    for(int site = 0; site < nst; ++site){

      real_t bc2;

      real_t u_0, u_1, u_2, u_3, u_4, u_5;
      real_t u_6, u_7, u_8, u_9, u10, u11;
      real_t u12, u13, u14, u15, u16, u17;
      real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
      real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
      real_t wt1r, wt1i, wt2r, wt2i;

      real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
      real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
      real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
      real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

#include "afopr_Wilson_vector2_xyz-inc.h"

#include "afopr_Wilson_vector2_t_chiral-inc.h"

      // aypx and write back to global memory
      v2[IDXV_SP_R(nst, 0,0,site)] = v1[IDXV_SP_R(nst, 0,0,site)] - kappa * v2_01;
      v2[IDXV_SP_I(nst, 0,0,site)] = v1[IDXV_SP_I(nst, 0,0,site)] - kappa * v2_11;
      v2[IDXV_SP_R(nst, 1,0,site)] = v1[IDXV_SP_R(nst, 1,0,site)] - kappa * v2_21;
      v2[IDXV_SP_I(nst, 1,0,site)] = v1[IDXV_SP_I(nst, 1,0,site)] - kappa * v2_31;
      v2[IDXV_SP_R(nst, 2,0,site)] = v1[IDXV_SP_R(nst, 2,0,site)] - kappa * v2_41;
      v2[IDXV_SP_I(nst, 2,0,site)] = v1[IDXV_SP_I(nst, 2,0,site)] - kappa * v2_51;

      v2[IDXV_SP_R(nst, 0,1,site)] = v1[IDXV_SP_R(nst, 0,1,site)] - kappa * v2_02;
      v2[IDXV_SP_I(nst, 0,1,site)] = v1[IDXV_SP_I(nst, 0,1,site)] - kappa * v2_12;
      v2[IDXV_SP_R(nst, 1,1,site)] = v1[IDXV_SP_R(nst, 1,1,site)] - kappa * v2_22;
      v2[IDXV_SP_I(nst, 1,1,site)] = v1[IDXV_SP_I(nst, 1,1,site)] - kappa * v2_32;
      v2[IDXV_SP_R(nst, 2,1,site)] = v1[IDXV_SP_R(nst, 2,1,site)] - kappa * v2_42;
      v2[IDXV_SP_I(nst, 2,1,site)] = v1[IDXV_SP_I(nst, 2,1,site)] - kappa * v2_52;

      v2[IDXV_SP_R(nst, 0,2,site)] = v1[IDXV_SP_R(nst, 0,2,site)] - kappa * v2_03;
      v2[IDXV_SP_I(nst, 0,2,site)] = v1[IDXV_SP_I(nst, 0,2,site)] - kappa * v2_13;
      v2[IDXV_SP_R(nst, 1,2,site)] = v1[IDXV_SP_R(nst, 1,2,site)] - kappa * v2_23;
      v2[IDXV_SP_I(nst, 1,2,site)] = v1[IDXV_SP_I(nst, 1,2,site)] - kappa * v2_33;
      v2[IDXV_SP_R(nst, 2,2,site)] = v1[IDXV_SP_R(nst, 2,2,site)] - kappa * v2_43;
      v2[IDXV_SP_I(nst, 2,2,site)] = v1[IDXV_SP_I(nst, 2,2,site)] - kappa * v2_53;

      v2[IDXV_SP_R(nst, 0,3,site)] = v1[IDXV_SP_R(nst, 0,3,site)] - kappa * v2_04;
      v2[IDXV_SP_I(nst, 0,3,site)] = v1[IDXV_SP_I(nst, 0,3,site)] - kappa * v2_14;
      v2[IDXV_SP_R(nst, 1,3,site)] = v1[IDXV_SP_R(nst, 1,3,site)] - kappa * v2_24;
      v2[IDXV_SP_I(nst, 1,3,site)] = v1[IDXV_SP_I(nst, 1,3,site)] - kappa * v2_34;
      v2[IDXV_SP_R(nst, 2,3,site)] = v1[IDXV_SP_R(nst, 2,3,site)] - kappa * v2_44;
      v2[IDXV_SP_I(nst, 2,3,site)] = v1[IDXV_SP_I(nst, 2,3,site)] - kappa * v2_54;

   }

  }

}

//====================================================================
void mult_wilson_1_dirac(
                    real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                    real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                    real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                    real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                    real_t *RESTRICT u, real_t *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;


  if(do_comm[0] > 0){

    int idir = 0;
    int Nyzt = Ny * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nyzt);

    //    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    for(int iyzt = is; iyzt < ns; ++iyzt){
      int ix  = 0;
      int ist = ix + Nx * iyzt;
      real_t bc2 = bc[0];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        real_t vt1[2], vt2[2];
        vt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
        vt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
        vt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
        vt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
        buf_xp[IDXV_2SP_R(Nyzt, ic, 0, iyzt)] = bc2 * vt1[0];
        buf_xp[IDXV_2SP_I(Nyzt, ic, 0, iyzt)] = bc2 * vt1[1];
        buf_xp[IDXV_2SP_R(Nyzt, ic, 1, iyzt)] = bc2 * vt2[0];
        buf_xp[IDXV_2SP_I(Nyzt, ic, 1, iyzt)] = bc2 * vt2[1];
      }
    }

    //    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    for(int iyzt = is; iyzt < ns; ++iyzt){
      int ix  = Nx-1;
      int ist = ix + Nx * iyzt;

      real_t vt1[NVC], vt2[NVC], ut[NVC];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
        vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
        vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
        vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
      }

      for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, idir)];
          ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, idir)];
        }

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        buf_xm[IDXV_2SP_R(Nyzt, ic, 0, iyzt)] = wt1r;
        buf_xm[IDXV_2SP_I(Nyzt, ic, 0, iyzt)] = wt1i;
        buf_xm[IDXV_2SP_R(Nyzt, ic, 1, iyzt)] = wt2r;
        buf_xm[IDXV_2SP_I(Nyzt, ic, 1, iyzt)] = wt2i;

      }

    }

 } // do_comm[0]


 // idir = 1;
  if(do_comm[1] > 0){

    int idir = 1;
    int Nzt = Nz * Nt;
    int Nxzt = Nx * Nzt;

    int ith, nth, is, ns;
    set_threadtask_mult(ith, nth, is, ns, Nxzt);

    //    for(int izt = 0; izt < Nzt; ++izt){
    //     for(int ix = 0; ix < Nx; ++ix){
    for(int ixzt = is; ixzt < ns; ++ixzt){
     {
       int ix  = ixzt % Nx;
       int izt = ixzt/Nx;

       int iy   = 0;
       int ist  = ix + Nx * (iy + Ny * izt);
       int ixzt = ix + Nx * izt;
       real_t bc2 = bc[1];

#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){
         real_t vt1[2], vt2[2];
         vt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
         vt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
         vt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
         vt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
         buf_yp[IDXV_2SP_R(Nxzt, ic, 0, ixzt)] = bc2 * vt1[0];
         buf_yp[IDXV_2SP_I(Nxzt, ic, 0, ixzt)] = bc2 * vt1[1];
         buf_yp[IDXV_2SP_R(Nxzt, ic, 1, ixzt)] = bc2 * vt2[0];
         buf_yp[IDXV_2SP_I(Nxzt, ic, 1, ixzt)] = bc2 * vt2[1];
       }

     }
    }

    //    for(int izt = 0; izt < Nzt; ++izt){
    //     for(int ix = 0; ix < Nx; ++ix){
    for(int ixzt = is; ixzt < ns; ++ixzt){
     {
       int ix  = ixzt % Nx;
       int izt = ixzt/Nx;

       int iy   = Ny-1;
       int ist  = ix + Nx * (iy + Ny*izt);
       int ixzt = ix + Nx * izt;

       real_t vt1[NVC], vt2[NVC];
#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){
         int icr = 2*ic;
         int ici = 2*ic + 1;
         vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
         vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
         vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
         vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
       }

#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){

         real_t ut[NVC];

#pragma _NEC unroll(NC)
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

         buf_ym[IDXV_2SP_R(Nxzt, ic, 0, ixzt)] = wt1[0];
         buf_ym[IDXV_2SP_I(Nxzt, ic, 0, ixzt)] = wt1[1];
         buf_ym[IDXV_2SP_R(Nxzt, ic, 1, ixzt)] = wt2[0];
         buf_ym[IDXV_2SP_I(Nxzt, ic, 1, ixzt)] = wt2[1];

       }

     }
    }

  } // do_comm[1]

  //  idir = 2;
  if(do_comm[2] > 0){

    int idir = 2;
    int Nxy  = Nx * Ny;
    int Nxyt = Nxy * Nt;

    for(int it = 0; it < Nt; ++it){
     for(int ixy = 0; ixy < Nxy; ++ixy){
       int iz   = 0;
       int ist  = ixy + Nxy * (iz + Nz * it);
       int ixyt = ixy + Nxy * it;
       real_t bc2 = bc[2];

#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){
         real_t wt1[2], wt2[2];
         wt1[0] = v1[IDXV_SP_R(Nst, ic, 0, ist)] - v1[IDXV_SP_I(Nst, ic, 2, ist)];
         wt1[1] = v1[IDXV_SP_I(Nst, ic, 0, ist)] + v1[IDXV_SP_R(Nst, ic, 2, ist)];
         wt2[0] = v1[IDXV_SP_R(Nst, ic, 1, ist)] + v1[IDXV_SP_I(Nst, ic, 3, ist)];
         wt2[1] = v1[IDXV_SP_I(Nst, ic, 1, ist)] - v1[IDXV_SP_R(Nst, ic, 3, ist)];
         buf_zp[IDXV_2SP_R(Nxyt, ic, 0, ixyt)] = bc2 * wt1[0];
         buf_zp[IDXV_2SP_I(Nxyt, ic, 0, ixyt)] = bc2 * wt1[1];
         buf_zp[IDXV_2SP_R(Nxyt, ic, 1, ixyt)] = bc2 * wt2[0];
         buf_zp[IDXV_2SP_I(Nxyt, ic, 1, ixyt)] = bc2 * wt2[1];
       }

     }
    }

    for(int it = 0; it < Nt; ++it){
     for(int ixy = 0; ixy < Nxy; ++ixy){
       int iz = Nz-1;
       int ist  = ixy + Nxy * (iz + Nz * it);
       int ixyt = ixy + Nxy * it;

       real_t vt1[NVC], vt2[NVC];

#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){
         int icr = 2*ic;
         int ici = 2*ic + 1;
         vt1[icr] = v1[IDXV_SP_R(Nst, ic, 0, ist)] + v1[IDXV_SP_I(Nst, ic, 2, ist)];
         vt1[ici] = v1[IDXV_SP_I(Nst, ic, 0, ist)] - v1[IDXV_SP_R(Nst, ic, 2, ist)];
         vt2[icr] = v1[IDXV_SP_R(Nst, ic, 1, ist)] - v1[IDXV_SP_I(Nst, ic, 3, ist)];
         vt2[ici] = v1[IDXV_SP_I(Nst, ic, 1, ist)] + v1[IDXV_SP_R(Nst, ic, 3, ist)];
       }

#pragma _NEC unroll(NC)
       for(int ic = 0; ic < NC; ++ic){

         real_t ut[NVC];
#pragma _NEC unroll(NC)
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
         buf_zm[IDXV_2SP_R(Nxyt, ic, 0, ixyt)] = wt1[0];
         buf_zm[IDXV_2SP_I(Nxyt, ic, 0, ixyt)] = wt1[1];
         buf_zm[IDXV_2SP_R(Nxyt, ic, 1, ixyt)] = wt2[0];
         buf_zm[IDXV_2SP_I(Nxyt, ic, 1, ixyt)] = wt2[1];

       }

     }
    }

  } // do_comm[2]


  //  idir = 3;
  if(do_comm[3] > 0){

    int idir = 3;
    int Nxyz = Nx * Ny * Nz;

    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it = 0;
      int ist  = ixyz + Nxyz * it;
      real_t bc2 = 2.0 * bc[3];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        buf_tp[IDXV_2SP_R(Nxyz, ic, 0, ixyz)] = bc2 * v1[IDXV_SP_R(Nst, ic, 2, ist)];
        buf_tp[IDXV_2SP_I(Nxyz, ic, 0, ixyz)] = bc2 * v1[IDXV_SP_I(Nst, ic, 2, ist)];
        buf_tp[IDXV_2SP_R(Nxyz, ic, 1, ixyz)] = bc2 * v1[IDXV_SP_R(Nst, ic, 3, ist)];
        buf_tp[IDXV_2SP_I(Nxyz, ic, 1, ixyz)] = bc2 * v1[IDXV_SP_I(Nst, ic, 3, ist)];
      }
    }

    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it   = Nt-1;
      int ist  = ixyz + Nxyz * it;

      real_t vt1[NVC], vt2[NVC];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        int icr = 2*ic;
        int ici = 2*ic + 1;
        vt1[icr] = 2.0 * v1[IDXV_SP_R(Nst, ic, 0, ist)];
        vt1[ici] = 2.0 * v1[IDXV_SP_I(Nst, ic, 0, ist)];
        vt2[icr] = 2.0 * v1[IDXV_SP_R(Nst, ic, 1, ist)];
        vt2[ici] = 2.0 * v1[IDXV_SP_I(Nst, ic, 1, ist)];
      }

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

        real_t ut[NVC];
#pragma _NEC unroll(NC)
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
        buf_tm[IDXV_2SP_R(Nxyz, ic, 0, ixyz)] = wt1[0];
        buf_tm[IDXV_2SP_I(Nxyz, ic, 0, ixyz)] = wt1[1];
        buf_tm[IDXV_2SP_R(Nxyz, ic, 1, ixyz)] = wt2[0];
        buf_tm[IDXV_2SP_I(Nxyz, ic, 1, ixyz)] = wt2[1];
      }

    }

  } // do_comm[3]

}

//====================================================================
void mult_wilson_1_chiral(
                    real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                    real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                    real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                    real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                    real_t *RESTRICT u, real_t *RESTRICT v1,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * Ny * Nz * Nt;
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;


#include "afopr_Wilson_vector2_1xyz-inc.h"


 //  idir = 3;
 if(do_comm[3] > 0){

   int Nxyz = Nx * Ny * Nz;
   int Nbf = Nx * Ny * Nz;

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;
    real_t bc2 = bc[3];
    for(int ivc = 0; ivc < NVC; ++ivc){
      buf_tp[IDXV_2SP(Nbf, ivc, 0, ixyz)]
       = bc2 * (v1[IDXV_SP(Nst, ivc, 0, ist)] + v1[IDXV_SP(Nst, ivc, 2, ist)]);
      buf_tp[IDXV_2SP(Nbf, ivc, 1, ixyz)]
       = bc2 * (v1[IDXV_SP(Nst, ivc, 1, ist)] + v1[IDXV_SP(Nst, ivc, 3, ist)]);
    }
  }

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;
    int Nbf = Nx * Ny * Nz;
    int idir = 3;

    real_t vt1[NVC], vt2[NVC];
    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = v1[IDXV_SP(Nst, ivc, 0, ist)] - v1[IDXV_SP(Nst, ivc, 2, ist)];
      vt2[ivc] = v1[IDXV_SP(Nst, ivc, 1, ist)] - v1[IDXV_SP(Nst, ivc, 3, ist)];
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
      buf_tm[IDXV_2SP_R(Nbf, ic, 0, ixyz)] = wt1[0];
      buf_tm[IDXV_2SP_I(Nbf, ic, 0, ixyz)] = wt1[1];
      buf_tm[IDXV_2SP_R(Nbf, ic, 1, ixyz)] = wt2[0];
      buf_tm[IDXV_2SP_I(Nbf, ic, 1, ixyz)] = wt2[1];

    }

  }

 } // do_comm[3]


}

//====================================================================
void mult_wilson_2_dirac(real_t *RESTRICT v2, real_t *RESTRICT u, 
                    real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                    real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                    real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                    real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                    real_t kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;


  if(do_comm[0] > 0){
    int idir = 0;
    int Nyzt  = Ny * Nz * Nt;

    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix = Nx-1;
      int ist = ix + Nx * iyzt;

      real_t vt1[NVC], vt2[NVC], ut[NVC];

#pragma _NEC unroll(NVC)
      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_xp[IDXV_2SP(Nyzt, ivc, 0, iyzt)];
        vt2[ivc] = buf_xp[IDXV_2SP(Nyzt, ivc, 1, iyzt)];
      }

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, idir)];
          ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, idir)];
        }

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2[IDXV_SP_R(Nst, ic, 0, ist)] += -kappa * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] += -kappa * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] += -kappa * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] += -kappa * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] += -kappa * wt2i;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  kappa * wt2r;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] += -kappa * wt1i;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  kappa * wt1r;

      }

    }

    for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
      int ix = 0;
      int ist = ix + Nx * iyzt;
      real_t bc2 = - kappa * bc[0];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = buf_xm[IDXV_2SP_R(Nyzt, ic, 0, iyzt)];
        wt1i = buf_xm[IDXV_2SP_I(Nyzt, ic, 0, iyzt)];
        wt2r = buf_xm[IDXV_2SP_R(Nyzt, ic, 1, iyzt)];
        wt2i = buf_xm[IDXV_2SP_I(Nyzt, ic, 1, iyzt)];

        v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] += -bc2 * wt2i;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt2r;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] += -bc2 * wt1i;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  bc2 * wt1r;

      }

    }

  } // do_comm[0]                                                              


  if(do_comm[1] > 0){
    int idir = 1;
    int Nzt  = Nz * Nt;
    int Nxzt = Nx * Nzt;

    for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
      int iy = Ny-1;
      int ix = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist = ix + Nx *(iy + Ny * izt);

      real_t vt1[NVC], vt2[NVC], ut[NVC];

#pragma _NEC unroll(NVC)
      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_yp[IDXV_2SP(Nxzt, ivc, 0, ixzt)];
        vt2[ivc] = buf_yp[IDXV_2SP(Nxzt, ivc, 1, ixzt)];
      }

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, idir)];
          ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, idir)];
        }

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2[IDXV_SP_R(Nst, ic, 0, ist)] += -kappa * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] += -kappa * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] += -kappa * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] += -kappa * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  kappa * wt2r;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  kappa * wt2i;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] += -kappa * wt1r;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] += -kappa * wt1i;
      }

    }

    for(int ixzt = 0; ixzt < Nxzt; ++ixzt){
      int iy = 0;
      int ix = ixzt % Nx;
      int izt = ixzt/Nx;
      int ist = ix + Nx *(iy + Ny * izt);
      real_t bc2 = -kappa * bc[1];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = buf_ym[IDXV_2SP_R(Nxzt, ic, 0, ixzt)];
        wt1i = buf_ym[IDXV_2SP_I(Nxzt, ic, 0, ixzt)];
        wt2r = buf_ym[IDXV_2SP_R(Nxzt, ic, 1, ixzt)];
        wt2i = buf_ym[IDXV_2SP_I(Nxzt, ic, 1, ixzt)];

        v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  bc2 * wt2r;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt2i;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] += -bc2 * wt1r;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] += -bc2 * wt1i;
      }

    }

  } // do_comm[1]


  if(do_comm[2] > 0){

    int idir = 2;
    int Nxy  = Nx * Ny;
    int Nxyt = Nxy * Nt;

    for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
      int iz = Nz-1;
      int ixy = ixyt % Nxy;
      int it = ixyt/Nxy;
      int ist = ixy + Nxy *(iz + Nz * it);

      real_t vt1[NVC], vt2[NVC], ut[NVC];

#pragma _NEC unroll(NVC)
      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_zp[IDXV_2SP(Nxyt, ivc, 0, ixyt)];
        vt2[ivc] = buf_zp[IDXV_2SP(Nxyt, ivc, 1, ixyt)];
      }

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, idir)];
          ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, idir)];
        }

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2[IDXV_SP_R(Nst, ic, 0, ist)] += -kappa * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] += -kappa * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] += -kappa * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] += -kappa * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] += -kappa * wt1i;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  kappa * wt1r;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  kappa * wt2i;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] += -kappa * wt2r;

      }

    }

    for(int ixyt = 0; ixyt < Nxyt; ++ixyt){
      int iz = 0;
      int ixy = ixyt % Nxy;
      int it = ixyt/Nxy;
      int ist = ixy + Nxy *(iz + Nz * it);
      real_t bc2 = -kappa * bc[2];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){
        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = buf_zm[IDXV_2SP_R(Nxyt, ic, 0, ixyt)];
        wt1i = buf_zm[IDXV_2SP_I(Nxyt, ic, 0, ixyt)];
        wt2r = buf_zm[IDXV_2SP_R(Nxyt, ic, 1, ixyt)];
        wt2i = buf_zm[IDXV_2SP_I(Nxyt, ic, 1, ixyt)];

        v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2i;
        v2[IDXV_SP_R(Nst, ic, 2, ist)] += -bc2 * wt1i;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt1r;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  bc2 * wt2i;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] += -bc2 * wt2r;
      }

    }

  } // do_comm[2]


  if(do_comm[3] > 0){

    int idir = 3;
    int Nxyz = Nx * Ny * Nz;

    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it = Nt-1;
      int ist = ixyz + Nxyz * it;

      real_t vt1[NVC], vt2[NVC], ut[NVC];

#pragma _NEC unroll(NVC)
      for(int ivc = 0; ivc < NVC; ++ivc){
        vt1[ivc] = buf_tp[IDXV_2SP(Nxyz, ivc, 0, ixyz)];
        vt2[ivc] = buf_tp[IDXV_2SP(Nxyz, ivc, 1, ixyz)];
      }

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

#pragma _NEC unroll(NC)
        for(int ic2 = 0; ic2 < NC; ++ic2){
          ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, idir)];
          ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, idir)];
        }

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt1i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
        wt2r = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
        wt2i = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                        vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

        v2[IDXV_SP_R(Nst, ic, 2, ist)] += -kappa * wt1r;
        v2[IDXV_SP_I(Nst, ic, 2, ist)] += -kappa * wt1i;
        v2[IDXV_SP_R(Nst, ic, 3, ist)] += -kappa * wt2r;
        v2[IDXV_SP_I(Nst, ic, 3, ist)] += -kappa * wt2i;

      }

    }

    for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
      int it = 0;
      int ist = ixyz + Nxyz * it;
      real_t bc2 = -kappa * bc[3];

#pragma _NEC unroll(NC)
      for(int ic = 0; ic < NC; ++ic){

        real_t wt1r, wt1i, wt2r, wt2i;
        wt1r = buf_tm[IDXV_2SP_R(Nxyz, ic, 0, ixyz)];
        wt1i = buf_tm[IDXV_2SP_I(Nxyz, ic, 0, ixyz)];
        wt2r = buf_tm[IDXV_2SP_R(Nxyz, ic, 1, ixyz)];
        wt2i = buf_tm[IDXV_2SP_I(Nxyz, ic, 1, ixyz)];

        v2[IDXV_SP_R(Nst, ic, 0, ist)] += bc2 * wt1r;
        v2[IDXV_SP_I(Nst, ic, 0, ist)] += bc2 * wt1i;
        v2[IDXV_SP_R(Nst, ic, 1, ist)] += bc2 * wt2r;
        v2[IDXV_SP_I(Nst, ic, 1, ist)] += bc2 * wt2i;
      }

    }

  } // do_comm[3]


}

//==================================================================== 
void mult_wilson_2_chiral(real_t *RESTRICT v2, real_t *RESTRICT u,
                    real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                    real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                    real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                    real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                    real_t kappa,
                    int *Nsize, int *bc, int *do_comm, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  int size    = NVC * ND * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = NVC * 2 * Ny * Nz * Nt;
  int size_by = NVC * 2 * Nx * Nz * Nt;
  int size_bz = NVC * 2 * Nx * Ny * Nt;
  int size_bt = NVC * 2 * Nx * Ny * Nz;

  for(int ist = 0; ist < Nst; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx*Ny)) % Nz;
    int it = ist/(Nx*Ny*Nz);

    real_t v2L[NVC*ND];
    for(int id = 0; id < ND; ++id){
      for(int ivc = 0; ivc < NVC; ++ivc){
        v2L[ivc + NVC * id] = 0.0;
      }
    }

#include "afopr_Wilson_vector2_2xyz-inc.h"

    // idir = 3
    if(do_comm[3] > 0){

      if(it == Nt-1){
        int ixyz = ix + Nx * (iy + Ny * iz);
        int Nbf = Nx * Ny * Nz;
        int idir = 3;

        real_t vt1[NVC], vt2[NVC], ut[NVC];
        real_t wt1[2], wt2[2];

        for(int ivc = 0; ivc < NVC; ++ivc){
          vt1[ivc] = buf_tp[IDXV_2SP(Nbf, ivc, 0, ixyz)];
          vt2[ivc] = buf_tp[IDXV_2SP(Nbf, ivc, 1, ixyz)];
        }

        for(int ic = 0; ic < NC; ++ic){

          for(int ic2 = 0; ic2 < NC; ++ic2){
            ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, idir)];
            ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, idir)];
          }

          wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
          wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
          wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

          v2L[ic*2   + ID1] += wt1[0];
          v2L[ic*2+1 + ID1] += wt1[1];
          v2L[ic*2   + ID2] += wt2[0];
          v2L[ic*2+1 + ID2] += wt2[1];
          v2L[ic*2   + ID3] += wt1[0];
          v2L[ic*2+1 + ID3] += wt1[1];
          v2L[ic*2   + ID4] += wt2[0];
          v2L[ic*2+1 + ID4] += wt2[1];
	}
      }

      if(it == 0){
        int ixyz = ix + Nx * (iy + Ny * iz);
        int Nbf = Nx * Ny * Nz;
        real_t bc2 = bc[3];
        real_t wt1[2], wt2[2];

        for(int ic = 0; ic < NC; ++ic){
          wt1[0] = bc2 * buf_tm[IDXV_2SP_R(Nbf, ic, 0, ixyz)];
          wt1[1] = bc2 * buf_tm[IDXV_2SP_I(Nbf, ic, 0, ixyz)];
          wt2[0] = bc2 * buf_tm[IDXV_2SP_R(Nbf, ic, 1, ixyz)];
          wt2[1] = bc2 * buf_tm[IDXV_2SP_I(Nbf, ic, 1, ixyz)];
          v2L[ic*2   + ID1] += wt1[0];
          v2L[ic*2+1 + ID1] += wt1[1];
          v2L[ic*2   + ID2] += wt2[0];
          v2L[ic*2+1 + ID2] += wt2[1];
          v2L[ic*2   + ID3] -= wt1[0];
          v2L[ic*2+1 + ID3] -= wt1[1];
          v2L[ic*2   + ID4] -= wt2[0];
          v2L[ic*2+1 + ID4] -= wt2[1];
	}
      }

    }

    //axpy
    for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        v2[IDXV_SP_R(Nst, ic, id, ist)] += -kappa * v2L[ic*2   + NVC*id];
        v2[IDXV_SP_I(Nst, ic, id, ist)] += -kappa * v2L[ic*2+1 + NVC*id];
      }
    }

  }

}

//============================================================END=====
