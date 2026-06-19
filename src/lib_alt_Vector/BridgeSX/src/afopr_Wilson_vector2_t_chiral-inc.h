/*!
      @file    afopr_Wilson_vector2_t_chiral.inc
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/


      //mult_tp
      idir = 3;
      nn = (it + 1) % nt;

      int ixyz = site % (nx*ny*nz);

      isn = ixyz + nn*nx*ny*nz;
      isg = ixyz + it*nx*ny*nz;

      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] + v1[IDXV_SP_R(nst, 0,2,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] + v1[IDXV_SP_I(nst, 0,2,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] + v1[IDXV_SP_R(nst, 1,2,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] + v1[IDXV_SP_I(nst, 1,2,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] + v1[IDXV_SP_R(nst, 2,2,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] + v1[IDXV_SP_I(nst, 2,2,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] + v1[IDXV_SP_R(nst, 0,3,isn)] ;
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] + v1[IDXV_SP_I(nst, 0,3,isn)] ;
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] + v1[IDXV_SP_R(nst, 1,3,isn)] ;
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] + v1[IDXV_SP_I(nst, 1,3,isn)] ;
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] + v1[IDXV_SP_R(nst, 2,3,isn)] ;
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] + v1[IDXV_SP_I(nst, 2,3,isn)] ;

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
      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
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
      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;
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
      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      v2_43 +=  bc2 * wt1r;
      v2_53 +=  bc2 * wt1i;
      v2_44 +=  bc2 * wt2r;
      v2_54 +=  bc2 * wt2i;


      //mult_tm
      nn = (it + nt - 1) % nt;

      isn = ixyz + nn*nx*ny*nz;
      isg = ixyz + nn*nx*ny*nz;


//qqqqq
      vt1_0 = v1[IDXV_SP_R(nst, 0,0,isn)] - v1[IDXV_SP_R(nst, 0,2,isn)];
      vt1_1 = v1[IDXV_SP_I(nst, 0,0,isn)] - v1[IDXV_SP_I(nst, 0,2,isn)];
      vt1_2 = v1[IDXV_SP_R(nst, 1,0,isn)] - v1[IDXV_SP_R(nst, 1,2,isn)];
      vt1_3 = v1[IDXV_SP_I(nst, 1,0,isn)] - v1[IDXV_SP_I(nst, 1,2,isn)];
      vt1_4 = v1[IDXV_SP_R(nst, 2,0,isn)] - v1[IDXV_SP_R(nst, 2,2,isn)];
      vt1_5 = v1[IDXV_SP_I(nst, 2,0,isn)] - v1[IDXV_SP_I(nst, 2,2,isn)];

      vt2_0 = v1[IDXV_SP_R(nst, 0,1,isn)] - v1[IDXV_SP_R(nst, 0,3,isn)] ;
      vt2_1 = v1[IDXV_SP_I(nst, 0,1,isn)] - v1[IDXV_SP_I(nst, 0,3,isn)] ;
      vt2_2 = v1[IDXV_SP_R(nst, 1,1,isn)] - v1[IDXV_SP_R(nst, 1,3,isn)] ;
      vt2_3 = v1[IDXV_SP_I(nst, 1,1,isn)] - v1[IDXV_SP_I(nst, 1,3,isn)] ;
      vt2_4 = v1[IDXV_SP_R(nst, 2,1,isn)] - v1[IDXV_SP_R(nst, 2,3,isn)] ;
      vt2_5 = v1[IDXV_SP_I(nst, 2,1,isn)] - v1[IDXV_SP_I(nst, 2,3,isn)] ;

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
      v2_03 += -bc2 * wt1r;
      v2_13 += -bc2 * wt1i;
      v2_04 += -bc2 * wt2r;
      v2_14 += -bc2 * wt2i;

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
      v2_23 += -bc2 * wt1r;
      v2_33 += -bc2 * wt1i;
      v2_24 += -bc2 * wt2r;
      v2_34 += -bc2 * wt2i;

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
      v2_43 += -bc2 * wt1r;
      v2_53 += -bc2 * wt1i;
      v2_44 += -bc2 * wt2r;
      v2_54 += -bc2 * wt2i;


//============================================================END=====
