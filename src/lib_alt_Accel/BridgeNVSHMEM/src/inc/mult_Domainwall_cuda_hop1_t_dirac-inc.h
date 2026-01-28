/*!
      @file    mult_Domainwall_cuda_hopb-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-05-08 13:51:53 #$
      @version $LastChangedRevision: 2606 $
*/

      int nx = Nx;
      int ny = Ny;
      int nz = Nz;
      int nt = Nt;
      int nst = nx * ny * nz * Nt;
      idir = 0;

{
      // mult_xp
      if(( ix < nx - 1 ) || (do_comm_x == 0)){

      int ix2  = (ix + 1) % nx; 
      int isn  = ix2  + nx  * iyzt;
      int isg  = site + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];

      u_0 = u_up[IDX2_G_R(0,0,isg)];
      u_1 = u_up[IDX2_G_I(0,0,isg)];
      u_2 = u_up[IDX2_G_R(1,0,isg)];
      u_3 = u_up[IDX2_G_I(1,0,isg)];
      u_4 = u_up[IDX2_G_R(2,0,isg)];
      u_5 = u_up[IDX2_G_I(2,0,isg)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;

      real_t bc2 = ( ix == nx - 1 ) ? bc_x : 1.0;

      v2_01 +=   bc2 * wt1r;
      v2_11 +=   bc2 * wt1i;
      v2_02 +=   bc2 * wt2r;
      v2_12 +=   bc2 * wt2i;
      v2_03 +=   bc2 * wt2i;
      v2_13 +=  -bc2 * wt2r;
      v2_04 +=   bc2 * wt1i;
      v2_14 +=  -bc2 * wt1r;

      u_6 = u_up[IDX2_G_R(0,1,isg)];
      u_7 = u_up[IDX2_G_I(0,1,isg)];
      u_8 = u_up[IDX2_G_R(1,1,isg)];
      u_9 = u_up[IDX2_G_I(1,1,isg)];
      u10 = u_up[IDX2_G_R(2,1,isg)];
      u11 = u_up[IDX2_G_I(2,1,isg)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 1;
      v2_21 +=   bc2 * wt1r;
      v2_31 +=   bc2 * wt1i;
      v2_22 +=   bc2 * wt2r;
      v2_32 +=   bc2 * wt2i;
      v2_23 +=   bc2 * wt2i;
      v2_33 +=  -bc2 * wt2r;
      v2_24 +=   bc2 * wt1i;
      v2_34 +=  -bc2 * wt1r;

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_up[IDX2_G_R(0,2,isg)];
      u13 = u_up[IDX2_G_I(0,2,isg)];
      u14 = u_up[IDX2_G_R(1,2,isg)];
      u15 = u_up[IDX2_G_I(1,2,isg)];
      u16 = u_up[IDX2_G_R(2,2,isg)];
      u17 = u_up[IDX2_G_I(2,2,isg)];
#endif

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 2;
      v2_41 +=   bc2 * wt1r;
      v2_51 +=   bc2 * wt1i;
      v2_42 +=   bc2 * wt2r;
      v2_52 +=   bc2 * wt2i;
      v2_43 +=   bc2 * wt2i;
      v2_53 +=  -bc2 * wt2r;
      v2_44 +=   bc2 * wt1i;
      v2_54 +=  -bc2 * wt1r;
      } 

      // mult_xm
      if(( ix > 0 ) || (do_comm_x == 0)){

      int ix2 = (ix - 1 + nx) % nx;
      int isn = ix2 + nx  * iyzt; //nei
      int isg = isn + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];

      u_0 = u_dn[IDX2_G_R(0,0,isg)];
      u_1 = u_dn[IDX2_G_I(0,0,isg)];
      u_2 = u_dn[IDX2_G_R(0,1,isg)];
      u_3 = u_dn[IDX2_G_I(0,1,isg)];
      u_4 = u_dn[IDX2_G_R(0,2,isg)];
      u_5 = u_dn[IDX2_G_I(0,2,isg)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      // ic = 0;
      real_t bc2 = ( ix == 0 ) ? bc_x : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt2i;
      v2_13 += +bc2 * wt2r;
      v2_04 += -bc2 * wt1i;
      v2_14 += +bc2 * wt1r;

      u_6 = u_dn[IDX2_G_R(1,0,isg)];
      u_7 = u_dn[IDX2_G_I(1,0,isg)];
      u_8 = u_dn[IDX2_G_R(1,1,isg)];
      u_9 = u_dn[IDX2_G_I(1,1,isg)];
      u10 = u_dn[IDX2_G_R(1,2,isg)];
      u11 = u_dn[IDX2_G_I(1,2,isg)];

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

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_dn[IDX2_G_R(2,0,isg)];
      u13 = u_dn[IDX2_G_I(2,0,isg)];
      u14 = u_dn[IDX2_G_R(2,1,isg)];
      u15 = u_dn[IDX2_G_I(2,1,isg)];
      u16 = u_dn[IDX2_G_R(2,2,isg)];
      u17 = u_dn[IDX2_G_I(2,2,isg)];
#endif

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
      }

      // mult_yp
      idir = 1;

      if(( iy < ny -1 ) || (do_comm_y == 0)){

      int nn = (iy + 1) % ny;
      int isn = ix + nx * (nn + izt*ny);
      int isg = site + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];

      u_0 = u_up[IDX2_G_R(0,0,isg)];
      u_1 = u_up[IDX2_G_I(0,0,isg)];
      u_2 = u_up[IDX2_G_R(1,0,isg)];
      u_3 = u_up[IDX2_G_I(1,0,isg)];
      u_4 = u_up[IDX2_G_R(2,0,isg)];
      u_5 = u_up[IDX2_G_I(2,0,isg)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //  ic = 0;
      real_t bc2 = ( iy == ny - 1 ) ? bc_y : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt2r;
      v2_13 += -bc2 * wt2i;
      v2_04 +=  bc2 * wt1r;
      v2_14 +=  bc2 * wt1i;

      u_6 = u_up[IDX2_G_R(0,1,isg)];
      u_7 = u_up[IDX2_G_I(0,1,isg)];
      u_8 = u_up[IDX2_G_R(1,1,isg)];
      u_9 = u_up[IDX2_G_I(1,1,isg)];
      u10 = u_up[IDX2_G_R(2,1,isg)];
      u11 = u_up[IDX2_G_I(2,1,isg)];

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

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_up[IDX2_G_R(0,2,isg)];
      u13 = u_up[IDX2_G_I(0,2,isg)];
      u14 = u_up[IDX2_G_R(1,2,isg)];
      u15 = u_up[IDX2_G_I(1,2,isg)];
      u16 = u_up[IDX2_G_R(2,2,isg)];
      u17 = u_up[IDX2_G_I(2,2,isg)];
#endif

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
      }

      //mult_ym
      if(( iy > 0 ) || (do_comm_y == 0)){

      int nn = (iy + ny - 1) % ny;
      int isn = ix + nx* (nn + izt*ny);
      int isg = isn + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];

      u_0 = u_dn[IDX2_G_R(0,0,isg)];
      u_1 = u_dn[IDX2_G_I(0,0,isg)];
      u_2 = u_dn[IDX2_G_R(0,1,isg)];
      u_3 = u_dn[IDX2_G_I(0,1,isg)];
      u_4 = u_dn[IDX2_G_R(0,2,isg)];
      u_5 = u_dn[IDX2_G_I(0,2,isg)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //    ic = 0;
      real_t bc2 = ( iy == 0 ) ? bc_y : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 +=  bc2 * wt2r;
      v2_13 +=  bc2 * wt2i;
      v2_04 += -bc2 * wt1r;
      v2_14 += -bc2 * wt1i;

      u_6 = u_dn[IDX2_G_R(1,0,isg)];
      u_7 = u_dn[IDX2_G_I(1,0,isg)];
      u_8 = u_dn[IDX2_G_R(1,1,isg)];
      u_9 = u_dn[IDX2_G_I(1,1,isg)];
      u10 = u_dn[IDX2_G_R(1,2,isg)];
      u11 = u_dn[IDX2_G_I(1,2,isg)];

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

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_dn[IDX2_G_R(2,0,isg)];
      u13 = u_dn[IDX2_G_I(2,0,isg)];
      u14 = u_dn[IDX2_G_R(2,1,isg)];
      u15 = u_dn[IDX2_G_I(2,1,isg)];
      u16 = u_dn[IDX2_G_R(2,2,isg)];
      u17 = u_dn[IDX2_G_I(2,2,isg)];
#endif

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
      }

      //mult_zp
      idir = 2;

      if(( iz < nz -1 ) || (do_comm_z == 0)){

      int nn = (iz + 1) % nz;
      int isn = ixy  + nn*nx*ny + it*nx*ny*nz;
      int isg = site + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];

      u_0 = u_up[IDX2_G_R(0,0,isg)];
      u_1 = u_up[IDX2_G_I(0,0,isg)];
      u_2 = u_up[IDX2_G_R(1,0,isg)];
      u_3 = u_up[IDX2_G_I(1,0,isg)];
      u_4 = u_up[IDX2_G_R(2,0,isg)];
      u_5 = u_up[IDX2_G_I(2,0,isg)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //   ic = 0;
      real_t bc2 = ( iz == nz - 1 ) ? bc_z : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 +=  bc2 * wt1i;
      v2_13 += -bc2 * wt1r;
      v2_04 += -bc2 * wt2i;
      v2_14 +=  bc2 * wt2r;

      u_6 = u_up[IDX2_G_R(0,1,isg)];
      u_7 = u_up[IDX2_G_I(0,1,isg)];
      u_8 = u_up[IDX2_G_R(1,1,isg)];
      u_9 = u_up[IDX2_G_I(1,1,isg)];
      u10 = u_up[IDX2_G_R(2,1,isg)];
      u11 = u_up[IDX2_G_I(2,1,isg)];

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

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_up[IDX2_G_R(0,2,isg)];
      u13 = u_up[IDX2_G_I(0,2,isg)];
      u14 = u_up[IDX2_G_R(1,2,isg)];
      u15 = u_up[IDX2_G_I(1,2,isg)];
      u16 = u_up[IDX2_G_R(2,2,isg)];
      u17 = u_up[IDX2_G_I(2,2,isg)];
#endif

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
      }


      //mult_zm
      if(( iz > 0 ) || (do_comm_z == 0)){

      int nn = (iz + nz - 1) % nz;
      int isn = ixy + nn*nx*ny + it*nx*ny*nz;
      int isg = isn + nst * idir;

      vt1_0 = v1[IDX2_SP_5D_R(0,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt1_1 = v1[IDX2_SP_5D_I(0,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt1_2 = v1[IDX2_SP_5D_R(1,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt1_3 = v1[IDX2_SP_5D_I(1,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt1_4 = v1[IDX2_SP_5D_R(2,0,is,Ns,isn)] + v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];
      vt1_5 = v1[IDX2_SP_5D_I(2,0,is,Ns,isn)] - v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];

      vt2_0 = v1[IDX2_SP_5D_R(0,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt2_1 = v1[IDX2_SP_5D_I(0,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt2_2 = v1[IDX2_SP_5D_R(1,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt2_3 = v1[IDX2_SP_5D_I(1,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt2_4 = v1[IDX2_SP_5D_R(2,1,is,Ns,isn)] - v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];
      vt2_5 = v1[IDX2_SP_5D_I(2,1,is,Ns,isn)] + v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];

      u_0 = u_dn[IDX2_G_R(0,0,isg)];
      u_1 = u_dn[IDX2_G_I(0,0,isg)];
      u_2 = u_dn[IDX2_G_R(0,1,isg)];
      u_3 = u_dn[IDX2_G_I(0,1,isg)];
      u_4 = u_dn[IDX2_G_R(0,2,isg)];
      u_5 = u_dn[IDX2_G_I(0,2,isg)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      //    ic = 0;
      real_t bc2 = ( iz == 0 ) ? bc_z : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;
      v2_03 += -bc2 * wt1i;
      v2_13 +=  bc2 * wt1r;
      v2_04 +=  bc2 * wt2i;
      v2_14 += -bc2 * wt2r;

      u_6 = u_dn[IDX2_G_R(1,0,isg)];
      u_7 = u_dn[IDX2_G_I(1,0,isg)];
      u_8 = u_dn[IDX2_G_R(1,1,isg)];
      u_9 = u_dn[IDX2_G_I(1,1,isg)];
      u10 = u_dn[IDX2_G_R(1,2,isg)];
      u11 = u_dn[IDX2_G_I(1,2,isg)];

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

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_dn[IDX2_G_R(2,0,isg)];
      u13 = u_dn[IDX2_G_I(2,0,isg)];
      u14 = u_dn[IDX2_G_R(2,1,isg)];
      u15 = u_dn[IDX2_G_I(2,1,isg)];
      u16 = u_dn[IDX2_G_R(2,2,isg)];
      u17 = u_dn[IDX2_G_I(2,2,isg)];
#endif

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
      }

      idir = 3;

      if(( it < nt -1 ) || (do_comm_t == 0)){

      int nn = (it + 1) % Nt;
      int ixyz = site % (Nx*Ny*Nz);
      int isn  = ixyz + nn*Nx*Ny*Nz;
      int isg  = site + Nst*idir;

      vt1_0 = 2.0 * v1[IDX2_SP_5D_R(0,2,is,Ns,isn)];
      vt1_1 = 2.0 * v1[IDX2_SP_5D_I(0,2,is,Ns,isn)];
      vt1_2 = 2.0 * v1[IDX2_SP_5D_R(1,2,is,Ns,isn)];
      vt1_3 = 2.0 * v1[IDX2_SP_5D_I(1,2,is,Ns,isn)];
      vt1_4 = 2.0 * v1[IDX2_SP_5D_R(2,2,is,Ns,isn)];
      vt1_5 = 2.0 * v1[IDX2_SP_5D_I(2,2,is,Ns,isn)];

      vt2_0 = 2.0 * v1[IDX2_SP_5D_R(0,3,is,Ns,isn)];
      vt2_1 = 2.0 * v1[IDX2_SP_5D_I(0,3,is,Ns,isn)];
      vt2_2 = 2.0 * v1[IDX2_SP_5D_R(1,3,is,Ns,isn)];
      vt2_3 = 2.0 * v1[IDX2_SP_5D_I(1,3,is,Ns,isn)];
      vt2_4 = 2.0 * v1[IDX2_SP_5D_R(2,3,is,Ns,isn)];
      vt2_5 = 2.0 * v1[IDX2_SP_5D_I(2,3,is,Ns,isn)];

      u_0 = u_up[IDX2_G_R(0,0,isg)];
      u_1 = u_up[IDX2_G_I(0,0,isg)];
      u_2 = u_up[IDX2_G_R(1,0,isg)];
      u_3 = u_up[IDX2_G_I(1,0,isg)];
      u_4 = u_up[IDX2_G_R(2,0,isg)];
      u_5 = u_up[IDX2_G_I(2,0,isg)];

      wt1r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_0, u_1, u_2, u_3, u_4, u_5,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      real_t bc2 = ( it == nt -1 ) ? bc_t : 1.0;

      v2_03 +=  bc2 * wt1r;
      v2_13 +=  bc2 * wt1i;
      v2_04 +=  bc2 * wt2r;
      v2_14 +=  bc2 * wt2i;

      u_6 = u_up[IDX2_G_R(0,1,isg)];
      u_7 = u_up[IDX2_G_I(0,1,isg)];
      u_8 = u_up[IDX2_G_R(1,1,isg)];
      u_9 = u_up[IDX2_G_I(1,1,isg)];
      u10 = u_up[IDX2_G_R(2,1,isg)];
      u11 = u_up[IDX2_G_I(2,1,isg)];

      wt1r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u_6, u_7, u_8, u_9, u10, u11,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      v2_23 +=  bc2 * wt1r;
      v2_33 +=  bc2 * wt1i;
      v2_24 +=  bc2 * wt2r;
      v2_34 +=  bc2 * wt2i;

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_up[IDX2_G_R(0,2,isg)];
      u13 = u_up[IDX2_G_I(0,2,isg)];
      u14 = u_up[IDX2_G_R(1,2,isg)];
      u15 = u_up[IDX2_G_I(1,2,isg)];
      u16 = u_up[IDX2_G_R(2,2,isg)];
      u17 = u_up[IDX2_G_I(2,2,isg)];
#endif

      wt1r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GXr(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GXi(u12, u13, u14, u15, u16, u17,
                      vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      v2_43 +=  bc2 * wt1r;
      v2_53 +=  bc2 * wt1i;
      v2_44 +=  bc2 * wt2r;
      v2_54 +=  bc2 * wt2i;
      }

      //mult_tm
      if(( it > 0 ) || (do_comm_t == 0)){

      int nn = (it + nt - 1) % nt;
      int isn = ixyz + nn*nx*ny*nz;
      int isg = ixyz + nn*nx*ny*nz + nst*idir;

      vt1_0 = 2.0 * v1[IDX2_SP_5D_R(0,0,is,Ns,isn)];
      vt1_1 = 2.0 * v1[IDX2_SP_5D_I(0,0,is,Ns,isn)];
      vt1_2 = 2.0 * v1[IDX2_SP_5D_R(1,0,is,Ns,isn)];
      vt1_3 = 2.0 * v1[IDX2_SP_5D_I(1,0,is,Ns,isn)];
      vt1_4 = 2.0 * v1[IDX2_SP_5D_R(2,0,is,Ns,isn)];
      vt1_5 = 2.0 * v1[IDX2_SP_5D_I(2,0,is,Ns,isn)];

      vt2_0 = 2.0 * v1[IDX2_SP_5D_R(0,1,is,Ns,isn)];
      vt2_1 = 2.0 * v1[IDX2_SP_5D_I(0,1,is,Ns,isn)];
      vt2_2 = 2.0 * v1[IDX2_SP_5D_R(1,1,is,Ns,isn)];
      vt2_3 = 2.0 * v1[IDX2_SP_5D_I(1,1,is,Ns,isn)];
      vt2_4 = 2.0 * v1[IDX2_SP_5D_R(2,1,is,Ns,isn)];
      vt2_5 = 2.0 * v1[IDX2_SP_5D_I(2,1,is,Ns,isn)];

      u_0 = u_dn[IDX2_G_R(0,0,isg)];
      u_1 = u_dn[IDX2_G_I(0,0,isg)];
      u_2 = u_dn[IDX2_G_R(0,1,isg)];
      u_3 = u_dn[IDX2_G_I(0,1,isg)];
      u_4 = u_dn[IDX2_G_R(0,2,isg)];
      u_5 = u_dn[IDX2_G_I(0,2,isg)];

      wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      real_t bc2 = ( it == 0 ) ? bc_t : 1.0;

      v2_01 +=  bc2 * wt1r;
      v2_11 +=  bc2 * wt1i;
      v2_02 +=  bc2 * wt2r;
      v2_12 +=  bc2 * wt2i;

      u_6 = u_dn[IDX2_G_R(1,0,isg)];
      u_7 = u_dn[IDX2_G_I(1,0,isg)];
      u_8 = u_dn[IDX2_G_R(1,1,isg)];
      u_9 = u_dn[IDX2_G_I(1,1,isg)];
      u10 = u_dn[IDX2_G_R(1,2,isg)];
      u11 = u_dn[IDX2_G_I(1,2,isg)];

      wt1r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u_6, u_7, u_8, u_9, u10, u11,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      v2_21 +=  bc2 * wt1r;
      v2_31 +=  bc2 * wt1i;
      v2_22 +=  bc2 * wt2r;
      v2_32 +=  bc2 * wt2i;

#ifdef SU3_3RD_ROW_RECONST
      u12 = EXT_IMG_R(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u13 = EXT_IMG_I(u_2, u_3, u_4, u_5, u_8, u_9, u10, u11);
      u14 = EXT_IMG_R(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u15 = EXT_IMG_I(u_4, u_5, u_0, u_1, u10, u11, u_6, u_7);
      u16 = EXT_IMG_R(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
      u17 = EXT_IMG_I(u_0, u_1, u_2, u_3, u_6, u_7, u_8, u_9);
#else
      u12 = u_dn[IDX2_G_R(2,0,isg)];
      u13 = u_dn[IDX2_G_I(2,0,isg)];
      u14 = u_dn[IDX2_G_R(2,1,isg)];
      u15 = u_dn[IDX2_G_I(2,1,isg)];
      u16 = u_dn[IDX2_G_R(2,2,isg)];
      u17 = u_dn[IDX2_G_I(2,2,isg)];
#endif

      wt1r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt1i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
      wt2r = MULT_GDXr(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
      wt2i = MULT_GDXi(u12, u13, u14, u15, u16, u17,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      v2_41 +=  bc2 * wt1r;
      v2_51 +=  bc2 * wt1i;
      v2_42 +=  bc2 * wt2r;
      v2_52 +=  bc2 * wt2i;
      }
}
//============================================================END=====