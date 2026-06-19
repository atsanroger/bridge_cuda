/*!
      @file    mult_Domainwall_cuda_hop1-dirac-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-05-08 13:51:53 #$
      @version $LastChangedRevision: 2606 $
*/

int nx = Nx;
int ny = Ny;
int nz = Nz;
int nt = Nt;
int nst = nx * ny * nz * nt;

idir = 0;
if( do_comm_x > 0){

	int ix2  = 0;
	int isn  = ix2 + nx * iyzt;
	real_t bc2 = bc_x;

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

	v2_01 += bc2 * vt1_0;
	v2_11 += bc2 * vt1_1;
	v2_02 += bc2 * vt1_2;
	v2_12 += bc2 * vt1_3;
	v2_03 += bc2 * vt1_4;
	v2_13 += bc2 * vt1_5;

	v2_21 += bc2 * vt2_0;
	v2_31 += bc2 * vt2_1;
	v2_22 += bc2 * vt2_2;
	v2_32 += bc2 * vt2_3;
	v2_23 += bc2 * vt2_4;
	v2_33 += bc2 * vt2_5;

	buf_xp[IDX2_SP_5D_R(0,0,is,Ns,isn)] = v2_01;
	buf_xp[IDX2_SP_5D_I(0,0,is,Ns,isn)] = v2_11;
	buf_xp[IDX2_SP_5D_R(1,0,is,Ns,isn)] = v2_21;
	buf_xp[IDX2_SP_5D_I(1,0,is,Ns,isn)] = v2_31;
	buf_xp[IDX2_SP_5D_R(2,0,is,Ns,isn)] = v2_41;
	buf_xp[IDX2_SP_5D_I(2,0,is,Ns,isn)] = v2_51;

	buf_xp[IDX2_SP_5D_R(0,1,is,Ns,isn)] = v2_02;
	buf_xp[IDX2_SP_5D_I(0,1,is,Ns,isn)] = v2_12;
	buf_xp[IDX2_SP_5D_R(1,1,is,Ns,isn)] = v2_22;
	buf_xp[IDX2_SP_5D_I(1,1,is,Ns,isn)] = v2_32;
	buf_xp[IDX2_SP_5D_R(2,1,is,Ns,isn)] = v2_42;
	buf_xp[IDX2_SP_5D_I(2,1,is,Ns,isn)] = v2_52;

    // minus direction of x
    int ix2 = nx  - 1 ;
    int isn = ix2 + nx  * iyzt; //nei
    int isg = isn + nst * idir;

//  Load_u
    u_0 = u_dn[IDX2_G_R(0,0,isg)];
    u_1 = u_dn[IDX2_G_I(0,0,isg)];
    u_2 = u_dn[IDX2_G_R(0,1,isg)];
    u_3 = u_dn[IDX2_G_I(0,1,isg)];
    u_4 = u_dn[IDX2_G_R(0,2,isg)];
    u_5 = u_dn[IDX2_G_I(0,2,isg)];
    u_6 = u_dn[IDX2_G_R(1,0,isg)];
    u_7 = u_dn[IDX2_G_I(1,0,isg)];
    u_8 = u_dn[IDX2_G_R(1,1,isg)];
    u_9 = u_dn[IDX2_G_I(1,1,isg)];
    u10 = u_dn[IDX2_G_R(1,2,isg)];
    u11 = u_dn[IDX2_G_I(1,2,isg)];
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




    wt1r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt1i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);
    wt2r = MULT_GDXr(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);
    wt2i = MULT_GDXi(u_0, u_1, u_2, u_3, u_4, u_5,
                       vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

    // ic = 0;
    real_t bc2 = bc_x;

    v2_01 +=  bc2 * wt1r;
    v2_11 +=  bc2 * wt1i;
    v2_02 +=  bc2 * wt2r;
    v2_12 +=  bc2 * wt2i;
    v2_03 += -bc2 * wt2i;
    v2_13 += +bc2 * wt2r;
    v2_04 += -bc2 * wt1i;
    v2_14 += +bc2 * wt1r;


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

    buf_xm[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
    buf_xm[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
    buf_xm[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
    buf_xm[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
    buf_xm[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
    buf_xm[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

    buf_xm[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
    buf_xm[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
    buf_xm[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
    buf_xm[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
    buf_xm[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
    buf_xm[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

    buf_xm[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
    buf_xm[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
    buf_xm[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
    buf_xm[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
    buf_xm[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
    buf_xm[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

    buf_xm[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
    buf_xm[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
    buf_xm[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
    buf_xm[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
    buf_xm[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
    buf_xm[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;


}


dir = 1;
if( do_comm_y > 0){

//printf()

}

dir = 2;
if( do_comm_z > 0){


}


dir = 3;
if( do_comm_t == 0){

}
