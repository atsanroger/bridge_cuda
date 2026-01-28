/*!
      @file    mult_Domainwall_cuda_flag-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2024-05-08 13:51:53 #$
      @version $LastChangedRevision: 2606 $
*/

{
        if( flag == 0 ){

            v2_01 = 0.0;
            v2_11 = 0.0;
            v2_21 = 0.0;
            v2_31 = 0.0;
            v2_41 = 0.0;
            v2_51 = 0.0;
            v2_02 = 0.0;
            v2_12 = 0.0;
            v2_22 = 0.0;
            v2_32 = 0.0;
            v2_42 = 0.0;
            v2_52 = 0.0;
            v2_03 = 0.0;
            v2_13 = 0.0;
            v2_23 = 0.0;
            v2_33 = 0.0;
            v2_43 = 0.0;
            v2_53 = 0.0;
            v2_04 = 0.0;
            v2_14 = 0.0;
            v2_24 = 0.0;
            v2_34 = 0.0;
            v2_44 = 0.0;
            v2_54 = 0.0;

        }else if ( flag == 1 ){

            v2_01 = vp[IDX2_SP_5D_R(0,0,is,Ns,site)];
            v2_11 = vp[IDX2_SP_5D_I(0,0,is,Ns,site)];
            v2_21 = vp[IDX2_SP_5D_R(1,0,is,Ns,site)];
            v2_31 = vp[IDX2_SP_5D_I(1,0,is,Ns,site)];
            v2_41 = vp[IDX2_SP_5D_R(2,0,is,Ns,site)];
            v2_51 = vp[IDX2_SP_5D_I(2,0,is,Ns,site)];
            v2_02 = vp[IDX2_SP_5D_R(0,1,is,Ns,site)];
            v2_12 = vp[IDX2_SP_5D_I(0,1,is,Ns,site)];
            v2_22 = vp[IDX2_SP_5D_R(1,1,is,Ns,site)];
            v2_32 = vp[IDX2_SP_5D_I(1,1,is,Ns,site)];
            v2_42 = vp[IDX2_SP_5D_R(2,1,is,Ns,site)];
            v2_52 = vp[IDX2_SP_5D_I(2,1,is,Ns,site)];
            v2_03 = vp[IDX2_SP_5D_R(0,2,is,Ns,site)];
            v2_13 = vp[IDX2_SP_5D_I(0,2,is,Ns,site)];
            v2_23 = vp[IDX2_SP_5D_R(1,2,is,Ns,site)];
            v2_33 = vp[IDX2_SP_5D_I(1,2,is,Ns,site)];
            v2_43 = vp[IDX2_SP_5D_R(2,2,is,Ns,site)];
            v2_53 = vp[IDX2_SP_5D_I(2,2,is,Ns,site)];
            v2_04 = vp[IDX2_SP_5D_R(0,3,is,Ns,site)];
            v2_14 = vp[IDX2_SP_5D_I(0,3,is,Ns,site)];
            v2_24 = vp[IDX2_SP_5D_R(1,3,is,Ns,site)];
            v2_34 = vp[IDX2_SP_5D_I(1,3,is,Ns,site)];
            v2_44 = vp[IDX2_SP_5D_R(2,3,is,Ns,site)];
            v2_54 = vp[IDX2_SP_5D_I(2,3,is,Ns,site)];

        }else{
            printf("Flag is missing");
        }
        __syncthreads();
}