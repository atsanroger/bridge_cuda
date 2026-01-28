/*!
      @file    mult_Clover_openacc2_csw_dirac-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

      v2_01 = -kappa * v2_01;
      v2_11 = -kappa * v2_11;
      v2_21 = -kappa * v2_21;
      v2_31 = -kappa * v2_31;
      v2_41 = -kappa * v2_41;
      v2_51 = -kappa * v2_51;

      v2_02 = -kappa * v2_02;
      v2_12 = -kappa * v2_12;
      v2_22 = -kappa * v2_22;
      v2_32 = -kappa * v2_32;
      v2_42 = -kappa * v2_42;
      v2_52 = -kappa * v2_52;

      v2_03 = -kappa * v2_03;
      v2_13 = -kappa * v2_13;
      v2_23 = -kappa * v2_23;
      v2_33 = -kappa * v2_33;
      v2_43 = -kappa * v2_43;
      v2_53 = -kappa * v2_53;

      v2_04 = -kappa * v2_04;
      v2_14 = -kappa * v2_14;
      v2_24 = -kappa * v2_24;
      v2_34 = -kappa * v2_34;
      v2_44 = -kappa * v2_44;
      v2_54 = -kappa * v2_54;

      for(int id = 0; id < ND2; ++id){
        int id2 = (id + ND2) % ND;

        vt1_0 = v1[IDX2_SP_R(0, id, site)];
        vt1_1 = v1[IDX2_SP_I(0, id, site)];
        vt1_2 = v1[IDX2_SP_R(1, id, site)];
        vt1_3 = v1[IDX2_SP_I(1, id, site)];
        vt1_4 = v1[IDX2_SP_R(2, id, site)];
        vt1_5 = v1[IDX2_SP_I(2, id, site)];

        vt2_0 = v1[IDX2_SP_R(0, id2, site)];
        vt2_1 = v1[IDX2_SP_I(0, id2, site)];
        vt2_2 = v1[IDX2_SP_R(1, id2, site)];
        vt2_3 = v1[IDX2_SP_I(1, id2, site)];
        vt2_4 = v1[IDX2_SP_R(2, id2, site)];
        vt2_5 = v1[IDX2_SP_I(2, id2, site)];

        int icst1 = site + Nst * (id  + ND * 0);
        int icst2 = site + Nst * (id2 + ND * 0);
        int icst3 = site + Nst * (id  + ND * 1);
        int icst4 = site + Nst * (id2 + ND * 1);

        int ic2 = 0;
        u_0 = ct[IDX2_G_R(0, ic2, icst1)];
        u_1 = ct[IDX2_G_I(0, ic2, icst1)];
        u_2 = ct[IDX2_G_R(1, ic2, icst1)];
        u_3 = ct[IDX2_G_I(1, ic2, icst1)];
        u_4 = ct[IDX2_G_R(2, ic2, icst1)];
        u_5 = ct[IDX2_G_I(2, ic2, icst1)];

        v2_01 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_11 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_03 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_13 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst2)];
        u_7 = ct[IDX2_G_I(0, ic2, icst2)];
        u_8 = ct[IDX2_G_R(1, ic2, icst2)];
        u_9 = ct[IDX2_G_I(1, ic2, icst2)];
        u10 = ct[IDX2_G_R(2, ic2, icst2)];
        u11 = ct[IDX2_G_I(2, ic2, icst2)];

        v2_03 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_13 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_01 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_11 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_0 = ct[IDX2_G_R(0, ic2, icst3)];
        u_1 = ct[IDX2_G_I(0, ic2, icst3)];
        u_2 = ct[IDX2_G_R(1, ic2, icst3)];
        u_3 = ct[IDX2_G_I(1, ic2, icst3)];
        u_4 = ct[IDX2_G_R(2, ic2, icst3)];
        u_5 = ct[IDX2_G_I(2, ic2, icst3)];

        v2_02 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_12 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_04 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_14 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst4)];
        u_7 = ct[IDX2_G_I(0, ic2, icst4)];
        u_8 = ct[IDX2_G_R(1, ic2, icst4)];
        u_9 = ct[IDX2_G_I(1, ic2, icst4)];
        u10 = ct[IDX2_G_R(2, ic2, icst4)];
        u11 = ct[IDX2_G_I(2, ic2, icst4)];

        v2_04 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_14 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_02 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_12 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);


        ic2 = 1;

        u_0 = ct[IDX2_G_R(0, ic2, icst1)];
        u_1 = ct[IDX2_G_I(0, ic2, icst1)];
        u_2 = ct[IDX2_G_R(1, ic2, icst1)];
        u_3 = ct[IDX2_G_I(1, ic2, icst1)];
        u_4 = ct[IDX2_G_R(2, ic2, icst1)];
        u_5 = ct[IDX2_G_I(2, ic2, icst1)];

        v2_21 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_31 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_23 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_33 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst2)];
        u_7 = ct[IDX2_G_I(0, ic2, icst2)];
        u_8 = ct[IDX2_G_R(1, ic2, icst2)];
        u_9 = ct[IDX2_G_I(1, ic2, icst2)];
        u10 = ct[IDX2_G_R(2, ic2, icst2)];
        u11 = ct[IDX2_G_I(2, ic2, icst2)];

        v2_23 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_33 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_21 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_31 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_0 = ct[IDX2_G_R(0, ic2, icst3)];
        u_1 = ct[IDX2_G_I(0, ic2, icst3)];
        u_2 = ct[IDX2_G_R(1, ic2, icst3)];
        u_3 = ct[IDX2_G_I(1, ic2, icst3)];
        u_4 = ct[IDX2_G_R(2, ic2, icst3)];
        u_5 = ct[IDX2_G_I(2, ic2, icst3)];

        v2_22 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_32 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_24 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_34 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst4)];
        u_7 = ct[IDX2_G_I(0, ic2, icst4)];
        u_8 = ct[IDX2_G_R(1, ic2, icst4)];
        u_9 = ct[IDX2_G_I(1, ic2, icst4)];
        u10 = ct[IDX2_G_R(2, ic2, icst4)];
        u11 = ct[IDX2_G_I(2, ic2, icst4)];

        v2_24 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_34 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_22 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_32 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);


        ic2 = 2;

        u_0 = ct[IDX2_G_R(0, ic2, icst1)];
        u_1 = ct[IDX2_G_I(0, ic2, icst1)];
        u_2 = ct[IDX2_G_R(1, ic2, icst1)];
        u_3 = ct[IDX2_G_I(1, ic2, icst1)];
        u_4 = ct[IDX2_G_R(2, ic2, icst1)];
        u_5 = ct[IDX2_G_I(2, ic2, icst1)];

        v2_41 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_51 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_43 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_53 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst2)];
        u_7 = ct[IDX2_G_I(0, ic2, icst2)];
        u_8 = ct[IDX2_G_R(1, ic2, icst2)];
        u_9 = ct[IDX2_G_I(1, ic2, icst2)];
        u10 = ct[IDX2_G_R(2, ic2, icst2)];
        u11 = ct[IDX2_G_I(2, ic2, icst2)];

        v2_43 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_53 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_41 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_51 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_0 = ct[IDX2_G_R(0, ic2, icst3)];
        u_1 = ct[IDX2_G_I(0, ic2, icst3)];
        u_2 = ct[IDX2_G_R(1, ic2, icst3)];
        u_3 = ct[IDX2_G_I(1, ic2, icst3)];
        u_4 = ct[IDX2_G_R(2, ic2, icst3)];
        u_5 = ct[IDX2_G_I(2, ic2, icst3)];

        v2_42 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_52 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_44 += MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_54 += MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        u_6 = ct[IDX2_G_R(0, ic2, icst4)];
        u_7 = ct[IDX2_G_I(0, ic2, icst4)];
        u_8 = ct[IDX2_G_R(1, ic2, icst4)];
        u_9 = ct[IDX2_G_I(1, ic2, icst4)];
        u10 = ct[IDX2_G_R(2, ic2, icst4)];
        u11 = ct[IDX2_G_I(2, ic2, icst4)];

        v2_44 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_54 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

        v2_42 += MULT_UV_R(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

        v2_52 += MULT_UV_I(u_6, u_7, u_8, u_9, u10, u11,
                           vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

      }

//============================================================END=====
