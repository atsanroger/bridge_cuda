/*!
      @file    mult_Clover_openacc2_csw_chiral-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2574 $
*/

      real_t wt[NVC*ND];
      for(int icd = 0; icd < NVC*ND; ++icd){
        wt[icd] = 0.0;
      }

      for(int jd = 0; jd < ND2; ++jd){
        for(int id = 0; id < ND2; ++id){
          int id2 = id + ND2;
          int icst  = site + Nst * (id  + ND2*jd);

          vt1_0 = v1[IDX2_SP_R(0, id, site)];
          vt1_1 = v1[IDX2_SP_I(0, id, site)];
          vt1_2 = v1[IDX2_SP_R(1, id, site)];
          vt1_3 = v1[IDX2_SP_I(1, id, site)];
          vt1_4 = v1[IDX2_SP_R(2, id, site)];
          vt1_5 = v1[IDX2_SP_I(2, id, site)];

          for(int ic2 = 0; ic2 < NC; ++ic2){
            u_0 = ct[IDX2_G_R(0, ic2, icst)];
            u_1 = ct[IDX2_G_I(0, ic2, icst)];
            u_2 = ct[IDX2_G_R(1, ic2, icst)];
            u_3 = ct[IDX2_G_I(1, ic2, icst)];
            u_4 = ct[IDX2_G_R(2, ic2, icst)];
            u_5 = ct[IDX2_G_I(2, ic2, icst)];

            wt1r = MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                             vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

            wt1i = MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                             vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5);

            int jd2 = jd + ND2;
            wt[2*ic2   + NVC*jd]  += wt1r;
            wt[2*ic2+1 + NVC*jd]  += wt1i;
          }

        }
      }

      for(int jd = 0; jd < ND2; ++jd){
        for(int id = 0; id < ND2; ++id){
          int id2 = id + ND2;
          int icst  = site + Nst * (id  + ND2*jd + ND);

          vt2_0 = v1[IDX2_SP_R(0, id2, site)];
          vt2_1 = v1[IDX2_SP_I(0, id2, site)];
          vt2_2 = v1[IDX2_SP_R(1, id2, site)];
          vt2_3 = v1[IDX2_SP_I(1, id2, site)];
          vt2_4 = v1[IDX2_SP_R(2, id2, site)];
          vt2_5 = v1[IDX2_SP_I(2, id2, site)];

          for(int ic2 = 0; ic2 < NC; ++ic2){
            u_0 = ct[IDX2_G_R(0, ic2, icst)];
            u_1 = ct[IDX2_G_I(0, ic2, icst)];
            u_2 = ct[IDX2_G_R(1, ic2, icst)];
            u_3 = ct[IDX2_G_I(1, ic2, icst)];
            u_4 = ct[IDX2_G_R(2, ic2, icst)];
            u_5 = ct[IDX2_G_I(2, ic2, icst)];

            wt2r = MULT_UV_R(u_0, u_1, u_2, u_3, u_4, u_5,
                             vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

            wt2i = MULT_UV_I(u_0, u_1, u_2, u_3, u_4, u_5,
                             vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5);

            int jd2 = jd + ND2;
            wt[2*ic2   + NVC*jd2] += wt2r;
            wt[2*ic2+1 + NVC*jd2] += wt2i;
          }

        }
      }

      v2_01 = -kappa * v2_01 + wt[0 + NVC*0];
      v2_11 = -kappa * v2_11 + wt[1 + NVC*0];
      v2_21 = -kappa * v2_21 + wt[2 + NVC*0];
      v2_31 = -kappa * v2_31 + wt[3 + NVC*0];
      v2_41 = -kappa * v2_41 + wt[4 + NVC*0];
      v2_51 = -kappa * v2_51 + wt[5 + NVC*0];

      v2_02 = -kappa * v2_02 + wt[0 + NVC*1];
      v2_12 = -kappa * v2_12 + wt[1 + NVC*1];
      v2_22 = -kappa * v2_22 + wt[2 + NVC*1];
      v2_32 = -kappa * v2_32 + wt[3 + NVC*1];
      v2_42 = -kappa * v2_42 + wt[4 + NVC*1];
      v2_52 = -kappa * v2_52 + wt[5 + NVC*1];

      v2_03 = -kappa * v2_03 + wt[0 + NVC*2];
      v2_13 = -kappa * v2_13 + wt[1 + NVC*2];
      v2_23 = -kappa * v2_23 + wt[2 + NVC*2];
      v2_33 = -kappa * v2_33 + wt[3 + NVC*2];
      v2_43 = -kappa * v2_43 + wt[4 + NVC*2];
      v2_53 = -kappa * v2_53 + wt[5 + NVC*2];

      v2_04 = -kappa * v2_04 + wt[0 + NVC*3];
      v2_14 = -kappa * v2_14 + wt[1 + NVC*3];
      v2_24 = -kappa * v2_24 + wt[2 + NVC*3];
      v2_34 = -kappa * v2_34 + wt[3 + NVC*3];
      v2_44 = -kappa * v2_44 + wt[4 + NVC*3];
      v2_54 = -kappa * v2_54 + wt[5 + NVC*3];

//============================================================END=====
