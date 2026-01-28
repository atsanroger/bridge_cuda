/*!
      @file    mult_Staggered_openacc2R_uvdn-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

      vt0 = v1[ IDX2_1SP_R(0, nei) ];
      vt1 = v1[ IDX2_1SP_I(0, nei) ];
      vt2 = v1[ IDX2_1SP_R(1, nei) ];
      vt3 = v1[ IDX2_1SP_I(1, nei) ];
      vt4 = v1[ IDX2_1SP_R(2, nei) ];
      vt5 = v1[ IDX2_1SP_I(2, nei) ];

      // ic = 0
      ut0 =  u_dn[IDX2_G_R(0, 0, neig)];
      ut1 = -u_dn[IDX2_G_I(0, 0, neig)];
      ut2 =  u_dn[IDX2_G_R(0, 1, neig)];
      ut3 = -u_dn[IDX2_G_I(0, 1, neig)];
      ut4 =  u_dn[IDX2_G_R(0, 2, neig)];
      ut5 = -u_dn[IDX2_G_I(0, 2, neig)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt0 += -bc2 * wtr;
      xt1 += -bc2 * wti;

      // ic = 1
      ut0 =  u_dn[IDX2_G_R(1, 0, neig)];
      ut1 = -u_dn[IDX2_G_I(1, 0, neig)];
      ut2 =  u_dn[IDX2_G_R(1, 1, neig)];
      ut3 = -u_dn[IDX2_G_I(1, 1, neig)];
      ut4 =  u_dn[IDX2_G_R(1, 2, neig)];
      ut5 = -u_dn[IDX2_G_I(1, 2, neig)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt2 += -bc2 * wtr;
      xt3 += -bc2 * wti;

      // ic = 2
      ut0 =  u_dn[IDX2_G_R(2, 0, neig)];
      ut1 = -u_dn[IDX2_G_I(2, 0, neig)];
      ut2 =  u_dn[IDX2_G_R(2, 1, neig)];
      ut3 = -u_dn[IDX2_G_I(2, 1, neig)];
      ut4 =  u_dn[IDX2_G_R(2, 2, neig)];
      ut5 = -u_dn[IDX2_G_I(2, 2, neig)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt4 += -bc2 * wtr;
      xt5 += -bc2 * wti;

//============================================================END=====
