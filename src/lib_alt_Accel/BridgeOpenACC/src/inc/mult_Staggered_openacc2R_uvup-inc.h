/*!
      @file    mult_Staggered_openacc2R_uvup-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

      vt0 = bc2 * v1[ IDX2_1SP_R(0, nei) ];
      vt1 = bc2 * v1[ IDX2_1SP_I(0, nei) ];
      vt2 = bc2 * v1[ IDX2_1SP_R(1, nei) ];
      vt3 = bc2 * v1[ IDX2_1SP_I(1, nei) ];
      vt4 = bc2 * v1[ IDX2_1SP_R(2, nei) ];
      vt5 = bc2 * v1[ IDX2_1SP_I(2, nei) ];

      // ic = 0
      ut0 = u_up[IDX2_G_R(0, 0, istg)];
      ut1 = u_up[IDX2_G_I(0, 0, istg)];
      ut2 = u_up[IDX2_G_R(1, 0, istg)];
      ut3 = u_up[IDX2_G_I(1, 0, istg)];
      ut4 = u_up[IDX2_G_R(2, 0, istg)];
      ut5 = u_up[IDX2_G_I(2, 0, istg)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt0 += wtr;
      xt1 += wti;

      // ic = 1
      ut0 = u_up[IDX2_G_R(0, 1, istg)];
      ut1 = u_up[IDX2_G_I(0, 1, istg)];
      ut2 = u_up[IDX2_G_R(1, 1, istg)];
      ut3 = u_up[IDX2_G_I(1, 1, istg)];
      ut4 = u_up[IDX2_G_R(2, 1, istg)];
      ut5 = u_up[IDX2_G_I(2, 1, istg)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt2 += wtr;
      xt3 += wti;

      // ic = 2
      ut0 = u_up[IDX2_G_R(0, 2, istg)];
      ut1 = u_up[IDX2_G_I(0, 2, istg)];
      ut2 = u_up[IDX2_G_R(1, 2, istg)];
      ut3 = u_up[IDX2_G_I(1, 2, istg)];
      ut4 = u_up[IDX2_G_R(2, 2, istg)];
      ut5 = u_up[IDX2_G_I(2, 2, istg)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt4 += wtr;
      xt5 += wti;

//============================================================END=====
