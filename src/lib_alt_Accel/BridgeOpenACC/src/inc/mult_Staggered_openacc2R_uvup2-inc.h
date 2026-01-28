/*!
      @file    mult_Staggered_openacc2R_uvup2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

      real_t vt0, vt1, vt2, vt3, vt4, vt5;
      real_t ut0, ut1, ut2, ut3, ut4, ut5;
      real_t wtr, wti;

      vt0 = buf[IDX2_1SP_R(0, ibuf)];
      vt1 = buf[IDX2_1SP_I(0, ibuf)];
      vt2 = buf[IDX2_1SP_R(1, ibuf)];
      vt3 = buf[IDX2_1SP_I(1, ibuf)];
      vt4 = buf[IDX2_1SP_R(2, ibuf)];
      vt5 = buf[IDX2_1SP_I(2, ibuf)];

      // ic = 0
      ut0 = u_up[IDX2_G_R(0, 0, istu)];
      ut1 = u_up[IDX2_G_I(0, 0, istu)];
      ut2 = u_up[IDX2_G_R(1, 0, istu)];
      ut3 = u_up[IDX2_G_I(1, 0, istu)];
      ut4 = u_up[IDX2_G_R(2, 0, istu)];
      ut5 = u_up[IDX2_G_I(2, 0, istu)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt0 += wtr;
      xt1 += wti;

      // ic = 1
      ut0 = u_up[IDX2_G_R(0, 1, istu)];
      ut1 = u_up[IDX2_G_I(0, 1, istu)];
      ut2 = u_up[IDX2_G_R(1, 1, istu)];
      ut3 = u_up[IDX2_G_I(1, 1, istu)];
      ut4 = u_up[IDX2_G_R(2, 1, istu)];
      ut5 = u_up[IDX2_G_I(2, 1, istu)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt2 += wtr;
      xt3 += wti;

      // ic = 2
      ut0 = u_up[IDX2_G_R(0, 2, istu)];
      ut1 = u_up[IDX2_G_I(0, 2, istu)];
      ut2 = u_up[IDX2_G_R(1, 2, istu)];
      ut3 = u_up[IDX2_G_I(1, 2, istu)];
      ut4 = u_up[IDX2_G_R(2, 2, istu)];
      ut5 = u_up[IDX2_G_I(2, 2, istu)];

      wtr = MULT_UV_R(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);
      wti = MULT_UV_I(ut0, ut1, ut2, ut3, ut4, ut5,
                      vt0, vt1, vt2, vt3, vt4, vt5);

      xt4 += wtr;
      xt5 += wti;

//============================================================END=====
