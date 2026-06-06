/*!
      @file    mult_Wilson_t_dirac_cuda_qdw-inc.h
      @brief   QDW version of mult_Wilson_t_dirac_cuda-inc.h
*/

      // mult_tp (T plus)
      idir = 3;
      nn = (it + 1) % Nt;

      int ixyz = site % (Nx * Ny * Nz);

      isn = ixyz + (Nx * Ny * Nz) * nn;
      isg = ixyz + (Nx * Ny * Nz) * it + Nst_pad * idir;

      // vt1_0 = 2.0 * v1[IDX2_SP_R(0,2,isn)];
      // vt1_1 = 2.0 * v1[IDX2_SP_I(0,2,isn)];
      // Original: 2 * spinor projection P_t+
      // (1 - gamma_4) psi.
      // gamma_4 = diag(1, 1, -1, -1).
      // 1 - g4 = diag(0, 0, 2, 2).
      // s0 -> 0, s1 -> 0, s2 -> 2*s2, s3 -> 2*s3.
      
      // Variables vs0..vs3, vt1/2_c0..c2, u_val, tmp_prod, wt1/2_c0..c2,
      // bc2, tmp_scaled, v2_c*_s* are already declared in xyz include.

      // Color 0
      vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)];
      QDW_SCAL(vt1_c0, 2.0, vs2); // s2 -> 2*s2
      
      vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)];
      QDW_SCAL(vt2_c0, 2.0, vs3); // s3 -> 2*s3
      
      // Color 1
      vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)];
      QDW_SCAL(vt1_c1, 2.0, vs2); 
      vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)];
      QDW_SCAL(vt2_c1, 2.0, vs3); 
      
      // Color 2
      vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)];
      QDW_SCAL(vt1_c2, 2.0, vs2);
      vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)];
      QDW_SCAL(vt2_c2, 2.0, vs3);

      // Color 0 Result
      u_val.x = u_up[IDX2_G_R(0,0,isg)]; u_val.y = u_up[IDX2_G_I(0,0,isg)];
      wt1_c0 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c0 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(0,1,isg)]; u_val.y = u_up[IDX2_G_I(0,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      u_val.x = u_up[IDX2_G_R(0,2,isg)]; u_val.y = u_up[IDX2_G_I(0,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      
      // Color 1 Result
      u_val.x = u_up[IDX2_G_R(1,0,isg)]; u_val.y = u_up[IDX2_G_I(1,0,isg)];
      wt1_c1 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c1 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(1,1,isg)]; u_val.y = u_up[IDX2_G_I(1,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(1,2,isg)]; u_val.y = u_up[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      
      // Color 2 Result
      u_val.x = u_up[IDX2_G_R(2,0,isg)]; u_val.y = u_up[IDX2_G_I(2,0,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(2,1,isg)]; u_val.y = u_up[IDX2_G_I(2,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      u_val.x = u_up[IDX2_G_R(2,2,isg)]; u_val.y = u_up[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      
      // Accumulate T_P
      bc2 = 1.0;
      if(it == Nt-1) bc2 = bc_t;
      
      // Color 0
      // Proj TP maps to s2, s3.
      // wt1 -> s2. wt2 -> s3.
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled);
      
      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);
      
      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);
      
      
      // mult_tm (T minus)
      nn = (it + Nt - 1) % Nt;
      isn = ixyz + (Nx * Ny * Nz) * nn;
      isg = isn + Nst_pad * idir;
      
      // (1 + gamma_4) = diag(2, 2, 0, 0).
      // s0 -> 2*s0, s1 -> 2*s1.
      
      // Color 0
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)];
      PROJ_T_P(vt1_c0, vs0); // Use same macro (2*a)
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)];
      PROJ_T_P(vt2_c0, vs1);
      
      // Color 1
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)];
      PROJ_T_P(vt1_c1, vs0);
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)];
      PROJ_T_P(vt2_c1, vs1);
      
      // Color 2
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)];
      PROJ_T_P(vt1_c2, vs0);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)];
      PROJ_T_P(vt2_c2, vs1);
      
      // Gauge Mul (u_dn Dagger)
      // Row 0
      u_val.x = u_dn[IDX2_G_R(0,0,isg)]; u_val.y = -u_dn[IDX2_G_I(0,0,isg)];
      wt1_c0 = qdw_mult_uc(u_val, vt1_c0); wt2_c0 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,0,isg)]; u_val.y = -u_dn[IDX2_G_I(1,0,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,0,isg)]; u_val.y = -u_dn[IDX2_G_I(2,0,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      
      // Row 1
      u_val.x = u_dn[IDX2_G_R(0,1,isg)]; u_val.y = -u_dn[IDX2_G_I(0,1,isg)];
      wt1_c1 = qdw_mult_uc(u_val, vt1_c0); wt2_c1 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,1,isg)]; u_val.y = -u_dn[IDX2_G_I(1,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,1,isg)]; u_val.y = -u_dn[IDX2_G_I(2,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      
      // Row 2
      u_val.x = u_dn[IDX2_G_R(0,2,isg)]; u_val.y = -u_dn[IDX2_G_I(0,2,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0); wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,2,isg)]; u_val.y = -u_dn[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,2,isg)]; u_val.y = -u_dn[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      
      // Accumulate T_M
      bc2 = 1.0;
      if(it == 0) bc2 = bc_t;
      
      // Color 0 
      // wt1 -> s0. wt2 -> s1.
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      
      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      
      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
