/*!
      @file    mult_Wilson_t_chiral_cuda_qdw-inc.h
      @brief   QDW version of mult_Wilson_t_chiral_cuda-inc.h
*/

      // mult_tp
      idir = 3;
      nn = (it + 1) % Nt;

      int ixyz = site % (Nx * Ny * Nz);

      isn = ixyz + (Nx * Ny * Nz) * nn;
      isg = ixyz + (Nx * Ny * Nz) * it + Nst_pad * idir;

      // Variables already declared in xyz include.

      // vt1 = s0 + s2
      // vt2 = s1 + s3
      
      // Color 0
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)];
      QDW_ADD(vt1_c0, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)];
      QDW_ADD(vt2_c0, vs1, vs3);
      
      // Color 1
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)];
      QDW_ADD(vt1_c1, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)];
      QDW_ADD(vt2_c1, vs1, vs3);
      
      // Color 2
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)];
      QDW_ADD(vt1_c2, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)];
      QDW_ADD(vt2_c2, vs1, vs3);
      
      // Gauge Mul (Row 0)
      u_val.x = u_up[IDX2_G_R(0,0,isg)]; u_val.y = u_up[IDX2_G_I(0,0,isg)];
      wt1_c0 = qdw_mult_uc(u_val, vt1_c0); wt2_c0 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(0,1,isg)]; u_val.y = u_up[IDX2_G_I(0,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      u_val.x = u_up[IDX2_G_R(0,2,isg)]; u_val.y = u_up[IDX2_G_I(0,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      
      // Gauge Mul Row 1
      u_val.x = u_up[IDX2_G_R(1,0,isg)]; u_val.y = u_up[IDX2_G_I(1,0,isg)];
      wt1_c1 = qdw_mult_uc(u_val, vt1_c0); wt2_c1 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(1,1,isg)]; u_val.y = u_up[IDX2_G_I(1,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(1,2,isg)]; u_val.y = u_up[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);

      // Gauge Mul Row 2
      u_val.x = u_up[IDX2_G_R(2,0,isg)]; u_val.y = u_up[IDX2_G_I(2,0,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0); wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(2,1,isg)]; u_val.y = u_up[IDX2_G_I(2,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      u_val.x = u_up[IDX2_G_R(2,2,isg)]; u_val.y = u_up[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);

      // Accumulate
      bc2 = 1.0;
      if(it == Nt-1) bc2 = bc_t;
      
      // Chiral logic:
      // s0 += wt1
      // s2 += wt1
      // s1 += wt2
      // s3 += wt2
      
      // Color 0
      QDW_SCAL(tmp_scaled, bc2, wt1_c0);
      QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c0);
      QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled);
      
      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1);
      QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c1);
      QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);
      
      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2);
      QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c2);
      QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);
      
      
      // mult_tm
      nn = (it + Nt - 1) % Nt;
      isn = ixyz + (Nx * Ny * Nz) * nn;
      isg = isn + Nst_pad * idir;
      
      // vt1 = s0 - s2
      // vt2 = s1 - s3
      
      // Color 0
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)];
      QDW_SUB(vt1_c0, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)];
      QDW_SUB(vt2_c0, vs1, vs3);
      
      // Color 1
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)];
      QDW_SUB(vt1_c1, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)];
      QDW_SUB(vt2_c1, vs1, vs3);
      
      // Color 2
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)];
      QDW_SUB(vt1_c2, vs0, vs2);
      
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)];
      QDW_SUB(vt2_c2, vs1, vs3);
      
      // Gauge Mul (u_dn Dagger)
      // Same logic as before
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
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      // Row 2
      u_val.x = u_dn[IDX2_G_R(0,2,isg)]; u_val.y = -u_dn[IDX2_G_I(0,2,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0); wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,2,isg)]; u_val.y = -u_dn[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,2,isg)]; u_val.y = -u_dn[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, tmp_prod);
      
      // Accumulate
      bc2 = 1.0;
      if(it == 0) bc2 = bc_t;
      
      // Chiral logic:
      // s0 += wt1
      // s2 -= wt1
      // s1 += wt2
      // s3 -= wt2
      
      // Color 0
      QDW_SCAL(tmp_scaled, bc2, wt1_c0);
      QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SUB(v2_c0_s2, v2_c0_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c0);
      QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      QDW_SUB(v2_c0_s3, v2_c0_s3, tmp_scaled);
      
      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1);
      QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SUB(v2_c1_s2, v2_c1_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c1);
      QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      QDW_SUB(v2_c1_s3, v2_c1_s3, tmp_scaled);
      
      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2);
      QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SUB(v2_c2_s2, v2_c2_s2, tmp_scaled);
      
      QDW_SCAL(tmp_scaled, bc2, wt2_c2);
      QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      QDW_SUB(v2_c2_s3, v2_c2_s3, tmp_scaled);
