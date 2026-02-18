/*!
      @file    mult_Wilson_xyz_cuda_qdw-inc.h
      @brief   QDW version of mult_Wilson_xyz_cuda-inc.h
      @author  Hideo Matsufuru (matufuru)
*/

      // =================================================================
      // X-Direction
      // =================================================================

      // mult_xp
      int idir = 0;

      int iyzt = site/Nx;
      int nn   = (ix+1) % Nx;

      int isn = nn + Nx * iyzt;
      int isg = ix + Nx * iyzt + Nst_pad * idir;

      double4 vs0, vs1, vs2, vs3;
      double4 vt1_c0, vt1_c1, vt1_c2;
      double4 vt2_c0, vt2_c1, vt2_c2;
      
      // Helper: res = a + i b
      #define PROJ_P(res, a, b) \
        res.x = a.x - b.y; \
        res.y = a.y + b.x; \
        res.z = a.z - b.w; \
        res.w = a.w + b.z;
        
      // Helper: res = a - i b
      #define PROJ_M(res, a, b) \
        res.x = a.x + b.y; \
        res.y = a.y - b.x; \
        res.z = a.z + b.w; \
        res.w = a.w - b.z;

      // Color 0
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)];
      PROJ_P(vt1_c0, vs0, vs3); 
      
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)];
      PROJ_P(vt2_c0, vs1, vs2); 
      
      // Color 1
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)];
      PROJ_P(vt1_c1, vs0, vs3); 
      
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)];
      PROJ_P(vt2_c1, vs1, vs2);
      
      // Color 2
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)];
      PROJ_P(vt1_c2, vs0, vs3);
      
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)];
      PROJ_P(vt2_c2, vs1, vs2);

      // Gauge Multiplication
      double2 u_val;
      double4 tmp_prod;
      double4 wt1_c0, wt2_c0, wt1_c1, wt2_c1, wt1_c2, wt2_c2;
      double bc2;
      double4 tmp_scaled;

      // --- Result Color 0 ---
      u_val.x = u_up[IDX2_G_R(0,0,isg)]; u_val.y = u_up[IDX2_G_I(0,0,isg)];
      wt1_c0 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c0 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(0,1,isg)]; u_val.y = u_up[IDX2_G_I(0,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      u_val.x = u_up[IDX2_G_R(0,2,isg)]; u_val.y = u_up[IDX2_G_I(0,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);

      // --- Result Color 1 ---
      u_val.x = u_up[IDX2_G_R(1,0,isg)]; u_val.y = u_up[IDX2_G_I(1,0,isg)];
      wt1_c1 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c1 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(1,1,isg)]; u_val.y = u_up[IDX2_G_I(1,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(1,2,isg)]; u_val.y = u_up[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);

      // --- Result Color 2 ---
      u_val.x = u_up[IDX2_G_R(2,0,isg)]; u_val.y = u_up[IDX2_G_I(2,0,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_up[IDX2_G_R(2,1,isg)]; u_val.y = u_up[IDX2_G_I(2,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(2,2,isg)]; u_val.y = u_up[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulation
      bc2 = 1.0;
      if(ix == Nx-1) bc2 = bc_x;

      // Color 0
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled);

      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);


      // mult_xm
      nn = (ix+Nx-1) % Nx;
      isn = nn + Nx * iyzt;
      isg = nn + Nx * iyzt + Nst_pad * idir;

      // Color 0
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)];
      PROJ_M(vt1_c0, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)];
      PROJ_M(vt2_c0, vs1, vs2); 
      // Color 1
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)];
      PROJ_M(vt1_c1, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)];
      PROJ_M(vt2_c1, vs1, vs2);
      // Color 2
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)];
      vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)];
      PROJ_M(vt1_c2, vs0, vs3);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)];
      vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)];
      PROJ_M(vt2_c2, vs1, vs2);

      // Gauge Mul (u_dn Dagger)
      // --- Result Color 0 ---
      u_val.x = u_dn[IDX2_G_R(0,0,isg)]; u_val.y = -u_dn[IDX2_G_I(0,0,isg)];
      wt1_c0 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c0 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,0,isg)]; u_val.y = -u_dn[IDX2_G_I(1,0,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,0,isg)]; u_val.y = -u_dn[IDX2_G_I(2,0,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, tmp_prod);

      // --- Result Color 1 ---
      u_val.x = u_dn[IDX2_G_R(0,1,isg)]; u_val.y = -u_dn[IDX2_G_I(0,1,isg)];
      wt1_c1 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c1 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,1,isg)]; u_val.y = -u_dn[IDX2_G_I(1,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,1,isg)]; u_val.y = -u_dn[IDX2_G_I(2,1,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, tmp_prod);

      // --- Result Color 2 ---
      u_val.x = u_dn[IDX2_G_R(0,2,isg)]; u_val.y = -u_dn[IDX2_G_I(0,2,isg)];
      wt1_c2 = qdw_mult_uc(u_val, vt1_c0);
      wt2_c2 = qdw_mult_uc(u_val, vt2_c0);
      u_val.x = u_dn[IDX2_G_R(1,2,isg)]; u_val.y = -u_dn[IDX2_G_I(1,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,2,isg)]; u_val.y = -u_dn[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulation
      bc2 = 1.0;
      if(ix == 0) bc2 = bc_x;

      // Color 0
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled);

      // Color 1
      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      // Color 2
      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);


      // =================================================================
      // Y-Direction
      // =================================================================
      
      // mult_yp
      nn = (iy + 1) % Ny;
      idir = 1;
      isn = ix + nn*Nx + izt*Nx*Ny;
      isg = ix + iy*Nx + izt*Nx*Ny + Nst_pad * idir;

      // Projection YP (1-sy)
      // sy = (0 0 0 -i), (0 0 i 0), (0 -i 0 0), (i 0 0 0)
      // vs0 + i vs3
      // vs1 - i vs2
      // vs2 + i vs1
      // vs3 - i vs0
      // NOTE: PROJ_P(a,b) is a+ib.
      // s0+is3 -> PROJ_P(t1, s0, s3)
      // s1-is2 -> PROJ_M(t2, s1, s2)
      // BUT output needs to be 2-spinor (t1, t2).
      // t1 = s0+is3.  t2 = s1-is2.
      // Reconstruction:
      // s0' = (1+sy) s0 ??
      // (1-sy)psi:
      // s0 + i s3. -> t1
      // s1 - i s2. -> t2
      // s2 + i s1. -> i(s1 - i s2) = i t2.
      // s3 - i s0. -> -i(s0 + i s3) = -i t1.
      
      // So s0=t1, s1=t2, s2=i*t2, s3=-i*t1.
      
      // Projection:
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)]; PROJ_P(vt1_c0, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)]; PROJ_M(vt2_c0, vs1, vs2); 
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)]; PROJ_P(vt1_c1, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)]; PROJ_M(vt2_c1, vs1, vs2);
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)]; PROJ_P(vt1_c2, vs0, vs3);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)]; PROJ_M(vt2_c2, vs1, vs2);

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
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(2,2,isg)]; u_val.y = u_up[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulate YP
      bc2 = 1.0; if(iy == Ny-1) bc2 = bc_y;
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled); // +i s1 -> +i t2
      MULT_MI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled); // -i s0 -> -i t1
      
      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);

      // mult_ym
      nn = (iy + Ny - 1) % Ny;
      isn = ix + nn*Nx + izt*Nx*Ny;
      isg = ix + nn*Nx + izt*Nx*Ny + Nst_pad * idir;

      // Projection YM (1+sy)
      // s0 - i s3. -> t1 (PROJ_M)
      // s1 + i s2. -> t2 (PROJ_P)
      // s2 - i s1 = -i(s1+is2) = -i t2. 
      // s3 + i s0 = i(s0-is3) = i t1.

      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)]; PROJ_M(vt1_c0, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)]; PROJ_P(vt2_c0, vs1, vs2); 
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)]; PROJ_M(vt1_c1, vs0, vs3); 
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)]; PROJ_P(vt2_c1, vs1, vs2);
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)]; PROJ_M(vt1_c2, vs0, vs3);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)]; PROJ_P(vt2_c2, vs1, vs2);

      // Gauge Mul (u_dn Dagger)
      // (Repeat gauge mul logic for 3 rows) - Copy paste from Xm but use loop/macro? 
      // I will write it explicitly.
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
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,2,isg)]; u_val.y = -u_dn[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulate YM
      bc2 = 1.0; if(iy == 0) bc2 = bc_y;
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled); // -i t2
      MULT_PI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled); // +i t1

      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);

      // =================================================================
      // Z-Direction
      // =================================================================
      
      // mult_zp
      idir = 2;
      nn = (iz + 1) % Nz;
      isn = ixy + Nx*Ny * (nn + it * Nz);
      isg = ixy + Nx*Ny * (iz + it * Nz) + Nst_pad * idir;

      // Projection ZP (1-gz)
      // gz = (0 0 -i 0), (0 0 0 i), (i 0 0 0), (0 -i 0 0)
      // s0 - (-i s2) = s0 + i s2. -> PROJ_P(t1, s0, s2).
      // s1 - (i s3) = s1 - i s3. -> PROJ_M(t2, s1, s3).
      // s2 - (i s0) = s2 - i s0 = -i(s0+is2) = -i t1.
      // s3 - (-i s1) = s3 + i s1 = i(s1-is3) = i t2.
      
      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)]; PROJ_P(vt1_c0, vs0, vs2); 
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)]; PROJ_M(vt2_c0, vs1, vs3); 
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)]; PROJ_P(vt1_c1, vs0, vs2); 
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)]; PROJ_M(vt2_c1, vs1, vs3);
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)]; PROJ_P(vt1_c2, vs0, vs2);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)]; PROJ_M(vt2_c2, vs1, vs3);

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
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_up[IDX2_G_R(2,2,isg)]; u_val.y = u_up[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulate ZP
      bc2 = 1.0; if(iz == Nz-1) bc2 = bc_z;
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled); // -i t1
      MULT_PI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled); // +i t2

      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_MI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_PI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);

      // mult_zm
      nn = (iz + Nz - 1) % Nz;
      isn = ixy + Nx*Ny * (nn + it * Nz);
      isg = ixy + Nx*Ny * (nn + it * Nz) + Nst_pad * idir;

      // Projection ZM (1+gz)
      // s0 + (-i s2) = s0 - i s2 -> PROJ_M(t1, s0, s2).
      // s1 + (i s3) = s1 + i s3 -> PROJ_P(t2, s1, s3).
      // s2 + (i s0) = i(s0 - i s2) = i t1.
      // s3 + (-i s1) = -i(s1+i s3) = -i t2.

      vs0 = ((double4*)v1)[IDX2_QDW(0, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(0, 2, isn)]; PROJ_M(vt1_c0, vs0, vs2); 
      vs1 = ((double4*)v1)[IDX2_QDW(0, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(0, 3, isn)]; PROJ_P(vt2_c0, vs1, vs3); 
      vs0 = ((double4*)v1)[IDX2_QDW(1, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(1, 2, isn)]; PROJ_M(vt1_c1, vs0, vs2); 
      vs1 = ((double4*)v1)[IDX2_QDW(1, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(1, 3, isn)]; PROJ_P(vt2_c1, vs1, vs3);
      vs0 = ((double4*)v1)[IDX2_QDW(2, 0, isn)]; vs2 = ((double4*)v1)[IDX2_QDW(2, 2, isn)]; PROJ_M(vt1_c2, vs0, vs2);
      vs1 = ((double4*)v1)[IDX2_QDW(2, 1, isn)]; vs3 = ((double4*)v1)[IDX2_QDW(2, 3, isn)]; PROJ_P(vt2_c2, vs1, vs3);

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
      tmp_prod = qdw_mult_uc(u_val, vt2_c1); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);
      u_val.x = u_dn[IDX2_G_R(2,2,isg)]; u_val.y = -u_dn[IDX2_G_I(2,2,isg)];
      tmp_prod = qdw_mult_uc(u_val, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, tmp_prod);
      tmp_prod = qdw_mult_uc(u_val, vt2_c2); QDW_ADD(wt2_c2, wt2_c1, tmp_prod);

      // Accumulate ZM
      bc2 = 1.0; if(iz == 0) bc2 = bc_z;
      QDW_SCAL(tmp_scaled, bc2, wt1_c0); QDW_ADD(v2_c0_s0, v2_c0_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c0); QDW_ADD(v2_c0_s1, v2_c0_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s2, v2_c0_s2, tmp_scaled); // +i t1
      MULT_MI(tmp_prod, wt2_c0); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c0_s3, v2_c0_s3, tmp_scaled); // -i t2

      QDW_SCAL(tmp_scaled, bc2, wt1_c1); QDW_ADD(v2_c1_s0, v2_c1_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c1); QDW_ADD(v2_c1_s1, v2_c1_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s2, v2_c1_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c1); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c1_s3, v2_c1_s3, tmp_scaled);

      QDW_SCAL(tmp_scaled, bc2, wt1_c2); QDW_ADD(v2_c2_s0, v2_c2_s0, tmp_scaled);
      QDW_SCAL(tmp_scaled, bc2, wt2_c2); QDW_ADD(v2_c2_s1, v2_c2_s1, tmp_scaled);
      MULT_PI(tmp_prod, wt1_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s2, v2_c2_s2, tmp_scaled);
      MULT_MI(tmp_prod, wt2_c2); QDW_SCAL(tmp_scaled, bc2, tmp_prod); QDW_ADD(v2_c2_s3, v2_c2_s3, tmp_scaled);
