/*!
      @file    mult_Domainwall_5din_eo_cuda_qdw-inc.h
      @brief   QDW Even-Odd Domain Wall fermion kernels (ee/eo 5dir + EO hopb)
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED

#include <type_traits>
#include "../inline/constant_memory_inline.h"

// All macros were undef'd by mult_Domainwall_5din_cuda_qdw-inc.h — redefine here.

// BLAS-compatible flat layout: NC*ND*Ns groups of real4 per NWP block.
// 4th arg is Ns (not Nst_pad) — matches the BLAS 4*IDX2(nin/4, in4, site) formula.
#define IDX_DWF_QDW(ic, id, is5, Ns_, site_) \
    IDX2(NC * ND * (Ns_), (ic) + NC * ((id) + ND * (is5)), (site_))

#define DWF_PROJ_P(res, a, b) \
    dw_add((a).x,(a).z,-(b).y,-(b).w,(res).x,(res).z); \
    dw_add((a).y,(a).w, (b).x, (b).z,(res).y,(res).w);

#define DWF_PROJ_M(res, a, b) \
    dw_add((a).x,(a).z, (b).y, (b).w,(res).x,(res).z); \
    dw_add((a).y,(a).w,-(b).x,-(b).z,(res).y,(res).w);

#define DWF_MULT_MI(res, a) \
    (res).x= (a).y; (res).y=-(a).x; (res).z= (a).w; (res).w=-(a).z;

#define DWF_MULT_PI(res, a) \
    (res).x=-(a).y; (res).y= (a).x; (res).z=-(a).w; (res).w= (a).z;

// Y-direction (gamma_2) helpers: combine two complex DD numbers with REAL
// +/-1, and reconstruct with REAL +/-1 (gamma_2 is real, unlike X/Z which use +/-i).
#define DWF_PROJ_RP(res, a, b) \
    dw_add((a).x,(a).z, (b).x, (b).z,(res).x,(res).z); \
    dw_add((a).y,(a).w, (b).y, (b).w,(res).y,(res).w);

#define DWF_PROJ_RM(res, a, b) \
    dw_add((a).x,(a).z,-(b).x,-(b).z,(res).x,(res).z); \
    dw_add((a).y,(a).w,-(b).y,-(b).w,(res).y,(res).w);

#define DWF_MULT_R1(res, a) \
    (res).x= (a).x; (res).y= (a).y; (res).z= (a).z; (res).w= (a).w;

#define DWF_MULT_RM1(res, a) \
    (res).x=-(a).x; (res).y=-(a).y; (res).z=-(a).z; (res).w=-(a).w;

#define DWF_PROJ_2(res, a) \
    dw_scal((real_t)2.0,(a).x,(a).z,(res).x,(res).z); \
    dw_scal((real_t)2.0,(a).y,(a).w,(res).y,(res).w);

// NOTE: gauge is contracted column-major (output[c] = sum_j U(j,c) vt[j]) to
// match the corelib FP hopb convention (forward link, plain U).
#define DWF_GMUL_FWD(u_ptr, isg_) \
{ \
    real2 _u; real4 _t; \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_uc(_u,vt1_c0); wt2_c0=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    wt1_c1=qdw_mult_uc(_u,vt1_c0); wt2_c1=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    wt1_c2=qdw_mult_uc(_u,vt1_c0); wt2_c2=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

// NOTE: gauge contracted row-major + conjugate (output[c] = sum_j conj(U(c,j)) vt[j])
// to match the corelib FP hopb convention (backward link, U^dagger).
#define DWF_GMUL_BCK(u_ptr, isg_) \
{ \
    real2 _u; real4 _t; \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_uc(_u,vt1_c0); wt2_c0=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    wt1_c1=qdw_mult_uc(_u,vt1_c0); wt2_c1=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    wt1_c2=qdw_mult_uc(_u,vt1_c0); wt2_c2=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

// Float-float (extended-precision) gauge variants: link = {hi from u_ptr,
// lo from ul_ptr}, multiplied into the DD vector via qdw_mult_cc.
//
// SU(3) 3rd-column reconstruction (RECON template arg, YAML su3_reconstruction):
// skip reading column 2 of each link from DRAM; rebuild U(:,2)=conj(U(:,0)xU(:,1))
// in dual-word. RECON is the kernel template param (constant-folded), so when off
// these branches dead-strip and the fast path is byte-identical.
#define DWF_LOAD_FF(_u, u_ptr, ul_ptr, a, b, isg_) \
    (_u).x=(u_ptr)[IDX2_G_R(a,b,(isg_))]; (_u).y=(u_ptr)[IDX2_G_I(a,b,(isg_))]; \
    (_u).z=(ul_ptr)[IDX2_G_R(a,b,(isg_))]; (_u).w=(ul_ptr)[IDX2_G_I(a,b,(isg_))];
// U(i,2)=conj(U(p,0)U(q,1)-U(q,0)U(p,1)), p=(i+1)%3, q=(i+2)%3. conj -> negate y,w.
#define DWF_RECON_COL2_FF(u_ptr, ul_ptr, isg_, rc0, rc1, rc2) \
{ \
    real4 _ra,_rb,_rc,_rd,_p1,_p2,_df; \
    DWF_LOAD_FF(_ra,u_ptr,ul_ptr,1,0,isg_) DWF_LOAD_FF(_rb,u_ptr,ul_ptr,2,1,isg_) \
    DWF_LOAD_FF(_rc,u_ptr,ul_ptr,2,0,isg_) DWF_LOAD_FF(_rd,u_ptr,ul_ptr,1,1,isg_) \
    _p1=qdw_mult_cc(_ra,_rb); _p2=qdw_mult_cc(_rc,_rd); QDW_SUB(_df,_p1,_p2); \
    (rc0).x=_df.x; (rc0).y=-_df.y; (rc0).z=_df.z; (rc0).w=-_df.w; \
    DWF_LOAD_FF(_ra,u_ptr,ul_ptr,2,0,isg_) DWF_LOAD_FF(_rb,u_ptr,ul_ptr,0,1,isg_) \
    DWF_LOAD_FF(_rc,u_ptr,ul_ptr,0,0,isg_) DWF_LOAD_FF(_rd,u_ptr,ul_ptr,2,1,isg_) \
    _p1=qdw_mult_cc(_ra,_rb); _p2=qdw_mult_cc(_rc,_rd); QDW_SUB(_df,_p1,_p2); \
    (rc1).x=_df.x; (rc1).y=-_df.y; (rc1).z=_df.z; (rc1).w=-_df.w; \
    DWF_LOAD_FF(_ra,u_ptr,ul_ptr,0,0,isg_) DWF_LOAD_FF(_rb,u_ptr,ul_ptr,1,1,isg_) \
    DWF_LOAD_FF(_rc,u_ptr,ul_ptr,1,0,isg_) DWF_LOAD_FF(_rd,u_ptr,ul_ptr,0,1,isg_) \
    _p1=qdw_mult_cc(_ra,_rb); _p2=qdw_mult_cc(_rc,_rd); QDW_SUB(_df,_p1,_p2); \
    (rc2).x=_df.x; (rc2).y=-_df.y; (rc2).z=_df.z; (rc2).w=-_df.w; \
}
#define DWF_GMUL_FWD_FF(u_ptr, ul_ptr, isg_) \
{ \
    real4 _u; real4 _t; \
    real4 _rc0,_rc1,_rc2; \
    if (RECON) { DWF_RECON_COL2_FF(u_ptr, ul_ptr, isg_, _rc0, _rc1, _rc2) } \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,0,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_cc(_u,vt1_c0); wt2_c0=qdw_mult_cc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,0,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(1,0,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,0,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(2,0,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,1,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(0,1,(isg_))]; \
    wt1_c1=qdw_mult_cc(_u,vt1_c0); wt2_c1=qdw_mult_cc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,1,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,1,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    if (RECON) { _u=_rc0; } else { _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,2,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(0,2,(isg_))]; } \
    wt1_c2=qdw_mult_cc(_u,vt1_c0); wt2_c2=qdw_mult_cc(_u,vt2_c0); \
    if (RECON) { _u=_rc1; } else { _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,2,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(1,2,(isg_))]; } \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    if (RECON) { _u=_rc2; } else { _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,2,(isg_))]; _u.w=(ul_ptr)[IDX2_G_I(2,2,(isg_))]; } \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

#define DWF_GMUL_BCK_FF(u_ptr, ul_ptr, isg_) \
{ \
    real4 _u; real4 _t; \
    real4 _rc0,_rc1,_rc2; \
    if (RECON) { DWF_RECON_COL2_FF(u_ptr, ul_ptr, isg_, _rc0, _rc1, _rc2) } \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,0,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_cc(_u,vt1_c0); wt2_c0=qdw_mult_cc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,1,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(0,1,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    if (RECON) { _u.x=_rc0.x; _u.y=-_rc0.y; _u.z=_rc0.z; _u.w=-_rc0.w; } else { _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(0,2,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(0,2,(isg_))]; } \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,0,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(1,0,(isg_))]; \
    wt1_c1=qdw_mult_cc(_u,vt1_c0); wt2_c1=qdw_mult_cc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,1,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    if (RECON) { _u.x=_rc1.x; _u.y=-_rc1.y; _u.z=_rc1.z; _u.w=-_rc1.w; } else { _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(1,2,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(1,2,(isg_))]; } \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,0,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,0,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(2,0,(isg_))]; \
    wt1_c2=qdw_mult_cc(_u,vt1_c0); wt2_c2=qdw_mult_cc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,1,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,1,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_cc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_cc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    if (RECON) { _u.x=_rc2.x; _u.y=-_rc2.y; _u.z=_rc2.z; _u.w=-_rc2.w; } else { _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,2,(isg_))]; _u.z=(ul_ptr)[IDX2_G_R(2,2,(isg_))]; _u.w=-(ul_ptr)[IDX2_G_I(2,2,(isg_))]; } \
    _t=qdw_mult_cc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_cc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

#define DWF_ACCUM_4(v0,v1,v2,v3, W1,W2, bc_, OP2,OP3) \
{ \
    real4 _ts,_tp; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
    OP2(_tp,W2); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v2,v2,_ts); \
    OP3(_tp,W1); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v3,v3,_ts); \
}

#define DWF_ACCUM_4_SW(v0,v1,v2,v3, W1,W2, bc_, OP2,OP3) \
{ \
    real4 _ts,_tp; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
    OP2(_tp,W1); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v2,v2,_ts); \
    OP3(_tp,W2); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v3,v3,_ts); \
}

#define DWF_ACCUM_TP(v2,v3, W1,W2, bc_) \
{ \
    real4 _ts; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v2,v2,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v3,v3,_ts); \
}
#define DWF_ACCUM_TM(v0,v1, W1,W2, bc_) \
{ \
    real4 _ts; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
}

#define DWF_LOAD_PROJ(wp_,is_,Np_,isn_,id_a,id_b,PROJ) \
{ \
    real4 _sa,_sb; \
    _sa=(wp_)[IDX_DWF_QDW(0,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(0,id_b,is_,Np_,isn_)]; PROJ(vt1_c0,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(1,id_b,is_,Np_,isn_)]; PROJ(vt1_c1,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(2,id_b,is_,Np_,isn_)]; PROJ(vt1_c2,_sa,_sb); \
}
#define DWF_LOAD_PROJ2(wp_,is_,Np_,isn_,id_a,id_b,PROJ) \
{ \
    real4 _sa,_sb; \
    _sa=(wp_)[IDX_DWF_QDW(0,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(0,id_b,is_,Np_,isn_)]; PROJ(vt2_c0,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(1,id_b,is_,Np_,isn_)]; PROJ(vt2_c1,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(2,id_b,is_,Np_,isn_)]; PROJ(vt2_c2,_sa,_sb); \
}
#define DWF_LOAD_PROJ_T(wp_,is_,Np_,isn_,id_a,id_b) \
{ \
    real4 _sa; \
    _sa=(wp_)[IDX_DWF_QDW(0,id_a,is_,Np_,isn_)]; DWF_PROJ_2(vt1_c0,_sa); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_a,is_,Np_,isn_)]; DWF_PROJ_2(vt1_c1,_sa); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_a,is_,Np_,isn_)]; DWF_PROJ_2(vt1_c2,_sa); \
    _sa=(wp_)[IDX_DWF_QDW(0,id_b,is_,Np_,isn_)]; DWF_PROJ_2(vt2_c0,_sa); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_b,is_,Np_,isn_)]; DWF_PROJ_2(vt2_c1,_sa); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_b,is_,Np_,isn_)]; DWF_PROJ_2(vt2_c2,_sa); \
}

// Helper: compute QDW of fac1*a + fac2*b
#define DWF_SCAL_ADD2(res, fac1_, a_, fac2_, b_) \
{ \
    real4 _t1, _t2; \
    QDW_SCAL(_t1, fac1_, a_); QDW_SCAL(_t2, fac2_, b_); QDW_ADD(res, _t1, _t2); \
}

// DD-paired scalar version: fac1 = (f1h, f1l), fac2 = (f2h, f2l). All FP32-only.
#define DWF_SCAL_ADD2_DD(res, f1h_, f1l_, a_, f2h_, f2l_, b_) \
{ \
    real4 _t1, _t2; \
    QDW_SCAL_DD(_t1, f1h_, f1l_, a_); QDW_SCAL_DD(_t2, f2h_, f2l_, b_); QDW_ADD(res, _t1, _t2); \
}

//====================================================================
// ee 5-direction (diagonal block)  vp = B_eff * wp + C_hop * vt
// (single-T variant — used by EXT=false and the double build.)
//====================================================================
__global__
void mult_domainwall_5din_ee_5dir_dirac_qdw_kernel_fp(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq, real_t M0, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();

    for (int is = 0; is < Ns; ++is) {
        const real_t FF1 = b_con[is] * (real_t(4.0) - M0) + real_t(1.0);
        const real_t FF2 = c_con[is] * (real_t(4.0) - M0) - real_t(1.0);

        // ---- s+1 hopping ----
        const int is_up = (is + 1) % Ns;
        const real_t Fup = (is == Ns-1) ? real_t(-0.5) * mq * FF2 : real_t(0.5) * FF2 * alpha;
        real4 tmp;

#define LOAD_SUB_SHop(ic) \
        real4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        real4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL(vt_##ic##_s0,Fup,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL(vt_##ic##_s1,Fup,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_SHop(0) LOAD_SUB_SHop(1) LOAD_SUB_SHop(2)
#undef LOAD_SUB_SHop

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const real_t Fdn = (is == 0) ? real_t(-0.5) * mq * FF2 : real_t(0.5) * FF2 * alpha;
        real4 sum02, sum13;

#define ADD_SUM_SHop(ic) \
        { \
        real4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_ADD(sum02,wd_s0,wd_s2); QDW_SCAL(tmp,Fdn,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,wd_s1,wd_s3); QDW_SCAL(tmp,Fdn,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SUM_SHop(0) ADD_SUM_SHop(1) ADD_SUM_SHop(2)
#undef ADD_SUM_SHop

        // ---- local diagonal ----
#define DIAG_AND_STORE(ic) \
        { \
        real4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0) { \
            real_t f1=FF1*real_t(0.5)*(real_t(1.0)+alpha), f2=FF1*real_t(0.5)*(real_t(-1.0)+alpha); \
            DWF_SCAL_ADD2(wt_s0,f1,w_s0,f2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,f1,w_s1,f2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,f1,w_s2,f2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,f1,w_s3,f2,w_s1); \
        } else if (is == Ns-1) { \
            real_t f1=FF1*real_t(0.5)*(real_t(1.0)+alpha), f2=FF1*real_t(0.5)*(real_t(1.0)-alpha); \
            DWF_SCAL_ADD2(wt_s0,f1,w_s0,f2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,f1,w_s1,f2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,f1,w_s2,f2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,f1,w_s3,f2,w_s1); \
        } else { \
            real_t fa=FF1*alpha; \
            QDW_SCAL(wt_s0,fa,w_s0); QDW_SCAL(wt_s1,fa,w_s1); \
            QDW_SCAL(wt_s2,fa,w_s2); QDW_SCAL(wt_s3,fa,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)]=wt_s0; \
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)]=wt_s1; \
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)]=wt_s2; \
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)]=wt_s3; \
        }

        DIAG_AND_STORE(0) DIAG_AND_STORE(1) DIAG_AND_STORE(2)
#undef DIAG_AND_STORE
    }
}

// Forward declaration (definition follows the host wrapper below).
__global__
void mult_domainwall_5din_ee_5dir_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, real_t M0_h, real_t M0_l,
    int Ns, real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad);

void mult_domainwall_5din_ee_5dir_dirac_qdw(
    real_t *vp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real4 *vp_dev = (real4 *)dev_ptr(vp);
    real4 *wp_dev = (real4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    if (ext && sizeof(real_t) == sizeof(float)) {
      const real_t mq_h = (real_t)mq;
      const real_t mq_l = (real_t)(mq - (double)mq_h);
      const real_t M0_h = (real_t)M0;
      const real_t M0_l = (real_t)(M0 - (double)M0_h);
      const real_t a_h  = (real_t)alpha;
      const real_t a_l  = (real_t)(alpha - (double)a_h);
      mult_domainwall_5din_ee_5dir_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
          vp_dev, wp_dev, mq_h, mq_l, M0_h, M0_l, Ns, a_h, a_l, Nst, Nst_pad);
    } else {
      mult_domainwall_5din_ee_5dir_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
          vp_dev, wp_dev, (real_t)mq, (real_t)M0, Ns, (real_t)alpha, Nst, Nst_pad);
    }
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// ee 5-direction (FP32-only DD variant — float build, extended_precision=true).
// All coefficients (b, c, mq, M0, alpha) are float pairs; arithmetic uses
// dw_add/dw_mul/dw_scal/QDW_SCAL_DD — no FP64 ALU.
//====================================================================
__global__
void mult_domainwall_5din_ee_5dir_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, real_t M0_h, real_t M0_l,
    int Ns, real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_h = ConstantMemoryTraits<real_t>::b();
    const real_t *b_l = ConstantMemoryTraits<real_t>::b_lo();
    const real_t *c_h = ConstantMemoryTraits<real_t>::c();
    const real_t *c_l = ConstantMemoryTraits<real_t>::c_lo();

    // (4 - M0): DD
    real_t fmM0_h, fmM0_l;
    dw_add(real_t(4.0), real_t(0.0), -M0_h, -M0_l, fmM0_h, fmM0_l);
    // -0.5 * mq: DD
    real_t mhmq_h, mhmq_l;
    dw_scal(real_t(-0.5), mq_h, mq_l, mhmq_h, mhmq_l);

    for (int is = 0; is < Ns; ++is) {
        // FF1 = b[is]*(4-M0) + 1
        real_t bm_h, bm_l;
        dw_mul(b_h[is], b_l[is], fmM0_h, fmM0_l, bm_h, bm_l);
        real_t FF1_h, FF1_l;
        dw_add(bm_h, bm_l, real_t(1.0), real_t(0.0), FF1_h, FF1_l);
        // FF2 = c[is]*(4-M0) - 1
        real_t cm_h, cm_l;
        dw_mul(c_h[is], c_l[is], fmM0_h, fmM0_l, cm_h, cm_l);
        real_t FF2_h, FF2_l;
        dw_add(cm_h, cm_l, real_t(-1.0), real_t(0.0), FF2_h, FF2_l);

        // ---- s+1 hopping ----
        const int is_up = (is + 1) % Ns;
        real_t Fup_h, Fup_l;
        if (is == Ns-1) {
            dw_mul(mhmq_h, mhmq_l, FF2_h, FF2_l, Fup_h, Fup_l);
        } else {
            real_t hFF2_h, hFF2_l;
            dw_scal(real_t(0.5), FF2_h, FF2_l, hFF2_h, hFF2_l);
            dw_mul(hFF2_h, hFF2_l, alpha_h, alpha_l, Fup_h, Fup_l);
        }
        real4 tmp;

#define LOAD_SUB_SHop_FF(ic) \
        real4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        real4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL_DD(vt_##ic##_s0,Fup_h,Fup_l,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL_DD(vt_##ic##_s1,Fup_h,Fup_l,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_SHop_FF(0) LOAD_SUB_SHop_FF(1) LOAD_SUB_SHop_FF(2)
#undef LOAD_SUB_SHop_FF

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn_h, Fdn_l;
        if (is == 0) {
            dw_mul(mhmq_h, mhmq_l, FF2_h, FF2_l, Fdn_h, Fdn_l);
        } else {
            real_t hFF2_h, hFF2_l;
            dw_scal(real_t(0.5), FF2_h, FF2_l, hFF2_h, hFF2_l);
            dw_mul(hFF2_h, hFF2_l, alpha_h, alpha_l, Fdn_h, Fdn_l);
        }
        real4 sum02, sum13;

#define ADD_SUM_SHop_FF(ic) \
        { \
        real4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_ADD(sum02,wd_s0,wd_s2); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,wd_s1,wd_s3); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SUM_SHop_FF(0) ADD_SUM_SHop_FF(1) ADD_SUM_SHop_FF(2)
#undef ADD_SUM_SHop_FF

        // ---- local diagonal: precompute alpha-derived factors ----
        real_t f1_h=0, f1_l=0, f2_h=0, f2_l=0, fa_h=0, fa_l=0;
        if (is == 0) {
            real_t opa_h, opa_l, mpa_h, mpa_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0),  real_t(0.0), alpha_h, alpha_l, opa_h, opa_l);
            dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, mpa_h, mpa_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), mpa_h, mpa_l, hma_h, hma_l);
            dw_mul(FF1_h, FF1_l, hpa_h, hpa_l, f1_h, f1_l);
            dw_mul(FF1_h, FF1_l, hma_h, hma_l, f2_h, f2_l);
        } else if (is == Ns-1) {
            real_t opa_h, opa_l, oma_h, oma_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0), real_t(0.0),  alpha_h,  alpha_l, opa_h, opa_l);
            dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, oma_h, oma_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), oma_h, oma_l, hma_h, hma_l);
            dw_mul(FF1_h, FF1_l, hpa_h, hpa_l, f1_h, f1_l);
            dw_mul(FF1_h, FF1_l, hma_h, hma_l, f2_h, f2_l);
        } else {
            dw_mul(FF1_h, FF1_l, alpha_h, alpha_l, fa_h, fa_l);
        }

#define DIAG_AND_STORE_FF(ic) \
        { \
        real4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0 || is == Ns-1) { \
            DWF_SCAL_ADD2_DD(wt_s0,f1_h,f1_l,w_s0,f2_h,f2_l,w_s2); \
            DWF_SCAL_ADD2_DD(wt_s1,f1_h,f1_l,w_s1,f2_h,f2_l,w_s3); \
            DWF_SCAL_ADD2_DD(wt_s2,f1_h,f1_l,w_s2,f2_h,f2_l,w_s0); \
            DWF_SCAL_ADD2_DD(wt_s3,f1_h,f1_l,w_s3,f2_h,f2_l,w_s1); \
        } else { \
            QDW_SCAL_DD(wt_s0,fa_h,fa_l,w_s0); QDW_SCAL_DD(wt_s1,fa_h,fa_l,w_s1); \
            QDW_SCAL_DD(wt_s2,fa_h,fa_l,w_s2); QDW_SCAL_DD(wt_s3,fa_h,fa_l,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)]=wt_s0; \
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)]=wt_s1; \
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)]=wt_s2; \
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)]=wt_s3; \
        }

        DIAG_AND_STORE_FF(0) DIAG_AND_STORE_FF(1) DIAG_AND_STORE_FF(2)
#undef DIAG_AND_STORE_FF
    }
}

//====================================================================
// eo 5-direction (off-diagonal block)  yp = -0.5*(b*w_diag + c*vt)
// (single-T variant.)
//====================================================================
__global__
void mult_domainwall_5din_eo_5dir_dirac_qdw_kernel_fp(
    real4 * __restrict__ yp, real4 * __restrict__ wp,
    real_t mq, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();

    for (int is = 0; is < Ns; ++is) {
        real4 tmp;

        // ---- s+1 hopping (C-scaled) ----
        const int is_up = (is + 1) % Ns;
        const real_t Fup = (is == Ns-1) ? real_t(-0.5) * mq * c_con[is] : real_t(0.5) * c_con[is] * alpha;

#define LOAD_SUB_EO(ic) \
        real4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        real4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL(vt_##ic##_s0,Fup,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL(vt_##ic##_s1,Fup,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_EO(0) LOAD_SUB_EO(1) LOAD_SUB_EO(2)
#undef LOAD_SUB_EO

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const real_t Fdn = (is == 0) ? real_t(-0.5) * mq * c_con[is] : real_t(0.5) * c_con[is] * alpha;
        real4 sum02, sum13;

#define ADD_SUM_EO(ic) \
        { \
        real4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_ADD(sum02,wd_s0,wd_s2); QDW_SCAL(tmp,Fdn,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,wd_s1,wd_s3); QDW_SCAL(tmp,Fdn,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SUM_EO(0) ADD_SUM_EO(1) ADD_SUM_EO(2)
#undef ADD_SUM_EO

        // ---- local B-scaled diagonal ----
#define DIAG_AND_STORE_EO(ic) \
        { \
        real4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0) { \
            real_t b1=b_con[is]*real_t(0.5)*(real_t(1.0)+alpha), b2=b_con[is]*real_t(0.5)*(real_t(-1.0)+alpha); \
            DWF_SCAL_ADD2(wt_s0,b1,w_s0,b2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,b1,w_s1,b2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,b1,w_s2,b2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,b1,w_s3,b2,w_s1); \
        } else if (is == Ns-1) { \
            real_t b1=b_con[is]*real_t(0.5)*(real_t(1.0)+alpha), b2=b_con[is]*real_t(0.5)*(real_t(1.0)-alpha); \
            DWF_SCAL_ADD2(wt_s0,b1,w_s0,b2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,b1,w_s1,b2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,b1,w_s2,b2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,b1,w_s3,b2,w_s1); \
        } else { \
            real_t bb=b_con[is]*alpha; \
            QDW_SCAL(wt_s0,bb,w_s0); QDW_SCAL(wt_s1,bb,w_s1); \
            QDW_SCAL(wt_s2,bb,w_s2); QDW_SCAL(wt_s3,bb,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        QDW_SCAL(wt_s0,real_t(-0.5),wt_s0); \
        QDW_SCAL(wt_s1,real_t(-0.5),wt_s1); \
        QDW_SCAL(wt_s2,real_t(-0.5),wt_s2); \
        QDW_SCAL(wt_s3,real_t(-0.5),wt_s3); \
        yp[IDX_DWF_QDW(ic,0,is,Ns,site)]=wt_s0; \
        yp[IDX_DWF_QDW(ic,1,is,Ns,site)]=wt_s1; \
        yp[IDX_DWF_QDW(ic,2,is,Ns,site)]=wt_s2; \
        yp[IDX_DWF_QDW(ic,3,is,Ns,site)]=wt_s3; \
        }

        DIAG_AND_STORE_EO(0) DIAG_AND_STORE_EO(1) DIAG_AND_STORE_EO(2)
#undef DIAG_AND_STORE_EO
    }
}

// Forward declaration (definition follows the host wrapper below).
__global__
void mult_domainwall_5din_eo_5dir_dirac_qdw_kernel_ff(
    real4 * __restrict__ yp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, int Ns,
    real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad);

void mult_domainwall_5din_eo_5dir_dirac_qdw(
    real_t *yp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real4 *yp_dev = (real4 *)dev_ptr(yp);
    real4 *wp_dev = (real4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    if (ext && sizeof(real_t) == sizeof(float)) {
      const real_t mq_h = (real_t)mq;
      const real_t mq_l = (real_t)(mq - (double)mq_h);
      const real_t a_h  = (real_t)alpha;
      const real_t a_l  = (real_t)(alpha - (double)a_h);
      mult_domainwall_5din_eo_5dir_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
          yp_dev, wp_dev, mq_h, mq_l, Ns, a_h, a_l, Nst, Nst_pad);
    } else {
      mult_domainwall_5din_eo_5dir_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
          yp_dev, wp_dev, (real_t)mq, Ns, (real_t)alpha, Nst, Nst_pad);
    }
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-direction (FP32-only DD variant). yp = -0.5*(b*w_diag + c*vt)
//====================================================================
__global__
void mult_domainwall_5din_eo_5dir_dirac_qdw_kernel_ff(
    real4 * __restrict__ yp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, int Ns,
    real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_h = ConstantMemoryTraits<real_t>::b();
    const real_t *b_l = ConstantMemoryTraits<real_t>::b_lo();
    const real_t *c_h = ConstantMemoryTraits<real_t>::c();
    const real_t *c_l = ConstantMemoryTraits<real_t>::c_lo();

    // -0.5 * mq: DD
    real_t mhmq_h, mhmq_l;
    dw_scal(real_t(-0.5), mq_h, mq_l, mhmq_h, mhmq_l);

    for (int is = 0; is < Ns; ++is) {
        real4 tmp;

        // ---- s+1 hopping (C-scaled) ----
        const int is_up = (is + 1) % Ns;
        real_t Fup_h, Fup_l;
        if (is == Ns-1) {
            dw_mul(mhmq_h, mhmq_l, c_h[is], c_l[is], Fup_h, Fup_l);
        } else {
            real_t hc_h, hc_l;
            dw_scal(real_t(0.5), c_h[is], c_l[is], hc_h, hc_l);
            dw_mul(hc_h, hc_l, alpha_h, alpha_l, Fup_h, Fup_l);
        }

#define LOAD_SUB_EO_FF(ic) \
        real4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        real4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL_DD(vt_##ic##_s0,Fup_h,Fup_l,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL_DD(vt_##ic##_s1,Fup_h,Fup_l,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_EO_FF(0) LOAD_SUB_EO_FF(1) LOAD_SUB_EO_FF(2)
#undef LOAD_SUB_EO_FF

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn_h, Fdn_l;
        if (is == 0) {
            dw_mul(mhmq_h, mhmq_l, c_h[is], c_l[is], Fdn_h, Fdn_l);
        } else {
            real_t hc_h, hc_l;
            dw_scal(real_t(0.5), c_h[is], c_l[is], hc_h, hc_l);
            dw_mul(hc_h, hc_l, alpha_h, alpha_l, Fdn_h, Fdn_l);
        }
        real4 sum02, sum13;

#define ADD_SUM_EO_FF(ic) \
        { \
        real4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_ADD(sum02,wd_s0,wd_s2); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,wd_s1,wd_s3); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SUM_EO_FF(0) ADD_SUM_EO_FF(1) ADD_SUM_EO_FF(2)
#undef ADD_SUM_EO_FF

        // ---- local B-scaled diagonal: precompute alpha-derived factors ----
        real_t b1_h=0, b1_l=0, b2_h=0, b2_l=0, bb_h=0, bb_l=0;
        if (is == 0) {
            real_t opa_h, opa_l, mpa_h, mpa_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0),  real_t(0.0), alpha_h, alpha_l, opa_h, opa_l);
            dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, mpa_h, mpa_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), mpa_h, mpa_l, hma_h, hma_l);
            dw_mul(b_h[is], b_l[is], hpa_h, hpa_l, b1_h, b1_l);
            dw_mul(b_h[is], b_l[is], hma_h, hma_l, b2_h, b2_l);
        } else if (is == Ns-1) {
            real_t opa_h, opa_l, oma_h, oma_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0), real_t(0.0),  alpha_h,  alpha_l, opa_h, opa_l);
            dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, oma_h, oma_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), oma_h, oma_l, hma_h, hma_l);
            dw_mul(b_h[is], b_l[is], hpa_h, hpa_l, b1_h, b1_l);
            dw_mul(b_h[is], b_l[is], hma_h, hma_l, b2_h, b2_l);
        } else {
            dw_mul(b_h[is], b_l[is], alpha_h, alpha_l, bb_h, bb_l);
        }

#define DIAG_AND_STORE_EO_FF(ic) \
        { \
        real4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0 || is == Ns-1) { \
            DWF_SCAL_ADD2_DD(wt_s0,b1_h,b1_l,w_s0,b2_h,b2_l,w_s2); \
            DWF_SCAL_ADD2_DD(wt_s1,b1_h,b1_l,w_s1,b2_h,b2_l,w_s3); \
            DWF_SCAL_ADD2_DD(wt_s2,b1_h,b1_l,w_s2,b2_h,b2_l,w_s0); \
            DWF_SCAL_ADD2_DD(wt_s3,b1_h,b1_l,w_s3,b2_h,b2_l,w_s1); \
        } else { \
            QDW_SCAL_DD(wt_s0,bb_h,bb_l,w_s0); QDW_SCAL_DD(wt_s1,bb_h,bb_l,w_s1); \
            QDW_SCAL_DD(wt_s2,bb_h,bb_l,w_s2); QDW_SCAL_DD(wt_s3,bb_h,bb_l,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        QDW_SCAL(wt_s0,real_t(-0.5),wt_s0); \
        QDW_SCAL(wt_s1,real_t(-0.5),wt_s1); \
        QDW_SCAL(wt_s2,real_t(-0.5),wt_s2); \
        QDW_SCAL(wt_s3,real_t(-0.5),wt_s3); \
        yp[IDX_DWF_QDW(ic,0,is,Ns,site)]=wt_s0; \
        yp[IDX_DWF_QDW(ic,1,is,Ns,site)]=wt_s1; \
        yp[IDX_DWF_QDW(ic,2,is,Ns,site)]=wt_s2; \
        yp[IDX_DWF_QDW(ic,3,is,Ns,site)]=wt_s3; \
        }

        DIAG_AND_STORE_EO_FF(0) DIAG_AND_STORE_EO_FF(1) DIAG_AND_STORE_EO_FF(2)
#undef DIAG_AND_STORE_EO_FF
    }
}

//====================================================================
// ee 5-dirdag (adjoint diagonal block) — single-T variant.
//====================================================================
__global__
void mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel_fp(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq, real_t M0, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();

    for (int is = 0; is < Ns; ++is) {
        const real_t B1 = b_con[is] * (real_t(4.0) - M0) + real_t(1.0);
        real4 tmp;

        // ---- local diagonal (self-adjoint) ----
#define DIAG_DAG(ic) \
        real4 w_##ic##_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_##ic##_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_##ic##_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_##ic##_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0) { \
            real_t f1=B1*real_t(0.5)*(real_t(1.0)+alpha), f2=B1*real_t(0.5)*(real_t(-1.0)+alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,f1,w_##ic##_s0,f2,w_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s1,f1,w_##ic##_s1,f2,w_##ic##_s3); \
            DWF_SCAL_ADD2(vt_##ic##_s2,f1,w_##ic##_s2,f2,w_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s3,f1,w_##ic##_s3,f2,w_##ic##_s1); \
        } else if (is == Ns-1) { \
            real_t f1=B1*real_t(0.5)*(real_t(1.0)+alpha), f2=B1*real_t(0.5)*(real_t(1.0)-alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,f1,w_##ic##_s0,f2,w_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s1,f1,w_##ic##_s1,f2,w_##ic##_s3); \
            DWF_SCAL_ADD2(vt_##ic##_s2,f1,w_##ic##_s2,f2,w_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s3,f1,w_##ic##_s3,f2,w_##ic##_s1); \
        } else { \
            real_t fa=B1*alpha; \
            QDW_SCAL(vt_##ic##_s0,fa,w_##ic##_s0); QDW_SCAL(vt_##ic##_s1,fa,w_##ic##_s1); \
            QDW_SCAL(vt_##ic##_s2,fa,w_##ic##_s2); QDW_SCAL(vt_##ic##_s3,fa,w_##ic##_s3); \
        }

        DIAG_DAG(0) DIAG_DAG(1) DIAG_DAG(2)
#undef DIAG_DAG

        // ---- s+1 contribution (adjoint: symmetric) ----
        const int is_up = (is + 1) % Ns;
        const real_t C2up = c_con[is_up] * (real_t(4.0) - M0) - real_t(1.0);
        const real_t Fup = (is == Ns-1) ? real_t(-0.5) * mq * C2up * real_t(0.5) : real_t(0.5) * C2up * alpha * real_t(0.5);
        real4 sum02, sum13;

#define ADD_SYM_UP(ic) \
        { \
        real4 xu_s0=wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 xu_s1=wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 xu_s2=wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 xu_s3=wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        QDW_ADD(sum02,xu_s0,xu_s2); QDW_SCAL(tmp,Fup,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,xu_s1,xu_s3); QDW_SCAL(tmp,Fup,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SYM_UP(0) ADD_SYM_UP(1) ADD_SYM_UP(2)
#undef ADD_SYM_UP

        // ---- s-1 contribution (adjoint: antisymmetric) ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const real_t C2dn = c_con[is_dn] * (real_t(4.0) - M0) - real_t(1.0);
        const real_t Fdn = (is == 0) ? real_t(-0.5) * mq * C2dn * real_t(0.5) : real_t(0.5) * C2dn * alpha * real_t(0.5);

#define ADD_ANTISYM_DN(ic) \
        { \
        real4 xd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 xd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 xd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 xd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_SUB(sum02,xd_s0,xd_s2); QDW_SCAL(tmp,Fdn,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,xd_s1,xd_s3); QDW_SCAL(tmp,Fdn,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_ANTISYM_DN(0) ADD_ANTISYM_DN(1) ADD_ANTISYM_DN(2)
#undef ADD_ANTISYM_DN

        vp[IDX_DWF_QDW(0,0,is,Ns,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Ns,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Ns,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Ns,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Ns,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Ns,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Ns,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Ns,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Ns,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Ns,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Ns,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Ns,site)]=vt_2_s3;
    }
}

// Forward declaration (definition follows the host wrapper below).
__global__
void mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, real_t M0_h, real_t M0_l,
    int Ns, real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad);

void mult_domainwall_5din_ee_5dirdag_dirac_qdw(
    real_t *vp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real4 *vp_dev = (real4 *)dev_ptr(vp);
    real4 *wp_dev = (real4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    if (ext && sizeof(real_t) == sizeof(float)) {
      const real_t mq_h = (real_t)mq;
      const real_t mq_l = (real_t)(mq - (double)mq_h);
      const real_t M0_h = (real_t)M0;
      const real_t M0_l = (real_t)(M0 - (double)M0_h);
      const real_t a_h  = (real_t)alpha;
      const real_t a_l  = (real_t)(alpha - (double)a_h);
      mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
          vp_dev, wp_dev, mq_h, mq_l, M0_h, M0_l, Ns, a_h, a_l, Nst, Nst_pad);
    } else {
      mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
          vp_dev, wp_dev, (real_t)mq, (real_t)M0, Ns, (real_t)alpha, Nst, Nst_pad);
    }
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// ee 5-dirdag (FP32-only DD variant).
//====================================================================
__global__
void mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ wp,
    real_t mq_h, real_t mq_l, real_t M0_h, real_t M0_l,
    int Ns, real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_h = ConstantMemoryTraits<real_t>::b();
    const real_t *b_l = ConstantMemoryTraits<real_t>::b_lo();
    const real_t *c_h = ConstantMemoryTraits<real_t>::c();
    const real_t *c_l = ConstantMemoryTraits<real_t>::c_lo();

    // (4 - M0): DD
    real_t fmM0_h, fmM0_l;
    dw_add(real_t(4.0), real_t(0.0), -M0_h, -M0_l, fmM0_h, fmM0_l);
    // -0.5 * mq: DD
    real_t mhmq_h, mhmq_l;
    dw_scal(real_t(-0.5), mq_h, mq_l, mhmq_h, mhmq_l);

    for (int is = 0; is < Ns; ++is) {
        // B1 = b[is]*(4-M0) + 1
        real_t bm_h, bm_l;
        dw_mul(b_h[is], b_l[is], fmM0_h, fmM0_l, bm_h, bm_l);
        real_t B1_h, B1_l;
        dw_add(bm_h, bm_l, real_t(1.0), real_t(0.0), B1_h, B1_l);
        real4 tmp;

        // ---- local diagonal: precompute alpha-derived factors ----
        real_t f1_h=0, f1_l=0, f2_h=0, f2_l=0, fa_h=0, fa_l=0;
        if (is == 0) {
            real_t opa_h, opa_l, mpa_h, mpa_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0),  real_t(0.0), alpha_h, alpha_l, opa_h, opa_l);
            dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, mpa_h, mpa_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), mpa_h, mpa_l, hma_h, hma_l);
            dw_mul(B1_h, B1_l, hpa_h, hpa_l, f1_h, f1_l);
            dw_mul(B1_h, B1_l, hma_h, hma_l, f2_h, f2_l);
        } else if (is == Ns-1) {
            real_t opa_h, opa_l, oma_h, oma_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0), real_t(0.0),  alpha_h,  alpha_l, opa_h, opa_l);
            dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, oma_h, oma_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), oma_h, oma_l, hma_h, hma_l);
            dw_mul(B1_h, B1_l, hpa_h, hpa_l, f1_h, f1_l);
            dw_mul(B1_h, B1_l, hma_h, hma_l, f2_h, f2_l);
        } else {
            dw_mul(B1_h, B1_l, alpha_h, alpha_l, fa_h, fa_l);
        }

#define DIAG_DAG_FF(ic) \
        real4 w_##ic##_s0=wp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 w_##ic##_s1=wp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 w_##ic##_s2=wp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 w_##ic##_s3=wp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0 || is == Ns-1) { \
            DWF_SCAL_ADD2_DD(vt_##ic##_s0,f1_h,f1_l,w_##ic##_s0,f2_h,f2_l,w_##ic##_s2); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s1,f1_h,f1_l,w_##ic##_s1,f2_h,f2_l,w_##ic##_s3); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s2,f1_h,f1_l,w_##ic##_s2,f2_h,f2_l,w_##ic##_s0); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s3,f1_h,f1_l,w_##ic##_s3,f2_h,f2_l,w_##ic##_s1); \
        } else { \
            QDW_SCAL_DD(vt_##ic##_s0,fa_h,fa_l,w_##ic##_s0); QDW_SCAL_DD(vt_##ic##_s1,fa_h,fa_l,w_##ic##_s1); \
            QDW_SCAL_DD(vt_##ic##_s2,fa_h,fa_l,w_##ic##_s2); QDW_SCAL_DD(vt_##ic##_s3,fa_h,fa_l,w_##ic##_s3); \
        }

        DIAG_DAG_FF(0) DIAG_DAG_FF(1) DIAG_DAG_FF(2)
#undef DIAG_DAG_FF

        // ---- s+1 contribution (adjoint: symmetric) ----
        const int is_up = (is + 1) % Ns;
        // C2up = c[is_up]*(4-M0) - 1
        real_t cmup_h, cmup_l;
        dw_mul(c_h[is_up], c_l[is_up], fmM0_h, fmM0_l, cmup_h, cmup_l);
        real_t C2up_h, C2up_l;
        dw_add(cmup_h, cmup_l, real_t(-1.0), real_t(0.0), C2up_h, C2up_l);
        // Fup = (is==Ns-1) ? -0.5*mq*C2up*0.5 : 0.5*C2up*alpha*0.5
        real_t Fup_h, Fup_l;
        if (is == Ns-1) {
            real_t t_h, t_l;
            dw_mul(mhmq_h, mhmq_l, C2up_h, C2up_l, t_h, t_l);
            dw_scal(real_t(0.5), t_h, t_l, Fup_h, Fup_l);
        } else {
            real_t qC_h, qC_l;
            dw_scal(real_t(0.25), C2up_h, C2up_l, qC_h, qC_l);
            dw_mul(qC_h, qC_l, alpha_h, alpha_l, Fup_h, Fup_l);
        }
        real4 sum02, sum13;

#define ADD_SYM_UP_FF(ic) \
        { \
        real4 xu_s0=wp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 xu_s1=wp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 xu_s2=wp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 xu_s3=wp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        QDW_ADD(sum02,xu_s0,xu_s2); QDW_SCAL_DD(tmp,Fup_h,Fup_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,xu_s1,xu_s3); QDW_SCAL_DD(tmp,Fup_h,Fup_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SYM_UP_FF(0) ADD_SYM_UP_FF(1) ADD_SYM_UP_FF(2)
#undef ADD_SYM_UP_FF

        // ---- s-1 contribution (adjoint: antisymmetric) ----
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t cmdn_h, cmdn_l;
        dw_mul(c_h[is_dn], c_l[is_dn], fmM0_h, fmM0_l, cmdn_h, cmdn_l);
        real_t C2dn_h, C2dn_l;
        dw_add(cmdn_h, cmdn_l, real_t(-1.0), real_t(0.0), C2dn_h, C2dn_l);
        real_t Fdn_h, Fdn_l;
        if (is == 0) {
            real_t t_h, t_l;
            dw_mul(mhmq_h, mhmq_l, C2dn_h, C2dn_l, t_h, t_l);
            dw_scal(real_t(0.5), t_h, t_l, Fdn_h, Fdn_l);
        } else {
            real_t qC_h, qC_l;
            dw_scal(real_t(0.25), C2dn_h, C2dn_l, qC_h, qC_l);
            dw_mul(qC_h, qC_l, alpha_h, alpha_l, Fdn_h, Fdn_l);
        }

#define ADD_ANTISYM_DN_FF(ic) \
        { \
        real4 xd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 xd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 xd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 xd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_SUB(sum02,xd_s0,xd_s2); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,xd_s1,xd_s3); QDW_SCAL_DD(tmp,Fdn_h,Fdn_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_ANTISYM_DN_FF(0) ADD_ANTISYM_DN_FF(1) ADD_ANTISYM_DN_FF(2)
#undef ADD_ANTISYM_DN_FF

        vp[IDX_DWF_QDW(0,0,is,Ns,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Ns,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Ns,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Ns,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Ns,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Ns,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Ns,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Ns,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Ns,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Ns,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Ns,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Ns,site)]=vt_2_s3;
    }
}

//====================================================================
// eo 5-dirdag (adjoint off-diagonal block)  vp = (M_eo)† * yp
// (single-T variant.)
//====================================================================
__global__
void mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel_fp(
    real4 * __restrict__ vp, real4 * __restrict__ yp,
    real_t mq, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();

    for (int is = 0; is < Ns; ++is) {
        real4 tmp;

        // ---- local adjoint (adjoint of -0.5*b[is]*B_alpha) ----
#define LOCAL_DAG_EO(ic) \
        real4 y_##ic##_s0=yp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 y_##ic##_s1=yp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 y_##ic##_s2=yp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 y_##ic##_s3=yp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0) { \
            real_t b1=real_t(-0.5)*b_con[is]*real_t(0.5)*(real_t(1.0)+alpha); \
            real_t b2=real_t(-0.5)*b_con[is]*real_t(0.5)*(real_t(-1.0)+alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,b1,y_##ic##_s2,b2,y_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s1,b1,y_##ic##_s3,b2,y_##ic##_s1); \
            DWF_SCAL_ADD2(vt_##ic##_s2,b1,y_##ic##_s0,b2,y_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s3,b1,y_##ic##_s1,b2,y_##ic##_s3); \
        } else if (is == Ns-1) { \
            real_t b1=real_t(-0.5)*b_con[is]*real_t(0.5)*(real_t(1.0)+alpha); \
            real_t b2=real_t(-0.5)*b_con[is]*real_t(0.5)*(real_t(1.0)-alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,b1,y_##ic##_s2,b2,y_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s1,b1,y_##ic##_s3,b2,y_##ic##_s1); \
            DWF_SCAL_ADD2(vt_##ic##_s2,b1,y_##ic##_s0,b2,y_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s3,b1,y_##ic##_s1,b2,y_##ic##_s3); \
        } else { \
            real_t bb=real_t(-0.5)*b_con[is]*alpha; \
            QDW_SCAL(vt_##ic##_s0,bb,y_##ic##_s2); \
            QDW_SCAL(vt_##ic##_s1,bb,y_##ic##_s3); \
            QDW_SCAL(vt_##ic##_s2,bb,y_##ic##_s0); \
            QDW_SCAL(vt_##ic##_s3,bb,y_##ic##_s1); \
        }

        LOCAL_DAG_EO(0) LOCAL_DAG_EO(1) LOCAL_DAG_EO(2)
#undef LOCAL_DAG_EO

        // ---- s+1 contribution: adjoint of c-hop, symmetric ----
        const int is_up = (is + 1) % Ns;
        const real_t Fup_d = (is == Ns-1)
            ? real_t(-0.5) * c_con[is_up] * real_t(-0.5) * mq
            :  real_t(-0.5) * c_con[is_up] * real_t(0.5) * alpha;
        real4 sum02, sum13;

#define ADD_DAG_UP(ic) \
        { \
        real4 yu_s0=yp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 yu_s1=yp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 yu_s2=yp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 yu_s3=yp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        QDW_ADD(sum02,yu_s0,yu_s2); QDW_SCAL(tmp,Fup_d,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,yu_s1,yu_s3); QDW_SCAL(tmp,Fup_d,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_UP(0) ADD_DAG_UP(1) ADD_DAG_UP(2)
#undef ADD_DAG_UP

        // ---- s-1 contribution: adjoint of c-hop, antisymmetric ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const real_t Fdn_d = (is == 0)
            ? real_t(-0.5) * c_con[is_dn] * real_t(-0.5) * mq
            :  real_t(-0.5) * c_con[is_dn] * real_t(0.5) * alpha;

#define ADD_DAG_DN(ic) \
        { \
        real4 yd_s0=yp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 yd_s1=yp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 yd_s2=yp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 yd_s3=yp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_SUB(sum02,yd_s2,yd_s0); QDW_SCAL(tmp,Fdn_d,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,yd_s3,yd_s1); QDW_SCAL(tmp,Fdn_d,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_DN(0) ADD_DAG_DN(1) ADD_DAG_DN(2)
#undef ADD_DAG_DN

        vp[IDX_DWF_QDW(0,0,is,Ns,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Ns,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Ns,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Ns,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Ns,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Ns,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Ns,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Ns,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Ns,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Ns,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Ns,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Ns,site)]=vt_2_s3;
    }
}

// Forward declaration (definition follows the host wrapper below).
__global__
void mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ yp,
    real_t mq_h, real_t mq_l, int Ns,
    real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad);

void mult_domainwall_5din_eo_5dirdag_dirac_qdw(
    real_t *vp, real_t *yp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real4 *vp_dev = (real4 *)dev_ptr(vp);
    real4 *yp_dev = (real4 *)dev_ptr(yp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    if (ext && sizeof(real_t) == sizeof(float)) {
      const real_t mq_h = (real_t)mq;
      const real_t mq_l = (real_t)(mq - (double)mq_h);
      const real_t a_h  = (real_t)alpha;
      const real_t a_l  = (real_t)(alpha - (double)a_h);
      mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
          vp_dev, yp_dev, mq_h, mq_l, Ns, a_h, a_l, Nst, Nst_pad);
    } else {
      mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
          vp_dev, yp_dev, (real_t)mq, Ns, (real_t)alpha, Nst, Nst_pad);
    }
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-dirdag (FP32-only DD variant)  vp = (M_eo)† * yp
//====================================================================
__global__
void mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel_ff(
    real4 * __restrict__ vp, real4 * __restrict__ yp,
    real_t mq_h, real_t mq_l, int Ns,
    real_t alpha_h, real_t alpha_l,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_h = ConstantMemoryTraits<real_t>::b();
    const real_t *b_l = ConstantMemoryTraits<real_t>::b_lo();
    const real_t *c_h = ConstantMemoryTraits<real_t>::c();
    const real_t *c_l = ConstantMemoryTraits<real_t>::c_lo();

    // -0.5 * mq: DD
    real_t mhmq_h, mhmq_l;
    dw_scal(real_t(-0.5), mq_h, mq_l, mhmq_h, mhmq_l);

    for (int is = 0; is < Ns; ++is) {
        real4 tmp;

        // ---- local diagonal: precompute alpha-derived factors ----
        real_t b1_h=0, b1_l=0, b2_h=0, b2_l=0, bb_h=0, bb_l=0;
        // mb = -0.5 * b[is]
        real_t mb_h, mb_l;
        dw_scal(real_t(-0.5), b_h[is], b_l[is], mb_h, mb_l);
        if (is == 0) {
            real_t opa_h, opa_l, mpa_h, mpa_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0),  real_t(0.0), alpha_h, alpha_l, opa_h, opa_l);
            dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, mpa_h, mpa_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), mpa_h, mpa_l, hma_h, hma_l);
            dw_mul(mb_h, mb_l, hpa_h, hpa_l, b1_h, b1_l);
            dw_mul(mb_h, mb_l, hma_h, hma_l, b2_h, b2_l);
        } else if (is == Ns-1) {
            real_t opa_h, opa_l, oma_h, oma_l, hpa_h, hpa_l, hma_h, hma_l;
            dw_add(real_t(1.0), real_t(0.0),  alpha_h,  alpha_l, opa_h, opa_l);
            dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, oma_h, oma_l);
            dw_scal(real_t(0.5), opa_h, opa_l, hpa_h, hpa_l);
            dw_scal(real_t(0.5), oma_h, oma_l, hma_h, hma_l);
            dw_mul(mb_h, mb_l, hpa_h, hpa_l, b1_h, b1_l);
            dw_mul(mb_h, mb_l, hma_h, hma_l, b2_h, b2_l);
        } else {
            dw_mul(mb_h, mb_l, alpha_h, alpha_l, bb_h, bb_l);
        }

#define LOCAL_DAG_EO_FF(ic) \
        real4 y_##ic##_s0=yp[IDX_DWF_QDW(ic,0,is,Ns,site)]; \
        real4 y_##ic##_s1=yp[IDX_DWF_QDW(ic,1,is,Ns,site)]; \
        real4 y_##ic##_s2=yp[IDX_DWF_QDW(ic,2,is,Ns,site)]; \
        real4 y_##ic##_s3=yp[IDX_DWF_QDW(ic,3,is,Ns,site)]; \
        real4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0 || is == Ns-1) { \
            DWF_SCAL_ADD2_DD(vt_##ic##_s0,b1_h,b1_l,y_##ic##_s2,b2_h,b2_l,y_##ic##_s0); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s1,b1_h,b1_l,y_##ic##_s3,b2_h,b2_l,y_##ic##_s1); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s2,b1_h,b1_l,y_##ic##_s0,b2_h,b2_l,y_##ic##_s2); \
            DWF_SCAL_ADD2_DD(vt_##ic##_s3,b1_h,b1_l,y_##ic##_s1,b2_h,b2_l,y_##ic##_s3); \
        } else { \
            QDW_SCAL_DD(vt_##ic##_s0,bb_h,bb_l,y_##ic##_s2); \
            QDW_SCAL_DD(vt_##ic##_s1,bb_h,bb_l,y_##ic##_s3); \
            QDW_SCAL_DD(vt_##ic##_s2,bb_h,bb_l,y_##ic##_s0); \
            QDW_SCAL_DD(vt_##ic##_s3,bb_h,bb_l,y_##ic##_s1); \
        }

        LOCAL_DAG_EO_FF(0) LOCAL_DAG_EO_FF(1) LOCAL_DAG_EO_FF(2)
#undef LOCAL_DAG_EO_FF

        // ---- s+1 contribution: adjoint of c-hop, symmetric ----
        const int is_up = (is + 1) % Ns;
        // mhc_up = -0.5 * c[is_up]
        real_t mhc_up_h, mhc_up_l;
        dw_scal(real_t(-0.5), c_h[is_up], c_l[is_up], mhc_up_h, mhc_up_l);
        real_t Fup_d_h, Fup_d_l;
        if (is == Ns-1) {
            // Fup_d = mhc_up * -0.5 * mq = mhc_up * mhmq
            dw_mul(mhc_up_h, mhc_up_l, mhmq_h, mhmq_l, Fup_d_h, Fup_d_l);
        } else {
            // Fup_d = mhc_up * 0.5 * alpha
            real_t ha_h, ha_l;
            dw_scal(real_t(0.5), alpha_h, alpha_l, ha_h, ha_l);
            dw_mul(mhc_up_h, mhc_up_l, ha_h, ha_l, Fup_d_h, Fup_d_l);
        }
        real4 sum02, sum13;

#define ADD_DAG_UP_FF(ic) \
        { \
        real4 yu_s0=yp[IDX_DWF_QDW(ic,0,is_up,Ns,site)]; \
        real4 yu_s1=yp[IDX_DWF_QDW(ic,1,is_up,Ns,site)]; \
        real4 yu_s2=yp[IDX_DWF_QDW(ic,2,is_up,Ns,site)]; \
        real4 yu_s3=yp[IDX_DWF_QDW(ic,3,is_up,Ns,site)]; \
        QDW_ADD(sum02,yu_s0,yu_s2); QDW_SCAL_DD(tmp,Fup_d_h,Fup_d_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,yu_s1,yu_s3); QDW_SCAL_DD(tmp,Fup_d_h,Fup_d_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_UP_FF(0) ADD_DAG_UP_FF(1) ADD_DAG_UP_FF(2)
#undef ADD_DAG_UP_FF

        // ---- s-1 contribution: adjoint of c-hop, antisymmetric ----
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t mhc_dn_h, mhc_dn_l;
        dw_scal(real_t(-0.5), c_h[is_dn], c_l[is_dn], mhc_dn_h, mhc_dn_l);
        real_t Fdn_d_h, Fdn_d_l;
        if (is == 0) {
            dw_mul(mhc_dn_h, mhc_dn_l, mhmq_h, mhmq_l, Fdn_d_h, Fdn_d_l);
        } else {
            real_t ha_h, ha_l;
            dw_scal(real_t(0.5), alpha_h, alpha_l, ha_h, ha_l);
            dw_mul(mhc_dn_h, mhc_dn_l, ha_h, ha_l, Fdn_d_h, Fdn_d_l);
        }

#define ADD_DAG_DN_FF(ic) \
        { \
        real4 yd_s0=yp[IDX_DWF_QDW(ic,0,is_dn,Ns,site)]; \
        real4 yd_s1=yp[IDX_DWF_QDW(ic,1,is_dn,Ns,site)]; \
        real4 yd_s2=yp[IDX_DWF_QDW(ic,2,is_dn,Ns,site)]; \
        real4 yd_s3=yp[IDX_DWF_QDW(ic,3,is_dn,Ns,site)]; \
        QDW_SUB(sum02,yd_s2,yd_s0); QDW_SCAL_DD(tmp,Fdn_d_h,Fdn_d_l,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,yd_s3,yd_s1); QDW_SCAL_DD(tmp,Fdn_d_h,Fdn_d_l,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_DN_FF(0) ADD_DAG_DN_FF(1) ADD_DAG_DN_FF(2)
#undef ADD_DAG_DN_FF

        vp[IDX_DWF_QDW(0,0,is,Ns,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Ns,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Ns,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Ns,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Ns,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Ns,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Ns,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Ns,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Ns,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Ns,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Ns,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Ns,site)]=vt_2_s3;
    }
}

//====================================================================
// EO 4D bulk hopping (QDW, EO gauge layout, jgm5 support)
// Gauge stride: Nst*(ieo + 2*idir) for forward, Nst*(1-ieo + 2*idir) for backward
//====================================================================
// EXT=false: single-precision gauge link (current behavior).
// EXT=true : float-float / double-double (extended) gauge link from up + up_lo.
// RECON: SU(3) 3rd-column reconstruction (only used on the EXT=true / FF path).
template<bool EXT, bool RECON>
__global__
void mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel(
    real4 * __restrict__ vp, const real_t * __restrict__ up,
    const real_t * __restrict__ up_lo, const real4 * __restrict__ wp,
    int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int Nx, int Ny, int Nz, int Nt,
    int ieo, int jeo,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
    int Nst, int Nst_pad, int jgm5)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const int Nxy  = Nx * Ny;
    const int Nxyz = Nxy * Nz;

    const int ix   = site % Nx;
    const int iyzt = site / Nx;
    const int ixy  = site % Nxy;
    const int iy   = iyzt % Ny;
    const int izt  = site / Nxy;
    const int iz   = izt  % Nz;
    const int it   = izt  / Nz;
    const int ixyz = site % Nxyz;
    const int keo  = (jeo + iy + iz + it) % 2;

    // jgm5: swap id 0↔2, 1↔3 before projection
    const int d0 = (jgm5 == 0) ? 0 : 2;
    const int d1 = (jgm5 == 0) ? 1 : 3;
    const int d2 = (jgm5 == 0) ? 2 : 0;
    const int d3 = (jgm5 == 0) ? 3 : 1;

    const real_t *u_up = up;
    const real_t *u_dn = up;
    const real_t *u_up_lo = up_lo;
    const real_t *u_dn_lo = up_lo;

    for (int is = 0; is < Ns; ++is) {

        real4 v2_c0_s0={0,0,0,0}, v2_c0_s1={0,0,0,0},
                v2_c0_s2={0,0,0,0}, v2_c0_s3={0,0,0,0};
        real4 v2_c1_s0={0,0,0,0}, v2_c1_s1={0,0,0,0},
                v2_c1_s2={0,0,0,0}, v2_c1_s3={0,0,0,0};
        real4 v2_c2_s0={0,0,0,0}, v2_c2_s1={0,0,0,0},
                v2_c2_s2={0,0,0,0}, v2_c2_s3={0,0,0,0};

        real4 vt1_c0, vt1_c1, vt1_c2;
        real4 vt2_c0, vt2_c1, vt2_c2;
        real4 wt1_c0, wt1_c1, wt1_c2;
        real4 wt2_c0, wt2_c1, wt2_c2;
        real_t bc2;

        // X+
        {
            int isn = ((ix + keo) % Nx) + Nx * iyzt;
            int isg = site + Nst * (ieo + 2*0);
            bc2 = (ix == Nx-1 && keo == 1) ? (real_t)bc_x : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d3, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d2, DWF_PROJ_P)
            if constexpr (EXT) { DWF_GMUL_FWD_FF(u_up, u_up_lo, isg) } else { DWF_GMUL_FWD(u_up, isg) }
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_MI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_MI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_MI,DWF_MULT_MI)
        }

        // X-
        {
            int ix2 = (ix - 1 + keo + Nx) % Nx;
            int isn = ix2 + Nx * iyzt;
            int isg = isn + Nst * (1-ieo + 2*0);
            bc2 = (ix == 0 && keo == 0) ? (real_t)bc_x : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d3, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d2, DWF_PROJ_M)
            if constexpr (EXT) { DWF_GMUL_BCK_FF(u_dn, u_dn_lo, isg) } else { DWF_GMUL_BCK(u_dn, isg) }
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_PI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_PI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_PI,DWF_MULT_PI)
        }

        // Y+
        if ((iy < Ny-1) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy+1)%Ny) + Ny*izt);
            int isg = site + Nst * (ieo + 2*1);
            bc2 = (iy == Ny-1) ? (real_t)bc_y : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d3, DWF_PROJ_RP)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d2, DWF_PROJ_RM)
            if constexpr (EXT) { DWF_GMUL_FWD_FF(u_up, u_up_lo, isg) } else { DWF_GMUL_FWD(u_up, isg) }
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_RM1,DWF_MULT_R1)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_RM1,DWF_MULT_R1)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_RM1,DWF_MULT_R1)
        }

        // Y-
        if ((iy > 0) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy-1+Ny)%Ny) + Ny*izt);
            int isg = isn + Nst * (1-ieo + 2*1);
            bc2 = (iy == 0) ? (real_t)bc_y : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d3, DWF_PROJ_RM)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d2, DWF_PROJ_RP)
            if constexpr (EXT) { DWF_GMUL_BCK_FF(u_dn, u_dn_lo, isg) } else { DWF_GMUL_BCK(u_dn, isg) }
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_R1,DWF_MULT_RM1)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_R1,DWF_MULT_RM1)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_R1,DWF_MULT_RM1)
        }

        // Z+
        if ((iz < Nz-1) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz+1)%Nz) + Nz*it);
            int isg = site + Nst * (ieo + 2*2);
            bc2 = (iz == Nz-1) ? (real_t)bc_z : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d2, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d3, DWF_PROJ_M)
            if constexpr (EXT) { DWF_GMUL_FWD_FF(u_up, u_up_lo, isg) } else { DWF_GMUL_FWD(u_up, isg) }
            DWF_ACCUM_4_SW(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_MI,DWF_MULT_PI)
        }

        // Z-
        if ((iz > 0) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz-1+Nz)%Nz) + Nz*it);
            int isg = isn + Nst * (1-ieo + 2*2);
            bc2 = (iz == 0) ? (real_t)bc_z : 1.0;
            DWF_LOAD_PROJ (wp, is, Ns, isn, d0, d2, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Ns, isn, d1, d3, DWF_PROJ_P)
            if constexpr (EXT) { DWF_GMUL_BCK_FF(u_dn, u_dn_lo, isg) } else { DWF_GMUL_BCK(u_dn, isg) }
            DWF_ACCUM_4_SW(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_PI,DWF_MULT_MI)
        }

        // T+
        if ((it < Nt-1) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it+1)%Nt);
            int isg = site + Nst * (ieo + 2*3);
            bc2 = (it == Nt-1) ? (real_t)bc_t : 1.0;
            DWF_LOAD_PROJ_T(wp, is, Ns, isn, d2, d3)
            if constexpr (EXT) { DWF_GMUL_FWD_FF(u_up, u_up_lo, isg) } else { DWF_GMUL_FWD(u_up, isg) }
            DWF_ACCUM_TP(v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2)
            DWF_ACCUM_TP(v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2)
            DWF_ACCUM_TP(v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2)
        }

        // T-
        if ((it > 0) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it-1+Nt)%Nt);
            int isg = isn + Nst * (1-ieo + 2*3);
            bc2 = (it == 0) ? (real_t)bc_t : 1.0;
            DWF_LOAD_PROJ_T(wp, is, Ns, isn, d0, d1)
            if constexpr (EXT) { DWF_GMUL_BCK_FF(u_dn, u_dn_lo, isg) } else { DWF_GMUL_BCK(u_dn, isg) }
            DWF_ACCUM_TM(v2_c0_s0,v2_c0_s1, wt1_c0,wt2_c0, bc2)
            DWF_ACCUM_TM(v2_c1_s0,v2_c1_s1, wt1_c1,wt2_c1, bc2)
            DWF_ACCUM_TM(v2_c2_s0,v2_c2_s1, wt1_c2,wt2_c2, bc2)
        }

        vp[IDX_DWF_QDW(0,0,is,Ns,site)]=v2_c0_s0; vp[IDX_DWF_QDW(0,1,is,Ns,site)]=v2_c0_s1;
        vp[IDX_DWF_QDW(0,2,is,Ns,site)]=v2_c0_s2; vp[IDX_DWF_QDW(0,3,is,Ns,site)]=v2_c0_s3;
        vp[IDX_DWF_QDW(1,0,is,Ns,site)]=v2_c1_s0; vp[IDX_DWF_QDW(1,1,is,Ns,site)]=v2_c1_s1;
        vp[IDX_DWF_QDW(1,2,is,Ns,site)]=v2_c1_s2; vp[IDX_DWF_QDW(1,3,is,Ns,site)]=v2_c1_s3;
        vp[IDX_DWF_QDW(2,0,is,Ns,site)]=v2_c2_s0; vp[IDX_DWF_QDW(2,1,is,Ns,site)]=v2_c2_s1;
        vp[IDX_DWF_QDW(2,2,is,Ns,site)]=v2_c2_s2; vp[IDX_DWF_QDW(2,3,is,Ns,site)]=v2_c2_s3;
    }
}

// up_lo == nullptr -> single-precision gauge link (default).
// up_lo != nullptr -> extended (float-float / double-double) gauge link.
void mult_domainwall_5din_eo_hopb_qdw_dirac_5d(
    real_t *vp, real_t *up, real_t *up_lo, real_t *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo, int jgm5, bool recon)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst     = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real4 *vp_dev = (real4 *)dev_ptr(vp);
    real_t  *up_dev = (real_t  *)dev_ptr(up);
    real_t  *up_lo_dev = up_lo ? (real_t *)dev_ptr(up_lo) : nullptr;
    real4 *wp_dev = (real4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;

    if (up_lo_dev) {
      if (recon)
        mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel<true,true><<<gridSize, blockSize>>>(
          vp_dev, up_dev, up_lo_dev, wp_dev, Ns,
          bc[0], bc[1], bc[2], bc[3], Nx, Ny, Nz, Nt, ieo, jeo,
          do_comm[0], do_comm[1], do_comm[2], do_comm[3], Nst, Nst_pad, jgm5);
      else
        mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel<true,false><<<gridSize, blockSize>>>(
          vp_dev, up_dev, up_lo_dev, wp_dev, Ns,
          bc[0], bc[1], bc[2], bc[3], Nx, Ny, Nz, Nt, ieo, jeo,
          do_comm[0], do_comm[1], do_comm[2], do_comm[3], Nst, Nst_pad, jgm5);
    } else {
      // single-float gauge has no separate columns to reconstruct; RECON=false.
      mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel<false,false><<<gridSize, blockSize>>>(
        vp_dev, up_dev, nullptr, wp_dev, Ns,
        bc[0], bc[1], bc[2], bc[3], Nx, Ny, Nz, Nt, ieo, jeo,
        do_comm[0], do_comm[1], do_comm[2], do_comm[3], Nst, Nst_pad, jgm5);
    }

    CHECK(cudaDeviceSynchronize());
}

// ---- Undef all macros ----
#undef IDX_DWF_QDW
#undef DWF_PROJ_P
#undef DWF_PROJ_M
#undef DWF_MULT_MI
#undef DWF_MULT_PI
#undef DWF_PROJ_2
#undef DWF_GMUL_FWD
#undef DWF_GMUL_BCK
#undef DWF_GMUL_FWD_FF
#undef DWF_GMUL_BCK_FF
#undef DWF_ACCUM_4
#undef DWF_ACCUM_4_SW
#undef DWF_ACCUM_TP
#undef DWF_ACCUM_TM
#undef DWF_LOAD_PROJ
#undef DWF_LOAD_PROJ2
#undef DWF_LOAD_PROJ_T
#undef DWF_SCAL_ADD2
#undef DWF_SCAL_ADD2_DD

#endif // MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED
