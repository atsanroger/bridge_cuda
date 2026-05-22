/*!
      @file    mult_Domainwall_5din_cuda_qdw-inc.h
      @brief   QDW Domain Wall fermion kernels (5dir + hopb, dirac basis)
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC_QDW_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC_QDW_INCLUDED

// ---------------------------------------------------------------------------
// Local helper macros (mirror Wilson QDW definitions, no conflict since
// mult_Wilson_cuda_qdw-inc.h is not included in the DWF compilation unit)
// ---------------------------------------------------------------------------
#ifndef DWF_QDW_HELPERS_DEFINED
#define DWF_QDW_HELPERS_DEFINED

// res = a + i*b  (both QDW complex)
#define DWF_PROJ_P(res, a, b) \
    dw_add((a).x, (a).z, -(b).y, -(b).w, (res).x, (res).z); \
    dw_add((a).y, (a).w,  (b).x,  (b).z, (res).y, (res).w);

// res = a - i*b
#define DWF_PROJ_M(res, a, b) \
    dw_add((a).x, (a).z,  (b).y,  (b).w, (res).x, (res).z); \
    dw_add((a).y, (a).w, -(b).x, -(b).z, (res).y, (res).w);

// res = -i * a  (exact, just component reorder)
#define DWF_MULT_MI(res, a) \
    (res).x =  (a).y; (res).y = -(a).x; (res).z =  (a).w; (res).w = -(a).z;

// res = +i * a
#define DWF_MULT_PI(res, a) \
    (res).x = -(a).y; (res).y =  (a).x; (res).z = -(a).w; (res).w =  (a).z;

// res = 2 * a  (used for T-direction projection)
#define DWF_PROJ_2(res, a) \
    dw_scal(2.0, (a).x, (a).z, (res).x, (res).z); \
    dw_scal(2.0, (a).y, (a).w, (res).y, (res).w);

// QDW field index for a DWF spinor: ic=color, id=spin, is5=s-slice, site=4D site
// Works because Nst_pad is always a multiple of NWP.
#define IDX_DWF_QDW(ic, id, is5, Nst_pad_, site_) \
    IDX2_QDW((ic), (id), (is5) * (Nst_pad_) + (site_))

#endif // DWF_QDW_HELPERS_DEFINED

// ---------------------------------------------------------------------------
// Compact gauge-multiply macros
// Require: vt1_c{0,1,2}, vt2_c{0,1,2} already set; isg already computed.
// Produce: wt1_c0, wt2_c0, wt1_c1, wt2_c1, wt1_c2, wt2_c2.
// u_ptr is either u_up (forward) or u_dn (backward).
// For backward (U†): pass u_ptr = u_dn, conj_sign = -1.0
// ---------------------------------------------------------------------------

// Forward U multiply (3x3 matrix rows, no conjugate)
#define DWF_GMUL_FWD(u_ptr, isg_) \
{ \
    double2 _u; double4 _t; \
    _u.x = (u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0 = qdw_mult_uc(_u, vt1_c0); wt2_c0 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, _t); \
    _u.x = (u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, _t); \
    _u.x = (u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    wt1_c1 = qdw_mult_uc(_u, vt1_c0); wt2_c1 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, _t); \
    _u.x = (u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, _t); \
    _u.x = (u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    wt1_c2 = qdw_mult_uc(_u, vt1_c0); wt2_c2 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, _t); \
    _u.x = (u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y = (u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, _t); \
}

// Backward U† multiply (conjugate-transpose: U†_{r,c} = conj(U_{c,r}))
#define DWF_GMUL_BCK(u_ptr, isg_) \
{ \
    double2 _u; double4 _t; \
    _u.x = (u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0 = qdw_mult_uc(_u, vt1_c0); wt2_c0 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c0, wt1_c0, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c0, wt2_c0, _t); \
    _u.x = (u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c0, wt1_c0, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c0, wt2_c0, _t); \
    _u.x = (u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    wt1_c1 = qdw_mult_uc(_u, vt1_c0); wt2_c1 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c1, wt1_c1, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c1, wt2_c1, _t); \
    _u.x = (u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c1, wt1_c1, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c1, wt2_c1, _t); \
    _u.x = (u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    wt1_c2 = qdw_mult_uc(_u, vt1_c0); wt2_c2 = qdw_mult_uc(_u, vt2_c0); \
    _u.x = (u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c1); QDW_ADD(wt1_c2, wt1_c2, _t); \
    _t = qdw_mult_uc(_u, vt2_c1); QDW_ADD(wt2_c2, wt2_c2, _t); \
    _u.x = (u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y = -(u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t = qdw_mult_uc(_u, vt1_c2); QDW_ADD(wt1_c2, wt1_c2, _t); \
    _t = qdw_mult_uc(_u, vt2_c2); QDW_ADD(wt2_c2, wt2_c2, _t); \
}

// Accumulate for one color: v2_s{0,1,2,3} += bc2 * {wt1, wt2, op2(wt2), op3(wt1)}
// XP: op2=-i*wt2, op3=-i*wt1  (DWF_MULT_MI)
// XM: op2=+i*wt2, op3=+i*wt1  (DWF_MULT_PI)
// YP: op2=+i*wt2, op3=-i*wt1
// YM: op2=-i*wt2, op3=+i*wt1
// ZP: op2=-i*wt1, op3=+i*wt2
// ZM: op2=+i*wt1, op3=-i*wt2
// TP: only s2,s3 ← wt1,wt2
// TM: only s0,s1 ← wt1,wt2
#define DWF_ACCUM_4(v0,v1,v2,v3, W1,W2, bc_, OP2_MACRO_ARG, OP3_MACRO_ARG) \
{ \
    double4 _ts, _tp; \
    QDW_SCAL(_ts, bc_, W1); QDW_ADD(v0, v0, _ts); \
    QDW_SCAL(_ts, bc_, W2); QDW_ADD(v1, v1, _ts); \
    OP2_MACRO_ARG(_tp, W2); QDW_SCAL(_ts, bc_, _tp); QDW_ADD(v2, v2, _ts); \
    OP3_MACRO_ARG(_tp, W1); QDW_SCAL(_ts, bc_, _tp); QDW_ADD(v3, v3, _ts); \
}

// TP / TM accumulate to only 2 spin components
#define DWF_ACCUM_TP(v2,v3, W1,W2, bc_) \
{ \
    double4 _ts; \
    QDW_SCAL(_ts, bc_, W1); QDW_ADD(v2, v2, _ts); \
    QDW_SCAL(_ts, bc_, W2); QDW_ADD(v3, v3, _ts); \
}
#define DWF_ACCUM_TM(v0,v1, W1,W2, bc_) \
{ \
    double4 _ts; \
    QDW_SCAL(_ts, bc_, W1); QDW_ADD(v0, v0, _ts); \
    QDW_SCAL(_ts, bc_, W2); QDW_ADD(v1, v1, _ts); \
}

// Swapped version: OP2 applies to W1 (->s2), OP3 applies to W2 (->s3).
// Used for Z directions where the projection half-spinors are swapped
// relative to X/Y: ZP/ZM have s2 from t1, s3 from t2.
#define DWF_ACCUM_4_SW(v0,v1,v2,v3, W1,W2, bc_, OP2_ON_W1, OP3_ON_W2) \
{ \
    double4 _ts, _tp; \
    QDW_SCAL(_ts, bc_, W1); QDW_ADD(v0, v0, _ts); \
    QDW_SCAL(_ts, bc_, W2); QDW_ADD(v1, v1, _ts); \
    OP2_ON_W1(_tp, W1); QDW_SCAL(_ts, bc_, _tp); QDW_ADD(v2, v2, _ts); \
    OP3_ON_W2(_tp, W2); QDW_SCAL(_ts, bc_, _tp); QDW_ADD(v3, v3, _ts); \
}

// Load projection half-spinors from DWF QDW field
// Sets vt1_c{0,1,2}, vt2_c{0,1,2} using PROJ macro with ids id_a,id_b
#define DWF_LOAD_PROJ(wp_, is_, Np_, isn_, id_a, id_b, PROJ_MACRO) \
{ \
    double4 _sa, _sb; \
    _sa = (wp_)[IDX_DWF_QDW(0, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(0, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt1_c0, _sa, _sb); \
    _sa = (wp_)[IDX_DWF_QDW(1, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(1, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt1_c1, _sa, _sb); \
    _sa = (wp_)[IDX_DWF_QDW(2, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(2, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt1_c2, _sa, _sb); \
}
// Second 2-spinor projection
#define DWF_LOAD_PROJ2(wp_, is_, Np_, isn_, id_a, id_b, PROJ_MACRO) \
{ \
    double4 _sa, _sb; \
    _sa = (wp_)[IDX_DWF_QDW(0, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(0, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt2_c0, _sa, _sb); \
    _sa = (wp_)[IDX_DWF_QDW(1, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(1, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt2_c1, _sa, _sb); \
    _sa = (wp_)[IDX_DWF_QDW(2, id_a, is_, Np_, isn_)]; \
    _sb = (wp_)[IDX_DWF_QDW(2, id_b, is_, Np_, isn_)]; \
    PROJ_MACRO(vt2_c2, _sa, _sb); \
}

// T-direction projection: vt1 = 2*s_{id_a}, vt2 = 2*s_{id_b}
#define DWF_LOAD_PROJ_T(wp_, is_, Np_, isn_, id_a, id_b) \
{ \
    double4 _sa; \
    _sa = (wp_)[IDX_DWF_QDW(0, id_a, is_, Np_, isn_)]; DWF_PROJ_2(vt1_c0, _sa); \
    _sa = (wp_)[IDX_DWF_QDW(1, id_a, is_, Np_, isn_)]; DWF_PROJ_2(vt1_c1, _sa); \
    _sa = (wp_)[IDX_DWF_QDW(2, id_a, is_, Np_, isn_)]; DWF_PROJ_2(vt1_c2, _sa); \
    _sa = (wp_)[IDX_DWF_QDW(0, id_b, is_, Np_, isn_)]; DWF_PROJ_2(vt2_c0, _sa); \
    _sa = (wp_)[IDX_DWF_QDW(1, id_b, is_, Np_, isn_)]; DWF_PROJ_2(vt2_c1, _sa); \
    _sa = (wp_)[IDX_DWF_QDW(2, id_b, is_, Np_, isn_)]; DWF_PROJ_2(vt2_c2, _sa); \
}


//====================================================================
// 5D s-direction hopping kernel (QDW)
// Computes:
//   vp[is] = B_is * wp[is] + C_is * vt[is]
//   yp[is] = -0.5 * (b[is] * wp[is] + c[is] * vt[is])
// where vt[is] is the s-direction hopping contribution.
//====================================================================
__global__
void mult_domainwall_5din_5dir_dirac_qdw_kernel(
    double4 *vp, double4 *yp, double4 *wp,
    double mq, double M0, int Ns, int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const double *b_con = const_b_double;
    const double *c_con = const_c_double;

    for (int is = 0; is < Ns; ++is) {

        // ----- s+1 hopping (up) -----
        const int is_up = (is + 1) % Ns;
        const double Fup = (is == Ns - 1) ? -0.5 * mq : 0.5;

        double4 wu_c0_s0 = wp[IDX_DWF_QDW(0, 0, is_up, Nst_pad, site)];
        double4 wu_c0_s1 = wp[IDX_DWF_QDW(0, 1, is_up, Nst_pad, site)];
        double4 wu_c0_s2 = wp[IDX_DWF_QDW(0, 2, is_up, Nst_pad, site)];
        double4 wu_c0_s3 = wp[IDX_DWF_QDW(0, 3, is_up, Nst_pad, site)];
        double4 wu_c1_s0 = wp[IDX_DWF_QDW(1, 0, is_up, Nst_pad, site)];
        double4 wu_c1_s1 = wp[IDX_DWF_QDW(1, 1, is_up, Nst_pad, site)];
        double4 wu_c1_s2 = wp[IDX_DWF_QDW(1, 2, is_up, Nst_pad, site)];
        double4 wu_c1_s3 = wp[IDX_DWF_QDW(1, 3, is_up, Nst_pad, site)];
        double4 wu_c2_s0 = wp[IDX_DWF_QDW(2, 0, is_up, Nst_pad, site)];
        double4 wu_c2_s1 = wp[IDX_DWF_QDW(2, 1, is_up, Nst_pad, site)];
        double4 wu_c2_s2 = wp[IDX_DWF_QDW(2, 2, is_up, Nst_pad, site)];
        double4 wu_c2_s3 = wp[IDX_DWF_QDW(2, 3, is_up, Nst_pad, site)];

        // vt[s0] = Fup*(wu[s0] - wu[s2]),  vt[s2] = -vt[s0]
        // vt[s1] = Fup*(wu[s1] - wu[s3]),  vt[s3] = -vt[s1]
        double4 vt_c0_s0, vt_c0_s1, vt_c0_s2, vt_c0_s3;
        double4 vt_c1_s0, vt_c1_s1, vt_c1_s2, vt_c1_s3;
        double4 vt_c2_s0, vt_c2_s1, vt_c2_s2, vt_c2_s3;
        double4 tmp;

        QDW_SUB(tmp, wu_c0_s0, wu_c0_s2); QDW_SCAL(vt_c0_s0, Fup, tmp);
        QDW_SUB(tmp, wu_c0_s1, wu_c0_s3); QDW_SCAL(vt_c0_s1, Fup, tmp);
        QDW_NEG(vt_c0_s2, vt_c0_s0);
        QDW_NEG(vt_c0_s3, vt_c0_s1);

        QDW_SUB(tmp, wu_c1_s0, wu_c1_s2); QDW_SCAL(vt_c1_s0, Fup, tmp);
        QDW_SUB(tmp, wu_c1_s1, wu_c1_s3); QDW_SCAL(vt_c1_s1, Fup, tmp);
        QDW_NEG(vt_c1_s2, vt_c1_s0);
        QDW_NEG(vt_c1_s3, vt_c1_s1);

        QDW_SUB(tmp, wu_c2_s0, wu_c2_s2); QDW_SCAL(vt_c2_s0, Fup, tmp);
        QDW_SUB(tmp, wu_c2_s1, wu_c2_s3); QDW_SCAL(vt_c2_s1, Fup, tmp);
        QDW_NEG(vt_c2_s2, vt_c2_s0);
        QDW_NEG(vt_c2_s3, vt_c2_s1);

        // ----- s-1 hopping (down) -----
        const int is_dn = (is - 1 + Ns) % Ns;
        const double Fdn = (is == 0) ? -0.5 * mq : 0.5;

        double4 wd_c0_s0 = wp[IDX_DWF_QDW(0, 0, is_dn, Nst_pad, site)];
        double4 wd_c0_s1 = wp[IDX_DWF_QDW(0, 1, is_dn, Nst_pad, site)];
        double4 wd_c0_s2 = wp[IDX_DWF_QDW(0, 2, is_dn, Nst_pad, site)];
        double4 wd_c0_s3 = wp[IDX_DWF_QDW(0, 3, is_dn, Nst_pad, site)];
        double4 wd_c1_s0 = wp[IDX_DWF_QDW(1, 0, is_dn, Nst_pad, site)];
        double4 wd_c1_s1 = wp[IDX_DWF_QDW(1, 1, is_dn, Nst_pad, site)];
        double4 wd_c1_s2 = wp[IDX_DWF_QDW(1, 2, is_dn, Nst_pad, site)];
        double4 wd_c1_s3 = wp[IDX_DWF_QDW(1, 3, is_dn, Nst_pad, site)];
        double4 wd_c2_s0 = wp[IDX_DWF_QDW(2, 0, is_dn, Nst_pad, site)];
        double4 wd_c2_s1 = wp[IDX_DWF_QDW(2, 1, is_dn, Nst_pad, site)];
        double4 wd_c2_s2 = wp[IDX_DWF_QDW(2, 2, is_dn, Nst_pad, site)];
        double4 wd_c2_s3 = wp[IDX_DWF_QDW(2, 3, is_dn, Nst_pad, site)];

        // vt[s0] += Fdn*(wd[s0]+wd[s2]),  vt[s2] += same  (s0 and s2 symmetric)
        // vt[s1] += Fdn*(wd[s1]+wd[s3]),  vt[s3] += same
        double4 sum02, sum13;

        QDW_ADD(sum02, wd_c0_s0, wd_c0_s2); QDW_SCAL(tmp, Fdn, sum02);
        QDW_ADD(vt_c0_s0, vt_c0_s0, tmp); QDW_ADD(vt_c0_s2, vt_c0_s2, tmp);
        QDW_ADD(sum13, wd_c0_s1, wd_c0_s3); QDW_SCAL(tmp, Fdn, sum13);
        QDW_ADD(vt_c0_s1, vt_c0_s1, tmp); QDW_ADD(vt_c0_s3, vt_c0_s3, tmp);

        QDW_ADD(sum02, wd_c1_s0, wd_c1_s2); QDW_SCAL(tmp, Fdn, sum02);
        QDW_ADD(vt_c1_s0, vt_c1_s0, tmp); QDW_ADD(vt_c1_s2, vt_c1_s2, tmp);
        QDW_ADD(sum13, wd_c1_s1, wd_c1_s3); QDW_SCAL(tmp, Fdn, sum13);
        QDW_ADD(vt_c1_s1, vt_c1_s1, tmp); QDW_ADD(vt_c1_s3, vt_c1_s3, tmp);

        QDW_ADD(sum02, wd_c2_s0, wd_c2_s2); QDW_SCAL(tmp, Fdn, sum02);
        QDW_ADD(vt_c2_s0, vt_c2_s0, tmp); QDW_ADD(vt_c2_s2, vt_c2_s2, tmp);
        QDW_ADD(sum13, wd_c2_s1, wd_c2_s3); QDW_SCAL(tmp, Fdn, sum13);
        QDW_ADD(vt_c2_s1, vt_c2_s1, tmp); QDW_ADD(vt_c2_s3, vt_c2_s3, tmp);

        // ----- Current slice w -----
        double4 w_c0_s0 = wp[IDX_DWF_QDW(0, 0, is, Nst_pad, site)];
        double4 w_c0_s1 = wp[IDX_DWF_QDW(0, 1, is, Nst_pad, site)];
        double4 w_c0_s2 = wp[IDX_DWF_QDW(0, 2, is, Nst_pad, site)];
        double4 w_c0_s3 = wp[IDX_DWF_QDW(0, 3, is, Nst_pad, site)];
        double4 w_c1_s0 = wp[IDX_DWF_QDW(1, 0, is, Nst_pad, site)];
        double4 w_c1_s1 = wp[IDX_DWF_QDW(1, 1, is, Nst_pad, site)];
        double4 w_c1_s2 = wp[IDX_DWF_QDW(1, 2, is, Nst_pad, site)];
        double4 w_c1_s3 = wp[IDX_DWF_QDW(1, 3, is, Nst_pad, site)];
        double4 w_c2_s0 = wp[IDX_DWF_QDW(2, 0, is, Nst_pad, site)];
        double4 w_c2_s1 = wp[IDX_DWF_QDW(2, 1, is, Nst_pad, site)];
        double4 w_c2_s2 = wp[IDX_DWF_QDW(2, 2, is, Nst_pad, site)];
        double4 w_c2_s3 = wp[IDX_DWF_QDW(2, 3, is, Nst_pad, site)];

        const double B_is = b_con[is] * (4.0 - M0) + 1.0;
        const double C_is = c_con[is] * (4.0 - M0) - 1.0;
        const double a_b  = -0.5 * b_con[is];
        const double a_c  = -0.5 * c_con[is];

        // vp[is] = B*w + C*vt
        double4 ov, oy;
#define STORE_COMP_VP_YP(ic, id, wv, vtv) \
        QDW_SCAL(ov, B_is, wv); QDW_SCAL(tmp, C_is, vtv); QDW_ADD(ov, ov, tmp); \
        QDW_SCAL(oy, a_b,  wv); QDW_SCAL(tmp, a_c,  vtv); QDW_ADD(oy, oy, tmp); \
        vp[IDX_DWF_QDW(ic, id, is, Nst_pad, site)] = ov; \
        yp[IDX_DWF_QDW(ic, id, is, Nst_pad, site)] = oy;

        STORE_COMP_VP_YP(0, 0, w_c0_s0, vt_c0_s0)
        STORE_COMP_VP_YP(0, 1, w_c0_s1, vt_c0_s1)
        STORE_COMP_VP_YP(0, 2, w_c0_s2, vt_c0_s2)
        STORE_COMP_VP_YP(0, 3, w_c0_s3, vt_c0_s3)
        STORE_COMP_VP_YP(1, 0, w_c1_s0, vt_c1_s0)
        STORE_COMP_VP_YP(1, 1, w_c1_s1, vt_c1_s1)
        STORE_COMP_VP_YP(1, 2, w_c1_s2, vt_c1_s2)
        STORE_COMP_VP_YP(1, 3, w_c1_s3, vt_c1_s3)
        STORE_COMP_VP_YP(2, 0, w_c2_s0, vt_c2_s0)
        STORE_COMP_VP_YP(2, 1, w_c2_s1, vt_c2_s1)
        STORE_COMP_VP_YP(2, 2, w_c2_s2, vt_c2_s2)
        STORE_COMP_VP_YP(2, 3, w_c2_s3, vt_c2_s3)
#undef STORE_COMP_VP_YP
    } // is loop
}


void mult_domainwall_5din_5dir_dirac_qdw(
    double *vp, double *yp, double *wp,
    double mq, double M0, int Ns, double *b, double *c, int *Nsize)
{
    int Nx  = Nsize[0];
    int Ny  = Nsize[1];
    int Nz  = Nsize[2];
    int Nt  = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    double4 *yp_dev = (double4 *)dev_ptr(yp);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;

    mult_domainwall_5din_5dir_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, yp_dev, wp_dev, mq, M0, Ns, Nst, Nst_pad);

    CHECK(cudaDeviceSynchronize());
}


//====================================================================
// 4D bulk hopping kernel (QDW, dirac basis)
// Adds the 8-direction Wilson hopping to vp[is] for each s-slice.
// flag=0: zero-initialize vp before adding; flag=1: accumulate.
//====================================================================
__global__
void mult_domainwall_5din_hopb_qdw_dirac_kernel(
    double4 *vp, real_t *up, double4 *wp, int Ns,
    int bc_x, int bc_y, int bc_z, int bc_t,
    int Nx, int Ny, int Nz, int Nt,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
    int flag, int Nst, int Nst_pad)
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

    real_t *u_up = up;
    real_t *u_dn = up;

    for (int is = 0; is < Ns; ++is) {

        // 12 accumulators (3 colors × 4 spins)
        double4 v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3;
        double4 v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3;
        double4 v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3;

        if (flag == 0) {
            double4 _z = {0.0, 0.0, 0.0, 0.0};
            v2_c0_s0=_z; v2_c0_s1=_z; v2_c0_s2=_z; v2_c0_s3=_z;
            v2_c1_s0=_z; v2_c1_s1=_z; v2_c1_s2=_z; v2_c1_s3=_z;
            v2_c2_s0=_z; v2_c2_s1=_z; v2_c2_s2=_z; v2_c2_s3=_z;
        } else {
            v2_c0_s0 = vp[IDX_DWF_QDW(0,0,is,Nst_pad,site)];
            v2_c0_s1 = vp[IDX_DWF_QDW(0,1,is,Nst_pad,site)];
            v2_c0_s2 = vp[IDX_DWF_QDW(0,2,is,Nst_pad,site)];
            v2_c0_s3 = vp[IDX_DWF_QDW(0,3,is,Nst_pad,site)];
            v2_c1_s0 = vp[IDX_DWF_QDW(1,0,is,Nst_pad,site)];
            v2_c1_s1 = vp[IDX_DWF_QDW(1,1,is,Nst_pad,site)];
            v2_c1_s2 = vp[IDX_DWF_QDW(1,2,is,Nst_pad,site)];
            v2_c1_s3 = vp[IDX_DWF_QDW(1,3,is,Nst_pad,site)];
            v2_c2_s0 = vp[IDX_DWF_QDW(2,0,is,Nst_pad,site)];
            v2_c2_s1 = vp[IDX_DWF_QDW(2,1,is,Nst_pad,site)];
            v2_c2_s2 = vp[IDX_DWF_QDW(2,2,is,Nst_pad,site)];
            v2_c2_s3 = vp[IDX_DWF_QDW(2,3,is,Nst_pad,site)];
        }

        // Projection/gauge-mul temporaries (reused across directions)
        double4 vt1_c0, vt1_c1, vt1_c2;
        double4 vt2_c0, vt2_c1, vt2_c2;
        double4 wt1_c0, wt1_c1, wt1_c2;
        double4 wt2_c0, wt2_c1, wt2_c2;
        double bc2;

        // =================================================================
        // X+ direction
        // =================================================================
        if ((ix < Nx - 1) || (do_comm_x == 0)) {
            int isn = ((ix + 1) % Nx) + Nx * iyzt;
            int isg = ix + Nx * iyzt + Nst_pad * 0;
            bc2 = (ix == Nx - 1) ? (double)bc_x : 1.0;

            // Projection XP: (1-gx), t1=s0+i*s3, t2=s1+i*s2
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 3, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 2, DWF_PROJ_P)
            DWF_GMUL_FWD(u_up, isg)

            // Reconstruct: s0=t1, s1=t2, s2=-i*t2, s3=-i*t1
            DWF_ACCUM_4(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                        wt1_c0, wt2_c0, bc2, DWF_MULT_MI, DWF_MULT_MI)
            DWF_ACCUM_4(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                        wt1_c1, wt2_c1, bc2, DWF_MULT_MI, DWF_MULT_MI)
            DWF_ACCUM_4(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                        wt1_c2, wt2_c2, bc2, DWF_MULT_MI, DWF_MULT_MI)
        }

        // =================================================================
        // X- direction
        // =================================================================
        if ((ix > 0) || (do_comm_x == 0)) {
            int isn = ((ix - 1 + Nx) % Nx) + Nx * iyzt;
            int isg = isn + Nst_pad * 0;
            bc2 = (ix == 0) ? (double)bc_x : 1.0;

            // Projection XM: (1+gx), t1=s0-i*s3, t2=s1-i*s2
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 3, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 2, DWF_PROJ_M)
            DWF_GMUL_BCK(u_dn, isg)

            // Reconstruct: s0=t1, s1=t2, s2=+i*t2, s3=+i*t1
            DWF_ACCUM_4(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                        wt1_c0, wt2_c0, bc2, DWF_MULT_PI, DWF_MULT_PI)
            DWF_ACCUM_4(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                        wt1_c1, wt2_c1, bc2, DWF_MULT_PI, DWF_MULT_PI)
            DWF_ACCUM_4(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                        wt1_c2, wt2_c2, bc2, DWF_MULT_PI, DWF_MULT_PI)
        }

        // =================================================================
        // Y+ direction
        // =================================================================
        if ((iy < Ny - 1) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy + 1) % Ny) + Ny * izt);
            int isg = ix + Nx * (iy + Ny * izt) + Nst_pad * 1;
            bc2 = (iy == Ny - 1) ? (double)bc_y : 1.0;

            // Projection YP: t1=s0+i*s3, t2=s1-i*s2
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 3, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 2, DWF_PROJ_M)
            DWF_GMUL_FWD(u_up, isg)

            // Reconstruct: s0=t1, s1=t2, s2=+i*t2, s3=-i*t1
            DWF_ACCUM_4(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                        wt1_c0, wt2_c0, bc2, DWF_MULT_PI, DWF_MULT_MI)
            DWF_ACCUM_4(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                        wt1_c1, wt2_c1, bc2, DWF_MULT_PI, DWF_MULT_MI)
            DWF_ACCUM_4(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                        wt1_c2, wt2_c2, bc2, DWF_MULT_PI, DWF_MULT_MI)
        }

        // =================================================================
        // Y- direction
        // =================================================================
        if ((iy > 0) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy - 1 + Ny) % Ny) + Ny * izt);
            int isg = isn + Nst_pad * 1;
            bc2 = (iy == 0) ? (double)bc_y : 1.0;

            // Projection YM: t1=s0-i*s3, t2=s1+i*s2
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 3, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 2, DWF_PROJ_P)
            DWF_GMUL_BCK(u_dn, isg)

            // Reconstruct: s0=t1, s1=t2, s2=-i*t2, s3=+i*t1
            DWF_ACCUM_4(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                        wt1_c0, wt2_c0, bc2, DWF_MULT_MI, DWF_MULT_PI)
            DWF_ACCUM_4(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                        wt1_c1, wt2_c1, bc2, DWF_MULT_MI, DWF_MULT_PI)
            DWF_ACCUM_4(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                        wt1_c2, wt2_c2, bc2, DWF_MULT_MI, DWF_MULT_PI)
        }

        // =================================================================
        // Z+ direction
        // =================================================================
        if ((iz < Nz - 1) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz + 1) % Nz) + Nz * it);
            int isg = ixy + Nxy * (iz + Nz * it) + Nst_pad * 2;
            bc2 = (iz == Nz - 1) ? (double)bc_z : 1.0;

            // Projection ZP: t1=s0+i*s2, t2=s1-i*s3
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 2, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 3, DWF_PROJ_M)
            DWF_GMUL_FWD(u_up, isg)

            // Reconstruct: s0=t1, s1=t2, s2=-i*t1, s3=+i*t2
            // (Z: s2 comes from t1, s3 from t2 — use _SW variant)
            DWF_ACCUM_4_SW(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                           wt1_c0, wt2_c0, bc2, DWF_MULT_MI, DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                           wt1_c1, wt2_c1, bc2, DWF_MULT_MI, DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                           wt1_c2, wt2_c2, bc2, DWF_MULT_MI, DWF_MULT_PI)
        }

        // =================================================================
        // Z- direction
        // =================================================================
        if ((iz > 0) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz - 1 + Nz) % Nz) + Nz * it);
            int isg = isn + Nst_pad * 2;
            bc2 = (iz == 0) ? (double)bc_z : 1.0;

            // Projection ZM: t1=s0-i*s2, t2=s1+i*s3
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, 0, 2, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, 1, 3, DWF_PROJ_P)
            DWF_GMUL_BCK(u_dn, isg)

            // Reconstruct: s0=t1, s1=t2, s2=+i*t1, s3=-i*t2
            DWF_ACCUM_4_SW(v2_c0_s0, v2_c0_s1, v2_c0_s2, v2_c0_s3,
                           wt1_c0, wt2_c0, bc2, DWF_MULT_PI, DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c1_s0, v2_c1_s1, v2_c1_s2, v2_c1_s3,
                           wt1_c1, wt2_c1, bc2, DWF_MULT_PI, DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c2_s0, v2_c2_s1, v2_c2_s2, v2_c2_s3,
                           wt1_c2, wt2_c2, bc2, DWF_MULT_PI, DWF_MULT_MI)
        }

        // =================================================================
        // T+ direction (dirac: gamma4=diag(1,1,-1,-1), 1-g4 picks s2,s3)
        // =================================================================
        if ((it < Nt - 1) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it + 1) % Nt);
            int isg = ixyz + Nxyz * it + Nst_pad * 3;
            bc2 = (it == Nt - 1) ? (double)bc_t : 1.0;

            // Projection TP: t1=2*s2, t2=2*s3
            DWF_LOAD_PROJ_T(wp, is, Nst_pad, isn, 2, 3)
            DWF_GMUL_FWD(u_up, isg)

            // Reconstruct: s2 += wt1, s3 += wt2
            DWF_ACCUM_TP(v2_c0_s2, v2_c0_s3, wt1_c0, wt2_c0, bc2)
            DWF_ACCUM_TP(v2_c1_s2, v2_c1_s3, wt1_c1, wt2_c1, bc2)
            DWF_ACCUM_TP(v2_c2_s2, v2_c2_s3, wt1_c2, wt2_c2, bc2)
        }

        // =================================================================
        // T- direction (1+g4 picks s0,s1)
        // =================================================================
        if ((it > 0) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it - 1 + Nt) % Nt);
            int isg = isn + Nst_pad * 3;
            bc2 = (it == 0) ? (double)bc_t : 1.0;

            // Projection TM: t1=2*s0, t2=2*s1
            DWF_LOAD_PROJ_T(wp, is, Nst_pad, isn, 0, 1)
            DWF_GMUL_BCK(u_dn, isg)

            // Reconstruct: s0 += wt1, s1 += wt2
            DWF_ACCUM_TM(v2_c0_s0, v2_c0_s1, wt1_c0, wt2_c0, bc2)
            DWF_ACCUM_TM(v2_c1_s0, v2_c1_s1, wt1_c1, wt2_c1, bc2)
            DWF_ACCUM_TM(v2_c2_s0, v2_c2_s1, wt1_c2, wt2_c2, bc2)
        }

        // Store
        vp[IDX_DWF_QDW(0,0,is,Nst_pad,site)] = v2_c0_s0;
        vp[IDX_DWF_QDW(0,1,is,Nst_pad,site)] = v2_c0_s1;
        vp[IDX_DWF_QDW(0,2,is,Nst_pad,site)] = v2_c0_s2;
        vp[IDX_DWF_QDW(0,3,is,Nst_pad,site)] = v2_c0_s3;
        vp[IDX_DWF_QDW(1,0,is,Nst_pad,site)] = v2_c1_s0;
        vp[IDX_DWF_QDW(1,1,is,Nst_pad,site)] = v2_c1_s1;
        vp[IDX_DWF_QDW(1,2,is,Nst_pad,site)] = v2_c1_s2;
        vp[IDX_DWF_QDW(1,3,is,Nst_pad,site)] = v2_c1_s3;
        vp[IDX_DWF_QDW(2,0,is,Nst_pad,site)] = v2_c2_s0;
        vp[IDX_DWF_QDW(2,1,is,Nst_pad,site)] = v2_c2_s1;
        vp[IDX_DWF_QDW(2,2,is,Nst_pad,site)] = v2_c2_s2;
        vp[IDX_DWF_QDW(2,3,is,Nst_pad,site)] = v2_c2_s3;

    } // is loop
}


void mult_domainwall_5din_hopb_qdw_dirac(
    double *vp, double *up, double *wp, int Ns,
    int *bc, int *Nsize, int *do_comm, int flag)
{
    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst     = Nx * Ny * Nz * Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    real_t  *up_dev = (real_t  *)dev_ptr(up);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;

    mult_domainwall_5din_hopb_qdw_dirac_kernel<<<gridSize, blockSize>>>(
        vp_dev, up_dev, wp_dev, Ns,
        bc[0], bc[1], bc[2], bc[3],
        Nx, Ny, Nz, Nt,
        do_comm[0], do_comm[1], do_comm[2], do_comm[3],
        flag, Nst, Nst_pad);

    CHECK(cudaDeviceSynchronize());
}

// Undefine local helper macros to avoid polluting other includes
#undef DWF_PROJ_P
#undef DWF_PROJ_M
#undef DWF_MULT_MI
#undef DWF_MULT_PI
#undef DWF_PROJ_2
#undef IDX_DWF_QDW
#undef DWF_GMUL_FWD
#undef DWF_GMUL_BCK
#undef DWF_ACCUM_4
#undef DWF_ACCUM_4_SW
#undef DWF_ACCUM_TP
#undef DWF_ACCUM_TM
#undef DWF_LOAD_PROJ
#undef DWF_LOAD_PROJ2
#undef DWF_LOAD_PROJ_T

#endif // MULT_DOMAINWALL_5DIN_ACC_QDW_INCLUDED
