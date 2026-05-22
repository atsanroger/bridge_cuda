/*!
      @file    mult_Domainwall_5din_eo_cuda_qdw-inc.h
      @brief   QDW Even-Odd Domain Wall fermion kernels (ee/eo 5dir + EO hopb)
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED

// All macros were undef'd by mult_Domainwall_5din_cuda_qdw-inc.h — redefine here.

#define IDX_DWF_QDW(ic, id, is5, Nst_pad_, site_) \
    IDX2_QDW((ic), (id), (is5) * (Nst_pad_) + (site_))

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

#define DWF_PROJ_2(res, a) \
    dw_scal(2.0,(a).x,(a).z,(res).x,(res).z); \
    dw_scal(2.0,(a).y,(a).w,(res).y,(res).w);

#define DWF_GMUL_FWD(u_ptr, isg_) \
{ \
    double2 _u; double4 _t; \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_uc(_u,vt1_c0); wt2_c0=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    wt1_c1=qdw_mult_uc(_u,vt1_c0); wt2_c1=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    wt1_c2=qdw_mult_uc(_u,vt1_c0); wt2_c2=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=(u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

#define DWF_GMUL_BCK(u_ptr, isg_) \
{ \
    double2 _u; double4 _t; \
    _u.x=(u_ptr)[IDX2_G_R(0,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    wt1_c0=qdw_mult_uc(_u,vt1_c0); wt2_c0=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,0,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,0,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,0,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c0,wt1_c0,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c0,wt2_c0,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,1,(isg_))]; \
    wt1_c1=qdw_mult_uc(_u,vt1_c0); wt2_c1=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,1,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,1,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c1,wt1_c1,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c1,wt2_c1,_t); \
    _u.x=(u_ptr)[IDX2_G_R(0,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(0,2,(isg_))]; \
    wt1_c2=qdw_mult_uc(_u,vt1_c0); wt2_c2=qdw_mult_uc(_u,vt2_c0); \
    _u.x=(u_ptr)[IDX2_G_R(1,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(1,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c1); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c1); QDW_ADD(wt2_c2,wt2_c2,_t); \
    _u.x=(u_ptr)[IDX2_G_R(2,2,(isg_))]; _u.y=-(u_ptr)[IDX2_G_I(2,2,(isg_))]; \
    _t=qdw_mult_uc(_u,vt1_c2); QDW_ADD(wt1_c2,wt1_c2,_t); \
    _t=qdw_mult_uc(_u,vt2_c2); QDW_ADD(wt2_c2,wt2_c2,_t); \
}

#define DWF_ACCUM_4(v0,v1,v2,v3, W1,W2, bc_, OP2,OP3) \
{ \
    double4 _ts,_tp; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
    OP2(_tp,W2); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v2,v2,_ts); \
    OP3(_tp,W1); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v3,v3,_ts); \
}

#define DWF_ACCUM_4_SW(v0,v1,v2,v3, W1,W2, bc_, OP2,OP3) \
{ \
    double4 _ts,_tp; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
    OP2(_tp,W1); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v2,v2,_ts); \
    OP3(_tp,W2); QDW_SCAL(_ts,bc_,_tp); QDW_ADD(v3,v3,_ts); \
}

#define DWF_ACCUM_TP(v2,v3, W1,W2, bc_) \
{ \
    double4 _ts; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v2,v2,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v3,v3,_ts); \
}
#define DWF_ACCUM_TM(v0,v1, W1,W2, bc_) \
{ \
    double4 _ts; \
    QDW_SCAL(_ts,bc_,W1); QDW_ADD(v0,v0,_ts); \
    QDW_SCAL(_ts,bc_,W2); QDW_ADD(v1,v1,_ts); \
}

#define DWF_LOAD_PROJ(wp_,is_,Np_,isn_,id_a,id_b,PROJ) \
{ \
    double4 _sa,_sb; \
    _sa=(wp_)[IDX_DWF_QDW(0,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(0,id_b,is_,Np_,isn_)]; PROJ(vt1_c0,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(1,id_b,is_,Np_,isn_)]; PROJ(vt1_c1,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(2,id_b,is_,Np_,isn_)]; PROJ(vt1_c2,_sa,_sb); \
}
#define DWF_LOAD_PROJ2(wp_,is_,Np_,isn_,id_a,id_b,PROJ) \
{ \
    double4 _sa,_sb; \
    _sa=(wp_)[IDX_DWF_QDW(0,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(0,id_b,is_,Np_,isn_)]; PROJ(vt2_c0,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(1,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(1,id_b,is_,Np_,isn_)]; PROJ(vt2_c1,_sa,_sb); \
    _sa=(wp_)[IDX_DWF_QDW(2,id_a,is_,Np_,isn_)]; \
    _sb=(wp_)[IDX_DWF_QDW(2,id_b,is_,Np_,isn_)]; PROJ(vt2_c2,_sa,_sb); \
}
#define DWF_LOAD_PROJ_T(wp_,is_,Np_,isn_,id_a,id_b) \
{ \
    double4 _sa; \
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
    double4 _t1, _t2; \
    QDW_SCAL(_t1, fac1_, a_); QDW_SCAL(_t2, fac2_, b_); QDW_ADD(res, _t1, _t2); \
}

//====================================================================
// ee 5-direction (diagonal block)  vp = B_eff * wp + C_hop * vt
//====================================================================
__global__
void mult_domainwall_5din_ee_5dir_dirac_qdw_kernel(
    double4 * __restrict__ vp, double4 * __restrict__ wp,
    double mq, double M0, int Ns, double alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const double *b_con = const_b;
    const double *c_con = const_c;

    for (int is = 0; is < Ns; ++is) {
        const double FF1 = b_con[is] * (4.0 - M0) + 1.0;
        const double FF2 = c_con[is] * (4.0 - M0) - 1.0;

        // ---- s+1 hopping ----
        const int is_up = (is + 1) % Ns;
        const double Fup = (is == Ns-1) ? -0.5 * mq * FF2 : 0.5 * FF2 * alpha;
        double4 tmp;

#define LOAD_SUB_SHop(ic) \
        double4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Nst_pad,site)]; \
        double4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL(vt_##ic##_s0,Fup,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL(vt_##ic##_s1,Fup,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_SHop(0) LOAD_SUB_SHop(1) LOAD_SUB_SHop(2)
#undef LOAD_SUB_SHop

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const double Fdn = (is == 0) ? -0.5 * mq * FF2 : 0.5 * FF2 * alpha;
        double4 sum02, sum13;

#define ADD_SUM_SHop(ic) \
        { \
        double4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Nst_pad,site)]; \
        double4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Nst_pad,site)]; \
        double4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Nst_pad,site)]; \
        double4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Nst_pad,site)]; \
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
        double4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]; \
        double4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]; \
        double4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]; \
        double4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]; \
        double4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0) { \
            double f1=FF1*0.5*(1.0+alpha), f2=FF1*0.5*(-1.0+alpha); \
            DWF_SCAL_ADD2(wt_s0,f1,w_s0,f2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,f1,w_s1,f2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,f1,w_s2,f2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,f1,w_s3,f2,w_s1); \
        } else if (is == Ns-1) { \
            double f1=FF1*0.5*(1.0+alpha), f2=FF1*0.5*(1.0-alpha); \
            DWF_SCAL_ADD2(wt_s0,f1,w_s0,f2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,f1,w_s1,f2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,f1,w_s2,f2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,f1,w_s3,f2,w_s1); \
        } else { \
            double fa=FF1*alpha; \
            QDW_SCAL(wt_s0,fa,w_s0); QDW_SCAL(wt_s1,fa,w_s1); \
            QDW_SCAL(wt_s2,fa,w_s2); QDW_SCAL(wt_s3,fa,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]=wt_s0; \
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]=wt_s1; \
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]=wt_s2; \
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]=wt_s3; \
        }

        DIAG_AND_STORE(0) DIAG_AND_STORE(1) DIAG_AND_STORE(2)
#undef DIAG_AND_STORE
    }
}

void mult_domainwall_5din_ee_5dir_dirac_qdw(
    double *vp, double *wp, double mq, double M0, int Ns,
    double *b, double *c, double alpha, int *Nsize)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    mult_domainwall_5din_ee_5dir_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, mq, M0, Ns, alpha, Nst, Nst_pad);
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-direction (off-diagonal block)  yp = -0.5*(b*w_diag + c*vt)
//====================================================================
__global__
void mult_domainwall_5din_eo_5dir_dirac_qdw_kernel(
    double4 * __restrict__ yp, double4 * __restrict__ wp,
    double mq, int Ns, double alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const double *b_con = const_b;
    const double *c_con = const_c;

    for (int is = 0; is < Ns; ++is) {
        double4 tmp;

        // ---- s+1 hopping (C-scaled) ----
        const int is_up = (is + 1) % Ns;
        const double Fup = (is == Ns-1) ? -0.5 * mq * c_con[is] : 0.5 * c_con[is] * alpha;

#define LOAD_SUB_EO(ic) \
        double4 wu_##ic##_s0 = wp[IDX_DWF_QDW(ic,0,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s1 = wp[IDX_DWF_QDW(ic,1,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s2 = wp[IDX_DWF_QDW(ic,2,is_up,Nst_pad,site)]; \
        double4 wu_##ic##_s3 = wp[IDX_DWF_QDW(ic,3,is_up,Nst_pad,site)]; \
        double4 vt_##ic##_s0, vt_##ic##_s1, vt_##ic##_s2, vt_##ic##_s3; \
        QDW_SUB(tmp,wu_##ic##_s0,wu_##ic##_s2); QDW_SCAL(vt_##ic##_s0,Fup,tmp); \
        QDW_SUB(tmp,wu_##ic##_s1,wu_##ic##_s3); QDW_SCAL(vt_##ic##_s1,Fup,tmp); \
        QDW_NEG(vt_##ic##_s2,vt_##ic##_s0); \
        QDW_NEG(vt_##ic##_s3,vt_##ic##_s1);

        LOAD_SUB_EO(0) LOAD_SUB_EO(1) LOAD_SUB_EO(2)
#undef LOAD_SUB_EO

        // ---- s-1 hopping ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const double Fdn = (is == 0) ? -0.5 * mq * c_con[is] : 0.5 * c_con[is] * alpha;
        double4 sum02, sum13;

#define ADD_SUM_EO(ic) \
        { \
        double4 wd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Nst_pad,site)]; \
        double4 wd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Nst_pad,site)]; \
        double4 wd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Nst_pad,site)]; \
        double4 wd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Nst_pad,site)]; \
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
        double4 w_s0=wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]; \
        double4 w_s1=wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]; \
        double4 w_s2=wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]; \
        double4 w_s3=wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]; \
        double4 wt_s0,wt_s1,wt_s2,wt_s3; \
        if (is == 0) { \
            double b1=b_con[is]*0.5*(1.0+alpha), b2=b_con[is]*0.5*(-1.0+alpha); \
            DWF_SCAL_ADD2(wt_s0,b1,w_s0,b2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,b1,w_s1,b2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,b1,w_s2,b2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,b1,w_s3,b2,w_s1); \
        } else if (is == Ns-1) { \
            double b1=b_con[is]*0.5*(1.0+alpha), b2=b_con[is]*0.5*(1.0-alpha); \
            DWF_SCAL_ADD2(wt_s0,b1,w_s0,b2,w_s2); \
            DWF_SCAL_ADD2(wt_s1,b1,w_s1,b2,w_s3); \
            DWF_SCAL_ADD2(wt_s2,b1,w_s2,b2,w_s0); \
            DWF_SCAL_ADD2(wt_s3,b1,w_s3,b2,w_s1); \
        } else { \
            double bb=b_con[is]*alpha; \
            QDW_SCAL(wt_s0,bb,w_s0); QDW_SCAL(wt_s1,bb,w_s1); \
            QDW_SCAL(wt_s2,bb,w_s2); QDW_SCAL(wt_s3,bb,w_s3); \
        } \
        QDW_ADD(wt_s0,wt_s0,vt_##ic##_s0); \
        QDW_ADD(wt_s1,wt_s1,vt_##ic##_s1); \
        QDW_ADD(wt_s2,wt_s2,vt_##ic##_s2); \
        QDW_ADD(wt_s3,wt_s3,vt_##ic##_s3); \
        QDW_SCAL(wt_s0,-0.5,wt_s0); \
        QDW_SCAL(wt_s1,-0.5,wt_s1); \
        QDW_SCAL(wt_s2,-0.5,wt_s2); \
        QDW_SCAL(wt_s3,-0.5,wt_s3); \
        yp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]=wt_s0; \
        yp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]=wt_s1; \
        yp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]=wt_s2; \
        yp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]=wt_s3; \
        }

        DIAG_AND_STORE_EO(0) DIAG_AND_STORE_EO(1) DIAG_AND_STORE_EO(2)
#undef DIAG_AND_STORE_EO
    }
}

void mult_domainwall_5din_eo_5dir_dirac_qdw(
    double *yp, double *wp, double mq, double M0, int Ns,
    double *b, double *c, double alpha, int *Nsize)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *yp_dev = (double4 *)dev_ptr(yp);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    mult_domainwall_5din_eo_5dir_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        yp_dev, wp_dev, mq, Ns, alpha, Nst, Nst_pad);
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// ee 5-dirdag (adjoint diagonal block)
//====================================================================
__global__
void mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel(
    double4 * __restrict__ vp, double4 * __restrict__ wp,
    double mq, double M0, int Ns, double alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const double *b_con = const_b;
    const double *c_con = const_c;

    for (int is = 0; is < Ns; ++is) {
        const double B1 = b_con[is] * (4.0 - M0) + 1.0;
        double4 tmp;

        // ---- local diagonal (self-adjoint) ----
#define DIAG_DAG(ic) \
        double4 w_##ic##_s0=wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]; \
        double4 w_##ic##_s1=wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]; \
        double4 w_##ic##_s2=wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]; \
        double4 w_##ic##_s3=wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]; \
        double4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0) { \
            double f1=B1*0.5*(1.0+alpha), f2=B1*0.5*(-1.0+alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,f1,w_##ic##_s0,f2,w_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s1,f1,w_##ic##_s1,f2,w_##ic##_s3); \
            DWF_SCAL_ADD2(vt_##ic##_s2,f1,w_##ic##_s2,f2,w_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s3,f1,w_##ic##_s3,f2,w_##ic##_s1); \
        } else if (is == Ns-1) { \
            double f1=B1*0.5*(1.0+alpha), f2=B1*0.5*(1.0-alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,f1,w_##ic##_s0,f2,w_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s1,f1,w_##ic##_s1,f2,w_##ic##_s3); \
            DWF_SCAL_ADD2(vt_##ic##_s2,f1,w_##ic##_s2,f2,w_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s3,f1,w_##ic##_s3,f2,w_##ic##_s1); \
        } else { \
            double fa=B1*alpha; \
            QDW_SCAL(vt_##ic##_s0,fa,w_##ic##_s0); QDW_SCAL(vt_##ic##_s1,fa,w_##ic##_s1); \
            QDW_SCAL(vt_##ic##_s2,fa,w_##ic##_s2); QDW_SCAL(vt_##ic##_s3,fa,w_##ic##_s3); \
        }

        DIAG_DAG(0) DIAG_DAG(1) DIAG_DAG(2)
#undef DIAG_DAG

        // ---- s+1 contribution (adjoint: symmetric) ----
        const int is_up = (is + 1) % Ns;
        const double C2up = c_con[is_up] * (4.0 - M0) - 1.0;
        const double Fup = (is == Ns-1) ? -0.5 * mq * C2up * 0.5 : 0.5 * C2up * alpha * 0.5;
        double4 sum02, sum13;

#define ADD_SYM_UP(ic) \
        { \
        double4 xu_s0=wp[IDX_DWF_QDW(ic,0,is_up,Nst_pad,site)]; \
        double4 xu_s1=wp[IDX_DWF_QDW(ic,1,is_up,Nst_pad,site)]; \
        double4 xu_s2=wp[IDX_DWF_QDW(ic,2,is_up,Nst_pad,site)]; \
        double4 xu_s3=wp[IDX_DWF_QDW(ic,3,is_up,Nst_pad,site)]; \
        QDW_ADD(sum02,xu_s0,xu_s2); QDW_SCAL(tmp,Fup,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,xu_s1,xu_s3); QDW_SCAL(tmp,Fup,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_SYM_UP(0) ADD_SYM_UP(1) ADD_SYM_UP(2)
#undef ADD_SYM_UP

        // ---- s-1 contribution (adjoint: antisymmetric) ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const double C2dn = c_con[is_dn] * (4.0 - M0) - 1.0;
        const double Fdn = (is == 0) ? -0.5 * mq * C2dn * 0.5 : 0.5 * C2dn * alpha * 0.5;

#define ADD_ANTISYM_DN(ic) \
        { \
        double4 xd_s0=wp[IDX_DWF_QDW(ic,0,is_dn,Nst_pad,site)]; \
        double4 xd_s1=wp[IDX_DWF_QDW(ic,1,is_dn,Nst_pad,site)]; \
        double4 xd_s2=wp[IDX_DWF_QDW(ic,2,is_dn,Nst_pad,site)]; \
        double4 xd_s3=wp[IDX_DWF_QDW(ic,3,is_dn,Nst_pad,site)]; \
        QDW_SUB(sum02,xd_s0,xd_s2); QDW_SCAL(tmp,Fdn,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,xd_s1,xd_s3); QDW_SCAL(tmp,Fdn,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_ANTISYM_DN(0) ADD_ANTISYM_DN(1) ADD_ANTISYM_DN(2)
#undef ADD_ANTISYM_DN

        vp[IDX_DWF_QDW(0,0,is,Nst_pad,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Nst_pad,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Nst_pad,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Nst_pad,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Nst_pad,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Nst_pad,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Nst_pad,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Nst_pad,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Nst_pad,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Nst_pad,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Nst_pad,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Nst_pad,site)]=vt_2_s3;
    }
}

void mult_domainwall_5din_ee_5dirdag_dirac_qdw(
    double *vp, double *wp, double mq, double M0, int Ns,
    double *b, double *c, double alpha, int *Nsize)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    mult_domainwall_5din_ee_5dirdag_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, mq, M0, Ns, alpha, Nst, Nst_pad);
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-dirdag (adjoint off-diagonal block)  vp = (M_eo)† * yp
//====================================================================
__global__
void mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel(
    double4 * __restrict__ vp, double4 * __restrict__ yp,
    double mq, int Ns, double alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const double *b_con = const_b;
    const double *c_con = const_c;

    for (int is = 0; is < Ns; ++is) {
        double4 tmp;

        // ---- local adjoint (adjoint of -0.5*b[is]*B_alpha) ----
#define LOCAL_DAG_EO(ic) \
        double4 y_##ic##_s0=yp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)]; \
        double4 y_##ic##_s1=yp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)]; \
        double4 y_##ic##_s2=yp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)]; \
        double4 y_##ic##_s3=yp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)]; \
        double4 vt_##ic##_s0,vt_##ic##_s1,vt_##ic##_s2,vt_##ic##_s3; \
        if (is == 0) { \
            double b1=-0.5*b_con[is]*0.5*(1.0+alpha); \
            double b2=-0.5*b_con[is]*0.5*(-1.0+alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,b1,y_##ic##_s2,b2,y_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s1,b1,y_##ic##_s3,b2,y_##ic##_s1); \
            DWF_SCAL_ADD2(vt_##ic##_s2,b1,y_##ic##_s0,b2,y_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s3,b1,y_##ic##_s1,b2,y_##ic##_s3); \
        } else if (is == Ns-1) { \
            double b1=-0.5*b_con[is]*0.5*(1.0+alpha); \
            double b2=-0.5*b_con[is]*0.5*(1.0-alpha); \
            DWF_SCAL_ADD2(vt_##ic##_s0,b1,y_##ic##_s2,b2,y_##ic##_s0); \
            DWF_SCAL_ADD2(vt_##ic##_s1,b1,y_##ic##_s3,b2,y_##ic##_s1); \
            DWF_SCAL_ADD2(vt_##ic##_s2,b1,y_##ic##_s0,b2,y_##ic##_s2); \
            DWF_SCAL_ADD2(vt_##ic##_s3,b1,y_##ic##_s1,b2,y_##ic##_s3); \
        } else { \
            double bb=-0.5*b_con[is]*alpha; \
            QDW_SCAL(vt_##ic##_s0,bb,y_##ic##_s2); \
            QDW_SCAL(vt_##ic##_s1,bb,y_##ic##_s3); \
            QDW_SCAL(vt_##ic##_s2,bb,y_##ic##_s0); \
            QDW_SCAL(vt_##ic##_s3,bb,y_##ic##_s1); \
        }

        LOCAL_DAG_EO(0) LOCAL_DAG_EO(1) LOCAL_DAG_EO(2)
#undef LOCAL_DAG_EO

        // ---- s+1 contribution: adjoint of c-hop, symmetric ----
        const int is_up = (is + 1) % Ns;
        const double Fup_d = (is == Ns-1)
            ? -0.5 * c_con[is_up] * (-0.5) * mq
            :  -0.5 * c_con[is_up] * 0.5 * alpha;
        double4 sum02, sum13;

#define ADD_DAG_UP(ic) \
        { \
        double4 yu_s0=yp[IDX_DWF_QDW(ic,0,is_up,Nst_pad,site)]; \
        double4 yu_s1=yp[IDX_DWF_QDW(ic,1,is_up,Nst_pad,site)]; \
        double4 yu_s2=yp[IDX_DWF_QDW(ic,2,is_up,Nst_pad,site)]; \
        double4 yu_s3=yp[IDX_DWF_QDW(ic,3,is_up,Nst_pad,site)]; \
        QDW_ADD(sum02,yu_s0,yu_s2); QDW_SCAL(tmp,Fup_d,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_ADD(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_ADD(sum13,yu_s1,yu_s3); QDW_SCAL(tmp,Fup_d,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_ADD(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_UP(0) ADD_DAG_UP(1) ADD_DAG_UP(2)
#undef ADD_DAG_UP

        // ---- s-1 contribution: adjoint of c-hop, antisymmetric ----
        const int is_dn = (is - 1 + Ns) % Ns;
        const double Fdn_d = (is == 0)
            ? -0.5 * c_con[is_dn] * (-0.5) * mq
            :  -0.5 * c_con[is_dn] * 0.5 * alpha;

#define ADD_DAG_DN(ic) \
        { \
        double4 yd_s0=yp[IDX_DWF_QDW(ic,0,is_dn,Nst_pad,site)]; \
        double4 yd_s1=yp[IDX_DWF_QDW(ic,1,is_dn,Nst_pad,site)]; \
        double4 yd_s2=yp[IDX_DWF_QDW(ic,2,is_dn,Nst_pad,site)]; \
        double4 yd_s3=yp[IDX_DWF_QDW(ic,3,is_dn,Nst_pad,site)]; \
        QDW_SUB(sum02,yd_s2,yd_s0); QDW_SCAL(tmp,Fdn_d,sum02); \
        QDW_ADD(vt_##ic##_s0,vt_##ic##_s0,tmp); QDW_SUB(vt_##ic##_s2,vt_##ic##_s2,tmp); \
        QDW_SUB(sum13,yd_s3,yd_s1); QDW_SCAL(tmp,Fdn_d,sum13); \
        QDW_ADD(vt_##ic##_s1,vt_##ic##_s1,tmp); QDW_SUB(vt_##ic##_s3,vt_##ic##_s3,tmp); \
        }

        ADD_DAG_DN(0) ADD_DAG_DN(1) ADD_DAG_DN(2)
#undef ADD_DAG_DN

        vp[IDX_DWF_QDW(0,0,is,Nst_pad,site)]=vt_0_s0; vp[IDX_DWF_QDW(0,1,is,Nst_pad,site)]=vt_0_s1;
        vp[IDX_DWF_QDW(0,2,is,Nst_pad,site)]=vt_0_s2; vp[IDX_DWF_QDW(0,3,is,Nst_pad,site)]=vt_0_s3;
        vp[IDX_DWF_QDW(1,0,is,Nst_pad,site)]=vt_1_s0; vp[IDX_DWF_QDW(1,1,is,Nst_pad,site)]=vt_1_s1;
        vp[IDX_DWF_QDW(1,2,is,Nst_pad,site)]=vt_1_s2; vp[IDX_DWF_QDW(1,3,is,Nst_pad,site)]=vt_1_s3;
        vp[IDX_DWF_QDW(2,0,is,Nst_pad,site)]=vt_2_s0; vp[IDX_DWF_QDW(2,1,is,Nst_pad,site)]=vt_2_s1;
        vp[IDX_DWF_QDW(2,2,is,Nst_pad,site)]=vt_2_s2; vp[IDX_DWF_QDW(2,3,is,Nst_pad,site)]=vt_2_s3;
    }
}

void mult_domainwall_5din_eo_5dirdag_dirac_qdw(
    double *vp, double *yp, double mq, double M0, int Ns,
    double *b, double *c, double alpha, int *Nsize)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    double4 *yp_dev = (double4 *)dev_ptr(yp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    mult_domainwall_5din_eo_5dirdag_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, yp_dev, mq, Ns, alpha, Nst, Nst_pad);
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// EO 4D bulk hopping (QDW, EO gauge layout, jgm5 support)
// Gauge stride: Nst*(ieo + 2*idir) for forward, Nst*(1-ieo + 2*idir) for backward
//====================================================================
__global__
void mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel(
    double4 * __restrict__ vp, real_t * __restrict__ up, double4 * __restrict__ wp,
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

    real_t *u_up = up;
    real_t *u_dn = up;

    for (int is = 0; is < Ns; ++is) {

        double4 v2_c0_s0={0,0,0,0}, v2_c0_s1={0,0,0,0},
                v2_c0_s2={0,0,0,0}, v2_c0_s3={0,0,0,0};
        double4 v2_c1_s0={0,0,0,0}, v2_c1_s1={0,0,0,0},
                v2_c1_s2={0,0,0,0}, v2_c1_s3={0,0,0,0};
        double4 v2_c2_s0={0,0,0,0}, v2_c2_s1={0,0,0,0},
                v2_c2_s2={0,0,0,0}, v2_c2_s3={0,0,0,0};

        double4 vt1_c0, vt1_c1, vt1_c2;
        double4 vt2_c0, vt2_c1, vt2_c2;
        double4 wt1_c0, wt1_c1, wt1_c2;
        double4 wt2_c0, wt2_c1, wt2_c2;
        double bc2;

        // X+
        {
            int isn = ((ix + keo) % Nx) + Nx * iyzt;
            int isg = site + Nst * (ieo + 2*0);
            bc2 = (ix == Nx-1 && keo == 1) ? (double)bc_x : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d3, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d2, DWF_PROJ_P)
            DWF_GMUL_FWD(u_up, isg)
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_MI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_MI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_MI,DWF_MULT_MI)
        }

        // X-
        {
            int ix2 = (ix - 1 + keo + Nx) % Nx;
            int isn = ix2 + Nx * iyzt;
            int isg = isn + Nst * (1-ieo + 2*0);
            bc2 = (ix == 0 && keo == 0) ? (double)bc_x : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d3, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d2, DWF_PROJ_M)
            DWF_GMUL_BCK(u_dn, isg)
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_PI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_PI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_PI,DWF_MULT_PI)
        }

        // Y+
        if ((iy < Ny-1) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy+1)%Ny) + Ny*izt);
            int isg = site + Nst * (ieo + 2*1);
            bc2 = (iy == Ny-1) ? (double)bc_y : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d3, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d2, DWF_PROJ_M)
            DWF_GMUL_FWD(u_up, isg)
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_PI,DWF_MULT_MI)
        }

        // Y-
        if ((iy > 0) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy-1+Ny)%Ny) + Ny*izt);
            int isg = isn + Nst * (1-ieo + 2*1);
            bc2 = (iy == 0) ? (double)bc_y : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d3, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d2, DWF_PROJ_P)
            DWF_GMUL_BCK(u_dn, isg)
            DWF_ACCUM_4(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_MI,DWF_MULT_PI)
        }

        // Z+
        if ((iz < Nz-1) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz+1)%Nz) + Nz*it);
            int isg = site + Nst * (ieo + 2*2);
            bc2 = (iz == Nz-1) ? (double)bc_z : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d2, DWF_PROJ_P)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d3, DWF_PROJ_M)
            DWF_GMUL_FWD(u_up, isg)
            DWF_ACCUM_4_SW(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_MI,DWF_MULT_PI)
            DWF_ACCUM_4_SW(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_MI,DWF_MULT_PI)
        }

        // Z-
        if ((iz > 0) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz-1+Nz)%Nz) + Nz*it);
            int isg = isn + Nst * (1-ieo + 2*2);
            bc2 = (iz == 0) ? (double)bc_z : 1.0;
            DWF_LOAD_PROJ (wp, is, Nst_pad, isn, d0, d2, DWF_PROJ_M)
            DWF_LOAD_PROJ2(wp, is, Nst_pad, isn, d1, d3, DWF_PROJ_P)
            DWF_GMUL_BCK(u_dn, isg)
            DWF_ACCUM_4_SW(v2_c0_s0,v2_c0_s1,v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c1_s0,v2_c1_s1,v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2, DWF_MULT_PI,DWF_MULT_MI)
            DWF_ACCUM_4_SW(v2_c2_s0,v2_c2_s1,v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2, DWF_MULT_PI,DWF_MULT_MI)
        }

        // T+
        if ((it < Nt-1) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it+1)%Nt);
            int isg = site + Nst * (ieo + 2*3);
            bc2 = (it == Nt-1) ? (double)bc_t : 1.0;
            DWF_LOAD_PROJ_T(wp, is, Nst_pad, isn, d2, d3)
            DWF_GMUL_FWD(u_up, isg)
            DWF_ACCUM_TP(v2_c0_s2,v2_c0_s3, wt1_c0,wt2_c0, bc2)
            DWF_ACCUM_TP(v2_c1_s2,v2_c1_s3, wt1_c1,wt2_c1, bc2)
            DWF_ACCUM_TP(v2_c2_s2,v2_c2_s3, wt1_c2,wt2_c2, bc2)
        }

        // T-
        if ((it > 0) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it-1+Nt)%Nt);
            int isg = isn + Nst * (1-ieo + 2*3);
            bc2 = (it == 0) ? (double)bc_t : 1.0;
            DWF_LOAD_PROJ_T(wp, is, Nst_pad, isn, d0, d1)
            DWF_GMUL_BCK(u_dn, isg)
            DWF_ACCUM_TM(v2_c0_s0,v2_c0_s1, wt1_c0,wt2_c0, bc2)
            DWF_ACCUM_TM(v2_c1_s0,v2_c1_s1, wt1_c1,wt2_c1, bc2)
            DWF_ACCUM_TM(v2_c2_s0,v2_c2_s1, wt1_c2,wt2_c2, bc2)
        }

        vp[IDX_DWF_QDW(0,0,is,Nst_pad,site)]=v2_c0_s0; vp[IDX_DWF_QDW(0,1,is,Nst_pad,site)]=v2_c0_s1;
        vp[IDX_DWF_QDW(0,2,is,Nst_pad,site)]=v2_c0_s2; vp[IDX_DWF_QDW(0,3,is,Nst_pad,site)]=v2_c0_s3;
        vp[IDX_DWF_QDW(1,0,is,Nst_pad,site)]=v2_c1_s0; vp[IDX_DWF_QDW(1,1,is,Nst_pad,site)]=v2_c1_s1;
        vp[IDX_DWF_QDW(1,2,is,Nst_pad,site)]=v2_c1_s2; vp[IDX_DWF_QDW(1,3,is,Nst_pad,site)]=v2_c1_s3;
        vp[IDX_DWF_QDW(2,0,is,Nst_pad,site)]=v2_c2_s0; vp[IDX_DWF_QDW(2,1,is,Nst_pad,site)]=v2_c2_s1;
        vp[IDX_DWF_QDW(2,2,is,Nst_pad,site)]=v2_c2_s2; vp[IDX_DWF_QDW(2,3,is,Nst_pad,site)]=v2_c2_s3;
    }
}

void mult_domainwall_5din_eo_hopb_qdw_dirac_5d(
    double *vp, double *up, double *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo, int jgm5)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst     = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    double4 *vp_dev = (double4 *)dev_ptr(vp);
    real_t  *up_dev = (real_t  *)dev_ptr(up);
    double4 *wp_dev = (double4 *)dev_ptr(wp);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;

    mult_domainwall_5din_eo_hopb_qdw_dirac_5d_kernel<<<gridSize, blockSize>>>(
        vp_dev, up_dev, wp_dev, Ns,
        bc[0], bc[1], bc[2], bc[3],
        Nx, Ny, Nz, Nt,
        ieo, jeo,
        do_comm[0], do_comm[1], do_comm[2], do_comm[3],
        Nst, Nst_pad, jgm5);

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
#undef DWF_ACCUM_4
#undef DWF_ACCUM_4_SW
#undef DWF_ACCUM_TP
#undef DWF_ACCUM_TM
#undef DWF_LOAD_PROJ
#undef DWF_LOAD_PROJ2
#undef DWF_LOAD_PROJ_T
#undef DWF_SCAL_ADD2

#endif // MULT_DOMAINWALL_5DIN_EO_ACC_QDW_INCLUDED
