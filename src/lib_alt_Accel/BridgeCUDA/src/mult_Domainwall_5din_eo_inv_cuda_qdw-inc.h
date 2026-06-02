/*!
      @file    mult_Domainwall_5din_eo_inv_cuda_qdw-inc.h
      @brief   QDW LU/LUdag inverse kernels for 5D even-odd domain-wall fermion.
               QDW is real_t-precision only.
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED

#include <type_traits>
#include "../inline/constant_memory_inline.h"

// IDX_DWF_QDW: BLAS-compatible flat layout — matches 4*IDX2(nin/4, in4, site)+k in BLAS.
// 4th arg is Ns (not Nst_pad) — same flat ordering as afield_cuda-inc.h QDW kernels.
#define IDX_DWF_QDW(ic, id, is5, Ns_, site_) \
    IDX2(NC * ND * (Ns_), (ic) + NC * ((id) + ND * (is5)), (site_))

// ===================================================================
// LU_inv QDW kernel
// Thread → (ic, site): NC * Nst threads.
// Each thread performs forward (L^{-1}) and backward (D^{-1}*U^{-1})
// sweep over s-slices for one (ic, site) pair with QDW arithmetic.
// ===================================================================
// Single-precision (or matched-T) LU_inv kernel. Alpha is a single T scalar.
// All arithmetic at T precision (float-float when T=float in QDW vectors,
// double-double when T=double).
__global__ void mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel_fp(
    real4* vp, real4* wp, int Ns, int Nst, int Nst_pad, real_t alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
    const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
    const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm();

    real4 vt0, vt1, vt2, vt3;
    real4 xt0, xt1, xt2, xt3;
    real4 yt0, yt1, yt2, yt3;
    real4 tmp, sum, diff;

    // === Forward sweep: L^{-1} ===
    vt0 = wp[IDX_DWF_QDW(ic,0,0,Ns,site)];
    vt1 = wp[IDX_DWF_QDW(ic,1,0,Ns,site)];
    vt2 = wp[IDX_DWF_QDW(ic,2,0,Ns,site)];
    vt3 = wp[IDX_DWF_QDW(ic,3,0,Ns,site)];
    vp[IDX_DWF_QDW(ic,0,0,Ns,site)] = vt0;
    vp[IDX_DWF_QDW(ic,1,0,Ns,site)] = vt1;
    vp[IDX_DWF_QDW(ic,2,0,Ns,site)] = vt2;
    vp[IDX_DWF_QDW(ic,3,0,Ns,site)] = vt3;
    QDW_SCAL(yt0, e_con[0], vt0);
    QDW_SCAL(yt1, e_con[0], vt1);
    QDW_SCAL(yt2, e_con[0], vt2);
    QDW_SCAL(yt3, e_con[0], vt3);

    for (int is = 1; is < Ns-1; ++is) {
        const real_t a   = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
        const real_t eis = e_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SCAL(tmp, eis, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL(tmp, eis, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL(tmp, eis, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL(tmp, eis, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    {
        const int is     = Ns-1;
        const real_t a_l = real_t(0.5) * dm_con[is] * dpinv_con[is-1];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a_l, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a_l, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        QDW_SUB(diff, yt0, yt2); QDW_SCAL(tmp, real_t(0.5), diff);
        QDW_SUB(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_SUB(diff, yt1, yt3); QDW_SCAL(tmp, real_t(0.5), diff);
        QDW_SUB(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
    }

    // === Backward sweep: D^{-1} * U^{-1} ===
    {
        const int is = Ns-1;
        const real_t aa = dpinv_con[Ns-1];
        const real_t f1 = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t f2 = real_t(0.5) * (real_t(-1.0) + alpha);
        real4 p0, p1, p2, p3, t;
        p0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        p1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        p2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        p3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SCAL(vt0,f1,p0); QDW_SCAL(t,f2,p2); QDW_ADD(vt0,vt0,t); QDW_SCAL(vt0,aa,vt0);
        QDW_SCAL(vt1,f1,p1); QDW_SCAL(t,f2,p3); QDW_ADD(vt1,vt1,t); QDW_SCAL(vt1,aa,vt1);
        QDW_SCAL(vt2,f1,p2); QDW_SCAL(t,f2,p0); QDW_ADD(vt2,vt2,t); QDW_SCAL(vt2,aa,vt2);
        QDW_SCAL(vt3,f1,p3); QDW_SCAL(t,f2,p1); QDW_ADD(vt3,vt3,t); QDW_SCAL(vt3,aa,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_ADD(sum, vt0, vt2); QDW_SCAL(yt0, real_t(0.5), sum); yt2 = yt0;
        QDW_ADD(sum, vt1, vt3); QDW_SCAL(yt1, real_t(0.5), sum); yt3 = yt1;
    }

    for (int is = Ns-2; is >= 0; --is) {
        const real_t a     = real_t(0.5) * dm_con[is];
        const real_t f_is  = f_con[is];
        const real_t aa_is = dpinv_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SUB(diff, xt0, xt2); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL(tmp, f_is, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL(tmp, f_is, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL(tmp, f_is, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL(tmp, f_is, yt3); QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL(vt0, aa_is, vt0);
        QDW_SCAL(vt1, aa_is, vt1);
        QDW_SCAL(vt2, aa_is, vt2);
        QDW_SCAL(vt3, aa_is, vt3);
        if (is == 0) {
            const real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
            const real_t ff2 = real_t(0.5) * (real_t(1.0) - alpha);
            real4 r0, r1, r2, r3, t;
            QDW_SCAL(r0,ff1,vt0); QDW_SCAL(t,ff2,vt2); QDW_ADD(r0,r0,t);
            QDW_SCAL(r1,ff1,vt1); QDW_SCAL(t,ff2,vt3); QDW_ADD(r1,r1,t);
            QDW_SCAL(r2,ff1,vt2); QDW_SCAL(t,ff2,vt0); QDW_ADD(r2,r2,t);
            QDW_SCAL(r3,ff1,vt3); QDW_SCAL(t,ff2,vt1); QDW_ADD(r3,r3,t);
            vp[IDX_DWF_QDW(ic,0,0,Ns,site)] = r0;
            vp[IDX_DWF_QDW(ic,1,0,Ns,site)] = r1;
            vp[IDX_DWF_QDW(ic,2,0,Ns,site)] = r2;
            vp[IDX_DWF_QDW(ic,3,0,Ns,site)] = r3;
        } else {
            vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
            vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
            vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
            vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        }
    }
}

// Float-float (FP32-only DD) LU_inv kernel. Coefficients are (hi,lo) pairs in
// real_t (= float in float build); all scalar arithmetic uses dw_add/dw_mul
// /dw_scal/QDW_SCAL_DD which are pure FP32 fma chains. No FP64 ALU is used.
// Only meaningful in the float build (NS lo arrays are populated by
// initDomainwallConstantMemoryFF when extended_precision is on).
__global__ void mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel_ff(
    real4* vp, real4* wp, int Ns, int Nst, int Nst_pad,
    real_t alpha_h, real_t alpha_l)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_h     = ConstantMemoryTraits<real_t>::e();
    const real_t* e_l     = ConstantMemoryTraits<real_t>::e_lo();
    const real_t* f_h     = ConstantMemoryTraits<real_t>::f();
    const real_t* f_l     = ConstantMemoryTraits<real_t>::f_lo();
    const real_t* dpinv_h = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dpinv_l = ConstantMemoryTraits<real_t>::dpinv_lo();
    const real_t* dm_h    = ConstantMemoryTraits<real_t>::dm();
    const real_t* dm_l    = ConstantMemoryTraits<real_t>::dm_lo();

    real4 vt0, vt1, vt2, vt3;
    real4 xt0, xt1, xt2, xt3;
    real4 yt0, yt1, yt2, yt3;
    real4 tmp, sum, diff;

    // === Forward sweep: L^{-1} ===
    vt0 = wp[IDX_DWF_QDW(ic,0,0,Ns,site)];
    vt1 = wp[IDX_DWF_QDW(ic,1,0,Ns,site)];
    vt2 = wp[IDX_DWF_QDW(ic,2,0,Ns,site)];
    vt3 = wp[IDX_DWF_QDW(ic,3,0,Ns,site)];
    vp[IDX_DWF_QDW(ic,0,0,Ns,site)] = vt0;
    vp[IDX_DWF_QDW(ic,1,0,Ns,site)] = vt1;
    vp[IDX_DWF_QDW(ic,2,0,Ns,site)] = vt2;
    vp[IDX_DWF_QDW(ic,3,0,Ns,site)] = vt3;
    QDW_SCAL_DD(yt0, e_h[0], e_l[0], vt0);
    QDW_SCAL_DD(yt1, e_h[0], e_l[0], vt1);
    QDW_SCAL_DD(yt2, e_h[0], e_l[0], vt2);
    QDW_SCAL_DD(yt3, e_h[0], e_l[0], vt3);

    for (int is = 1; is < Ns-1; ++is) {
        // a = 0.5 * dm[is] * dpinv[is-1]  — all DD, FP32-only
        real_t a_h, a_l;
        dw_mul(dm_h[is], dm_l[is], dpinv_h[is-1], dpinv_l[is-1], a_h, a_l);
        dw_scal(real_t(0.5), a_h, a_l, a_h, a_l);
        const real_t eis_h = e_h[is], eis_l = e_l[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_ADD(sum, xt0, xt2); QDW_SCAL_DD(tmp, a_h, a_l, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL_DD(tmp, a_h, a_l, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SCAL_DD(tmp, eis_h, eis_l, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    {
        const int is = Ns-1;
        real_t al_h, al_l;
        dw_mul(dm_h[is], dm_l[is], dpinv_h[is-1], dpinv_l[is-1], al_h, al_l);
        dw_scal(real_t(0.5), al_h, al_l, al_h, al_l);
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_ADD(sum, xt0, xt2); QDW_SCAL_DD(tmp, al_h, al_l, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL_DD(tmp, al_h, al_l, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        QDW_SUB(diff, yt0, yt2); QDW_SCAL(tmp, real_t(0.5), diff);
        QDW_SUB(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_SUB(diff, yt1, yt3); QDW_SCAL(tmp, real_t(0.5), diff);
        QDW_SUB(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
    }

    // === Backward sweep: D^{-1} * U^{-1} ===
    {
        const int is = Ns-1;
        const real_t aa_h = dpinv_h[Ns-1], aa_l = dpinv_l[Ns-1];
        // f1 = 0.5 * (1 + alpha), f2 = 0.5 * (-1 + alpha)  — DD
        real_t one_plus_a_h, one_plus_a_l;
        dw_add(real_t(1.0), real_t(0.0), alpha_h, alpha_l, one_plus_a_h, one_plus_a_l);
        real_t f1_h, f1_l;
        dw_scal(real_t(0.5), one_plus_a_h, one_plus_a_l, f1_h, f1_l);
        real_t mone_plus_a_h, mone_plus_a_l;
        dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, mone_plus_a_h, mone_plus_a_l);
        real_t f2_h, f2_l;
        dw_scal(real_t(0.5), mone_plus_a_h, mone_plus_a_l, f2_h, f2_l);
        real4 p0, p1, p2, p3, t;
        p0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        p1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        p2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        p3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SCAL_DD(vt0,f1_h,f1_l,p0); QDW_SCAL_DD(t,f2_h,f2_l,p2); QDW_ADD(vt0,vt0,t); QDW_SCAL_DD(vt0,aa_h,aa_l,vt0);
        QDW_SCAL_DD(vt1,f1_h,f1_l,p1); QDW_SCAL_DD(t,f2_h,f2_l,p3); QDW_ADD(vt1,vt1,t); QDW_SCAL_DD(vt1,aa_h,aa_l,vt1);
        QDW_SCAL_DD(vt2,f1_h,f1_l,p2); QDW_SCAL_DD(t,f2_h,f2_l,p0); QDW_ADD(vt2,vt2,t); QDW_SCAL_DD(vt2,aa_h,aa_l,vt2);
        QDW_SCAL_DD(vt3,f1_h,f1_l,p3); QDW_SCAL_DD(t,f2_h,f2_l,p1); QDW_ADD(vt3,vt3,t); QDW_SCAL_DD(vt3,aa_h,aa_l,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_ADD(sum, vt0, vt2); QDW_SCAL(yt0, real_t(0.5), sum); yt2 = yt0;
        QDW_ADD(sum, vt1, vt3); QDW_SCAL(yt1, real_t(0.5), sum); yt3 = yt1;
    }

    for (int is = Ns-2; is >= 0; --is) {
        real_t a_h, a_l;
        dw_scal(real_t(0.5), dm_h[is], dm_l[is], a_h, a_l);
        real_t fis_h = f_h[is];
        real_t fis_l = f_l[is];
        real_t aais_h = dpinv_h[is], aais_l = dpinv_l[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SUB(diff, xt0, xt2); QDW_SCAL_DD(tmp, a_h, a_l, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL_DD(tmp, a_h, a_l, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, yt3); QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL_DD(vt0, aais_h, aais_l, vt0);
        QDW_SCAL_DD(vt1, aais_h, aais_l, vt1);
        QDW_SCAL_DD(vt2, aais_h, aais_l, vt2);
        QDW_SCAL_DD(vt3, aais_h, aais_l, vt3);
        if (is == 0) {
            real_t one_plus_a_h, one_plus_a_l;
            dw_add(real_t(1.0), real_t(0.0), alpha_h, alpha_l, one_plus_a_h, one_plus_a_l);
            real_t ff1_h, ff1_l;
            dw_scal(real_t(0.5), one_plus_a_h, one_plus_a_l, ff1_h, ff1_l);
            real_t one_minus_a_h, one_minus_a_l;
            dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, one_minus_a_h, one_minus_a_l);
            real_t ff2_h, ff2_l;
            dw_scal(real_t(0.5), one_minus_a_h, one_minus_a_l, ff2_h, ff2_l);
            real4 r0, r1, r2, r3, t;
            QDW_SCAL_DD(r0,ff1_h,ff1_l,vt0); QDW_SCAL_DD(t,ff2_h,ff2_l,vt2); QDW_ADD(r0,r0,t);
            QDW_SCAL_DD(r1,ff1_h,ff1_l,vt1); QDW_SCAL_DD(t,ff2_h,ff2_l,vt3); QDW_ADD(r1,r1,t);
            QDW_SCAL_DD(r2,ff1_h,ff1_l,vt2); QDW_SCAL_DD(t,ff2_h,ff2_l,vt0); QDW_ADD(r2,r2,t);
            QDW_SCAL_DD(r3,ff1_h,ff1_l,vt3); QDW_SCAL_DD(t,ff2_h,ff2_l,vt1); QDW_ADD(r3,r3,t);
            vp[IDX_DWF_QDW(ic,0,0,Ns,site)] = r0;
            vp[IDX_DWF_QDW(ic,1,0,Ns,site)] = r1;
            vp[IDX_DWF_QDW(ic,2,0,Ns,site)] = r2;
            vp[IDX_DWF_QDW(ic,3,0,Ns,site)] = r3;
        } else {
            vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
            vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
            vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
            vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        }
    }
}


// ===================================================================
// LUdag_inv QDW kernel (adjoint of LU_inv)
// ===================================================================
__global__ void mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel_fp(
    real4* vp, real4* wp, int Ns, int Nst, int Nst_pad, real_t alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
    const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
    const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm();

    real4 vt0, vt1, vt2, vt3;
    real4 xt0, xt1, xt2, xt3;
    real4 yt0, yt1, yt2, yt3;
    real4 tmp, sum, diff;

    // === Forward sweep: Udag^{-1} * D^{-1} ===

    // is = 0: alpha rotation + scaling by dpinv[0]
    {
        const int is     = 0;
        const real_t a0  = dpinv_con[0];
        const real_t f1  = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t f2  = real_t(0.5) * (real_t(1.0) - alpha);
        real4 w0, w1, w2, w3, t;
        w0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        w1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        w2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        w3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SCAL(vt0,f1,w0); QDW_SCAL(t,f2,w2); QDW_ADD(vt0,vt0,t); QDW_SCAL(vt0,a0,vt0);
        QDW_SCAL(vt1,f1,w1); QDW_SCAL(t,f2,w3); QDW_ADD(vt1,vt1,t); QDW_SCAL(vt1,a0,vt1);
        QDW_SCAL(vt2,f1,w2); QDW_SCAL(t,f2,w0); QDW_ADD(vt2,vt2,t); QDW_SCAL(vt2,a0,vt2);
        QDW_SCAL(vt3,f1,w3); QDW_SCAL(t,f2,w1); QDW_ADD(vt3,vt3,t); QDW_SCAL(vt3,a0,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        // yt = f[0] * vt
        QDW_SCAL(yt0, f_con[0], vt0);
        QDW_SCAL(yt1, f_con[0], vt1);
        QDW_SCAL(yt2, f_con[0], vt2);
        QDW_SCAL(yt3, f_con[0], vt3);
    }

    // is = 1 .. Ns-2
    for (int is = 1; is < Ns-1; ++is) {
        const real_t a    = real_t(0.5) * dm_con[is-1];
        const real_t aa   = dpinv_con[is];
        const real_t f_is = f_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        // vt[0] += a*(xt[0]-xt[2])  (antisymmetric)
        QDW_SUB(diff, xt0, xt2); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        // scale by dpinv[is]
        QDW_SCAL(vt0, aa, vt0);
        QDW_SCAL(vt1, aa, vt1);
        QDW_SCAL(vt2, aa, vt2);
        QDW_SCAL(vt3, aa, vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SCAL(tmp, f_is, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL(tmp, f_is, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL(tmp, f_is, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL(tmp, f_is, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    // is = Ns-1
    {
        const int is     = Ns-1;
        const real_t a_l = real_t(0.5) * dm_con[is-1];
        const real_t aa  = dpinv_con[is];
        const real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t ff2 = real_t(0.5) * (real_t(-1.0) + alpha);
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        // vt[0] += a*(xt[0]-xt[2])  (antisymmetric)
        QDW_SUB(diff, xt0, xt2); QDW_SCAL(tmp, a_l, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL(tmp, a_l, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        // -0.5*(yt[s0]+yt[s2]) symmetric yt correction
        QDW_ADD(sum, yt0, yt2); QDW_SCAL(tmp, real_t(0.5), sum);
        QDW_SUB(vt0, vt0, tmp); QDW_SUB(vt2, vt2, tmp);
        QDW_ADD(sum, yt1, yt3); QDW_SCAL(tmp, real_t(0.5), sum);
        QDW_SUB(vt1, vt1, tmp); QDW_SUB(vt3, vt3, tmp);
        // scale by dpinv[Ns-1]
        QDW_SCAL(vt0, aa, vt0);
        QDW_SCAL(vt1, aa, vt1);
        QDW_SCAL(vt2, aa, vt2);
        QDW_SCAL(vt3, aa, vt3);
        // alpha rotation at top boundary
        real4 r0, r1, r2, r3, t;
        QDW_SCAL(r0,ff1,vt0); QDW_SCAL(t,ff2,vt2); QDW_ADD(r0,r0,t);
        QDW_SCAL(r1,ff1,vt1); QDW_SCAL(t,ff2,vt3); QDW_ADD(r1,r1,t);
        QDW_SCAL(r2,ff1,vt2); QDW_SCAL(t,ff2,vt0); QDW_ADD(r2,r2,t);
        QDW_SCAL(r3,ff1,vt3); QDW_SCAL(t,ff2,vt1); QDW_ADD(r3,r3,t);
        vt0=r0; vt1=r1; vt2=r2; vt3=r3;
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        // yt[id] = 0.5*(vt[id]-vt[id_mirror])  (antisymmetric)
        QDW_SUB(diff, vt0, vt2); QDW_SCAL(yt0,  real_t(0.5), diff);
        QDW_NEG(yt2, yt0);
        QDW_SUB(diff, vt1, vt3); QDW_SCAL(yt1,  real_t(0.5), diff);
        QDW_NEG(yt3, yt1);
    }
    // Udag^{-1} * D^{-1} completed

    // === Backward sweep: Ldag^{-1} ===
    for (int is = Ns-2; is >= 0; --is) {
        const real_t a    = real_t(0.5) * dm_con[is+1] * dpinv_con[is];
        const real_t e_is = e_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        // vt[0] += a*(xt[0]+xt[2]) - e_is*yt[0]  (symmetric hopping)
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        QDW_SCAL(tmp, e_is, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL(tmp, e_is, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL(tmp, e_is, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL(tmp, e_is, yt3); QDW_SUB(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
    }
}


// ===================================================================
// LUdag_inv: FP32-only DD variant (float build, extended_precision=true).
// Coefficients e/f/dpinv/dm and alpha are float pairs (h,l); all arithmetic
// uses dw_add/dw_mul/dw_scal/QDW_SCAL_DD — pure FP32 fma chain.
// ===================================================================
__global__ void mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel_ff(
    real4* vp, real4* wp, int Ns, int Nst, int Nst_pad,
    real_t alpha_h, real_t alpha_l)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_h     = ConstantMemoryTraits<real_t>::e();
    const real_t* e_l     = ConstantMemoryTraits<real_t>::e_lo();
    const real_t* f_h     = ConstantMemoryTraits<real_t>::f();
    const real_t* f_l     = ConstantMemoryTraits<real_t>::f_lo();
    const real_t* dpinv_h = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dpinv_l = ConstantMemoryTraits<real_t>::dpinv_lo();
    const real_t* dm_h    = ConstantMemoryTraits<real_t>::dm();
    const real_t* dm_l    = ConstantMemoryTraits<real_t>::dm_lo();

    real4 vt0, vt1, vt2, vt3;
    real4 xt0, xt1, xt2, xt3;
    real4 yt0, yt1, yt2, yt3;
    real4 tmp, sum, diff;

    // === Forward sweep: Udag^{-1} * D^{-1} ===
    // is = 0: alpha rotation + scaling by dpinv[0]
    {
        const int is = 0;
        const real_t a0_h = dpinv_h[0], a0_l = dpinv_l[0];
        real_t one_plus_a_h, one_plus_a_l;
        dw_add(real_t(1.0), real_t(0.0), alpha_h, alpha_l, one_plus_a_h, one_plus_a_l);
        real_t one_minus_a_h, one_minus_a_l;
        dw_add(real_t(1.0), real_t(0.0), -alpha_h, -alpha_l, one_minus_a_h, one_minus_a_l);
        real_t f1_h, f1_l, f2_h, f2_l;
        dw_scal(real_t(0.5), one_plus_a_h, one_plus_a_l, f1_h, f1_l);
        dw_scal(real_t(0.5), one_minus_a_h, one_minus_a_l, f2_h, f2_l);
        real4 w0, w1, w2, w3, t;
        w0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        w1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        w2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        w3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SCAL_DD(vt0,f1_h,f1_l,w0); QDW_SCAL_DD(t,f2_h,f2_l,w2); QDW_ADD(vt0,vt0,t); QDW_SCAL_DD(vt0,a0_h,a0_l,vt0);
        QDW_SCAL_DD(vt1,f1_h,f1_l,w1); QDW_SCAL_DD(t,f2_h,f2_l,w3); QDW_ADD(vt1,vt1,t); QDW_SCAL_DD(vt1,a0_h,a0_l,vt1);
        QDW_SCAL_DD(vt2,f1_h,f1_l,w2); QDW_SCAL_DD(t,f2_h,f2_l,w0); QDW_ADD(vt2,vt2,t); QDW_SCAL_DD(vt2,a0_h,a0_l,vt2);
        QDW_SCAL_DD(vt3,f1_h,f1_l,w3); QDW_SCAL_DD(t,f2_h,f2_l,w1); QDW_ADD(vt3,vt3,t); QDW_SCAL_DD(vt3,a0_h,a0_l,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SCAL_DD(yt0, f_h[0], f_l[0], vt0);
        QDW_SCAL_DD(yt1, f_h[0], f_l[0], vt1);
        QDW_SCAL_DD(yt2, f_h[0], f_l[0], vt2);
        QDW_SCAL_DD(yt3, f_h[0], f_l[0], vt3);
    }

    for (int is = 1; is < Ns-1; ++is) {
        real_t a_h, a_l;
        dw_scal(real_t(0.5), dm_h[is-1], dm_l[is-1], a_h, a_l);
        const real_t aa_h = dpinv_h[is], aa_l = dpinv_l[is];
        const real_t fis_h = f_h[is], fis_l = f_l[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SUB(diff, xt0, xt2); QDW_SCAL_DD(tmp, a_h, a_l, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL_DD(tmp, a_h, a_l, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL_DD(vt0, aa_h, aa_l, vt0);
        QDW_SCAL_DD(vt1, aa_h, aa_l, vt1);
        QDW_SCAL_DD(vt2, aa_h, aa_l, vt2);
        QDW_SCAL_DD(vt3, aa_h, aa_l, vt3);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SCAL_DD(tmp, fis_h, fis_l, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL_DD(tmp, fis_h, fis_l, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    {
        const int is = Ns-1;
        real_t al_h, al_l;
        dw_scal(real_t(0.5), dm_h[is-1], dm_l[is-1], al_h, al_l);
        const real_t aa_h = dpinv_h[is], aa_l = dpinv_l[is];
        real_t one_plus_a_h, one_plus_a_l, m_one_plus_a_h, m_one_plus_a_l;
        dw_add(real_t(1.0),  real_t(0.0), alpha_h, alpha_l, one_plus_a_h, one_plus_a_l);
        dw_add(real_t(-1.0), real_t(0.0), alpha_h, alpha_l, m_one_plus_a_h, m_one_plus_a_l);
        real_t ff1_h, ff1_l, ff2_h, ff2_l;
        dw_scal(real_t(0.5), one_plus_a_h,   one_plus_a_l,   ff1_h, ff1_l);
        dw_scal(real_t(0.5), m_one_plus_a_h, m_one_plus_a_l, ff2_h, ff2_l);

        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_SUB(diff, xt0, xt2); QDW_SCAL_DD(tmp, al_h, al_l, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL_DD(tmp, al_h, al_l, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        QDW_ADD(sum, yt0, yt2); QDW_SCAL(tmp, real_t(0.5), sum);
        QDW_SUB(vt0, vt0, tmp); QDW_SUB(vt2, vt2, tmp);
        QDW_ADD(sum, yt1, yt3); QDW_SCAL(tmp, real_t(0.5), sum);
        QDW_SUB(vt1, vt1, tmp); QDW_SUB(vt3, vt3, tmp);
        QDW_SCAL_DD(vt0, aa_h, aa_l, vt0);
        QDW_SCAL_DD(vt1, aa_h, aa_l, vt1);
        QDW_SCAL_DD(vt2, aa_h, aa_l, vt2);
        QDW_SCAL_DD(vt3, aa_h, aa_l, vt3);
        real4 r0, r1, r2, r3, t;
        QDW_SCAL_DD(r0,ff1_h,ff1_l,vt0); QDW_SCAL_DD(t,ff2_h,ff2_l,vt2); QDW_ADD(r0,r0,t);
        QDW_SCAL_DD(r1,ff1_h,ff1_l,vt1); QDW_SCAL_DD(t,ff2_h,ff2_l,vt3); QDW_ADD(r1,r1,t);
        QDW_SCAL_DD(r2,ff1_h,ff1_l,vt2); QDW_SCAL_DD(t,ff2_h,ff2_l,vt0); QDW_ADD(r2,r2,t);
        QDW_SCAL_DD(r3,ff1_h,ff1_l,vt3); QDW_SCAL_DD(t,ff2_h,ff2_l,vt1); QDW_ADD(r3,r3,t);
        vt0=r0; vt1=r1; vt2=r2; vt3=r3;
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
        QDW_SUB(diff, vt0, vt2); QDW_SCAL(yt0,  real_t(0.5), diff);
        QDW_NEG(yt2, yt0);
        QDW_SUB(diff, vt1, vt3); QDW_SCAL(yt1,  real_t(0.5), diff);
        QDW_NEG(yt3, yt1);
    }

    // === Backward sweep: Ldag^{-1} ===
    for (int is = Ns-2; is >= 0; --is) {
        real_t a_h, a_l;
        dw_mul(dm_h[is+1], dm_l[is+1], dpinv_h[is], dpinv_l[is], a_h, a_l);
        dw_scal(real_t(0.5), a_h, a_l, a_h, a_l);
        const real_t eis_h = e_h[is], eis_l = e_l[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Ns,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Ns,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Ns,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Ns,site)];
        QDW_ADD(sum, xt0, xt2); QDW_SCAL_DD(tmp, a_h, a_l, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL_DD(tmp, a_h, a_l, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL_DD(tmp, eis_h, eis_l, yt3); QDW_SUB(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Ns,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Ns,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Ns,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Ns,site)] = vt3;
    }
}


// ===================================================================
// Launcher: LU_inv QDW
// ===================================================================
void mult_domainwall_5din_ee_LUinv_dirac_qdw(
    real_t* vp, real_t* wp, int Ns, int* Nsize, double alpha, bool ext)
{
    const int Nx  = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    real4* vp_dev = (real4*)dev_ptr(vp);
    real4* wp_dev = (real4*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;

    // ext=true is only meaningful for the float build (where it raises the
    // coefficient precision to ~48-bit via float-float DD math, FP32-only).
    // For real_t=double the FF path is silently disabled (no benefit).
    if (ext && sizeof(real_t) == sizeof(float)) {
        const real_t a_h = (real_t)alpha;
        const real_t a_l = (real_t)(alpha - (double)a_h);
        mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
            vp_dev, wp_dev, Ns, Nst, Nst_pad, a_h, a_l);
    } else {
        mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
            vp_dev, wp_dev, Ns, Nst, Nst_pad, (real_t)alpha);
    }

    CHECK(cudaDeviceSynchronize());
}


// ===================================================================
// Launcher: LUdag_inv QDW
// ===================================================================
void mult_domainwall_5din_ee_LUdaginv_dirac_qdw(
    real_t* vp, real_t* wp, int Ns, int* Nsize, double alpha, bool ext)
{
    const int Nx  = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    real4* vp_dev = (real4*)dev_ptr(vp);
    real4* wp_dev = (real4*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;

    if (ext && sizeof(real_t) == sizeof(float)) {
        const real_t a_h = (real_t)alpha;
        const real_t a_l = (real_t)(alpha - (double)a_h);
        mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel_ff<<<gridSize, blockSize>>>(
            vp_dev, wp_dev, Ns, Nst, Nst_pad, a_h, a_l);
    } else {
        mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel_fp<<<gridSize, blockSize>>>(
            vp_dev, wp_dev, Ns, Nst, Nst_pad, (real_t)alpha);
    }

    CHECK(cudaDeviceSynchronize());
}

#undef IDX_DWF_QDW

#endif // MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED
