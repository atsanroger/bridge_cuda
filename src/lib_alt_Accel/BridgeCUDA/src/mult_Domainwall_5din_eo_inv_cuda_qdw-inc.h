/*!
      @file    mult_Domainwall_5din_eo_inv_cuda_qdw-inc.h
      @brief   QDW LU/LUdag inverse kernels for 5D even-odd domain-wall fermion.
               QDW is double-precision only.
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED

#include "../inline/constant_memory_inline.h"

// IDX_DWF_QDW: earlier QDW EO include undefines this macro — redefine here.
#define IDX_DWF_QDW(ic, id, is5, Nst_pad_, site_) \
    IDX2_QDW((ic), (id), (is5) * (Nst_pad_) + (site_))

// ===================================================================
// LU_inv QDW kernel
// Thread → (ic, site): NC * Nst threads.
// Each thread performs forward (L^{-1}) and backward (D^{-1}*U^{-1})
// sweep over s-slices for one (ic, site) pair with QDW arithmetic.
// ===================================================================
__global__ void mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel(
    double4* vp, double4* wp, int Ns, int Nst, int Nst_pad, double alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const double* e_con     = ConstantMemoryTraits<double>::e();
    const double* f_con     = ConstantMemoryTraits<double>::f();
    const double* dpinv_con = ConstantMemoryTraits<double>::dpinv();
    const double* dm_con    = ConstantMemoryTraits<double>::dm();

    double4 vt0, vt1, vt2, vt3;
    double4 xt0, xt1, xt2, xt3;
    double4 yt0, yt1, yt2, yt3;
    double4 tmp, sum, diff;

    // === Forward sweep: L^{-1} ===

    // is = 0
    vt0 = wp[IDX_DWF_QDW(ic,0,0,Nst_pad,site)];
    vt1 = wp[IDX_DWF_QDW(ic,1,0,Nst_pad,site)];
    vt2 = wp[IDX_DWF_QDW(ic,2,0,Nst_pad,site)];
    vt3 = wp[IDX_DWF_QDW(ic,3,0,Nst_pad,site)];
    vp[IDX_DWF_QDW(ic,0,0,Nst_pad,site)] = vt0;
    vp[IDX_DWF_QDW(ic,1,0,Nst_pad,site)] = vt1;
    vp[IDX_DWF_QDW(ic,2,0,Nst_pad,site)] = vt2;
    vp[IDX_DWF_QDW(ic,3,0,Nst_pad,site)] = vt3;
    QDW_SCAL(yt0, e_con[0], vt0);
    QDW_SCAL(yt1, e_con[0], vt1);
    QDW_SCAL(yt2, e_con[0], vt2);
    QDW_SCAL(yt3, e_con[0], vt3);

    // is = 1 .. Ns-2
    for (int is = 1; is < Ns-1; ++is) {
        const double a   = 0.5 * dm_con[is] * dpinv_con[is-1];
        const double eis = e_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // vt[s0] += a*(xt[s0]+xt[s2]),  vt[s2] += same  (symmetric)
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        QDW_SCAL(tmp, eis, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL(tmp, eis, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL(tmp, eis, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL(tmp, eis, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    // is = Ns-1
    {
        const int is     = Ns-1;
        const double a_l = 0.5 * dm_con[is] * dpinv_con[is-1];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // a*(xt[s0]+xt[s2]) symmetric correction
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a_l, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a_l, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        // -0.5*(yt[s0]-yt[s2]) antisymmetric yt correction
        QDW_SUB(diff, yt0, yt2); QDW_SCAL(tmp, 0.5, diff);
        QDW_SUB(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_SUB(diff, yt1, yt3); QDW_SCAL(tmp, 0.5, diff);
        QDW_SUB(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
    }
    // L^{-1} completed

    // === Backward sweep: D^{-1} * U^{-1} ===

    // is = Ns-1: alpha rotation + scaling by dpinv[Ns-1]
    {
        const int is = Ns-1;
        const double aa = dpinv_con[Ns-1];
        const double f1 = 0.5 * (1.0 + alpha);
        const double f2 = 0.5 * (-1.0 + alpha);
        double4 p0, p1, p2, p3, t;
        p0 = vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        p1 = vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        p2 = vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        p3 = vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // vt[0] = aa*(f1*p0 + f2*p2)
        QDW_SCAL(vt0,f1,p0); QDW_SCAL(t,f2,p2); QDW_ADD(vt0,vt0,t); QDW_SCAL(vt0,aa,vt0);
        QDW_SCAL(vt1,f1,p1); QDW_SCAL(t,f2,p3); QDW_ADD(vt1,vt1,t); QDW_SCAL(vt1,aa,vt1);
        QDW_SCAL(vt2,f1,p2); QDW_SCAL(t,f2,p0); QDW_ADD(vt2,vt2,t); QDW_SCAL(vt2,aa,vt2);
        QDW_SCAL(vt3,f1,p3); QDW_SCAL(t,f2,p1); QDW_ADD(vt3,vt3,t); QDW_SCAL(vt3,aa,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        // yt[id] = 0.5*(vt[id]+vt[id_mirror])  — both mirror values equal
        QDW_ADD(sum, vt0, vt2); QDW_SCAL(yt0, 0.5, sum); yt2 = yt0;
        QDW_ADD(sum, vt1, vt3); QDW_SCAL(yt1, 0.5, sum); yt3 = yt1;
    }

    // is = Ns-2 .. 0
    for (int is = Ns-2; is >= 0; --is) {
        const double a    = 0.5 * dm_con[is];
        const double f_is = f_con[is];
        const double aa_is = dpinv_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // vt[0] += a*(xt[0]-xt[2])  (antisymmetric)
        QDW_SUB(diff, xt0, xt2); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL(tmp, a, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        // subtract f_is * yt
        QDW_SCAL(tmp, f_is, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL(tmp, f_is, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL(tmp, f_is, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL(tmp, f_is, yt3); QDW_SUB(vt3, vt3, tmp);
        // scale by dpinv[is]
        QDW_SCAL(vt0, aa_is, vt0);
        QDW_SCAL(vt1, aa_is, vt1);
        QDW_SCAL(vt2, aa_is, vt2);
        QDW_SCAL(vt3, aa_is, vt3);
        if (is == 0) {
            // alpha rotation at the boundary slice
            const double ff1 = 0.5 * (1.0 + alpha);
            const double ff2 = 0.5 * (1.0 - alpha);
            double4 r0, r1, r2, r3, t;
            QDW_SCAL(r0,ff1,vt0); QDW_SCAL(t,ff2,vt2); QDW_ADD(r0,r0,t);
            QDW_SCAL(r1,ff1,vt1); QDW_SCAL(t,ff2,vt3); QDW_ADD(r1,r1,t);
            QDW_SCAL(r2,ff1,vt2); QDW_SCAL(t,ff2,vt0); QDW_ADD(r2,r2,t);
            QDW_SCAL(r3,ff1,vt3); QDW_SCAL(t,ff2,vt1); QDW_ADD(r3,r3,t);
            vp[IDX_DWF_QDW(ic,0,0,Nst_pad,site)] = r0;
            vp[IDX_DWF_QDW(ic,1,0,Nst_pad,site)] = r1;
            vp[IDX_DWF_QDW(ic,2,0,Nst_pad,site)] = r2;
            vp[IDX_DWF_QDW(ic,3,0,Nst_pad,site)] = r3;
        } else {
            vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
            vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
            vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
            vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        }
    }
}


// ===================================================================
// LUdag_inv QDW kernel (adjoint of LU_inv)
// ===================================================================
__global__ void mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel(
    double4* vp, double4* wp, int Ns, int Nst, int Nst_pad, double alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const double* e_con     = ConstantMemoryTraits<double>::e();
    const double* f_con     = ConstantMemoryTraits<double>::f();
    const double* dpinv_con = ConstantMemoryTraits<double>::dpinv();
    const double* dm_con    = ConstantMemoryTraits<double>::dm();

    double4 vt0, vt1, vt2, vt3;
    double4 xt0, xt1, xt2, xt3;
    double4 yt0, yt1, yt2, yt3;
    double4 tmp, sum, diff;

    // === Forward sweep: Udag^{-1} * D^{-1} ===

    // is = 0: alpha rotation + scaling by dpinv[0]
    {
        const int is     = 0;
        const double a0  = dpinv_con[0];
        const double f1  = 0.5 * (1.0 + alpha);
        const double f2  = 0.5 * (1.0 - alpha);
        double4 w0, w1, w2, w3, t;
        w0 = wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        w1 = wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        w2 = wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        w3 = wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        QDW_SCAL(vt0,f1,w0); QDW_SCAL(t,f2,w2); QDW_ADD(vt0,vt0,t); QDW_SCAL(vt0,a0,vt0);
        QDW_SCAL(vt1,f1,w1); QDW_SCAL(t,f2,w3); QDW_ADD(vt1,vt1,t); QDW_SCAL(vt1,a0,vt1);
        QDW_SCAL(vt2,f1,w2); QDW_SCAL(t,f2,w0); QDW_ADD(vt2,vt2,t); QDW_SCAL(vt2,a0,vt2);
        QDW_SCAL(vt3,f1,w3); QDW_SCAL(t,f2,w1); QDW_ADD(vt3,vt3,t); QDW_SCAL(vt3,a0,vt3);
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        // yt = f[0] * vt
        QDW_SCAL(yt0, f_con[0], vt0);
        QDW_SCAL(yt1, f_con[0], vt1);
        QDW_SCAL(yt2, f_con[0], vt2);
        QDW_SCAL(yt3, f_con[0], vt3);
    }

    // is = 1 .. Ns-2
    for (int is = 1; is < Ns-1; ++is) {
        const double a    = 0.5 * dm_con[is-1];
        const double aa   = dpinv_con[is];
        const double f_is = f_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
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
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        QDW_SCAL(tmp, f_is, vt0); QDW_ADD(yt0, yt0, tmp);
        QDW_SCAL(tmp, f_is, vt1); QDW_ADD(yt1, yt1, tmp);
        QDW_SCAL(tmp, f_is, vt2); QDW_ADD(yt2, yt2, tmp);
        QDW_SCAL(tmp, f_is, vt3); QDW_ADD(yt3, yt3, tmp);
    }

    // is = Ns-1
    {
        const int is     = Ns-1;
        const double a_l = 0.5 * dm_con[is-1];
        const double aa  = dpinv_con[is];
        const double ff1 = 0.5 * (1.0 + alpha);
        const double ff2 = 0.5 * (-1.0 + alpha);
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = wp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = wp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = wp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = wp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // vt[0] += a*(xt[0]-xt[2])  (antisymmetric)
        QDW_SUB(diff, xt0, xt2); QDW_SCAL(tmp, a_l, diff);
        QDW_ADD(vt0, vt0, tmp);  QDW_SUB(vt2, vt2, tmp);
        QDW_SUB(diff, xt1, xt3); QDW_SCAL(tmp, a_l, diff);
        QDW_ADD(vt1, vt1, tmp);  QDW_SUB(vt3, vt3, tmp);
        // -0.5*(yt[s0]+yt[s2]) symmetric yt correction
        QDW_ADD(sum, yt0, yt2); QDW_SCAL(tmp, 0.5, sum);
        QDW_SUB(vt0, vt0, tmp); QDW_SUB(vt2, vt2, tmp);
        QDW_ADD(sum, yt1, yt3); QDW_SCAL(tmp, 0.5, sum);
        QDW_SUB(vt1, vt1, tmp); QDW_SUB(vt3, vt3, tmp);
        // scale by dpinv[Ns-1]
        QDW_SCAL(vt0, aa, vt0);
        QDW_SCAL(vt1, aa, vt1);
        QDW_SCAL(vt2, aa, vt2);
        QDW_SCAL(vt3, aa, vt3);
        // alpha rotation at top boundary
        double4 r0, r1, r2, r3, t;
        QDW_SCAL(r0,ff1,vt0); QDW_SCAL(t,ff2,vt2); QDW_ADD(r0,r0,t);
        QDW_SCAL(r1,ff1,vt1); QDW_SCAL(t,ff2,vt3); QDW_ADD(r1,r1,t);
        QDW_SCAL(r2,ff1,vt2); QDW_SCAL(t,ff2,vt0); QDW_ADD(r2,r2,t);
        QDW_SCAL(r3,ff1,vt3); QDW_SCAL(t,ff2,vt1); QDW_ADD(r3,r3,t);
        vt0=r0; vt1=r1; vt2=r2; vt3=r3;
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
        // yt[id] = 0.5*(vt[id]-vt[id_mirror])  (antisymmetric)
        QDW_SUB(diff, vt0, vt2); QDW_SCAL(yt0,  0.5, diff);
        QDW_NEG(yt2, yt0);
        QDW_SUB(diff, vt1, vt3); QDW_SCAL(yt1,  0.5, diff);
        QDW_NEG(yt3, yt1);
    }
    // Udag^{-1} * D^{-1} completed

    // === Backward sweep: Ldag^{-1} ===
    for (int is = Ns-2; is >= 0; --is) {
        const double a    = 0.5 * dm_con[is+1] * dpinv_con[is];
        const double e_is = e_con[is];
        xt0 = vt0; xt1 = vt1; xt2 = vt2; xt3 = vt3;
        vt0 = vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)];
        vt1 = vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)];
        vt2 = vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)];
        vt3 = vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)];
        // vt[0] += a*(xt[0]+xt[2]) - e_is*yt[0]  (symmetric hopping)
        QDW_ADD(sum, xt0, xt2); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt0, vt0, tmp); QDW_ADD(vt2, vt2, tmp);
        QDW_ADD(sum, xt1, xt3); QDW_SCAL(tmp, a, sum);
        QDW_ADD(vt1, vt1, tmp); QDW_ADD(vt3, vt3, tmp);
        QDW_SCAL(tmp, e_is, yt0); QDW_SUB(vt0, vt0, tmp);
        QDW_SCAL(tmp, e_is, yt1); QDW_SUB(vt1, vt1, tmp);
        QDW_SCAL(tmp, e_is, yt2); QDW_SUB(vt2, vt2, tmp);
        QDW_SCAL(tmp, e_is, yt3); QDW_SUB(vt3, vt3, tmp);
        vp[IDX_DWF_QDW(ic,0,is,Nst_pad,site)] = vt0;
        vp[IDX_DWF_QDW(ic,1,is,Nst_pad,site)] = vt1;
        vp[IDX_DWF_QDW(ic,2,is,Nst_pad,site)] = vt2;
        vp[IDX_DWF_QDW(ic,3,is,Nst_pad,site)] = vt3;
    }
}


// ===================================================================
// Launcher: LU_inv QDW
// ===================================================================
void mult_domainwall_5din_ee_LUinv_dirac_qdw(
    double* vp, double* wp, int Ns, int* Nsize, double alpha)
{
    const int Nx  = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    double4* vp_dev = (double4*)dev_ptr(vp);
    double4* wp_dev = (double4*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;

    mult_domainwall_5din_ee_LUinv_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst, Nst_pad, alpha);

    CHECK(cudaDeviceSynchronize());
}


// ===================================================================
// Launcher: LUdag_inv QDW
// ===================================================================
void mult_domainwall_5din_ee_LUdaginv_dirac_qdw(
    double* vp, double* wp, int Ns, int* Nsize, double alpha)
{
    const int Nx  = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    double4* vp_dev = (double4*)dev_ptr(vp);
    double4* wp_dev = (double4*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;

    mult_domainwall_5din_ee_LUdaginv_dirac_qdw_kernel<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst, Nst_pad, alpha);

    CHECK(cudaDeviceSynchronize());
}

#undef IDX_DWF_QDW

#endif // MULT_DOMAINWALL_5DIN_EO_INV_ACC_QDW_INCLUDED
