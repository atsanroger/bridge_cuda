/*!
      @file    mult_Domainwall_5din_eo_inv_cuda_qtw-inc.h
      @brief   QTW (triple-word) LU/LUdag inverse kernels for 5D EO Domain Wall.
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_INV_ACC_QTW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_INV_ACC_QTW_INCLUDED

#include <type_traits>
#include "../inline/constant_memory_inline.h"
// qtw_prof_timer.h included at global scope by the enclosing .cu.

#define IDX_DWF_QTW(ic, id, is5, Ns_, site_) \
    IDX2(NC * ND * (Ns_), (ic) + NC * ((id) + ND * (is5)), (site_))

// Load 6 reals into named scalars prefix##_{rh,ih,rm,im,rl,il}.
#define QTW_LOAD6(prefix, ptr, base) \
    real_t prefix##_rh = (ptr)[(base) + 0]; \
    real_t prefix##_ih = (ptr)[(base) + 1]; \
    real_t prefix##_rm = (ptr)[(base) + 2]; \
    real_t prefix##_im = (ptr)[(base) + 3]; \
    real_t prefix##_rl = (ptr)[(base) + 4]; \
    real_t prefix##_il = (ptr)[(base) + 5]

// Store renormalizes both triples: inner-loop adds (QTW_ADDC/QTW_SUBC) run
// sloppy, so the canonical 3-word form is restored once here, at the end of
// the accumulation chain. Local temps keep the caller's registers reusable.
#define QTW_STORE6(ptr, base, prefix) do { \
    real_t _sr_h = prefix##_rh, _sr_m = prefix##_rm, _sr_l = prefix##_rl; \
    real_t _si_h = prefix##_ih, _si_m = prefix##_im, _si_l = prefix##_il; \
    Renormalize3(_sr_h, _sr_m, _sr_l); \
    Renormalize3(_si_h, _si_m, _si_l); \
    (ptr)[(base) + 0] = _sr_h; \
    (ptr)[(base) + 1] = _si_h; \
    (ptr)[(base) + 2] = _sr_m; \
    (ptr)[(base) + 3] = _si_m; \
    (ptr)[(base) + 4] = _sr_l; \
    (ptr)[(base) + 5] = _si_l; \
} while (0)

#define QTW_DECL(name) \
    real_t name##_rh, name##_ih, name##_rm, name##_im, name##_rl, name##_il

#define QTW_COPY(dst, src) \
    dst##_rh = src##_rh; dst##_ih = src##_ih; \
    dst##_rm = src##_rm; dst##_im = src##_im; \
    dst##_rl = src##_rl; dst##_il = src##_il

// Vector ops on the 6-tuple complex. Adds run sloppy (no Renormalize3);
// QTW_STORE6 restores the canonical form at the end of the chain.
#define QTW_ADDC(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, b##_rh, b##_rm, b##_rl, res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, b##_ih, b##_im, b##_il, res##_ih, res##_im, res##_il); \
} while (0)

#define QTW_SUBC(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, -b##_rh, -b##_rm, -b##_rl, res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, -b##_ih, -b##_im, -b##_il, res##_ih, res##_im, res##_il); \
} while (0)

#define QTW_NEGC(res, a) do { \
    res##_rh = -a##_rh; res##_rm = -a##_rm; res##_rl = -a##_rl; \
    res##_ih = -a##_ih; res##_im = -a##_im; res##_il = -a##_il; \
} while (0)

#define QTW_SCALC(res, u, a) do { \
    tw_scal((u), a##_rh, a##_rm, a##_rl, res##_rh, res##_rm, res##_rl); \
    tw_scal((u), a##_ih, a##_im, a##_il, res##_ih, res##_im, res##_il); \
} while (0)

// TW scalar (ch,cm,cl) * TW complex a -> res. Used for the error-carrying 5D
// coefficients (e, f, dpinv, dm and products of them) so the LU operator runs
// at full triple-word precision instead of single-precision coefficients.
#define QTW_SCALC_TT(res, ch, cm, cl, a) do { \
    tw_mul_sloppy((ch),(cm),(cl), a##_rh, a##_rm, a##_rl, res##_rh, res##_rm, res##_rl); \
    tw_mul_sloppy((ch),(cm),(cl), a##_ih, a##_im, a##_il, res##_ih, res##_im, res##_il); \
} while (0)

// Build a 3-word coefficient from the (hi,mid,lo) constant-memory arrays.
#define QTW_COEF3(name, HI, MID, LO, idx) \
    const real_t name##_h = (HI)[idx]; \
    const real_t name##_m = (MID)[idx]; \
    const real_t name##_l = (LO)[idx]

// =================================================================
// LU_inv QDW kernel — single-T coefficients on TW field.
// Thread → (ic, site). Per thread performs forward (L^{-1}) and backward
// (D^{-1}*U^{-1}) sweep over s-slices for one (ic, site) pair in TW.
// =================================================================
__global__ void mult_domainwall_5din_ee_LUinv_dirac_qtw_kernel_fp(
    real_t* vp, real_t* wp, int Ns, int Nst, int Nst_pad, real_t alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
    const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
    const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm();
    // 3-word (mid/lo) words of the error-carrying coefficients.
    const real_t* e_mid_con     = ConstantMemoryTraits<real_t>::e_mid();
    const real_t* e_lo_con      = ConstantMemoryTraits<real_t>::e_lo_tw();
    const real_t* f_mid_con     = ConstantMemoryTraits<real_t>::f_mid();
    const real_t* f_lo_con      = ConstantMemoryTraits<real_t>::f_lo_tw();
    const real_t* dpinv_mid_con = ConstantMemoryTraits<real_t>::dpinv_mid();
    const real_t* dpinv_lo_con  = ConstantMemoryTraits<real_t>::dpinv_lo_tw();
    const real_t* dm_mid_con    = ConstantMemoryTraits<real_t>::dm_mid();
    const real_t* dm_lo_con     = ConstantMemoryTraits<real_t>::dm_lo_tw();

// a = 0.5 * dm[id_dm] * dpinv[id_dp], 3-word (tw_mul then *0.5 exact per word).
#define QTW_HALF_DMDP(name, id_dm, id_dp) \
    real_t name##_h, name##_m, name##_l; \
    do { tw_mul_sloppy(dm_con[id_dm], dm_mid_con[id_dm], dm_lo_con[id_dm], \
                       dpinv_con[id_dp], dpinv_mid_con[id_dp], dpinv_lo_con[id_dp], \
                       name##_h, name##_m, name##_l); \
         name##_h *= real_t(0.5); name##_m *= real_t(0.5); name##_l *= real_t(0.5); } while (0)
// a = 0.5 * dm[id], 3-word (0.5 exact per word).
#define QTW_HALF_DM(name, id) \
    const real_t name##_h = real_t(0.5) * dm_con[id]; \
    const real_t name##_m = real_t(0.5) * dm_mid_con[id]; \
    const real_t name##_l = real_t(0.5) * dm_lo_con[id]

    QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
    QTW_DECL(xt0); QTW_DECL(xt1); QTW_DECL(xt2); QTW_DECL(xt3);
    QTW_DECL(yt0); QTW_DECL(yt1); QTW_DECL(yt2); QTW_DECL(yt3);
    QTW_DECL(tmp); QTW_DECL(sum); QTW_DECL(diff);

    // ===== Forward sweep: L^{-1} =====
    {
        QTW_LOAD6(_v0, wp, 6 * IDX_DWF_QTW(ic, 0, 0, Ns, site));
        QTW_LOAD6(_v1, wp, 6 * IDX_DWF_QTW(ic, 1, 0, Ns, site));
        QTW_LOAD6(_v2, wp, 6 * IDX_DWF_QTW(ic, 2, 0, Ns, site));
        QTW_LOAD6(_v3, wp, 6 * IDX_DWF_QTW(ic, 3, 0, Ns, site));
        QTW_COPY(vt0, _v0); QTW_COPY(vt1, _v1); QTW_COPY(vt2, _v2); QTW_COPY(vt3, _v3);
    }
    QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, 0, Ns, site), vt0);
    QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, 0, Ns, site), vt1);
    QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, 0, Ns, site), vt2);
    QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, 0, Ns, site), vt3);
    QTW_SCALC_TT(yt0, e_con[0], e_mid_con[0], e_lo_con[0], vt0);
    QTW_SCALC_TT(yt1, e_con[0], e_mid_con[0], e_lo_con[0], vt1);
    QTW_SCALC_TT(yt2, e_con[0], e_mid_con[0], e_lo_con[0], vt2);
    QTW_SCALC_TT(yt3, e_con[0], e_mid_con[0], e_lo_con[0], vt3);

    for (int is = 1; is < Ns-1; ++is) {
        QTW_HALF_DMDP(a, is, is-1);
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_v0, wp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_v1, wp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_v2, wp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_v3, wp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _v0); QTW_COPY(vt1, _v1); QTW_COPY(vt2, _v2); QTW_COPY(vt3, _v3);
        QTW_ADDC(sum, xt0, xt2); QTW_SCALC_TT(tmp, a_h, a_m, a_l, sum);
        QTW_ADDC(vt0, vt0, tmp); QTW_ADDC(vt2, vt2, tmp);
        QTW_ADDC(sum, xt1, xt3); QTW_SCALC_TT(tmp, a_h, a_m, a_l, sum);
        QTW_ADDC(vt1, vt1, tmp); QTW_ADDC(vt3, vt3, tmp);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], vt0); QTW_ADDC(yt0, yt0, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], vt1); QTW_ADDC(yt1, yt1, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], vt2); QTW_ADDC(yt2, yt2, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], vt3); QTW_ADDC(yt3, yt3, tmp);
    }

    {
        const int is     = Ns-1;
        QTW_HALF_DMDP(a_l, is, is-1);
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_v0, wp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_v1, wp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_v2, wp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_v3, wp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _v0); QTW_COPY(vt1, _v1); QTW_COPY(vt2, _v2); QTW_COPY(vt3, _v3);
        QTW_ADDC(sum, xt0, xt2); QTW_SCALC_TT(tmp, a_l_h, a_l_m, a_l_l, sum);
        QTW_ADDC(vt0, vt0, tmp); QTW_ADDC(vt2, vt2, tmp);
        QTW_ADDC(sum, xt1, xt3); QTW_SCALC_TT(tmp, a_l_h, a_l_m, a_l_l, sum);
        QTW_ADDC(vt1, vt1, tmp); QTW_ADDC(vt3, vt3, tmp);
        QTW_SUBC(diff, yt0, yt2); QTW_SCALC(tmp, real_t(0.5), diff);
        QTW_SUBC(vt0, vt0, tmp); QTW_ADDC(vt2, vt2, tmp);
        QTW_SUBC(diff, yt1, yt3); QTW_SCALC(tmp, real_t(0.5), diff);
        QTW_SUBC(vt1, vt1, tmp); QTW_ADDC(vt3, vt3, tmp);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
    }

    // ===== Backward sweep: D^{-1} * U^{-1} =====
    {
        const int is = Ns-1;
        // aa = dpinv[Ns-1] applied as 3-word; f1,f2 are exact at alpha=1.0.
        const real_t aa_h = dpinv_con[Ns-1], aa_m = dpinv_mid_con[Ns-1], aa_l = dpinv_lo_con[Ns-1];
        const real_t f1 = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t f2 = real_t(0.5) * (real_t(-1.0) + alpha);
        QTW_LOAD6(_p0, vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_p1, vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_p2, vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_p3, vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_DECL(t);
        QTW_SCALC(vt0, f1, _p0); QTW_SCALC(t, f2, _p2); QTW_ADDC(vt0, vt0, t); QTW_SCALC_TT(vt0, aa_h, aa_m, aa_l, vt0);
        QTW_SCALC(vt1, f1, _p1); QTW_SCALC(t, f2, _p3); QTW_ADDC(vt1, vt1, t); QTW_SCALC_TT(vt1, aa_h, aa_m, aa_l, vt1);
        QTW_SCALC(vt2, f1, _p2); QTW_SCALC(t, f2, _p0); QTW_ADDC(vt2, vt2, t); QTW_SCALC_TT(vt2, aa_h, aa_m, aa_l, vt2);
        QTW_SCALC(vt3, f1, _p3); QTW_SCALC(t, f2, _p1); QTW_ADDC(vt3, vt3, t); QTW_SCALC_TT(vt3, aa_h, aa_m, aa_l, vt3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        QTW_ADDC(sum, vt0, vt2); QTW_SCALC(yt0, real_t(0.5), sum); QTW_COPY(yt2, yt0);
        QTW_ADDC(sum, vt1, vt3); QTW_SCALC(yt1, real_t(0.5), sum); QTW_COPY(yt3, yt1);
    }

    for (int is = Ns-2; is >= 0; --is) {
        QTW_HALF_DM(a, is);
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_v0, vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_v1, vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_v2, vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_v3, vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _v0); QTW_COPY(vt1, _v1); QTW_COPY(vt2, _v2); QTW_COPY(vt3, _v3);
        QTW_SUBC(diff, xt0, xt2); QTW_SCALC_TT(tmp, a_h, a_m, a_l, diff);
        QTW_ADDC(vt0, vt0, tmp);  QTW_SUBC(vt2, vt2, tmp);
        QTW_SUBC(diff, xt1, xt3); QTW_SCALC_TT(tmp, a_h, a_m, a_l, diff);
        QTW_ADDC(vt1, vt1, tmp);  QTW_SUBC(vt3, vt3, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], yt0); QTW_SUBC(vt0, vt0, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], yt1); QTW_SUBC(vt1, vt1, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], yt2); QTW_SUBC(vt2, vt2, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], yt3); QTW_SUBC(vt3, vt3, tmp);
        QTW_SCALC_TT(vt0, dpinv_con[is], dpinv_mid_con[is], dpinv_lo_con[is], vt0);
        QTW_SCALC_TT(vt1, dpinv_con[is], dpinv_mid_con[is], dpinv_lo_con[is], vt1);
        QTW_SCALC_TT(vt2, dpinv_con[is], dpinv_mid_con[is], dpinv_lo_con[is], vt2);
        QTW_SCALC_TT(vt3, dpinv_con[is], dpinv_mid_con[is], dpinv_lo_con[is], vt3);
        if (is == 0) {
            const real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
            const real_t ff2 = real_t(0.5) * (real_t(1.0) - alpha);
            QTW_DECL(r0); QTW_DECL(r1); QTW_DECL(r2); QTW_DECL(r3); QTW_DECL(t);
            QTW_SCALC(r0, ff1, vt0); QTW_SCALC(t, ff2, vt2); QTW_ADDC(r0, r0, t);
            QTW_SCALC(r1, ff1, vt1); QTW_SCALC(t, ff2, vt3); QTW_ADDC(r1, r1, t);
            QTW_SCALC(r2, ff1, vt2); QTW_SCALC(t, ff2, vt0); QTW_ADDC(r2, r2, t);
            QTW_SCALC(r3, ff1, vt3); QTW_SCALC(t, ff2, vt1); QTW_ADDC(r3, r3, t);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, 0, Ns, site), r0);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, 0, Ns, site), r1);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, 0, Ns, site), r2);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, 0, Ns, site), r3);
        } else {
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
            QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        }
    }
}

// =================================================================
// LUdag_inv QDW kernel — single-T coefficients on TW field.
// =================================================================
__global__ void mult_domainwall_5din_ee_LUdaginv_dirac_qtw_kernel_fp(
    real_t* vp, real_t* wp, int Ns, int Nst, int Nst_pad, real_t alpha)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int site = tid % Nst_pad;
    const int ic   = tid / Nst_pad;
    if (ic >= NC || site >= Nst) return;

    const real_t* e_con     = ConstantMemoryTraits<real_t>::e();
    const real_t* f_con     = ConstantMemoryTraits<real_t>::f();
    const real_t* dpinv_con = ConstantMemoryTraits<real_t>::dpinv();
    const real_t* dm_con    = ConstantMemoryTraits<real_t>::dm();
    const real_t* e_mid_con     = ConstantMemoryTraits<real_t>::e_mid();
    const real_t* e_lo_con      = ConstantMemoryTraits<real_t>::e_lo_tw();
    const real_t* f_mid_con     = ConstantMemoryTraits<real_t>::f_mid();
    const real_t* f_lo_con      = ConstantMemoryTraits<real_t>::f_lo_tw();
    const real_t* dpinv_mid_con = ConstantMemoryTraits<real_t>::dpinv_mid();
    const real_t* dpinv_lo_con  = ConstantMemoryTraits<real_t>::dpinv_lo_tw();
    const real_t* dm_mid_con    = ConstantMemoryTraits<real_t>::dm_mid();
    const real_t* dm_lo_con     = ConstantMemoryTraits<real_t>::dm_lo_tw();
    (void)e_mid_con; (void)e_lo_con;  // e unused in Udag^{-1} but kept for symmetry

    QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
    QTW_DECL(xt0); QTW_DECL(xt1); QTW_DECL(xt2); QTW_DECL(xt3);
    QTW_DECL(yt0); QTW_DECL(yt1); QTW_DECL(yt2); QTW_DECL(yt3);
    QTW_DECL(tmp); QTW_DECL(sum); QTW_DECL(diff);

    // ===== Forward sweep: Udag^{-1} * D^{-1} =====
    {
        const int is = 0;
        const real_t a0_h = dpinv_con[0], a0_m = dpinv_mid_con[0], a0_l = dpinv_lo_con[0];
        const real_t f1 = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t f2 = real_t(0.5) * (real_t(1.0) - alpha);
        QTW_LOAD6(_w0, wp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_w1, wp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_w2, wp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_w3, wp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_DECL(t);
        QTW_SCALC(vt0, f1, _w0); QTW_SCALC(t, f2, _w2); QTW_ADDC(vt0, vt0, t); QTW_SCALC_TT(vt0, a0_h, a0_m, a0_l, vt0);
        QTW_SCALC(vt1, f1, _w1); QTW_SCALC(t, f2, _w3); QTW_ADDC(vt1, vt1, t); QTW_SCALC_TT(vt1, a0_h, a0_m, a0_l, vt1);
        QTW_SCALC(vt2, f1, _w2); QTW_SCALC(t, f2, _w0); QTW_ADDC(vt2, vt2, t); QTW_SCALC_TT(vt2, a0_h, a0_m, a0_l, vt2);
        QTW_SCALC(vt3, f1, _w3); QTW_SCALC(t, f2, _w1); QTW_ADDC(vt3, vt3, t); QTW_SCALC_TT(vt3, a0_h, a0_m, a0_l, vt3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        QTW_SCALC_TT(yt0, f_con[0], f_mid_con[0], f_lo_con[0], vt0);
        QTW_SCALC_TT(yt1, f_con[0], f_mid_con[0], f_lo_con[0], vt1);
        QTW_SCALC_TT(yt2, f_con[0], f_mid_con[0], f_lo_con[0], vt2);
        QTW_SCALC_TT(yt3, f_con[0], f_mid_con[0], f_lo_con[0], vt3);
    }

    for (int is = 1; is < Ns-1; ++is) {
        QTW_HALF_DM(a, is-1);
        const real_t aa_h = dpinv_con[is], aa_m = dpinv_mid_con[is], aa_l = dpinv_lo_con[is];
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_w0, wp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_w1, wp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_w2, wp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_w3, wp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _w0); QTW_COPY(vt1, _w1); QTW_COPY(vt2, _w2); QTW_COPY(vt3, _w3);
        QTW_SUBC(diff, xt0, xt2); QTW_SCALC_TT(tmp, a_h, a_m, a_l, diff);
        QTW_ADDC(vt0, vt0, tmp);  QTW_SUBC(vt2, vt2, tmp);
        QTW_SUBC(diff, xt1, xt3); QTW_SCALC_TT(tmp, a_h, a_m, a_l, diff);
        QTW_ADDC(vt1, vt1, tmp);  QTW_SUBC(vt3, vt3, tmp);
        QTW_SCALC_TT(vt0, aa_h, aa_m, aa_l, vt0);
        QTW_SCALC_TT(vt1, aa_h, aa_m, aa_l, vt1);
        QTW_SCALC_TT(vt2, aa_h, aa_m, aa_l, vt2);
        QTW_SCALC_TT(vt3, aa_h, aa_m, aa_l, vt3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], vt0); QTW_ADDC(yt0, yt0, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], vt1); QTW_ADDC(yt1, yt1, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], vt2); QTW_ADDC(yt2, yt2, tmp);
        QTW_SCALC_TT(tmp, f_con[is], f_mid_con[is], f_lo_con[is], vt3); QTW_ADDC(yt3, yt3, tmp);
    }

    {
        const int is     = Ns-1;
        QTW_HALF_DM(a_l, is-1);
        const real_t aa_h = dpinv_con[is], aa_m = dpinv_mid_con[is], aa_l = dpinv_lo_con[is];
        const real_t ff1 = real_t(0.5) * (real_t(1.0) + alpha);
        const real_t ff2 = real_t(0.5) * (real_t(-1.0) + alpha);
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_w0, wp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_w1, wp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_w2, wp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_w3, wp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _w0); QTW_COPY(vt1, _w1); QTW_COPY(vt2, _w2); QTW_COPY(vt3, _w3);
        QTW_SUBC(diff, xt0, xt2); QTW_SCALC_TT(tmp, a_l_h, a_l_m, a_l_l, diff);
        QTW_ADDC(vt0, vt0, tmp);  QTW_SUBC(vt2, vt2, tmp);
        QTW_SUBC(diff, xt1, xt3); QTW_SCALC_TT(tmp, a_l_h, a_l_m, a_l_l, diff);
        QTW_ADDC(vt1, vt1, tmp);  QTW_SUBC(vt3, vt3, tmp);
        QTW_ADDC(sum, yt0, yt2); QTW_SCALC(tmp, real_t(0.5), sum);
        QTW_SUBC(vt0, vt0, tmp); QTW_SUBC(vt2, vt2, tmp);
        QTW_ADDC(sum, yt1, yt3); QTW_SCALC(tmp, real_t(0.5), sum);
        QTW_SUBC(vt1, vt1, tmp); QTW_SUBC(vt3, vt3, tmp);
        QTW_SCALC_TT(vt0, aa_h, aa_m, aa_l, vt0);
        QTW_SCALC_TT(vt1, aa_h, aa_m, aa_l, vt1);
        QTW_SCALC_TT(vt2, aa_h, aa_m, aa_l, vt2);
        QTW_SCALC_TT(vt3, aa_h, aa_m, aa_l, vt3);
        QTW_DECL(r0); QTW_DECL(r1); QTW_DECL(r2); QTW_DECL(r3); QTW_DECL(t);
        QTW_SCALC(r0, ff1, vt0); QTW_SCALC(t, ff2, vt2); QTW_ADDC(r0, r0, t);
        QTW_SCALC(r1, ff1, vt1); QTW_SCALC(t, ff2, vt3); QTW_ADDC(r1, r1, t);
        QTW_SCALC(r2, ff1, vt2); QTW_SCALC(t, ff2, vt0); QTW_ADDC(r2, r2, t);
        QTW_SCALC(r3, ff1, vt3); QTW_SCALC(t, ff2, vt1); QTW_ADDC(r3, r3, t);
        QTW_COPY(vt0, r0); QTW_COPY(vt1, r1); QTW_COPY(vt2, r2); QTW_COPY(vt3, r3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
        QTW_SUBC(diff, vt0, vt2); QTW_SCALC(yt0, real_t(0.5), diff);
        QTW_NEGC(yt2, yt0);
        QTW_SUBC(diff, vt1, vt3); QTW_SCALC(yt1, real_t(0.5), diff);
        QTW_NEGC(yt3, yt1);
    }

    // ===== Backward sweep: Ldag^{-1} =====
    for (int is = Ns-2; is >= 0; --is) {
        QTW_HALF_DMDP(a, is+1, is);
        QTW_COPY(xt0, vt0); QTW_COPY(xt1, vt1); QTW_COPY(xt2, vt2); QTW_COPY(xt3, vt3);
        QTW_LOAD6(_v0, vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site));
        QTW_LOAD6(_v1, vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site));
        QTW_LOAD6(_v2, vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site));
        QTW_LOAD6(_v3, vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site));
        QTW_COPY(vt0, _v0); QTW_COPY(vt1, _v1); QTW_COPY(vt2, _v2); QTW_COPY(vt3, _v3);
        QTW_ADDC(sum, xt0, xt2); QTW_SCALC_TT(tmp, a_h, a_m, a_l, sum);
        QTW_ADDC(vt0, vt0, tmp); QTW_ADDC(vt2, vt2, tmp);
        QTW_ADDC(sum, xt1, xt3); QTW_SCALC_TT(tmp, a_h, a_m, a_l, sum);
        QTW_ADDC(vt1, vt1, tmp); QTW_ADDC(vt3, vt3, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], yt0); QTW_SUBC(vt0, vt0, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], yt1); QTW_SUBC(vt1, vt1, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], yt2); QTW_SUBC(vt2, vt2, tmp);
        QTW_SCALC_TT(tmp, e_con[is], e_mid_con[is], e_lo_con[is], yt3); QTW_SUBC(vt3, vt3, tmp);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 0, is, Ns, site), vt0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 1, is, Ns, site), vt1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 2, is, Ns, site), vt2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(ic, 3, is, Ns, site), vt3);
    }
}

// Host launchers
void mult_domainwall_5din_ee_LUinv_dirac_qtw(
    real_t* vp, real_t* wp, int Ns, int* Nsize, double alpha, bool ext)
{
    const int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_ee_LUinv_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst, Nst_pad, (real_t)alpha);
    QTW_PROF_END("LUinv_qtw");
    CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_ee_LUdaginv_dirac_qtw(
    real_t* vp, real_t* wp, int Ns, int* Nsize, double alpha, bool ext)
{
    const int Nx = Nsize[0], Ny = Nsize[1], Nz = Nsize[2], Nt = Nsize[3];
    const int Nst = Nx * Ny * Nz * Nt;
    const int Nst_pad = ceil_nwp(Nst);

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    const int blockSize = VECTOR_LENGTH;
    const int gridSize  = (NC * Nst_pad + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_ee_LUdaginv_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, Ns, Nst, Nst_pad, (real_t)alpha);
    QTW_PROF_END("LUdaginv_qtw");
    CHECK(cudaDeviceSynchronize());
}

#undef IDX_DWF_QTW
#undef QTW_LOAD6
#undef QTW_STORE6
#undef QTW_DECL
#undef QTW_COPY
#undef QTW_ADDC
#undef QTW_SUBC
#undef QTW_NEGC
#undef QTW_SCALC

#endif // MULT_DOMAINWALL_5DIN_EO_INV_ACC_QTW_INCLUDED
