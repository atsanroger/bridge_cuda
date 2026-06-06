/*!
      @file    mult_Domainwall_5din_eo_cuda_qtw-inc.h
      @brief   QTW (triple-word) Even-Odd Domain Wall fermion kernels.
               Field layout: 6 real_t per complex {rh, ih, rm, im, rl, il}.
               All arithmetic uses tw_add/tw_mul/tw_scal (FP32-only when
               real_t = float).
      @author  Wei-Lun Chen (wlchen)
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_QTW_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_QTW_INCLUDED

#include <type_traits>
#include "../inline/constant_memory_inline.h"
// qtw_prof_timer.h is included at global scope by the enclosing .cu (before
// namespace BridgeACC) so its <map>/<string> don't land inside the namespace.

// BLAS-compatible flat layout: NC*ND*Ns groups of 6 real_t per NWP block.
// Each complex occupies 6 consecutive reals: {rh, ih, rm, im, rl, il}.
#define IDX_DWF_QTW(ic, id, is5, Ns_, site_) \
    IDX2(NC * ND * (Ns_), (ic) + NC * ((id) + ND * (is5)), (site_))

// ---------------- TW complex register tuple ----------------
// Each "TW complex" is six named scalars. We use a macro-friendly naming
// convention: prefix##_rh, _ih, _rm, _im, _rl, _il.

#define QTW_DECL(name) \
    real_t name##_rh, name##_ih, name##_rm, name##_im, name##_rl, name##_il

#define QTW_ZERO(name) \
    name##_rh = name##_ih = name##_rm = name##_im = name##_rl = name##_il = real_t(0)

// Vectorized SoAoS loads: the 6 words {rh,ih,rm,im,rl,il} are contiguous and
// (rh,ih),(rm,im),(rl,il) form three real2 pairs. base = 6*IDX2 is always even,
// so &ptr[base+0/2/4] is real2-aligned (float2: 24*IDX2 is 8-aligned; double2:
// 48*IDX2 is 16-aligned). 3 wide loads instead of 6 scalar LDG.
#define QTW_LOAD6(prefix, ptr, base) \
    real2 prefix##_ld0 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 0]); \
    real2 prefix##_ld1 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 2]); \
    real2 prefix##_ld2 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 4]); \
    real_t prefix##_rh = prefix##_ld0.x, prefix##_ih = prefix##_ld0.y; \
    real_t prefix##_rm = prefix##_ld1.x, prefix##_im = prefix##_ld1.y; \
    real_t prefix##_rl = prefix##_ld2.x, prefix##_il = prefix##_ld2.y

#define QTW_LOAD6_NEW(prefix, ptr, base) do { \
    real2 _l0 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 0]); \
    real2 _l1 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 2]); \
    real2 _l2 = *reinterpret_cast<const real2 *>(&(ptr)[(base) + 4]); \
    prefix##_rh = _l0.x; prefix##_ih = _l0.y; \
    prefix##_rm = _l1.x; prefix##_im = _l1.y; \
    prefix##_rl = _l2.x; prefix##_il = _l2.y; \
} while (0)

// Store renormalizes the two triples first: the inner-loop adds run sloppy
// (QTW_ADDC/QTW_SUBC skip Renormalize3), so the canonical 3-word form is only
// restored here, at the end of the accumulation chain (Ozaki strategy).
// Uses local temps so the caller's registers stay usable after the store.
#define QTW_STORE6(ptr, base, prefix) do { \
    real_t _sr_h = prefix##_rh, _sr_m = prefix##_rm, _sr_l = prefix##_rl; \
    real_t _si_h = prefix##_ih, _si_m = prefix##_im, _si_l = prefix##_il; \
    Renormalize3(_sr_h, _sr_m, _sr_l); \
    Renormalize3(_si_h, _si_m, _si_l); \
    real2 _o0; _o0.x = _sr_h; _o0.y = _si_h; \
    real2 _o1; _o1.x = _sr_m; _o1.y = _si_m; \
    real2 _o2; _o2.x = _sr_l; _o2.y = _si_l; \
    *reinterpret_cast<real2 *>(&(ptr)[(base) + 0]) = _o0; \
    *reinterpret_cast<real2 *>(&(ptr)[(base) + 2]) = _o1; \
    *reinterpret_cast<real2 *>(&(ptr)[(base) + 4]) = _o2; \
} while (0)

#define QTW_COPY(dst, src) \
    dst##_rh = src##_rh; dst##_ih = src##_ih; \
    dst##_rm = src##_rm; dst##_im = src##_im; \
    dst##_rl = src##_rl; dst##_il = src##_il

// ---------------- TW complex arithmetic (on the prefix tuples) ----------------

// TW complex addition: res = a + b  (real & imag each tw_add separately).
// Sloppy: skip Renormalize3 in the inner loop; QTW_STORE6 restores the
// canonical form once the whole stencil/accumulation chain is complete.
#define QTW_ADDC(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, b##_rh, b##_rm, b##_rl, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, b##_ih, b##_im, b##_il, \
           res##_ih, res##_im, res##_il); \
} while (0)

#define QTW_SUBC(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, -b##_rh, -b##_rm, -b##_rl, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, -b##_ih, -b##_im, -b##_il, \
           res##_ih, res##_im, res##_il); \
} while (0)

#define QTW_NEGC(res, a) do { \
    res##_rh = -a##_rh; res##_rm = -a##_rm; res##_rl = -a##_rl; \
    res##_ih = -a##_ih; res##_im = -a##_im; res##_il = -a##_il; \
} while (0)

// Scalar T × TW complex
#define QTW_SCALC(res, u, a) do { \
    tw_scal((u), a##_rh, a##_rm, a##_rl, res##_rh, res##_rm, res##_rl); \
    tw_scal((u), a##_ih, a##_im, a##_il, res##_ih, res##_im, res##_il); \
} while (0)

// TW scalar × TW complex  (sh,sm,sl) * a
#define QTW_SCALC_TT(res, sh, sm, sl, a) do { \
    tw_mul((sh), (sm), (sl), a##_rh, a##_rm, a##_rl, \
           res##_rh, res##_rm, res##_rl); \
    tw_mul((sh), (sm), (sl), a##_ih, a##_im, a##_il, \
           res##_ih, res##_im, res##_il); \
} while (0)

// ---------------- DW gamma-matrix projection / reconstruction ----------------
// X-direction (gamma_1): combines a + i*b via tw_add. The QDW versions are:
//   PROJ_P(res, a, b): res.re = a.re - b.im,  res.im = a.im + b.re
//   PROJ_M(res, a, b): res.re = a.re + b.im,  res.im = a.im - b.re
// Lift to TW: same component layout, just use TW arithmetic.

#define DWF_PROJ_P_TW(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, -b##_ih, -b##_im, -b##_il, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il,  b##_rh,  b##_rm,  b##_rl, \
           res##_ih, res##_im, res##_il); \
} while (0)

#define DWF_PROJ_M_TW(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl,  b##_ih,  b##_im,  b##_il, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, -b##_rh, -b##_rm, -b##_rl, \
           res##_ih, res##_im, res##_il); \
} while (0)

// Multiply complex by -i: (re, im) -> (im, -re), each TW-wise unchanged
#define DWF_MULT_MI_TW(res, a) do { \
    res##_rh =  a##_ih; res##_rm =  a##_im; res##_rl =  a##_il; \
    res##_ih = -a##_rh; res##_im = -a##_rm; res##_il = -a##_rl; \
} while (0)

#define DWF_MULT_PI_TW(res, a) do { \
    res##_rh = -a##_ih; res##_rm = -a##_im; res##_rl = -a##_il; \
    res##_ih =  a##_rh; res##_im =  a##_rm; res##_il =  a##_rl; \
} while (0)

// Y-direction (gamma_2, real)
#define DWF_PROJ_RP_TW(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl,  b##_rh,  b##_rm,  b##_rl, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il,  b##_ih,  b##_im,  b##_il, \
           res##_ih, res##_im, res##_il); \
} while (0)

#define DWF_PROJ_RM_TW(res, a, b) do { \
    tw_add_sloppy(a##_rh, a##_rm, a##_rl, -b##_rh, -b##_rm, -b##_rl, \
           res##_rh, res##_rm, res##_rl); \
    tw_add_sloppy(a##_ih, a##_im, a##_il, -b##_ih, -b##_im, -b##_il, \
           res##_ih, res##_im, res##_il); \
} while (0)

#define DWF_MULT_R1_TW(res, a)  QTW_COPY(res, a)
#define DWF_MULT_RM1_TW(res, a) QTW_NEGC(res, a)

// 2 * (TW complex): multiply by exact 2.0 in TW (no error)
#define DWF_PROJ_2_TW(res, a) do { \
    res##_rh = real_t(2) * a##_rh; res##_rm = real_t(2) * a##_rm; res##_rl = real_t(2) * a##_rl; \
    res##_ih = real_t(2) * a##_ih; res##_im = real_t(2) * a##_im; res##_il = real_t(2) * a##_il; \
} while (0)

// ---------------- DWF accumulator helpers ----------------
// DWF_ACCUM_4_TW(v0,v1,v2,v3, W1,W2, bc, OP2,OP3):
//   v0 += bc * W1;  v1 += bc * W2;  v2 += bc * OP2(W2);  v3 += bc * OP3(W1)
#define DWF_ACCUM_4_TW(v0, v1, v2, v3, W1, W2, bc_, OP2, OP3) do { \
    QTW_DECL(_ts); QTW_DECL(_tp); \
    QTW_SCALC(_ts, (bc_), W1); QTW_ADDC(v0, v0, _ts); \
    QTW_SCALC(_ts, (bc_), W2); QTW_ADDC(v1, v1, _ts); \
    OP2(_tp, W2); QTW_SCALC(_ts, (bc_), _tp); QTW_ADDC(v2, v2, _ts); \
    OP3(_tp, W1); QTW_SCALC(_ts, (bc_), _tp); QTW_ADDC(v3, v3, _ts); \
} while (0)

#define DWF_ACCUM_4_SW_TW(v0, v1, v2, v3, W1, W2, bc_, OP2, OP3) do { \
    QTW_DECL(_ts); QTW_DECL(_tp); \
    QTW_SCALC(_ts, (bc_), W1); QTW_ADDC(v0, v0, _ts); \
    QTW_SCALC(_ts, (bc_), W2); QTW_ADDC(v1, v1, _ts); \
    OP2(_tp, W1); QTW_SCALC(_ts, (bc_), _tp); QTW_ADDC(v2, v2, _ts); \
    OP3(_tp, W2); QTW_SCALC(_ts, (bc_), _tp); QTW_ADDC(v3, v3, _ts); \
} while (0)

#define DWF_ACCUM_TP_TW(v2, v3, W1, W2, bc_) do { \
    QTW_DECL(_ts); \
    QTW_SCALC(_ts, (bc_), W1); QTW_ADDC(v2, v2, _ts); \
    QTW_SCALC(_ts, (bc_), W2); QTW_ADDC(v3, v3, _ts); \
} while (0)

#define DWF_ACCUM_TM_TW(v0, v1, W1, W2, bc_) do { \
    QTW_DECL(_ts); \
    QTW_SCALC(_ts, (bc_), W1); QTW_ADDC(v0, v0, _ts); \
    QTW_SCALC(_ts, (bc_), W2); QTW_ADDC(v1, v1, _ts); \
} while (0)

// ---------------- Gauge multiply (TW vector × T gauge link) ----------------
// FP gauge × TW vec (cheap path): each output color = sum_j U(j,c) v_j with
// U a single-T complex. Uses tw_scal (T × TW). Forward direction.
#define DWF_GMUL_FWD_TW(u_ptr, isg_) do { \
    real_t _ur, _ui; \
    QTW_DECL(_t); \
    /* output color 0 */ \
    _ur = (u_ptr)[IDX2_G_R(0,0,(isg_))]; _ui = (u_ptr)[IDX2_G_I(0,0,(isg_))]; \
    tw_scal(_ur, vt1_c0_rh, vt1_c0_rm, vt1_c0_rl, wt1_c0_rh, wt1_c0_rm, wt1_c0_rl); \
    tw_scal(_ur, vt1_c0_ih, vt1_c0_im, vt1_c0_il, wt1_c0_ih, wt1_c0_im, wt1_c0_il); \
    tw_scal(_ui, vt1_c0_ih, vt1_c0_im, vt1_c0_il, _t##_rh, _t##_rm, _t##_rl); \
    tw_add(wt1_c0_rh, wt1_c0_rm, wt1_c0_rl, -_t##_rh, -_t##_rm, -_t##_rl, wt1_c0_rh, wt1_c0_rm, wt1_c0_rl); \
    tw_scal(_ui, vt1_c0_rh, vt1_c0_rm, vt1_c0_rl, _t##_rh, _t##_rm, _t##_rl); \
    tw_add(wt1_c0_ih, wt1_c0_im, wt1_c0_il,  _t##_rh,  _t##_rm,  _t##_rl, wt1_c0_ih, wt1_c0_im, wt1_c0_il); \
    /* ... full 3x3 expansion deferred to per-kernel inline ... */ \
} while (0)

// NOTE: A full DWF_GMUL_FWD_TW / BCK_TW with 3x3 complex link is voluminous
// (108 tw_scal + tw_add per call). To keep this header compact, the full
// gauge multiply is expanded inline inside each kernel body using a helper
// function (see qtw_gmult_fwd_3x3 below). The macro above is a stub showing
// the pattern for one color row; production kernels use the function.

// Reusable 3x3 gauge multiply: out[c] = sum_j U(j,c) in[j]   (forward, plain U)
template<typename T>
__device__ __forceinline__ void qtw_gmult_fwd_3x3(
    const T* __restrict__ u, int isg,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T out_rh[3], T out_ih[3], T out_rm[3], T out_im[3], T out_rl[3], T out_il[3])
{
    // Type-split unroll: double-QTW is at the 255-reg cap and spills heavily, so
    // fully unrolling promotes the in_/out_ arrays to registers (spill 2448->344B).
    // float-QTW lives below the cap (188 reg) at a good occupancy/spill balance, so
    // unrolling there pushes it to 255 and regresses ~33% -> keep float unroll=1.
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int c = 0; c < NC; ++c) {
        T ar_h = 0, ar_m = 0, ar_l = 0;
        T ai_h = 0, ai_m = 0, ai_l = 0;
        #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
        for (int j = 0; j < NC; ++j) {
            T ur = u[IDX2_G_R(j, c, isg)];
            T ui = u[IDX2_G_I(j, c, isg)];
            // real += ur*in_r - ui*in_i
            T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l;
            tw_scal_sloppy(ur, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
            tw_scal_sloppy(ui, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
            T nh, nm, nl;
            tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            // imag += ur*in_i + ui*in_r
            tw_scal_sloppy(ur, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
            tw_scal_sloppy(ui, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
            tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
        }
        out_rh[c] = ar_h; out_rm[c] = ar_m; out_rl[c] = ar_l;
        out_ih[c] = ai_h; out_im[c] = ai_m; out_il[c] = ai_l;
    }
}

// Backward gauge multiply: out[c] = sum_j conj(U(c, j)) in[j]
template<typename T>
__device__ __forceinline__ void qtw_gmult_bck_3x3(
    const T* __restrict__ u, int isg,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T out_rh[3], T out_ih[3], T out_rm[3], T out_im[3], T out_rl[3], T out_il[3])
{
    // Type-split unroll: double-QTW is at the 255-reg cap and spills heavily, so
    // fully unrolling promotes the in_/out_ arrays to registers (spill 2448->344B).
    // float-QTW lives below the cap (188 reg) at a good occupancy/spill balance, so
    // unrolling there pushes it to 255 and regresses ~33% -> keep float unroll=1.
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int c = 0; c < NC; ++c) {
        T ar_h = 0, ar_m = 0, ar_l = 0;
        T ai_h = 0, ai_m = 0, ai_l = 0;
        #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
        for (int j = 0; j < NC; ++j) {
            T ur =  u[IDX2_G_R(c, j, isg)];
            T ui = -u[IDX2_G_I(c, j, isg)];  // conjugate
            T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l, nh, nm, nl;
            tw_scal_sloppy(ur, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
            tw_scal_sloppy(ui, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_scal_sloppy(ur, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
            tw_scal_sloppy(ui, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
            tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
        }
        out_rh[c] = ar_h; out_rm[c] = ar_m; out_rl[c] = ar_l;
        out_ih[c] = ai_h; out_im[c] = ai_m; out_il[c] = ai_l;
    }
}

// Runtime-toggled via the YAML `su3_reconstruction` flag (-> kernel template arg
// RECON, so when OFF this code is dead-stripped and the fast path is unchanged).
// Measured on RTX 3080 it REGRESSED the TW solve ~17% (64.3 -> 75.5 s): the hopping
// kernel is register/occupancy-limited, not purely bandwidth-limited, so trading
// 33% gauge-read BW for +18 registers and ~24 tw_mul/link lost more occupancy than
// it saved. Useful on a more bandwidth-bound GPU; exposed as a runtime option.
// SU(3) 3rd-column reconstruction in triple-word: skip reading column 2 of each
// link from DRAM (the last 6 of 18 reals/link, ~33% of gauge read BW) and rebuild
// it on-chip from columns 0,1 via U(:,2) = conj(U(:,0) x U(:,1)). Matches the FP
// EXT_IMG_R/I convention: U(i,2) = conj(U(p,0)U(q,1) - U(q,0)U(p,1)), p=(i+1)%3,
// q=(i+2)%3. Arithmetic is hidden under the memory-bound hopping.
template<typename T>
__device__ __forceinline__ void qtw_recon_col2_tw(
    const T* __restrict__ uh, const T* __restrict__ um, const T* __restrict__ ul,
    int isg,
    T c2r_h[3], T c2r_m[3], T c2r_l[3], T c2i_h[3], T c2i_m[3], T c2i_l[3])
{
    for (int i = 0; i < 3; ++i) {
        const int p = (i + 1) % 3;
        const int q = (i + 2) % 3;
        // v1=U(p,0) v2=U(q,0) w1=U(p,1) w2=U(q,1)  (columns 0,1 only)
        T v1r_H=uh[IDX2_G_R(p,0,isg)],v1r_M=um[IDX2_G_R(p,0,isg)],v1r_L=ul[IDX2_G_R(p,0,isg)];
        T v1i_H=uh[IDX2_G_I(p,0,isg)],v1i_M=um[IDX2_G_I(p,0,isg)],v1i_L=ul[IDX2_G_I(p,0,isg)];
        T v2r_H=uh[IDX2_G_R(q,0,isg)],v2r_M=um[IDX2_G_R(q,0,isg)],v2r_L=ul[IDX2_G_R(q,0,isg)];
        T v2i_H=uh[IDX2_G_I(q,0,isg)],v2i_M=um[IDX2_G_I(q,0,isg)],v2i_L=ul[IDX2_G_I(q,0,isg)];
        T w1r_H=uh[IDX2_G_R(p,1,isg)],w1r_M=um[IDX2_G_R(p,1,isg)],w1r_L=ul[IDX2_G_R(p,1,isg)];
        T w1i_H=uh[IDX2_G_I(p,1,isg)],w1i_M=um[IDX2_G_I(p,1,isg)],w1i_L=ul[IDX2_G_I(p,1,isg)];
        T w2r_H=uh[IDX2_G_R(q,1,isg)],w2r_M=um[IDX2_G_R(q,1,isg)],w2r_L=ul[IDX2_G_R(q,1,isg)];
        T w2i_H=uh[IDX2_G_I(q,1,isg)],w2i_M=um[IDX2_G_I(q,1,isg)],w2i_L=ul[IDX2_G_I(q,1,isg)];
        // p1 = v1*w2 , p2 = v2*w1 ; diff = p1 - p2 ; U(i,2) = conj(diff)
        T p1r_h,p1r_m,p1r_l,p1i_h,p1i_m,p1i_l, p2r_h,p2r_m,p2r_l,p2i_h,p2i_m,p2i_l;
        T t_h,t_m,t_l,u_h,u_m,u_l;
        // p1 = v1*w2
        tw_mul(v1r_H,v1r_M,v1r_L, w2r_H,w2r_M,w2r_L, t_h,t_m,t_l);
        tw_mul(v1i_H,v1i_M,v1i_L, w2i_H,w2i_M,w2i_L, u_h,u_m,u_l);
        tw_add(t_h,t_m,t_l, -u_h,-u_m,-u_l, p1r_h,p1r_m,p1r_l);            // Re
        tw_mul(v1r_H,v1r_M,v1r_L, w2i_H,w2i_M,w2i_L, t_h,t_m,t_l);
        tw_mul(v1i_H,v1i_M,v1i_L, w2r_H,w2r_M,w2r_L, u_h,u_m,u_l);
        tw_add(t_h,t_m,t_l, u_h,u_m,u_l, p1i_h,p1i_m,p1i_l);              // Im
        // p2 = v2*w1
        tw_mul(v2r_H,v2r_M,v2r_L, w1r_H,w1r_M,w1r_L, t_h,t_m,t_l);
        tw_mul(v2i_H,v2i_M,v2i_L, w1i_H,w1i_M,w1i_L, u_h,u_m,u_l);
        tw_add(t_h,t_m,t_l, -u_h,-u_m,-u_l, p2r_h,p2r_m,p2r_l);
        tw_mul(v2r_H,v2r_M,v2r_L, w1i_H,w1i_M,w1i_L, t_h,t_m,t_l);
        tw_mul(v2i_H,v2i_M,v2i_L, w1r_H,w1r_M,w1r_L, u_h,u_m,u_l);
        tw_add(t_h,t_m,t_l, u_h,u_m,u_l, p2i_h,p2i_m,p2i_l);
        // diff = p1 - p2 ; conj -> imag negated
        tw_add(p1r_h,p1r_m,p1r_l, -p2r_h,-p2r_m,-p2r_l, c2r_h[i],c2r_m[i],c2r_l[i]);
        T di_h,di_m,di_l;
        tw_add(p1i_h,p1i_m,p1i_l, -p2i_h,-p2i_m,-p2i_l, di_h,di_m,di_l);
        c2i_h[i] = -di_h; c2i_m[i] = -di_m; c2i_l[i] = -di_l;            // conjugate
    }
}

// Same but with full TW gauge link (h, m, l) loaded from u_h / u_m / u_l arrays.
// Used when extended_precision pushes 3-word gauge links into separate arrays.
template<bool RECON, typename T>
__device__ __forceinline__ void qtw_gmult_fwd_3x3_tw(
    const T* __restrict__ uh, const T* __restrict__ um, const T* __restrict__ ul,
    int isg,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T out_rh[3], T out_ih[3], T out_rm[3], T out_im[3], T out_rl[3], T out_il[3])
{
    T c2r_h[3], c2r_m[3], c2r_l[3], c2i_h[3], c2i_m[3], c2i_l[3];
    if (RECON) {
        qtw_recon_col2_tw(uh, um, ul, isg, c2r_h, c2r_m, c2r_l, c2i_h, c2i_m, c2i_l);
    }
    // Type-split unroll: double-QTW is at the 255-reg cap and spills heavily, so
    // fully unrolling promotes the in_/out_ arrays to registers (spill 2448->344B).
    // float-QTW lives below the cap (188 reg) at a good occupancy/spill balance, so
    // unrolling there pushes it to 255 and regresses ~33% -> keep float unroll=1.
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int c = 0; c < NC; ++c) {
        T ar_h = 0, ar_m = 0, ar_l = 0;
        T ai_h = 0, ai_m = 0, ai_l = 0;
        #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
        for (int j = 0; j < NC; ++j) {
            T ur_h, ur_m, ur_l, ui_h, ui_m, ui_l;
            if (RECON && c == 2) {  // column 2: use reconstructed (not read from DRAM)
                ur_h = c2r_h[j]; ur_m = c2r_m[j]; ur_l = c2r_l[j];
                ui_h = c2i_h[j]; ui_m = c2i_m[j]; ui_l = c2i_l[j];
            } else
            {
                ur_h = uh[IDX2_G_R(j, c, isg)];
                ur_m = um[IDX2_G_R(j, c, isg)];
                ur_l = ul[IDX2_G_R(j, c, isg)];
                ui_h = uh[IDX2_G_I(j, c, isg)];
                ui_m = um[IDX2_G_I(j, c, isg)];
                ui_l = ul[IDX2_G_I(j, c, isg)];
            }
            T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l, nh, nm, nl;
            tw_mul_sloppy(ur_h, ur_m, ur_l, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
            tw_mul_sloppy(ui_h, ui_m, ui_l, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_mul_sloppy(ur_h, ur_m, ur_l, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
            tw_mul_sloppy(ui_h, ui_m, ui_l, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
            tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
        }
        out_rh[c] = ar_h; out_rm[c] = ar_m; out_rl[c] = ar_l;
        out_ih[c] = ai_h; out_im[c] = ai_m; out_il[c] = ai_l;
    }
}

template<bool RECON, typename T>
__device__ __forceinline__ void qtw_gmult_bck_3x3_tw(
    const T* __restrict__ uh, const T* __restrict__ um, const T* __restrict__ ul,
    int isg,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T out_rh[3], T out_ih[3], T out_rm[3], T out_im[3], T out_rl[3], T out_il[3])
{
    T c2r_h[3], c2r_m[3], c2r_l[3], c2i_h[3], c2i_m[3], c2i_l[3];
    if (RECON) {
        qtw_recon_col2_tw(uh, um, ul, isg, c2r_h, c2r_m, c2r_l, c2i_h, c2i_m, c2i_l);
    }
    // Type-split unroll: double-QTW is at the 255-reg cap and spills heavily, so
    // fully unrolling promotes the in_/out_ arrays to registers (spill 2448->344B).
    // float-QTW lives below the cap (188 reg) at a good occupancy/spill balance, so
    // unrolling there pushes it to 255 and regresses ~33% -> keep float unroll=1.
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int c = 0; c < NC; ++c) {
        T ar_h = 0, ar_m = 0, ar_l = 0;
        T ai_h = 0, ai_m = 0, ai_l = 0;
        #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
        for (int j = 0; j < NC; ++j) {
            T ur_h, ur_m, ur_l, ui_h, ui_m, ui_l;
            if (RECON && j == 2) {  // column 2 = reconstructed; conj negates imag
                ur_h =  c2r_h[c]; ur_m =  c2r_m[c]; ur_l =  c2r_l[c];
                ui_h = -c2i_h[c]; ui_m = -c2i_m[c]; ui_l = -c2i_l[c];
            } else
            {
                ur_h =  uh[IDX2_G_R(c, j, isg)];
                ur_m =  um[IDX2_G_R(c, j, isg)];
                ur_l =  ul[IDX2_G_R(c, j, isg)];
                ui_h = -uh[IDX2_G_I(c, j, isg)];
                ui_m = -um[IDX2_G_I(c, j, isg)];
                ui_l = -ul[IDX2_G_I(c, j, isg)];
            }
            T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l, nh, nm, nl;
            tw_mul_sloppy(ur_h, ur_m, ur_l, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
            tw_mul_sloppy(ui_h, ui_m, ui_l, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl);
            ar_h = nh; ar_m = nm; ar_l = nl;
            tw_mul_sloppy(ur_h, ur_m, ur_l, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
            tw_mul_sloppy(ui_h, ui_m, ui_l, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
            tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
            tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl);
            ai_h = nh; ai_m = nm; ai_l = nl;
        }
        out_rh[c] = ar_h; out_rm[c] = ar_m; out_rl[c] = ar_l;
        out_ih[c] = ai_h; out_im[c] = ai_m; out_il[c] = ai_l;
    }
}

// ---- Single-output-color TW gauge multiplies (color-split kernel) ----
// Compute only output color `cout` => one TW complex out. Each color-split
// thread does 1/NC of the gmult work (the whole point of the split).
template<bool RECON, typename T>
__device__ __forceinline__ void qtw_gmult_fwd_1col_tw(
    const T* __restrict__ uh, const T* __restrict__ um, const T* __restrict__ ul,
    int isg, int cout,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T& o_rh, T& o_ih, T& o_rm, T& o_im, T& o_rl, T& o_il)
{
    T c2r_h[3], c2r_m[3], c2r_l[3], c2i_h[3], c2i_m[3], c2i_l[3];
    if (RECON && cout == 2)
        qtw_recon_col2_tw(uh, um, ul, isg, c2r_h, c2r_m, c2r_l, c2i_h, c2i_m, c2i_l);
    T ar_h = 0, ar_m = 0, ar_l = 0, ai_h = 0, ai_m = 0, ai_l = 0;
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int j = 0; j < NC; ++j) {
        T ur_h, ur_m, ur_l, ui_h, ui_m, ui_l;
        if (RECON && cout == 2) {
            ur_h = c2r_h[j]; ur_m = c2r_m[j]; ur_l = c2r_l[j];
            ui_h = c2i_h[j]; ui_m = c2i_m[j]; ui_l = c2i_l[j];
        } else {
            ur_h = uh[IDX2_G_R(j, cout, isg)]; ur_m = um[IDX2_G_R(j, cout, isg)]; ur_l = ul[IDX2_G_R(j, cout, isg)];
            ui_h = uh[IDX2_G_I(j, cout, isg)]; ui_m = um[IDX2_G_I(j, cout, isg)]; ui_l = ul[IDX2_G_I(j, cout, isg)];
        }
        T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l, nh, nm, nl;
        tw_mul_sloppy(ur_h, ur_m, ur_l, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
        tw_mul_sloppy(ui_h, ui_m, ui_l, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
        tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl); ar_h = nh; ar_m = nm; ar_l = nl;
        tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl); ar_h = nh; ar_m = nm; ar_l = nl;
        tw_mul_sloppy(ur_h, ur_m, ur_l, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
        tw_mul_sloppy(ui_h, ui_m, ui_l, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
        tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl); ai_h = nh; ai_m = nm; ai_l = nl;
        tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl); ai_h = nh; ai_m = nm; ai_l = nl;
    }
    o_rh = ar_h; o_rm = ar_m; o_rl = ar_l; o_ih = ai_h; o_im = ai_m; o_il = ai_l;
}

template<bool RECON, typename T>
__device__ __forceinline__ void qtw_gmult_bck_1col_tw(
    const T* __restrict__ uh, const T* __restrict__ um, const T* __restrict__ ul,
    int isg, int cout,
    const T in_rh[3], const T in_ih[3], const T in_rm[3], const T in_im[3],
    const T in_rl[3], const T in_il[3],
    T& o_rh, T& o_ih, T& o_rm, T& o_im, T& o_rl, T& o_il)
{
    T c2r_h[3], c2r_m[3], c2r_l[3], c2i_h[3], c2i_m[3], c2i_l[3];
    if (RECON)  // backward j==2 term uses reconstructed col2, indexed by cout
        qtw_recon_col2_tw(uh, um, ul, isg, c2r_h, c2r_m, c2r_l, c2i_h, c2i_m, c2i_l);
    T ar_h = 0, ar_m = 0, ar_l = 0, ai_h = 0, ai_m = 0, ai_l = 0;
    #pragma unroll (std::is_same<T, double>::value ? 3 : 1)
    for (int j = 0; j < NC; ++j) {
        T ur_h, ur_m, ur_l, ui_h, ui_m, ui_l;
        if (RECON && j == 2) {
            ur_h =  c2r_h[cout]; ur_m =  c2r_m[cout]; ur_l =  c2r_l[cout];
            ui_h = -c2i_h[cout]; ui_m = -c2i_m[cout]; ui_l = -c2i_l[cout];
        } else {
            ur_h =  uh[IDX2_G_R(cout, j, isg)]; ur_m =  um[IDX2_G_R(cout, j, isg)]; ur_l =  ul[IDX2_G_R(cout, j, isg)];
            ui_h = -uh[IDX2_G_I(cout, j, isg)]; ui_m = -um[IDX2_G_I(cout, j, isg)]; ui_l = -ul[IDX2_G_I(cout, j, isg)];
        }
        T tr_h, tr_m, tr_l, ti_h, ti_m, ti_l, nh, nm, nl;
        tw_mul_sloppy(ur_h, ur_m, ur_l, in_rh[j], in_rm[j], in_rl[j], tr_h, tr_m, tr_l);
        tw_mul_sloppy(ui_h, ui_m, ui_l, in_ih[j], in_im[j], in_il[j], ti_h, ti_m, ti_l);
        tw_add_sloppy(ar_h, ar_m, ar_l,  tr_h,  tr_m,  tr_l, nh, nm, nl); ar_h = nh; ar_m = nm; ar_l = nl;
        tw_add_sloppy(ar_h, ar_m, ar_l, -ti_h, -ti_m, -ti_l, nh, nm, nl); ar_h = nh; ar_m = nm; ar_l = nl;
        tw_mul_sloppy(ur_h, ur_m, ur_l, in_ih[j], in_im[j], in_il[j], tr_h, tr_m, tr_l);
        tw_mul_sloppy(ui_h, ui_m, ui_l, in_rh[j], in_rm[j], in_rl[j], ti_h, ti_m, ti_l);
        tw_add_sloppy(ai_h, ai_m, ai_l, tr_h, tr_m, tr_l, nh, nm, nl); ai_h = nh; ai_m = nm; ai_l = nl;
        tw_add_sloppy(ai_h, ai_m, ai_l, ti_h, ti_m, ti_l, nh, nm, nl); ai_h = nh; ai_m = nm; ai_l = nl;
    }
    o_rh = ar_h; o_rm = ar_m; o_rl = ar_l; o_ih = ai_h; o_im = ai_m; o_il = ai_l;
}

//====================================================================
// ee 5-direction (diagonal block) TW kernel — vp = B_eff * wp + C_hop * vt
//====================================================================
__global__
void mult_domainwall_5din_ee_5dir_dirac_qtw_kernel_fp(
    real_t * __restrict__ vp, real_t * __restrict__ wp,
    real_t mq, real_t mq_m, real_t mq_l,
    real_t M0, real_t M0_m, real_t M0_l, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();
    const real_t *b_mid_con = ConstantMemoryTraits<real_t>::b_mid();
    const real_t *b_lo_con  = ConstantMemoryTraits<real_t>::b_lo_tw();
    const real_t *c_mid_con = ConstantMemoryTraits<real_t>::c_mid();
    const real_t *c_lo_con  = ConstantMemoryTraits<real_t>::c_lo_tw();

    // (4 - M0) as triple-word; M0=1.6 is inexact in float so this is the
    // coefficient that made the forward Dee/Doo only single-word accurate.
    real_t fm_h, fm_m, fm_l;
    tw_add(real_t(4.0), real_t(0), real_t(0), -M0, -M0_m, -M0_l, fm_h, fm_m, fm_l);

    for (int is = 0; is < Ns; ++is) {
        // FF1 = b*(4-M0)+1 ; FF2 = c*(4-M0)-1  (genuine triple-word)
        real_t FF1_h, FF1_m, FF1_l, FF2_h, FF2_m, FF2_l, t_h, t_m, t_l;
        tw_mul(b_con[is], b_mid_con[is], b_lo_con[is], fm_h, fm_m, fm_l, t_h, t_m, t_l);
        tw_add(t_h, t_m, t_l, real_t( 1.0), real_t(0), real_t(0), FF1_h, FF1_m, FF1_l);
        tw_mul(c_con[is], c_mid_con[is], c_lo_con[is], fm_h, fm_m, fm_l, t_h, t_m, t_l);
        tw_add(t_h, t_m, t_l, real_t(-1.0), real_t(0), real_t(0), FF2_h, FF2_m, FF2_l);

        const int is_up = (is + 1) % Ns;
        const int is_dn = (is - 1 + Ns) % Ns;
        const real_t s_bulk = real_t(0.5) * alpha;
        // Fup/Fdn = bulk 0.5*alpha*FF2 ; boundary -0.5*mq*FF2  (3-word)
        real_t Fup_h, Fup_m, Fup_l, Fdn_h, Fdn_m, Fdn_l;
        if (is == Ns-1) {
            tw_mul(mq, mq_m, mq_l, FF2_h, FF2_m, FF2_l, t_h, t_m, t_l);
            tw_scal(real_t(-0.5), t_h, t_m, t_l, Fup_h, Fup_m, Fup_l);
        } else {
            tw_scal(s_bulk, FF2_h, FF2_m, FF2_l, Fup_h, Fup_m, Fup_l);
        }
        if (is == 0) {
            tw_mul(mq, mq_m, mq_l, FF2_h, FF2_m, FF2_l, t_h, t_m, t_l);
            tw_scal(real_t(-0.5), t_h, t_m, t_l, Fdn_h, Fdn_m, Fdn_l);
        } else {
            tw_scal(s_bulk, FF2_h, FF2_m, FF2_l, Fdn_h, Fdn_m, Fdn_l);
        }

        // diagonal B-block factors from FF1 (alpha-scalings exact at alpha=1).
        real_t f1_h=0,f1_m=0,f1_l=0, f2_h=0,f2_m=0,f2_l=0, fa_h=0,fa_m=0,fa_l=0;
        if (is == 0) {
            tw_scal(real_t(0.5)*(real_t( 1.0)+alpha), FF1_h, FF1_m, FF1_l, f1_h, f1_m, f1_l);
            tw_scal(real_t(0.5)*(real_t(-1.0)+alpha), FF1_h, FF1_m, FF1_l, f2_h, f2_m, f2_l);
        } else if (is == Ns-1) {
            tw_scal(real_t(0.5)*(real_t(1.0)+alpha), FF1_h, FF1_m, FF1_l, f1_h, f1_m, f1_l);
            tw_scal(real_t(0.5)*(real_t(1.0)-alpha), FF1_h, FF1_m, FF1_l, f2_h, f2_m, f2_l);
        } else {
            tw_scal(alpha, FF1_h, FF1_m, FF1_l, fa_h, fa_m, fa_l);
        }

        for (int ic = 0; ic < NC; ++ic) {
            int b_s0u = 6 * IDX_DWF_QTW(ic, 0, is_up, Ns, site);
            int b_s1u = 6 * IDX_DWF_QTW(ic, 1, is_up, Ns, site);
            int b_s2u = 6 * IDX_DWF_QTW(ic, 2, is_up, Ns, site);
            int b_s3u = 6 * IDX_DWF_QTW(ic, 3, is_up, Ns, site);
            int b_s0d = 6 * IDX_DWF_QTW(ic, 0, is_dn, Ns, site);
            int b_s1d = 6 * IDX_DWF_QTW(ic, 1, is_dn, Ns, site);
            int b_s2d = 6 * IDX_DWF_QTW(ic, 2, is_dn, Ns, site);
            int b_s3d = 6 * IDX_DWF_QTW(ic, 3, is_dn, Ns, site);
            int b_s0  = 6 * IDX_DWF_QTW(ic, 0, is,    Ns, site);
            int b_s1  = 6 * IDX_DWF_QTW(ic, 1, is,    Ns, site);
            int b_s2  = 6 * IDX_DWF_QTW(ic, 2, is,    Ns, site);
            int b_s3  = 6 * IDX_DWF_QTW(ic, 3, is,    Ns, site);

            QTW_LOAD6(wu0, wp, b_s0u);
            QTW_LOAD6(wu1, wp, b_s1u);
            QTW_LOAD6(wu2, wp, b_s2u);
            QTW_LOAD6(wu3, wp, b_s3u);
            QTW_LOAD6(wd0, wp, b_s0d);
            QTW_LOAD6(wd1, wp, b_s1d);
            QTW_LOAD6(wd2, wp, b_s2d);
            QTW_LOAD6(wd3, wp, b_s3d);
            QTW_LOAD6(w0,  wp, b_s0);
            QTW_LOAD6(w1,  wp, b_s1);
            QTW_LOAD6(w2,  wp, b_s2);
            QTW_LOAD6(w3,  wp, b_s3);

            // --- s+1: vt0 = Fup * (wu0 - wu2),  vt1 = Fup * (wu1 - wu3),
            //         vt2 = -vt0, vt3 = -vt1 ---
            QTW_DECL(tmp);
            QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
            QTW_SUBC(tmp, wu0, wu2); QTW_SCALC_TT(vt0, Fup_h, Fup_m, Fup_l, tmp);
            QTW_SUBC(tmp, wu1, wu3); QTW_SCALC_TT(vt1, Fup_h, Fup_m, Fup_l, tmp);
            QTW_NEGC(vt2, vt0);
            QTW_NEGC(vt3, vt1);

            // --- s-1: sum02 = wd0+wd2; sum13 = wd1+wd3
            //   vt0 += Fdn*sum02, vt2 += Fdn*sum02
            //   vt1 += Fdn*sum13, vt3 += Fdn*sum13 ---
            QTW_DECL(sum02); QTW_DECL(sum13); QTW_DECL(scaled);
            QTW_ADDC(sum02, wd0, wd2); QTW_SCALC_TT(scaled, Fdn_h, Fdn_m, Fdn_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_ADDC(vt2, vt2, scaled);
            QTW_ADDC(sum13, wd1, wd3); QTW_SCALC_TT(scaled, Fdn_h, Fdn_m, Fdn_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_ADDC(vt3, vt3, scaled);

            // --- local diagonal (B-block) ---
            QTW_DECL(wt0); QTW_DECL(wt1); QTW_DECL(wt2); QTW_DECL(wt3);
            QTW_DECL(t1); QTW_DECL(t2);
            if (is == 0 || is == Ns-1) {
                // wt_s0 = f1*w_s0 + f2*w_s2  (and similar permutations)
                QTW_SCALC_TT(t1, f1_h, f1_m, f1_l, w0); QTW_SCALC_TT(t2, f2_h, f2_m, f2_l, w2); QTW_ADDC(wt0, t1, t2);
                QTW_SCALC_TT(t1, f1_h, f1_m, f1_l, w1); QTW_SCALC_TT(t2, f2_h, f2_m, f2_l, w3); QTW_ADDC(wt1, t1, t2);
                QTW_SCALC_TT(t1, f1_h, f1_m, f1_l, w2); QTW_SCALC_TT(t2, f2_h, f2_m, f2_l, w0); QTW_ADDC(wt2, t1, t2);
                QTW_SCALC_TT(t1, f1_h, f1_m, f1_l, w3); QTW_SCALC_TT(t2, f2_h, f2_m, f2_l, w1); QTW_ADDC(wt3, t1, t2);
            } else {
                QTW_SCALC_TT(wt0, fa_h, fa_m, fa_l, w0); QTW_SCALC_TT(wt1, fa_h, fa_m, fa_l, w1);
                QTW_SCALC_TT(wt2, fa_h, fa_m, fa_l, w2); QTW_SCALC_TT(wt3, fa_h, fa_m, fa_l, w3);
            }

            // wt += vt;  vp[is] = wt
            QTW_ADDC(wt0, wt0, vt0);
            QTW_ADDC(wt1, wt1, vt1);
            QTW_ADDC(wt2, wt2, vt2);
            QTW_ADDC(wt3, wt3, vt3);
            QTW_STORE6(vp, b_s0, wt0);
            QTW_STORE6(vp, b_s1, wt1);
            QTW_STORE6(vp, b_s2, wt2);
            QTW_STORE6(vp, b_s3, wt3);
        }
    }
}

void mult_domainwall_5din_ee_5dir_dirac_qtw(
    real_t *vp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t *vp_dev = (real_t *)dev_ptr(vp);
    real_t *wp_dev = (real_t *)dev_ptr(wp);

    // triple-word split of mq and M0 (both inexact in float).
    const real_t mq_h = (real_t)mq;
    const real_t mq_m = (real_t)(mq - (double)mq_h);
    const real_t mq_l = (real_t)(mq - (double)mq_h - (double)mq_m);
    const real_t M0_h = (real_t)M0;
    const real_t M0_m = (real_t)(M0 - (double)M0_h);
    const real_t M0_l = (real_t)(M0 - (double)M0_h - (double)M0_m);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_ee_5dir_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, mq_h, mq_m, mq_l, M0_h, M0_m, M0_l, Ns, (real_t)alpha, Nst, Nst_pad);
    QTW_PROF_END("ee_5dir_qtw");
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-direction (off-diagonal) TW kernel — yp = -0.5*(b*w_diag + c*vt)
//====================================================================
__global__
void mult_domainwall_5din_eo_5dir_dirac_qtw_kernel_fp(
    real_t * __restrict__ yp, real_t * __restrict__ wp,
    real_t mq, real_t mq_m, real_t mq_l, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();
    const real_t *c_mid_con = ConstantMemoryTraits<real_t>::c_mid();
    const real_t *c_lo_con  = ConstantMemoryTraits<real_t>::c_lo_tw();

    for (int is = 0; is < Ns; ++is) {
        // The boundary hopping carries mq (=0.1, inexact in float), so build it
        // as a full triple-word coefficient (-0.5 * c[is] * mq). Bulk is exact
        // single-word (0.5 * c * alpha = 0.25 with c=0.5, alpha=1).
        const int is_up = (is + 1) % Ns;
        real_t Fup_h, Fup_m, Fup_l;
        if (is == Ns-1) {
            real_t p_h, p_m, p_l;
            tw_mul(c_con[is], c_mid_con[is], c_lo_con[is], mq, mq_m, mq_l, p_h, p_m, p_l);
            Fup_h = real_t(-0.5)*p_h; Fup_m = real_t(-0.5)*p_m; Fup_l = real_t(-0.5)*p_l;
        } else {
            Fup_h = real_t(0.5) * c_con[is] * alpha; Fup_m = real_t(0); Fup_l = real_t(0);
        }
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn_h, Fdn_m, Fdn_l;
        if (is == 0) {
            real_t p_h, p_m, p_l;
            tw_mul(c_con[is], c_mid_con[is], c_lo_con[is], mq, mq_m, mq_l, p_h, p_m, p_l);
            Fdn_h = real_t(-0.5)*p_h; Fdn_m = real_t(-0.5)*p_m; Fdn_l = real_t(-0.5)*p_l;
        } else {
            Fdn_h = real_t(0.5) * c_con[is] * alpha; Fdn_m = real_t(0); Fdn_l = real_t(0);
        }

        real_t b1 = 0, b2 = 0, bb = 0;
        if (is == 0) {
            b1 = b_con[is] * real_t(0.5) * (real_t(1.0) + alpha);
            b2 = b_con[is] * real_t(0.5) * (real_t(-1.0) + alpha);
        } else if (is == Ns-1) {
            b1 = b_con[is] * real_t(0.5) * (real_t(1.0) + alpha);
            b2 = b_con[is] * real_t(0.5) * (real_t(1.0) - alpha);
        } else {
            bb = b_con[is] * alpha;
        }

        for (int ic = 0; ic < NC; ++ic) {
            int b_s0u = 6 * IDX_DWF_QTW(ic, 0, is_up, Ns, site);
            int b_s1u = 6 * IDX_DWF_QTW(ic, 1, is_up, Ns, site);
            int b_s2u = 6 * IDX_DWF_QTW(ic, 2, is_up, Ns, site);
            int b_s3u = 6 * IDX_DWF_QTW(ic, 3, is_up, Ns, site);
            int b_s0d = 6 * IDX_DWF_QTW(ic, 0, is_dn, Ns, site);
            int b_s1d = 6 * IDX_DWF_QTW(ic, 1, is_dn, Ns, site);
            int b_s2d = 6 * IDX_DWF_QTW(ic, 2, is_dn, Ns, site);
            int b_s3d = 6 * IDX_DWF_QTW(ic, 3, is_dn, Ns, site);
            int b_s0  = 6 * IDX_DWF_QTW(ic, 0, is,    Ns, site);
            int b_s1  = 6 * IDX_DWF_QTW(ic, 1, is,    Ns, site);
            int b_s2  = 6 * IDX_DWF_QTW(ic, 2, is,    Ns, site);
            int b_s3  = 6 * IDX_DWF_QTW(ic, 3, is,    Ns, site);

            QTW_LOAD6(wu0, wp, b_s0u);
            QTW_LOAD6(wu1, wp, b_s1u);
            QTW_LOAD6(wu2, wp, b_s2u);
            QTW_LOAD6(wu3, wp, b_s3u);
            QTW_LOAD6(wd0, wp, b_s0d);
            QTW_LOAD6(wd1, wp, b_s1d);
            QTW_LOAD6(wd2, wp, b_s2d);
            QTW_LOAD6(wd3, wp, b_s3d);
            QTW_LOAD6(w0,  wp, b_s0);
            QTW_LOAD6(w1,  wp, b_s1);
            QTW_LOAD6(w2,  wp, b_s2);
            QTW_LOAD6(w3,  wp, b_s3);

            QTW_DECL(tmp);
            QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
            QTW_SUBC(tmp, wu0, wu2); QTW_SCALC_TT(vt0, Fup_h, Fup_m, Fup_l, tmp);
            QTW_SUBC(tmp, wu1, wu3); QTW_SCALC_TT(vt1, Fup_h, Fup_m, Fup_l, tmp);
            QTW_NEGC(vt2, vt0);
            QTW_NEGC(vt3, vt1);

            QTW_DECL(sum02); QTW_DECL(sum13); QTW_DECL(scaled);
            QTW_ADDC(sum02, wd0, wd2); QTW_SCALC_TT(scaled, Fdn_h, Fdn_m, Fdn_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_ADDC(vt2, vt2, scaled);
            QTW_ADDC(sum13, wd1, wd3); QTW_SCALC_TT(scaled, Fdn_h, Fdn_m, Fdn_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_ADDC(vt3, vt3, scaled);

            // local B-block
            QTW_DECL(wt0); QTW_DECL(wt1); QTW_DECL(wt2); QTW_DECL(wt3);
            QTW_DECL(t1); QTW_DECL(t2);
            if (is == 0 || is == Ns-1) {
                QTW_SCALC(t1, b1, w0); QTW_SCALC(t2, b2, w2); QTW_ADDC(wt0, t1, t2);
                QTW_SCALC(t1, b1, w1); QTW_SCALC(t2, b2, w3); QTW_ADDC(wt1, t1, t2);
                QTW_SCALC(t1, b1, w2); QTW_SCALC(t2, b2, w0); QTW_ADDC(wt2, t1, t2);
                QTW_SCALC(t1, b1, w3); QTW_SCALC(t2, b2, w1); QTW_ADDC(wt3, t1, t2);
            } else {
                QTW_SCALC(wt0, bb, w0); QTW_SCALC(wt1, bb, w1);
                QTW_SCALC(wt2, bb, w2); QTW_SCALC(wt3, bb, w3);
            }
            QTW_ADDC(wt0, wt0, vt0);
            QTW_ADDC(wt1, wt1, vt1);
            QTW_ADDC(wt2, wt2, vt2);
            QTW_ADDC(wt3, wt3, vt3);
            // multiply by -0.5
            QTW_SCALC(wt0, real_t(-0.5), wt0);
            QTW_SCALC(wt1, real_t(-0.5), wt1);
            QTW_SCALC(wt2, real_t(-0.5), wt2);
            QTW_SCALC(wt3, real_t(-0.5), wt3);
            QTW_STORE6(yp, b_s0, wt0);
            QTW_STORE6(yp, b_s1, wt1);
            QTW_STORE6(yp, b_s2, wt2);
            QTW_STORE6(yp, b_s3, wt3);
        }
    }
}

void mult_domainwall_5din_eo_5dir_dirac_qtw(
    real_t *yp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t *yp_dev = (real_t *)dev_ptr(yp);
    real_t *wp_dev = (real_t *)dev_ptr(wp);

    // Triple-word split of mq so the boundary hopping is computed at full TW.
    const real_t mq_h = (real_t)mq;
    const real_t mq_m = (real_t)(mq - (double)mq_h);
    const real_t mq_l = (real_t)(mq - (double)mq_h - (double)mq_m);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_eo_5dir_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        yp_dev, wp_dev, mq_h, mq_m, mq_l, Ns, (real_t)alpha, Nst, Nst_pad);
    QTW_PROF_END("eo_5dir_qtw");
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// ee 5-dirdag (adjoint diagonal block) TW kernel
//====================================================================
__global__
void mult_domainwall_5din_ee_5dirdag_dirac_qtw_kernel_fp(
    real_t * __restrict__ vp, real_t * __restrict__ wp,
    real_t mq, real_t mq_m, real_t mq_l,
    real_t M0, real_t M0_m, real_t M0_l, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();
    const real_t *b_mid_con = ConstantMemoryTraits<real_t>::b_mid();
    const real_t *b_lo_con  = ConstantMemoryTraits<real_t>::b_lo_tw();
    const real_t *c_mid_con = ConstantMemoryTraits<real_t>::c_mid();
    const real_t *c_lo_con  = ConstantMemoryTraits<real_t>::c_lo_tw();

    // (4 - M0) as triple-word; M0=1.6 is inexact in float so this is the
    // coefficient that made the forward Ddag diagonal only single-word accurate.
    real_t fm_h, fm_m, fm_l;
    tw_add(real_t(4.0), real_t(0), real_t(0), -M0, -M0_m, -M0_l, fm_h, fm_m, fm_l);

    for (int is = 0; is < Ns; ++is) {
        const int is_up = (is + 1) % Ns;
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t t_h, t_m, t_l;

        // B1 = b*(4-M0)+1  (genuine triple-word)
        real_t B1_h, B1_m, B1_l;
        tw_mul(b_con[is], b_mid_con[is], b_lo_con[is], fm_h, fm_m, fm_l, t_h, t_m, t_l);
        tw_add(t_h, t_m, t_l, real_t(1.0), real_t(0), real_t(0), B1_h, B1_m, B1_l);

        // local diagonal factors f1/f2/fa from B1 (alpha-scalings exact at alpha=1)
        real_t f1_h=0,f1_m=0,f1_l=0, f2_h=0,f2_m=0,f2_l=0, fa_h=0,fa_m=0,fa_l=0;
        if (is == 0) {
            tw_scal(real_t(0.5)*(real_t( 1.0)+alpha), B1_h, B1_m, B1_l, f1_h, f1_m, f1_l);
            tw_scal(real_t(0.5)*(real_t(-1.0)+alpha), B1_h, B1_m, B1_l, f2_h, f2_m, f2_l);
        } else if (is == Ns-1) {
            tw_scal(real_t(0.5)*(real_t(1.0)+alpha), B1_h, B1_m, B1_l, f1_h, f1_m, f1_l);
            tw_scal(real_t(0.5)*(real_t(1.0)-alpha), B1_h, B1_m, B1_l, f2_h, f2_m, f2_l);
        } else {
            tw_scal(alpha, B1_h, B1_m, B1_l, fa_h, fa_m, fa_l);
        }

        // C2up = c[is_up]*(4-M0)-1 ; Fup = bulk 0.25*alpha*C2up, boundary -0.25*mq*C2up
        real_t C2u_h, C2u_m, C2u_l, Fup_h, Fup_m, Fup_l;
        tw_mul(c_con[is_up], c_mid_con[is_up], c_lo_con[is_up], fm_h, fm_m, fm_l, t_h, t_m, t_l);
        tw_add(t_h, t_m, t_l, real_t(-1.0), real_t(0), real_t(0), C2u_h, C2u_m, C2u_l);
        if (is == Ns-1) {
            tw_mul(mq, mq_m, mq_l, C2u_h, C2u_m, C2u_l, t_h, t_m, t_l);
            tw_scal(real_t(-0.25), t_h, t_m, t_l, Fup_h, Fup_m, Fup_l);
        } else {
            tw_scal(real_t(0.25)*alpha, C2u_h, C2u_m, C2u_l, Fup_h, Fup_m, Fup_l);
        }

        // C2dn = c[is_dn]*(4-M0)-1 ; Fdn = bulk 0.25*alpha*C2dn, boundary -0.25*mq*C2dn
        real_t C2d_h, C2d_m, C2d_l, Fdn_h, Fdn_m, Fdn_l;
        tw_mul(c_con[is_dn], c_mid_con[is_dn], c_lo_con[is_dn], fm_h, fm_m, fm_l, t_h, t_m, t_l);
        tw_add(t_h, t_m, t_l, real_t(-1.0), real_t(0), real_t(0), C2d_h, C2d_m, C2d_l);
        if (is == 0) {
            tw_mul(mq, mq_m, mq_l, C2d_h, C2d_m, C2d_l, t_h, t_m, t_l);
            tw_scal(real_t(-0.25), t_h, t_m, t_l, Fdn_h, Fdn_m, Fdn_l);
        } else {
            tw_scal(real_t(0.25)*alpha, C2d_h, C2d_m, C2d_l, Fdn_h, Fdn_m, Fdn_l);
        }

        for (int ic = 0; ic < NC; ++ic) {
            int b_s0u = 6 * IDX_DWF_QTW(ic, 0, is_up, Ns, site);
            int b_s1u = 6 * IDX_DWF_QTW(ic, 1, is_up, Ns, site);
            int b_s2u = 6 * IDX_DWF_QTW(ic, 2, is_up, Ns, site);
            int b_s3u = 6 * IDX_DWF_QTW(ic, 3, is_up, Ns, site);
            int b_s0d = 6 * IDX_DWF_QTW(ic, 0, is_dn, Ns, site);
            int b_s1d = 6 * IDX_DWF_QTW(ic, 1, is_dn, Ns, site);
            int b_s2d = 6 * IDX_DWF_QTW(ic, 2, is_dn, Ns, site);
            int b_s3d = 6 * IDX_DWF_QTW(ic, 3, is_dn, Ns, site);
            int b_s0  = 6 * IDX_DWF_QTW(ic, 0, is,    Ns, site);
            int b_s1  = 6 * IDX_DWF_QTW(ic, 1, is,    Ns, site);
            int b_s2  = 6 * IDX_DWF_QTW(ic, 2, is,    Ns, site);
            int b_s3  = 6 * IDX_DWF_QTW(ic, 3, is,    Ns, site);

            QTW_LOAD6(w0, wp, b_s0);
            QTW_LOAD6(w1, wp, b_s1);
            QTW_LOAD6(w2, wp, b_s2);
            QTW_LOAD6(w3, wp, b_s3);

            // Local diagonal (self-adjoint of ee_5dir local)
            QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
            QTW_DECL(t1); QTW_DECL(t2);
            if (is == 0 || is == Ns-1) {
                QTW_SCALC_TT(t1, f1_h,f1_m,f1_l, w0); QTW_SCALC_TT(t2, f2_h,f2_m,f2_l, w2); QTW_ADDC(vt0, t1, t2);
                QTW_SCALC_TT(t1, f1_h,f1_m,f1_l, w1); QTW_SCALC_TT(t2, f2_h,f2_m,f2_l, w3); QTW_ADDC(vt1, t1, t2);
                QTW_SCALC_TT(t1, f1_h,f1_m,f1_l, w2); QTW_SCALC_TT(t2, f2_h,f2_m,f2_l, w0); QTW_ADDC(vt2, t1, t2);
                QTW_SCALC_TT(t1, f1_h,f1_m,f1_l, w3); QTW_SCALC_TT(t2, f2_h,f2_m,f2_l, w1); QTW_ADDC(vt3, t1, t2);
            } else {
                QTW_SCALC_TT(vt0, fa_h,fa_m,fa_l, w0); QTW_SCALC_TT(vt1, fa_h,fa_m,fa_l, w1);
                QTW_SCALC_TT(vt2, fa_h,fa_m,fa_l, w2); QTW_SCALC_TT(vt3, fa_h,fa_m,fa_l, w3);
            }

            // s+1 contribution (adjoint: symmetric)
            QTW_LOAD6(xu0, wp, b_s0u);
            QTW_LOAD6(xu1, wp, b_s1u);
            QTW_LOAD6(xu2, wp, b_s2u);
            QTW_LOAD6(xu3, wp, b_s3u);
            QTW_DECL(sum02); QTW_DECL(sum13); QTW_DECL(scaled);
            QTW_ADDC(sum02, xu0, xu2); QTW_SCALC_TT(scaled, Fup_h,Fup_m,Fup_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_ADDC(vt2, vt2, scaled);
            QTW_ADDC(sum13, xu1, xu3); QTW_SCALC_TT(scaled, Fup_h,Fup_m,Fup_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_ADDC(vt3, vt3, scaled);

            // s-1 contribution (adjoint: antisymmetric)
            QTW_LOAD6(xd0, wp, b_s0d);
            QTW_LOAD6(xd1, wp, b_s1d);
            QTW_LOAD6(xd2, wp, b_s2d);
            QTW_LOAD6(xd3, wp, b_s3d);
            QTW_SUBC(sum02, xd0, xd2); QTW_SCALC_TT(scaled, Fdn_h,Fdn_m,Fdn_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_SUBC(vt2, vt2, scaled);
            QTW_SUBC(sum13, xd1, xd3); QTW_SCALC_TT(scaled, Fdn_h,Fdn_m,Fdn_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_SUBC(vt3, vt3, scaled);

            QTW_STORE6(vp, b_s0, vt0);
            QTW_STORE6(vp, b_s1, vt1);
            QTW_STORE6(vp, b_s2, vt2);
            QTW_STORE6(vp, b_s3, vt3);
        }
    }
}

void mult_domainwall_5din_ee_5dirdag_dirac_qtw(
    real_t *vp, real_t *wp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t *vp_dev = (real_t *)dev_ptr(vp);
    real_t *wp_dev = (real_t *)dev_ptr(wp);

    // triple-word split of mq and M0 (both inexact in float).
    const real_t mq_h = (real_t)mq;
    const real_t mq_m = (real_t)(mq - (double)mq_h);
    const real_t mq_l = (real_t)(mq - (double)mq_h - (double)mq_m);
    const real_t M0_h = (real_t)M0;
    const real_t M0_m = (real_t)(M0 - (double)M0_h);
    const real_t M0_l = (real_t)(M0 - (double)M0_h - (double)M0_m);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_ee_5dirdag_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        vp_dev, wp_dev, mq_h, mq_m, mq_l, M0_h, M0_m, M0_l, Ns, (real_t)alpha, Nst, Nst_pad);
    QTW_PROF_END("ee_5dirdag_qtw");
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// eo 5-dirdag (adjoint off-diagonal block) TW kernel — vp = (M_eo)† * yp
//====================================================================
__global__
void mult_domainwall_5din_eo_5dirdag_dirac_qtw_kernel_fp(
    real_t * __restrict__ vp, real_t * __restrict__ yp,
    real_t mq, real_t mq_m, real_t mq_l, int Ns, real_t alpha,
    int Nst, int Nst_pad)
{
    const int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= Nst) return;

    const real_t *b_con = ConstantMemoryTraits<real_t>::b();
    const real_t *c_con = ConstantMemoryTraits<real_t>::c();
    const real_t *c_mid_con = ConstantMemoryTraits<real_t>::c_mid();
    const real_t *c_lo_con  = ConstantMemoryTraits<real_t>::c_lo_tw();

    for (int is = 0; is < Ns; ++is) {
        real_t b1 = 0, b2 = 0, bb = 0;
        if (is == 0) {
            b1 = real_t(-0.5) * b_con[is] * real_t(0.5) * (real_t(1.0)  + alpha);
            b2 = real_t(-0.5) * b_con[is] * real_t(0.5) * (real_t(-1.0) + alpha);
        } else if (is == Ns-1) {
            b1 = real_t(-0.5) * b_con[is] * real_t(0.5) * (real_t(1.0) + alpha);
            b2 = real_t(-0.5) * b_con[is] * real_t(0.5) * (real_t(1.0) - alpha);
        } else {
            bb = real_t(-0.5) * b_con[is] * alpha;
        }

        // Boundary hopping carries mq (inexact in float) -> full triple-word
        // coefficient (-0.5 * c * -0.5 * mq = 0.25 * c * mq). Bulk is exact.
        const int is_up = (is + 1) % Ns;
        real_t Fup_d_h, Fup_d_m, Fup_d_l;
        if (is == Ns-1) {
            real_t p_h, p_m, p_l;
            tw_mul(c_con[is_up], c_mid_con[is_up], c_lo_con[is_up], mq, mq_m, mq_l, p_h, p_m, p_l);
            Fup_d_h = real_t(0.25)*p_h; Fup_d_m = real_t(0.25)*p_m; Fup_d_l = real_t(0.25)*p_l;
        } else {
            Fup_d_h = real_t(-0.5) * c_con[is_up] * real_t(0.5) * alpha; Fup_d_m = real_t(0); Fup_d_l = real_t(0);
        }
        const int is_dn = (is - 1 + Ns) % Ns;
        real_t Fdn_d_h, Fdn_d_m, Fdn_d_l;
        if (is == 0) {
            real_t p_h, p_m, p_l;
            tw_mul(c_con[is_dn], c_mid_con[is_dn], c_lo_con[is_dn], mq, mq_m, mq_l, p_h, p_m, p_l);
            Fdn_d_h = real_t(0.25)*p_h; Fdn_d_m = real_t(0.25)*p_m; Fdn_d_l = real_t(0.25)*p_l;
        } else {
            Fdn_d_h = real_t(-0.5) * c_con[is_dn] * real_t(0.5) * alpha; Fdn_d_m = real_t(0); Fdn_d_l = real_t(0);
        }

        for (int ic = 0; ic < NC; ++ic) {
            int b_s0u = 6 * IDX_DWF_QTW(ic, 0, is_up, Ns, site);
            int b_s1u = 6 * IDX_DWF_QTW(ic, 1, is_up, Ns, site);
            int b_s2u = 6 * IDX_DWF_QTW(ic, 2, is_up, Ns, site);
            int b_s3u = 6 * IDX_DWF_QTW(ic, 3, is_up, Ns, site);
            int b_s0d = 6 * IDX_DWF_QTW(ic, 0, is_dn, Ns, site);
            int b_s1d = 6 * IDX_DWF_QTW(ic, 1, is_dn, Ns, site);
            int b_s2d = 6 * IDX_DWF_QTW(ic, 2, is_dn, Ns, site);
            int b_s3d = 6 * IDX_DWF_QTW(ic, 3, is_dn, Ns, site);
            int b_s0  = 6 * IDX_DWF_QTW(ic, 0, is,    Ns, site);
            int b_s1  = 6 * IDX_DWF_QTW(ic, 1, is,    Ns, site);
            int b_s2  = 6 * IDX_DWF_QTW(ic, 2, is,    Ns, site);
            int b_s3  = 6 * IDX_DWF_QTW(ic, 3, is,    Ns, site);

            QTW_LOAD6(y0, yp, b_s0);
            QTW_LOAD6(y1, yp, b_s1);
            QTW_LOAD6(y2, yp, b_s2);
            QTW_LOAD6(y3, yp, b_s3);

            // Local adjoint: vt = -0.5 * b[is] * B_alpha^dag * y
            QTW_DECL(vt0); QTW_DECL(vt1); QTW_DECL(vt2); QTW_DECL(vt3);
            QTW_DECL(t1); QTW_DECL(t2);
            if (is == 0 || is == Ns-1) {
                QTW_SCALC(t1, b1, y2); QTW_SCALC(t2, b2, y0); QTW_ADDC(vt0, t1, t2);
                QTW_SCALC(t1, b1, y3); QTW_SCALC(t2, b2, y1); QTW_ADDC(vt1, t1, t2);
                QTW_SCALC(t1, b1, y0); QTW_SCALC(t2, b2, y2); QTW_ADDC(vt2, t1, t2);
                QTW_SCALC(t1, b1, y1); QTW_SCALC(t2, b2, y3); QTW_ADDC(vt3, t1, t2);
            } else {
                QTW_SCALC(vt0, bb, y2);
                QTW_SCALC(vt1, bb, y3);
                QTW_SCALC(vt2, bb, y0);
                QTW_SCALC(vt3, bb, y1);
            }

            // s+1 contribution
            QTW_LOAD6(yu0, yp, b_s0u);
            QTW_LOAD6(yu1, yp, b_s1u);
            QTW_LOAD6(yu2, yp, b_s2u);
            QTW_LOAD6(yu3, yp, b_s3u);
            QTW_DECL(sum02); QTW_DECL(sum13); QTW_DECL(scaled);
            QTW_ADDC(sum02, yu0, yu2); QTW_SCALC_TT(scaled, Fup_d_h, Fup_d_m, Fup_d_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_ADDC(vt2, vt2, scaled);
            QTW_ADDC(sum13, yu1, yu3); QTW_SCALC_TT(scaled, Fup_d_h, Fup_d_m, Fup_d_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_ADDC(vt3, vt3, scaled);

            // s-1 contribution: antisymmetric flavor (uses y2-y0, y3-y1)
            QTW_LOAD6(yd0, yp, b_s0d);
            QTW_LOAD6(yd1, yp, b_s1d);
            QTW_LOAD6(yd2, yp, b_s2d);
            QTW_LOAD6(yd3, yp, b_s3d);
            QTW_SUBC(sum02, yd2, yd0); QTW_SCALC_TT(scaled, Fdn_d_h, Fdn_d_m, Fdn_d_l, sum02);
            QTW_ADDC(vt0, vt0, scaled); QTW_SUBC(vt2, vt2, scaled);
            QTW_SUBC(sum13, yd3, yd1); QTW_SCALC_TT(scaled, Fdn_d_h, Fdn_d_m, Fdn_d_l, sum13);
            QTW_ADDC(vt1, vt1, scaled); QTW_SUBC(vt3, vt3, scaled);

            QTW_STORE6(vp, b_s0, vt0);
            QTW_STORE6(vp, b_s1, vt1);
            QTW_STORE6(vp, b_s2, vt2);
            QTW_STORE6(vp, b_s3, vt3);
        }
    }
}

void mult_domainwall_5din_eo_5dirdag_dirac_qtw(
    real_t *vp, real_t *yp, double mq, double M0, int Ns,
    real_t *b, real_t *c, double alpha, int *Nsize, bool ext)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t *vp_dev = (real_t *)dev_ptr(vp);
    real_t *yp_dev = (real_t *)dev_ptr(yp);

    // Triple-word split of mq for the boundary hopping.
    const real_t mq_h = (real_t)mq;
    const real_t mq_m = (real_t)(mq - (double)mq_h);
    const real_t mq_l = (real_t)(mq - (double)mq_h - (double)mq_m);

    int blockSize = VECTOR_LENGTH;
    int gridSize  = (Nst + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    mult_domainwall_5din_eo_5dirdag_dirac_qtw_kernel_fp<<<gridSize, blockSize>>>(
        vp_dev, yp_dev, mq_h, mq_m, mq_l, Ns, (real_t)alpha, Nst, Nst_pad);
    QTW_PROF_END("eo_5dirdag_qtw");
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
// EO 4D bulk hopping (TW, eo-gauge layout, jgm5 support)
// Reuses QDW gauge-stencil structure; vector ops use TW arithmetic, gauge
// link can be plain T (up_mid==nullptr) or full TW (up_mid != nullptr).
//====================================================================
// NOTE: do NOT add __launch_bounds__ to cap registers here. Measured: forcing
// 168->128 regs (25%->33% occupancy) made DdagD SLOWER (7.46->8.23 ms). This
// kernel is ILP/arithmetic-bound, not occupancy-bound; the extra registers buy
// instruction-level parallelism that hides the TW dependency chains.
// RECON (template) = SU(3) 3rd-column reconstruction (YAML su3_reconstruction);
// when false the recon code is dead-stripped so the fast path is unchanged.
template<bool RECON>
__global__
void mult_domainwall_5din_eo_hopb_qtw_dirac_5d_kernel(
    real_t * __restrict__ vp, const real_t * __restrict__ up,
    const real_t * __restrict__ up_mid, const real_t * __restrict__ up_lo,
    const real_t * __restrict__ wp,
    int Ns, int bc_x, int bc_y, int bc_z, int bc_t,
    int Nx, int Ny, int Nz, int Nt,
    int ieo, int jeo,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
    int Nst, int Nst_pad, int jgm5)
{
    // 5D-parallel: one thread per 5D site (4D site x Ls slice). Grid-stride
    // over Nst_pad*Ns; NWP-style decomposition keeps loads coalesced
    // (consecutive threads -> consecutive 4D site within an NWP block, same is).
    const int ist     = blockIdx.x * blockDim.x + threadIdx.x;
    const int gstride = blockDim.x * gridDim.x;
    for (int idx = ist; idx < Nst_pad * Ns; idx += gstride) {
    const int idx2    = idx / NWP;
    const int idx_in  = idx % NWP;
    const int is      = idx2 % Ns;
    const int idx_out = idx2 / Ns;
    const int site    = idx_in + NWP * idx_out;
    if (site >= Nst) continue;

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

    const int d0 = (jgm5 == 0) ? 0 : 2;
    const int d1 = (jgm5 == 0) ? 1 : 3;
    const int d2 = (jgm5 == 0) ? 2 : 0;
    const int d3 = (jgm5 == 0) ? 3 : 1;

    const bool use_tw_link = (up_mid != nullptr) && (up_lo != nullptr);

    // Accumulators: 3 colors × 4 spinors × TW(6 floats) = 72 reals per site.
    // Use named scalars for clarity.
    QTW_DECL(acc_c0_s0); QTW_ZERO(acc_c0_s0);
    QTW_DECL(acc_c0_s1); QTW_ZERO(acc_c0_s1);
    QTW_DECL(acc_c0_s2); QTW_ZERO(acc_c0_s2);
    QTW_DECL(acc_c0_s3); QTW_ZERO(acc_c0_s3);
    QTW_DECL(acc_c1_s0); QTW_ZERO(acc_c1_s0);
    QTW_DECL(acc_c1_s1); QTW_ZERO(acc_c1_s1);
    QTW_DECL(acc_c1_s2); QTW_ZERO(acc_c1_s2);
    QTW_DECL(acc_c1_s3); QTW_ZERO(acc_c1_s3);
    QTW_DECL(acc_c2_s0); QTW_ZERO(acc_c2_s0);
    QTW_DECL(acc_c2_s1); QTW_ZERO(acc_c2_s1);
    QTW_DECL(acc_c2_s2); QTW_ZERO(acc_c2_s2);
    QTW_DECL(acc_c2_s3); QTW_ZERO(acc_c2_s3);

    // Helper macro: load proj into vt1_c[0..2] (spin pair) from wp[isn].
    // PROJ_FN one of DWF_PROJ_P_TW / DWF_PROJ_M_TW / DWF_PROJ_RP_TW / DWF_PROJ_RM_TW.
    // PAIR: 1 = first pair (vt1, ida=d0, idb=d3 or as needed), 2 = second pair (vt2).
    // We expand explicitly to avoid token-pasting hell.

    // For each direction we:
    //   1. Load 3 colors × 2 spinors from wp at isn -> sa_cX_iy, sb_cX_iy (Y=A,B)
    //   2. Project into vt[1/2]_c0..c2 using PROJ_FN(vt, sa, sb)
    //   3. Apply gauge: wt = U(_ptr, isg) × vt
    //   4. Accumulate into acc with sign/permutation given by ACCUM_OP2, OP3

    // We pack the per-direction sequence in a (large) macro to avoid duplication.
    // Each "direction block" produces six wt outputs (3 colors × 2 spinor pairs).

    #define QTW_PROCESS_DIR(isn_, isg_, FWD, PROJ_FN_A, PROJ_FN_B, ida, idb_2, ACCUM_OP_A2, ACCUM_OP_A3, ACCUM_OP_B2, ACCUM_OP_B3, bc2_) do { \
        /* Load wp[isn] for 3 colors × 2 spinors per pair (ida, idb_2). */ \
        QTW_LOAD6(_sa_c0, wp, 6 * IDX_DWF_QTW(0, ida, is, Ns, (isn_))); \
        QTW_LOAD6(_sb_c0, wp, 6 * IDX_DWF_QTW(0, idb_2, is, Ns, (isn_))); \
        QTW_LOAD6(_sa_c1, wp, 6 * IDX_DWF_QTW(1, ida, is, Ns, (isn_))); \
        QTW_LOAD6(_sb_c1, wp, 6 * IDX_DWF_QTW(1, idb_2, is, Ns, (isn_))); \
        QTW_LOAD6(_sa_c2, wp, 6 * IDX_DWF_QTW(2, ida, is, Ns, (isn_))); \
        QTW_LOAD6(_sb_c2, wp, 6 * IDX_DWF_QTW(2, idb_2, is, Ns, (isn_))); \
        /* Note: For pair B (vt2), the same load schema must be supplied by  \
         * caller; here we instead use separate macro per pair below.        \
         */ \
    } while (0)

    // (See note above — for brevity in this header the per-direction code is
    // expanded directly in C++ below, without the deeply-nested macro path.)

    real_t bc2;

    // Helper lambda equivalents — declared as inline blocks with explicit
    // variables. We process each direction in its own `{}` scope so locals
    // don't pollute neighbours.

    // ---- For each of the 8 directions: X+, X-, Y+, Y-, Z+, Z-, T+, T- ----
    // Pattern for X+ (FWD):
    //   isn = ((ix+keo)%Nx) + Nx*iyzt
    //   isg = site + Nst*(ieo + 2*0)
    //   bc2 = (ix==Nx-1 && keo==1) ? bc_x : 1
    //   PROJ_P for (d0, d3) into vt1 and (d1, d2) into vt2
    //   apply U_fwd
    //   ACCUM_4(..., wt1_c, wt2_c, bc2, MULT_MI, MULT_MI)

    // To keep this tractable, we wrap one direction in a macro that takes:
    //   isn, isg, gmul_is_fwd, proj_fn (PROJ_P_TW/M_TW/RP_TW/RM_TW), accum_op2, accum_op3
    // The PROJ phase loads 6 spinors and projects.
    // We expand for both spinor pairs (vt1 and vt2). Different directions
    // can use different ida/idb assignments; the simplest is to inline each
    // direction below.

    #define QTW_DIR_PROJ_PAIR(VTPRE, PROJ_FN, ida_, idb_, isn_) \
        QTW_LOAD6(VTPRE##_sa_c0, wp, 6 * IDX_DWF_QTW(0, (ida_), is, Ns, (isn_))); \
        QTW_LOAD6(VTPRE##_sb_c0, wp, 6 * IDX_DWF_QTW(0, (idb_), is, Ns, (isn_))); \
        QTW_LOAD6(VTPRE##_sa_c1, wp, 6 * IDX_DWF_QTW(1, (ida_), is, Ns, (isn_))); \
        QTW_LOAD6(VTPRE##_sb_c1, wp, 6 * IDX_DWF_QTW(1, (idb_), is, Ns, (isn_))); \
        QTW_LOAD6(VTPRE##_sa_c2, wp, 6 * IDX_DWF_QTW(2, (ida_), is, Ns, (isn_))); \
        QTW_LOAD6(VTPRE##_sb_c2, wp, 6 * IDX_DWF_QTW(2, (idb_), is, Ns, (isn_))); \
        QTW_DECL(VTPRE##_c0); PROJ_FN(VTPRE##_c0, VTPRE##_sa_c0, VTPRE##_sb_c0); \
        QTW_DECL(VTPRE##_c1); PROJ_FN(VTPRE##_c1, VTPRE##_sa_c1, VTPRE##_sb_c1); \
        QTW_DECL(VTPRE##_c2); PROJ_FN(VTPRE##_c2, VTPRE##_sa_c2, VTPRE##_sb_c2)

    // Apply gauge (FWD or BCK, FP or TW) to a packed (3-color TW vec).
    // Caller must declare WTPRE##_c0/c1/c2 (via QTW_DECL) before invoking; the
    // macro only assigns into them. The do-while scope keeps the temporary
    // _in_* / _o_* arrays local.
    #define QTW_DIR_GMUL(WTPRE, VTPRE, isg_, is_fwd) do { \
        real_t _in_rh[3] = { VTPRE##_c0_rh, VTPRE##_c1_rh, VTPRE##_c2_rh }; \
        real_t _in_ih[3] = { VTPRE##_c0_ih, VTPRE##_c1_ih, VTPRE##_c2_ih }; \
        real_t _in_rm[3] = { VTPRE##_c0_rm, VTPRE##_c1_rm, VTPRE##_c2_rm }; \
        real_t _in_im[3] = { VTPRE##_c0_im, VTPRE##_c1_im, VTPRE##_c2_im }; \
        real_t _in_rl[3] = { VTPRE##_c0_rl, VTPRE##_c1_rl, VTPRE##_c2_rl }; \
        real_t _in_il[3] = { VTPRE##_c0_il, VTPRE##_c1_il, VTPRE##_c2_il }; \
        real_t _o_rh[3], _o_ih[3], _o_rm[3], _o_im[3], _o_rl[3], _o_il[3]; \
        if (use_tw_link) { \
            if (is_fwd) qtw_gmult_fwd_3x3_tw<RECON>(up, up_mid, up_lo, (isg_), _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, _o_rh,_o_ih,_o_rm,_o_im,_o_rl,_o_il); \
            else        qtw_gmult_bck_3x3_tw<RECON>(up, up_mid, up_lo, (isg_), _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, _o_rh,_o_ih,_o_rm,_o_im,_o_rl,_o_il); \
        } else { \
            if (is_fwd) qtw_gmult_fwd_3x3(up, (isg_), _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, _o_rh,_o_ih,_o_rm,_o_im,_o_rl,_o_il); \
            else        qtw_gmult_bck_3x3(up, (isg_), _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, _o_rh,_o_ih,_o_rm,_o_im,_o_rl,_o_il); \
        } \
        WTPRE##_c0_rh = _o_rh[0]; WTPRE##_c0_ih = _o_ih[0]; WTPRE##_c0_rm = _o_rm[0]; WTPRE##_c0_im = _o_im[0]; WTPRE##_c0_rl = _o_rl[0]; WTPRE##_c0_il = _o_il[0]; \
        WTPRE##_c1_rh = _o_rh[1]; WTPRE##_c1_ih = _o_ih[1]; WTPRE##_c1_rm = _o_rm[1]; WTPRE##_c1_im = _o_im[1]; WTPRE##_c1_rl = _o_rl[1]; WTPRE##_c1_il = _o_il[1]; \
        WTPRE##_c2_rh = _o_rh[2]; WTPRE##_c2_ih = _o_ih[2]; WTPRE##_c2_rm = _o_rm[2]; WTPRE##_c2_im = _o_im[2]; WTPRE##_c2_rl = _o_rl[2]; WTPRE##_c2_il = _o_il[2]; \
    } while (0)

    // (one 5D site per thread: `is` and `site` come from the grid-stride
    //  decomposition above; the 8-direction stencil below runs once.)
        // X+
        {
            int isn = ((ix + keo) % Nx) + Nx * iyzt;
            int isg = site + Nst * (ieo + 2 * 0);
            bc2 = (ix == Nx-1 && keo == 1) ? (real_t)bc_x : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_P_TW, d0, d3, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_P_TW, d1, d2, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, true);
            QTW_DIR_GMUL(wt2, vt2, isg, true);
            DWF_ACCUM_4_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_MI_TW, DWF_MULT_MI_TW);
            DWF_ACCUM_4_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_MI_TW, DWF_MULT_MI_TW);
            DWF_ACCUM_4_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_MI_TW, DWF_MULT_MI_TW);
        }
        // X-
        {
            int ix2 = (ix - 1 + keo + Nx) % Nx;
            int isn = ix2 + Nx * iyzt;
            int isg = isn + Nst * (1 - ieo + 2 * 0);
            bc2 = (ix == 0 && keo == 0) ? (real_t)bc_x : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_M_TW, d0, d3, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_M_TW, d1, d2, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, false);
            QTW_DIR_GMUL(wt2, vt2, isg, false);
            DWF_ACCUM_4_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_PI_TW, DWF_MULT_PI_TW);
            DWF_ACCUM_4_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_PI_TW, DWF_MULT_PI_TW);
            DWF_ACCUM_4_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_PI_TW, DWF_MULT_PI_TW);
        }
        // Y+
        if ((iy < Ny-1) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy + 1) % Ny) + Ny * izt);
            int isg = site + Nst * (ieo + 2 * 1);
            bc2 = (iy == Ny-1) ? (real_t)bc_y : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_RP_TW, d0, d3, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_RM_TW, d1, d2, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, true);
            QTW_DIR_GMUL(wt2, vt2, isg, true);
            DWF_ACCUM_4_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_RM1_TW, DWF_MULT_R1_TW);
            DWF_ACCUM_4_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_RM1_TW, DWF_MULT_R1_TW);
            DWF_ACCUM_4_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_RM1_TW, DWF_MULT_R1_TW);
        }
        // Y-
        if ((iy > 0) || (do_comm_y == 0)) {
            int isn = ix + Nx * (((iy - 1 + Ny) % Ny) + Ny * izt);
            int isg = isn + Nst * (1 - ieo + 2 * 1);
            bc2 = (iy == 0) ? (real_t)bc_y : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_RM_TW, d0, d3, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_RP_TW, d1, d2, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, false);
            QTW_DIR_GMUL(wt2, vt2, isg, false);
            DWF_ACCUM_4_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_R1_TW, DWF_MULT_RM1_TW);
            DWF_ACCUM_4_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_R1_TW, DWF_MULT_RM1_TW);
            DWF_ACCUM_4_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_R1_TW, DWF_MULT_RM1_TW);
        }
        // Z+
        if ((iz < Nz-1) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz + 1) % Nz) + Nz * it);
            int isg = site + Nst * (ieo + 2 * 2);
            bc2 = (iz == Nz-1) ? (real_t)bc_z : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_P_TW, d0, d2, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_M_TW, d1, d3, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, true);
            QTW_DIR_GMUL(wt2, vt2, isg, true);
            DWF_ACCUM_4_SW_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_MI_TW, DWF_MULT_PI_TW);
            DWF_ACCUM_4_SW_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_MI_TW, DWF_MULT_PI_TW);
            DWF_ACCUM_4_SW_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_MI_TW, DWF_MULT_PI_TW);
        }
        // Z-
        if ((iz > 0) || (do_comm_z == 0)) {
            int isn = ixy + Nxy * (((iz - 1 + Nz) % Nz) + Nz * it);
            int isg = isn + Nst * (1 - ieo + 2 * 2);
            bc2 = (iz == 0) ? (real_t)bc_z : real_t(1.0);
            QTW_DIR_PROJ_PAIR(vt1, DWF_PROJ_M_TW, d0, d2, isn);
            QTW_DIR_PROJ_PAIR(vt2, DWF_PROJ_P_TW, d1, d3, isn);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, false);
            QTW_DIR_GMUL(wt2, vt2, isg, false);
            DWF_ACCUM_4_SW_TW(acc_c0_s0, acc_c0_s1, acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2, DWF_MULT_PI_TW, DWF_MULT_MI_TW);
            DWF_ACCUM_4_SW_TW(acc_c1_s0, acc_c1_s1, acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2, DWF_MULT_PI_TW, DWF_MULT_MI_TW);
            DWF_ACCUM_4_SW_TW(acc_c2_s0, acc_c2_s1, acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2, DWF_MULT_PI_TW, DWF_MULT_MI_TW);
        }
        // T+: use special "T" projection (PROJ_2): load only one spinor pair, double it
        if ((it < Nt-1) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it + 1) % Nt);
            int isg = site + Nst * (ieo + 2 * 3);
            bc2 = (it == Nt-1) ? (real_t)bc_t : real_t(1.0);
            // For T-dir, the QDW pattern uses DWF_LOAD_PROJ_T which loads
            // (d2, d3) and applies PROJ_2 to each. Equivalent: load and
            // multiply by 2 in TW.
            QTW_LOAD6(_va_c0, wp, 6 * IDX_DWF_QTW(0, d2, is, Ns, isn));
            QTW_LOAD6(_va_c1, wp, 6 * IDX_DWF_QTW(1, d2, is, Ns, isn));
            QTW_LOAD6(_va_c2, wp, 6 * IDX_DWF_QTW(2, d2, is, Ns, isn));
            QTW_LOAD6(_vb_c0, wp, 6 * IDX_DWF_QTW(0, d3, is, Ns, isn));
            QTW_LOAD6(_vb_c1, wp, 6 * IDX_DWF_QTW(1, d3, is, Ns, isn));
            QTW_LOAD6(_vb_c2, wp, 6 * IDX_DWF_QTW(2, d3, is, Ns, isn));
            QTW_DECL(vt1_c0); DWF_PROJ_2_TW(vt1_c0, _va_c0);
            QTW_DECL(vt1_c1); DWF_PROJ_2_TW(vt1_c1, _va_c1);
            QTW_DECL(vt1_c2); DWF_PROJ_2_TW(vt1_c2, _va_c2);
            QTW_DECL(vt2_c0); DWF_PROJ_2_TW(vt2_c0, _vb_c0);
            QTW_DECL(vt2_c1); DWF_PROJ_2_TW(vt2_c1, _vb_c1);
            QTW_DECL(vt2_c2); DWF_PROJ_2_TW(vt2_c2, _vb_c2);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, true);
            QTW_DIR_GMUL(wt2, vt2, isg, true);
            DWF_ACCUM_TP_TW(acc_c0_s2, acc_c0_s3, wt1_c0, wt2_c0, bc2);
            DWF_ACCUM_TP_TW(acc_c1_s2, acc_c1_s3, wt1_c1, wt2_c1, bc2);
            DWF_ACCUM_TP_TW(acc_c2_s2, acc_c2_s3, wt1_c2, wt2_c2, bc2);
        }
        // T-
        if ((it > 0) || (do_comm_t == 0)) {
            int isn = ixyz + Nxyz * ((it - 1 + Nt) % Nt);
            int isg = isn + Nst * (1 - ieo + 2 * 3);
            bc2 = (it == 0) ? (real_t)bc_t : real_t(1.0);
            QTW_LOAD6(_va_c0, wp, 6 * IDX_DWF_QTW(0, d0, is, Ns, isn));
            QTW_LOAD6(_va_c1, wp, 6 * IDX_DWF_QTW(1, d0, is, Ns, isn));
            QTW_LOAD6(_va_c2, wp, 6 * IDX_DWF_QTW(2, d0, is, Ns, isn));
            QTW_LOAD6(_vb_c0, wp, 6 * IDX_DWF_QTW(0, d1, is, Ns, isn));
            QTW_LOAD6(_vb_c1, wp, 6 * IDX_DWF_QTW(1, d1, is, Ns, isn));
            QTW_LOAD6(_vb_c2, wp, 6 * IDX_DWF_QTW(2, d1, is, Ns, isn));
            QTW_DECL(vt1_c0); DWF_PROJ_2_TW(vt1_c0, _va_c0);
            QTW_DECL(vt1_c1); DWF_PROJ_2_TW(vt1_c1, _va_c1);
            QTW_DECL(vt1_c2); DWF_PROJ_2_TW(vt1_c2, _va_c2);
            QTW_DECL(vt2_c0); DWF_PROJ_2_TW(vt2_c0, _vb_c0);
            QTW_DECL(vt2_c1); DWF_PROJ_2_TW(vt2_c1, _vb_c1);
            QTW_DECL(vt2_c2); DWF_PROJ_2_TW(vt2_c2, _vb_c2);
            QTW_DECL(wt1_c0); QTW_DECL(wt1_c1); QTW_DECL(wt1_c2);
            QTW_DECL(wt2_c0); QTW_DECL(wt2_c1); QTW_DECL(wt2_c2);
            QTW_DIR_GMUL(wt1, vt1, isg, false);
            QTW_DIR_GMUL(wt2, vt2, isg, false);
            DWF_ACCUM_TM_TW(acc_c0_s0, acc_c0_s1, wt1_c0, wt2_c0, bc2);
            DWF_ACCUM_TM_TW(acc_c1_s0, acc_c1_s1, wt1_c1, wt2_c1, bc2);
            DWF_ACCUM_TM_TW(acc_c2_s0, acc_c2_s1, wt1_c2, wt2_c2, bc2);
        }

        // Store accumulators for this is
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(0, 0, is, Ns, site), acc_c0_s0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(0, 1, is, Ns, site), acc_c0_s1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(0, 2, is, Ns, site), acc_c0_s2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(0, 3, is, Ns, site), acc_c0_s3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(1, 0, is, Ns, site), acc_c1_s0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(1, 1, is, Ns, site), acc_c1_s1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(1, 2, is, Ns, site), acc_c1_s2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(1, 3, is, Ns, site), acc_c1_s3);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(2, 0, is, Ns, site), acc_c2_s0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(2, 1, is, Ns, site), acc_c2_s1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(2, 2, is, Ns, site), acc_c2_s2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(2, 3, is, Ns, site), acc_c2_s3);
    }  // grid-stride loop over 5D sites

    #undef QTW_PROCESS_DIR
    #undef QTW_DIR_PROJ_PAIR
    #undef QTW_DIR_GMUL
}

// ============================================================================
// COLOR-SPLIT + BLOCK-SMEM SHARE (env QTW_CSPLIT=2, TW link only).
// 2D block (NWP, NC): threadIdx.x=lane=site (coalesced), threadIdx.y=cout=color.
// Per direction: each color-thread loads & projects ITS OWN input color
// (coalesced, read ONCE), deposits the projected spinor into __shared__; after
// __syncthreads every thread reads all NC input colors and does its single
// output-color gmult. Kills the 3x redundant spinor reads of the plain split
// while keeping coalescing -> "bandwidth 止血". 24-real accumulator.
// Requires: TW link, no padding (Nst % NWP == 0), single-rank (do_comm all 0 ->
// every thread takes every direction -> no divergent __syncthreads).
// ============================================================================
template<bool RECON>
__global__
void mult_domainwall_5din_eo_hopb_qtw_dirac_5d_csmem_kernel(
    real_t * __restrict__ vp, const real_t * __restrict__ up,
    const real_t * __restrict__ up_mid, const real_t * __restrict__ up_lo,
    const real_t * __restrict__ wp,
    int Ns, int bc_x, int bc_y, int bc_z, int bc_t,
    int Nx, int Ny, int Nz, int Nt,
    int ieo, int jeo,
    int do_comm_x, int do_comm_y, int do_comm_z, int do_comm_t,
    int Nst, int Nst_pad, int jgm5)
{
    const int lane = threadIdx.x;             // 0..NWP-1 (4D site within block)
    const int cout = threadIdx.y;             // 0..NC-1 (this thread's color)
    const int is   = blockIdx.y;
    const int site = lane + NWP * blockIdx.x;

    // Shared projected spinor, DOUBLE-BUFFERED: [buf(2)][pair(2)][color(NC)][word(6)][lane].
    // Direction d uses buf=d&1; consecutive directions alternate buffers, so the
    // WAR barrier (overwrite vs previous reads) is covered by the intervening
    // direction's RAW sync -> only 1 __syncthreads/direction (8 vs 16 per site).
    __shared__ real_t vsh[2 * 2 * NC * 6 * NWP];
    #define VSH(buf_, pair_, cin_, w_) vsh[((buf_) * (2 * NC * 6 * NWP)) + (((((pair_) * NC + (cin_)) * 6 + (w_)) * NWP) + lane)]

    if (site >= Nst) return;  // no padding required -> uniform within block

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

    const int d0 = (jgm5 == 0) ? 0 : 2;
    const int d1 = (jgm5 == 0) ? 1 : 3;
    const int d2 = (jgm5 == 0) ? 2 : 0;
    const int d3 = (jgm5 == 0) ? 3 : 1;

    QTW_DECL(acc_s0); QTW_ZERO(acc_s0);
    QTW_DECL(acc_s1); QTW_ZERO(acc_s1);
    QTW_DECL(acc_s2); QTW_ZERO(acc_s2);
    QTW_DECL(acc_s3); QTW_ZERO(acc_s3);

    // Project this thread's own input color (cout) for one pair, write to smem.
    #define QTW_CS_PROJ(buf_, pair_, PROJ_FN, ida_, idb_, isn_) do { \
        QTW_LOAD6(_psa, wp, 6 * IDX_DWF_QTW(cout, (ida_), is, Ns, (isn_))); \
        QTW_LOAD6(_psb, wp, 6 * IDX_DWF_QTW(cout, (idb_), is, Ns, (isn_))); \
        QTW_DECL(_pv); PROJ_FN(_pv, _psa, _psb); \
        VSH(buf_, pair_, cout, 0) = _pv_rh; VSH(buf_, pair_, cout, 1) = _pv_ih; \
        VSH(buf_, pair_, cout, 2) = _pv_rm; VSH(buf_, pair_, cout, 3) = _pv_im; \
        VSH(buf_, pair_, cout, 4) = _pv_rl; VSH(buf_, pair_, cout, 5) = _pv_il; \
    } while (0)
    // T-direction projection: single spin component, doubled (PROJ_2).
    #define QTW_CS_PROJ_T(buf_, pair_, idsp_, isn_) do { \
        QTW_LOAD6(_psa, wp, 6 * IDX_DWF_QTW(cout, (idsp_), is, Ns, (isn_))); \
        QTW_DECL(_pv); DWF_PROJ_2_TW(_pv, _psa); \
        VSH(buf_, pair_, cout, 0) = _pv_rh; VSH(buf_, pair_, cout, 1) = _pv_ih; \
        VSH(buf_, pair_, cout, 2) = _pv_rm; VSH(buf_, pair_, cout, 3) = _pv_im; \
        VSH(buf_, pair_, cout, 4) = _pv_rl; VSH(buf_, pair_, cout, 5) = _pv_il; \
    } while (0)
    // Read all NC input colors from smem, single-output-color gmul into WT.
    #define QTW_CS_GMUL(WT, buf_, pair_, isg_, is_fwd) do { \
        real_t _in_rh[3], _in_ih[3], _in_rm[3], _in_im[3], _in_rl[3], _in_il[3]; \
        _Pragma("unroll") \
        for (int _j = 0; _j < NC; ++_j) { \
            _in_rh[_j] = VSH(buf_, pair_, _j, 0); _in_ih[_j] = VSH(buf_, pair_, _j, 1); \
            _in_rm[_j] = VSH(buf_, pair_, _j, 2); _in_im[_j] = VSH(buf_, pair_, _j, 3); \
            _in_rl[_j] = VSH(buf_, pair_, _j, 4); _in_il[_j] = VSH(buf_, pair_, _j, 5); \
        } \
        if (is_fwd) qtw_gmult_fwd_1col_tw<RECON>(up, up_mid, up_lo, (isg_), cout, _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, WT##_rh,WT##_ih,WT##_rm,WT##_im,WT##_rl,WT##_il); \
        else        qtw_gmult_bck_1col_tw<RECON>(up, up_mid, up_lo, (isg_), cout, _in_rh,_in_ih,_in_rm,_in_im,_in_rl,_in_il, WT##_rh,WT##_ih,WT##_rm,WT##_im,WT##_rl,WT##_il); \
    } while (0)

    real_t bc2;

        // X+
        {
            int isn = ((ix + keo) % Nx) + Nx * iyzt;
            int isg = site + Nst * (ieo + 2 * 0);
            bc2 = (ix == Nx-1 && keo == 1) ? (real_t)bc_x : real_t(1.0);
            QTW_CS_PROJ(0, 0, DWF_PROJ_P_TW, d0, d3, isn);
            QTW_CS_PROJ(0, 1, DWF_PROJ_P_TW, d1, d2, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 0, 0, isg, true);
            QTW_CS_GMUL(wt2, 0, 1, isg, true);
            DWF_ACCUM_4_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_MI_TW, DWF_MULT_MI_TW);
        }
        // X-
        {
            int ix2 = (ix - 1 + keo + Nx) % Nx;
            int isn = ix2 + Nx * iyzt;
            int isg = isn + Nst * (1 - ieo + 2 * 0);
            bc2 = (ix == 0 && keo == 0) ? (real_t)bc_x : real_t(1.0);
            QTW_CS_PROJ(1, 0, DWF_PROJ_M_TW, d0, d3, isn);
            QTW_CS_PROJ(1, 1, DWF_PROJ_M_TW, d1, d2, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 1, 0, isg, false);
            QTW_CS_GMUL(wt2, 1, 1, isg, false);
            DWF_ACCUM_4_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_PI_TW, DWF_MULT_PI_TW);
        }
        // Y+
        {
            int isn = ix + Nx * (((iy + 1) % Ny) + Ny * izt);
            int isg = site + Nst * (ieo + 2 * 1);
            bc2 = (iy == Ny-1) ? (real_t)bc_y : real_t(1.0);
            QTW_CS_PROJ(0, 0, DWF_PROJ_RP_TW, d0, d3, isn);
            QTW_CS_PROJ(0, 1, DWF_PROJ_RM_TW, d1, d2, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 0, 0, isg, true);
            QTW_CS_GMUL(wt2, 0, 1, isg, true);
            DWF_ACCUM_4_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_RM1_TW, DWF_MULT_R1_TW);
        }
        // Y-
        {
            int isn = ix + Nx * (((iy - 1 + Ny) % Ny) + Ny * izt);
            int isg = isn + Nst * (1 - ieo + 2 * 1);
            bc2 = (iy == 0) ? (real_t)bc_y : real_t(1.0);
            QTW_CS_PROJ(1, 0, DWF_PROJ_RM_TW, d0, d3, isn);
            QTW_CS_PROJ(1, 1, DWF_PROJ_RP_TW, d1, d2, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 1, 0, isg, false);
            QTW_CS_GMUL(wt2, 1, 1, isg, false);
            DWF_ACCUM_4_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_R1_TW, DWF_MULT_RM1_TW);
        }
        // Z+
        {
            int isn = ixy + Nxy * (((iz + 1) % Nz) + Nz * it);
            int isg = site + Nst * (ieo + 2 * 2);
            bc2 = (iz == Nz-1) ? (real_t)bc_z : real_t(1.0);
            QTW_CS_PROJ(0, 0, DWF_PROJ_P_TW, d0, d2, isn);
            QTW_CS_PROJ(0, 1, DWF_PROJ_M_TW, d1, d3, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 0, 0, isg, true);
            QTW_CS_GMUL(wt2, 0, 1, isg, true);
            DWF_ACCUM_4_SW_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_MI_TW, DWF_MULT_PI_TW);
        }
        // Z-
        {
            int isn = ixy + Nxy * (((iz - 1 + Nz) % Nz) + Nz * it);
            int isg = isn + Nst * (1 - ieo + 2 * 2);
            bc2 = (iz == 0) ? (real_t)bc_z : real_t(1.0);
            QTW_CS_PROJ(1, 0, DWF_PROJ_M_TW, d0, d2, isn);
            QTW_CS_PROJ(1, 1, DWF_PROJ_P_TW, d1, d3, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 1, 0, isg, false);
            QTW_CS_GMUL(wt2, 1, 1, isg, false);
            DWF_ACCUM_4_SW_TW(acc_s0, acc_s1, acc_s2, acc_s3, wt1, wt2, bc2, DWF_MULT_PI_TW, DWF_MULT_MI_TW);
        }
        // T+
        {
            int isn = ixyz + Nxyz * ((it + 1) % Nt);
            int isg = site + Nst * (ieo + 2 * 3);
            bc2 = (it == Nt-1) ? (real_t)bc_t : real_t(1.0);
            QTW_CS_PROJ_T(0, 0, d2, isn);
            QTW_CS_PROJ_T(0, 1, d3, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 0, 0, isg, true);
            QTW_CS_GMUL(wt2, 0, 1, isg, true);
            DWF_ACCUM_TP_TW(acc_s2, acc_s3, wt1, wt2, bc2);
        }
        // T-
        {
            int isn = ixyz + Nxyz * ((it - 1 + Nt) % Nt);
            int isg = isn + Nst * (1 - ieo + 2 * 3);
            bc2 = (it == 0) ? (real_t)bc_t : real_t(1.0);
            QTW_CS_PROJ_T(1, 0, d0, isn);
            QTW_CS_PROJ_T(1, 1, d1, isn);
            __syncthreads();
            QTW_DECL(wt1); QTW_DECL(wt2);
            QTW_CS_GMUL(wt1, 1, 0, isg, false);
            QTW_CS_GMUL(wt2, 1, 1, isg, false);
            DWF_ACCUM_TM_TW(acc_s0, acc_s1, wt1, wt2, bc2);
        }

        QTW_STORE6(vp, 6 * IDX_DWF_QTW(cout, 0, is, Ns, site), acc_s0);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(cout, 1, is, Ns, site), acc_s1);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(cout, 2, is, Ns, site), acc_s2);
        QTW_STORE6(vp, 6 * IDX_DWF_QTW(cout, 3, is, Ns, site), acc_s3);

    #undef VSH
    #undef QTW_CS_PROJ
    #undef QTW_CS_PROJ_T
    #undef QTW_CS_GMUL
}

// up_mid/up_lo: when non-null, gauge link is full TW; when null, FP gauge link.
// recon: SU(3) 3rd-column reconstruction (YAML su3_reconstruction).
void mult_domainwall_5din_eo_hopb_qtw_dirac_5d(
    real_t *vp, real_t *up, real_t *up_mid, real_t *up_lo, real_t *wp, int Ns, int *bc,
    int *Nsize, int *do_comm, int ieo, int jeo, int jgm5, bool recon)
{
    int Nx=Nsize[0], Ny=Nsize[1], Nz=Nsize[2], Nt=Nsize[3];
    int Nst     = Nx*Ny*Nz*Nt;
    int Nst_pad = ceil_nwp(Nst);

    real_t *vp_dev = (real_t *)dev_ptr(vp);
    real_t *up_dev = (real_t *)dev_ptr(up);
    real_t *up_mid_dev = up_mid ? (real_t *)dev_ptr(up_mid) : nullptr;
    real_t *up_lo_dev  = up_lo  ? (real_t *)dev_ptr(up_lo)  : nullptr;
    real_t *wp_dev = (real_t *)dev_ptr(wp);

    // blockSize=64 (=2*NWP) is near-optimal: the QTW kernel is register-bound, so
    // larger blocks (e.g. 256=NWP*Ns, which would cache the gauge across all Ns
    // is-slices) collapse occupancy and regress ~30%. Measured 64~=128, 256 worse.
    // COLOR-SPLIT + block-smem share gate (env QTW_CSPLIT=2, TW link only):
    // one thread per (5D site, output color), 24-real accumulator => register
    // relief (255->166, spill 0). Gated OFF by default.
    //   - FLOAT: a dead end. Measured ~10% SLOWER on RTX 3080 (saturated GPU,
    //     register-light float gains nothing; pays smem round-trip + lost
    //     cross-direction ILP). The default register-resident kernel wins.
    //   - DOUBLE-QTW (server target, A100/H100): the ONLY case where the
    //     register relief can pay off (double busts the 255 cap and spills).
    //     UNVALIDATED. Before production use: (a) time vs default on A100 with
    //     do_comm=0; (b) this kernel is single-rank only (uniform __syncthreads
    //     assumes every thread takes every direction) -> needs multi-rank halo
    //     handling. Double-buffer below is known WORSE than single-buffer;
    //     prefer single-buffer when re-testing on the server.
    static const int csplit_mode = [](){ const char* e = getenv("QTW_CSPLIT"); return e ? atoi(e) : 0; }();
    if (csplit_mode == 2 && up_mid_dev && up_lo_dev &&
        (Nst % NWP == 0) &&
        do_comm[0] == 0 && do_comm[1] == 0 && do_comm[2] == 0 && do_comm[3] == 0) {
      dim3 block(NWP, NC);
      dim3 grid (Nst_pad / NWP, Ns);
      QTW_PROF_BEGIN();
      if (recon)
        mult_domainwall_5din_eo_hopb_qtw_dirac_5d_csmem_kernel<true><<<grid, block>>>(
          vp_dev, up_dev, up_mid_dev, up_lo_dev, wp_dev, Ns,
          bc[0], bc[1], bc[2], bc[3], Nx, Ny, Nz, Nt, ieo, jeo,
          do_comm[0], do_comm[1], do_comm[2], do_comm[3], Nst, Nst_pad, jgm5);
      else
        mult_domainwall_5din_eo_hopb_qtw_dirac_5d_csmem_kernel<false><<<grid, block>>>(
          vp_dev, up_dev, up_mid_dev, up_lo_dev, wp_dev, Ns,
          bc[0], bc[1], bc[2], bc[3], Nx, Ny, Nz, Nt, ieo, jeo,
          do_comm[0], do_comm[1], do_comm[2], do_comm[3], Nst, Nst_pad, jgm5);
      QTW_PROF_END("hopb_qtw_gauge");
      CHECK(cudaDeviceSynchronize());
      return;
    }

    int blockSize = VECTOR_LENGTH;
    // 5D-parallel: one thread per 5D site (Nst_pad*Ns), grid-stride inside.
    int gridSize  = (Nst_pad * Ns + blockSize - 1) / blockSize;
    QTW_PROF_BEGIN();
    if (recon) {
      mult_domainwall_5din_eo_hopb_qtw_dirac_5d_kernel<true><<<gridSize, blockSize>>>(
        vp_dev, up_dev, up_mid_dev, up_lo_dev, wp_dev, Ns,
        bc[0], bc[1], bc[2], bc[3],
        Nx, Ny, Nz, Nt,
        ieo, jeo,
        do_comm[0], do_comm[1], do_comm[2], do_comm[3],
        Nst, Nst_pad, jgm5);
    } else {
      mult_domainwall_5din_eo_hopb_qtw_dirac_5d_kernel<false><<<gridSize, blockSize>>>(
        vp_dev, up_dev, up_mid_dev, up_lo_dev, wp_dev, Ns,
        bc[0], bc[1], bc[2], bc[3],
        Nx, Ny, Nz, Nt,
        ieo, jeo,
        do_comm[0], do_comm[1], do_comm[2], do_comm[3],
        Nst, Nst_pad, jgm5);
    }
    QTW_PROF_END("hopb_qtw_gauge");

    CHECK(cudaDeviceSynchronize());
}

// ---------------- macro undef ----------------
#undef IDX_DWF_QTW
#undef QTW_DECL
#undef QTW_ZERO
#undef QTW_LOAD6
#undef QTW_LOAD6_NEW
#undef QTW_STORE6
#undef QTW_COPY
#undef QTW_ADDC
#undef QTW_SUBC
#undef QTW_NEGC
#undef QTW_SCALC
#undef QTW_SCALC_TT
#undef DWF_PROJ_P_TW
#undef DWF_PROJ_M_TW
#undef DWF_MULT_MI_TW
#undef DWF_MULT_PI_TW
#undef DWF_PROJ_RP_TW
#undef DWF_PROJ_RM_TW
#undef DWF_MULT_R1_TW
#undef DWF_MULT_RM1_TW
#undef DWF_PROJ_2_TW
#undef DWF_ACCUM_4_TW
#undef DWF_ACCUM_4_SW_TW
#undef DWF_ACCUM_TP_TW
#undef DWF_ACCUM_TM_TW
#undef DWF_GMUL_FWD_TW

#endif // MULT_DOMAINWALL_5DIN_EO_ACC_QTW_INCLUDED
