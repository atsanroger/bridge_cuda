/*!
      @file    bridgeACC_AField.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#ifndef BRIDGEACC_AFIELD_INCLUDED
#define BRIDGEACC_AFIELD_INCLUDED

#include <type_traits>

namespace BridgeACC {

#ifdef __CUDACC__
// QDW (Quasi Double-Word) Helper Macros & Functions
// Data Layout Option C: double4 {r_hi, i_hi, r_lo, i_lo}

#define QDW_R_HI(val) (val.x)
#define QDW_I_HI(val) (val.y)
#define QDW_R_LO(val) (val.z)
#define QDW_I_LO(val) (val.w)

template<typename T>
__device__ __forceinline__ void TwoSum(T a, T b, T &s, T &e) {
    s = a + b;
    T v = s - a;
    e = (a - (s - v)) + (b - v);
}

template<typename T>
__device__ __forceinline__ void TwoProd(T a, T b, T &p, T &e) {
    p = a * b;
    e = fma(a, b, -p);
}

// DW Addition: d = a + b 
// Algorithm by Knuth (requires |a.hi| >= |b.hi|, but usually TwoSum is safer without branching)
// However, since a and b are both DW, the accurate formula is:
// s1, e1 = TwoSum(a.hi, b.hi)
// s2 = a.lo + b.lo
// s3 = e1 + s2
// d.hi, d.lo = TwoSum(s1, s3)
template<typename T>
__device__ __forceinline__ void dw_add(T ah, T al, T bh, T bl, T &dh, T &dl) {
    T s1, e1, d_hi, d_lo;
    TwoSum(ah, bh, s1, e1);
    T s2 = al + bl;
    T s3 = e1 + s2;
    TwoSum(s1, s3, d_hi, d_lo);
    dh = d_hi;
    dl = d_lo;
}

// QDW Add: res = a + b
#define QDW_ADD(res, a, b) \
    dw_add(QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(b), QDW_R_LO(b), QDW_R_HI(res), QDW_R_LO(res)); \
    dw_add(QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(b), QDW_I_LO(b), QDW_I_HI(res), QDW_I_LO(res));

// QDW Subtract: res = a - b
#define QDW_SUB(res, a, b) \
    dw_add(QDW_R_HI(a), QDW_R_LO(a), -QDW_R_HI(b), -QDW_R_LO(b), QDW_R_HI(res), QDW_R_LO(res)); \
    dw_add(QDW_I_HI(a), QDW_I_LO(a), -QDW_I_HI(b), -QDW_I_LO(b), QDW_I_HI(res), QDW_I_LO(res));

// QDW Negate: res = -a
#define QDW_NEG(res, a) \
    QDW_R_HI(res) = -QDW_R_HI(a); \
    QDW_I_HI(res) = -QDW_I_HI(a); \
    QDW_R_LO(res) = -QDW_R_LO(a); \
    QDW_I_LO(res) = -QDW_I_LO(a);


// QDW Scalar Multiply: res = u * a  where u is scalar T, a is QDW pair (vh,vl).
// All arithmetic in T (no precision promotion): float-float when T=float,
// double-double when T=double. Single rounding on the cross-term via fma.
template<typename T>
__device__ __forceinline__ void dw_scal(T u, T vh, T vl, T &dh, T &dl) {
    T s1, e1;
    TwoProd(u, vh, s1, e1);
    T s3 = fma(u, vl, e1);
    TwoSum(s1, s3, dh, dl);
}

#define QDW_SCAL(res, u, a) \
    dw_scal(u, QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(res), QDW_R_LO(res)); \
    dw_scal(u, QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(res), QDW_I_LO(res));

// DD scalar × QDW complex: (sh,sl) is a T-pair scalar (float-float when
// T=float). Uses dw_mul on each component — entirely T arithmetic, no
// promotion to a wider type. This is the float-only path that lets a
// float-base QDW carry coefficients at ~48-bit precision.
#define QDW_SCAL_DD(res, sh, sl, a) \
    dw_mul((sh), (sl), QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(res), QDW_R_LO(res)); \
    dw_mul((sh), (sl), QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(res), QDW_I_LO(res));

// QDW Multiply: res = u * v (Complex * Complex QDW)
// u is real2 (r, i), v is real4 (QDW); works for double (double-double) and
// float (float-float) base via template deduction. S = scalar (double/float).
template<typename T2, typename T4>
__device__ __forceinline__ T4 qdw_mult_uc(T2 u, T4 v) {
    typedef decltype(v.x) S;
    T4 res;
    S p_r_hi, e_r_hi;
    S p_i_hi, e_i_hi;

    // Real part: (ur*vr - ui*vi)
    S p1, e1, p2, e2;
    TwoProd(u.x, v.x, p1, e1); // ur * vr_hi
    TwoProd(u.y, v.y, p2, e2); // ui * vi_hi
    TwoSum(p1, -p2, p_r_hi, e_r_hi); // p1 - p2


    // Low part approximation — fma'd to keep one rounding instead of three.
    S p_r_lo = fma(u.x, v.z, -u.y * v.w);

    res.x = p_r_hi;
    res.z = p_r_lo + e1 - e2 + e_r_hi;

    // Imaginary part: (ur*vi + ui*vr)
    TwoProd(u.x, v.y, p1, e1); // ur * vi_hi
    TwoProd(u.y, v.x, p2, e2); // ui * vr_hi
    TwoSum(p1, p2, p_i_hi, e_i_hi); // p1 + p2

    S p_i_lo = fma(u.x, v.w, u.y * v.z);

    res.y = p_i_hi;
    res.w = p_i_lo + e1 + e2 + e_i_hi;

    return res;
}

// DW scalar multiply: (ch,cl) = (ah,al) * (bh,bl)   (double-double / float-float)
template<typename T>
__device__ __forceinline__ void dw_mul(T ah, T al, T bh, T bl, T &ch, T &cl) {
    T p, e;
    TwoProd(ah, bh, p, e);
    // Two fmas keep the cross-term additions single-rounded; the al*bl term
    // is below precision and dropped.
    e = fma(ah, bl, e);
    e = fma(al, bh, e);
    TwoSum(p, e, ch, cl);
}

// QDW Multiply: res = u * v (Complex DD * Complex DD). Both u and v are real4
// (QDW {r_hi, i_hi, r_lo, i_lo}). Used for the float-float / double-double gauge
// link multiply. S = scalar type (double/float).
template<typename T4>
__device__ __forceinline__ T4 qdw_mult_cc(T4 u, T4 v) {
    typedef decltype(v.x) S;
    T4 res;
    S ar_h, ar_l, ai_h, ai_l;   // partial products
    // Real part: ur*vr - ui*vi
    dw_mul(u.x, u.z, v.x, v.z, ar_h, ar_l);  // ur * vr
    dw_mul(u.y, u.w, v.y, v.w, ai_h, ai_l);  // ui * vi
    dw_add(ar_h, ar_l, -ai_h, -ai_l, res.x, res.z);
    // Imag part: ur*vi + ui*vr
    dw_mul(u.x, u.z, v.y, v.w, ar_h, ar_l);  // ur * vi
    dw_mul(u.y, u.w, v.x, v.z, ai_h, ai_l);  // ui * vr
    dw_add(ar_h, ar_l, ai_h, ai_l, res.y, res.w);
    return res;
}

#define QDW_LOAD(var, ptr, idx) var = ptr[idx]
#define QDW_STORE(ptr, idx, var) ptr[idx] = var


// ========================================================================
// QTW (Quasi Triple-Word) Helpers
// Hida-Bailey-Li QD-library style: every op renormalizes to (h, m, l) with
// |m| <= ulp(h)/2 and |l| <= ulp(m)/2. When T = float the entire chain stays
// on FP32 fma ALUs, simulating ~72-bit mantissa.
// ========================================================================

// Per-complex layout: 6 floats {r_hi, i_hi, r_mid, i_mid, r_lo, i_lo}.
// We store as two real4: lane0 = {r_hi, i_hi, r_mid, i_mid}, lane1 = {r_lo, i_lo, 0, 0}.
// Accessors abstract the layout — kernels work with explicit scalar tuples.
//
// Naming: hi = .x, mid = .y, lo = .z components of a real3-like view; for
// complex, we keep (re_hi, im_hi, re_mid, im_mid, re_lo, im_lo) in flat order.
//
// Single TW real:  (h, m, l)
// Single TW complex: (rh, ih, rm, im, rl, il)

// ---- Renormalize 3-word: ensure |m| <= ulp(h)/2, |l| <= ulp(m)/2 ----
// Algorithm: two successive TwoSum passes (Knuth-style; sufficient for the
// numerical-stability proofs in Hida-Bailey-Li for TQ-arithmetic).
template<typename T>
__device__ __forceinline__ void Renormalize3(T &h, T &m, T &l)
{
    T s, e1;
    TwoSum(m, l, s, e1);     // m+l -> s + e1
    TwoSum(h, s, h, m);      // h+s -> h + m   (e1 absorbed below)
    TwoSum(m, e1, m, l);     // distribute e1
}

// ---- ThreeSum: (a, b, c) -> renormalized (h, m, l) ----
template<typename T>
__device__ __forceinline__ void ThreeSum(T a, T b, T c, T &h, T &m, T &l)
{
    T t1, t2, t3, t4;
    TwoSum(a, b, t1, t2);    // a+b
    TwoSum(t1, c, h, t3);    // sum hi
    TwoSum(t2, t3, m, l);    // distribute residuals
}

// ---- ThreeProd: a*b -> exact (p_h, p_m, p_l) in T ----
// Uses two fma to split the product into hi/lo, then nothing remains (a*b is
// exactly two T words). Stored in 3-word form with l = 0; consumers feed it
// into tw_add chains that absorb the zero properly.
template<typename T>
__device__ __forceinline__ void ThreeProd(T a, T b, T &h, T &m, T &l)
{
    h = a * b;
    m = fma(a, b, -h);
    l = T(0);
}

// ---- tw_add (sloppy): TW + TW -> TW without final renormalization.
// Inner-loop primitive. Caller must renormalize the chain at the end
// (e.g. after a stencil accumulation) via Renormalize3. ~22 flops.
template<typename T>
__device__ __forceinline__ void tw_add_sloppy(T ah, T am, T al,
                                              T bh, T bm, T bl,
                                              T &dh, T &dm, T &dl)
{
    // Hi: TwoSum to preserve top precision (6 ops)
    T sh, eh;
    TwoSum(ah, bh, sh, eh);
    // Mid: TwoSum, capture mid error (6 ops)
    T sm, em;
    TwoSum(am, bm, sm, em);
    // Fold hi-error into mid with a TwoSum so the rounding error (~eps^2*value)
    // is captured, not dropped — this is what makes the result genuinely
    // triple-word (~eps^3) rather than dual-word (~eps^2). A plain `sm + eh`
    // here silently loses ~eps^2 and caps the whole TW operator at ~double.
    T mid, ef;
    TwoSum(sm, eh, mid, ef);   // mid = sm+eh exactly = mid + ef
    T lo  = (al + bl) + em + ef; // residuals collected at ~eps^3
    // One pass of Knuth split to keep mid/lo separated (6 ops)
    TwoSum(mid, lo, dm, dl);
    dh = sh;
    // Skip Renormalize3 — caller does it after the full chain.
}

// ---- tw_add (precise): TW + TW -> TW renormalized.
// ~45 flops. Use only when single-shot precision is needed.
template<typename T>
__device__ __forceinline__ void tw_add(T ah, T am, T al,
                                       T bh, T bm, T bl,
                                       T &dh, T &dm, T &dl)
{
    tw_add_sloppy(ah, am, al, bh, bm, bl, dh, dm, dl);
    Renormalize3(dh, dm, dl);
}

// ---- tw_mul: TW * TW -> TW (4 fma, then 2 renormalizations) ----
// We keep the leading three terms of the Cauchy product:
//   ah*bh + ah*bm + am*bh   (+ residuals folded in)
template<typename T>
__device__ __forceinline__ void tw_mul(T ah, T am, T al,
                                       T bh, T bm, T bl,
                                       T &dh, T &dm, T &dl)
{
    // Genuine triple-word product (~eps^3). The leading three Cauchy terms
    // (ah*bh, ah*bm, am*bh) are formed with error-free TwoProd so their
    // rounding (~eps^2*value) is captured in the low word rather than dropped.
    // Computing the mid via a plain fma chain (the old code) silently lost
    // ~eps^2 and capped the product at dual-word (~double) accuracy.
    T p_h, p_l;
    TwoProd(ah, bh, p_h, p_l);            // ah*bh = p_h + p_l (exact)
    T q1_h, q1_l, q2_h, q2_l;
    TwoProd(ah, bm, q1_h, q1_l);          // ah*bm (exact)
    TwoProd(am, bh, q2_h, q2_l);          // am*bh (exact)
    // mid = p_l + q1_h + q2_h, capturing the two adds exactly.
    T s1, e1;  TwoSum(q1_h, q2_h, s1, e1);
    T m_h, e2; TwoSum(p_l, s1, m_h, e2);
    // lo: all remaining eps^2 terms + captured carries (rounded at ~eps^3).
    T lo = q1_l + q2_l + e1 + e2
         + fma(ah, bl, fma(am, bm, fma(al, bh, T(0))));
    ThreeSum(p_h, m_h, lo, dh, dm, dl);
}

// ---- tw_scal: T scalar * TW -> TW (cheaper than tw_mul) ----
template<typename T>
__device__ __forceinline__ void tw_scal(T u,
                                        T vh, T vm, T vl,
                                        T &dh, T &dm, T &dl)
{
    // Genuine triple-word scale (~eps^3). u*vh and u*vm are formed with TwoProd
    // so the rounding of the mid word is captured in the low word instead of
    // being dropped by a plain fma (which capped this at dual-word accuracy).
    T a_h, a_l;  TwoProd(u, vh, a_h, a_l);   // u*vh (exact)
    T b_h, b_l;  TwoProd(u, vm, b_h, b_l);   // u*vm (exact)
    T c_h = u * vl;                          // eps^2 term
    T m1, e1;    TwoSum(a_l, b_h, m1, e1);   // mid
    T l1 = (b_l + c_h) + e1;                 // lo (eps^2 terms, ~eps^3 rounding)
    ThreeSum(a_h, m1, l1, dh, dm, dl);
}

// ---- ThreeSum (sloppy): value-preserving, looser normalization (2 TwoSum). ----
// Skips ThreeSum's final residual-distribution pass. The value h+m+l equals
// a+b+c exactly (both TwoSum are error-free); the three words may overlap
// slightly. Intended for inner-loop chains that are renormalized later (e.g.
// the gauge-mult accumulation, renormalized at QTW_STORE6). ~12 vs ~18 ops.
template<typename T>
__device__ __forceinline__ void ThreeSum_sloppy(T a, T b, T c, T &h, T &m, T &l)
{
    T e;
    TwoSum(a, b, h, e);     // h = a+b, e = exact error
    TwoSum(e, c, m, l);     // m = e+c, l = exact error
}

// ---- tw_mul (sloppy): now delegates to the genuine triple-word product. ----
// The previous fma-chain form lost ~eps^2 (dual-word only); since the gauge
// multiply feeds these into the operator, that capped the whole solve at
// ~double precision. Use the correct product so the gauge runs at true ~eps^3.
template<typename T>
__device__ __forceinline__ void tw_mul_sloppy(T ah, T am, T al,
                                              T bh, T bm, T bl,
                                              T &dh, T &dm, T &dl)
{
    tw_mul(ah, am, al, bh, bm, bl, dh, dm, dl);
}

// ---- tw_scal (sloppy): delegates to the genuine triple-word scale. ----
template<typename T>
__device__ __forceinline__ void tw_scal_sloppy(T u,
                                               T vh, T vm, T vl,
                                               T &dh, T &dm, T &dl)
{
    tw_scal(u, vh, vm, vl, dh, dm, dl);
}

// ---- Macros for TW complex stored as 6 floats per complex ----
// We use two real4 packs per complex: tw0 = {rh, ih, rm, im}, tw1 = {rl, il, _, _}.
// QTW_R_HI/MID/LO and QTW_I_HI/MID/LO retrieve the six components.
#define QTW_R_HI(t)  ((t).first.x)
#define QTW_I_HI(t)  ((t).first.y)
#define QTW_R_MID(t) ((t).first.z)
#define QTW_I_MID(t) ((t).first.w)
#define QTW_R_LO(t)  ((t).second.x)
#define QTW_I_LO(t)  ((t).second.y)

// Aggregate "TW complex" container: a pair of real4 holding the six components.
// We avoid a custom struct (no constructors / triviality across translation
// units) by using a simple C++ struct of two real4. CUDA inlines this fully.
template<typename T4>
struct qtw_cplx_t { T4 first; T4 second; };

// QTW complex add: res = a + b
#define QTW_ADD(res, a, b) do { \
    tw_add(QTW_R_HI(a), QTW_R_MID(a), QTW_R_LO(a), \
           QTW_R_HI(b), QTW_R_MID(b), QTW_R_LO(b), \
           QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res)); \
    tw_add(QTW_I_HI(a), QTW_I_MID(a), QTW_I_LO(a), \
           QTW_I_HI(b), QTW_I_MID(b), QTW_I_LO(b), \
           QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res)); \
} while (0)

// QTW complex subtract: res = a - b  (negate b in place using -h, -m, -l)
#define QTW_SUB(res, a, b) do { \
    tw_add(QTW_R_HI(a), QTW_R_MID(a), QTW_R_LO(a), \
           -QTW_R_HI(b), -QTW_R_MID(b), -QTW_R_LO(b), \
           QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res)); \
    tw_add(QTW_I_HI(a), QTW_I_MID(a), QTW_I_LO(a), \
           -QTW_I_HI(b), -QTW_I_MID(b), -QTW_I_LO(b), \
           QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res)); \
} while (0)

// QTW negate
#define QTW_NEG(res, a) do { \
    QTW_R_HI(res)  = -QTW_R_HI(a);  QTW_I_HI(res)  = -QTW_I_HI(a); \
    QTW_R_MID(res) = -QTW_R_MID(a); QTW_I_MID(res) = -QTW_I_MID(a); \
    QTW_R_LO(res)  = -QTW_R_LO(a);  QTW_I_LO(res)  = -QTW_I_LO(a); \
} while (0)

// QTW scalar multiply: res = u * a  (u single T)
#define QTW_SCAL(res, u, a) do { \
    tw_scal(u, QTW_R_HI(a), QTW_R_MID(a), QTW_R_LO(a), \
            QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res)); \
    tw_scal(u, QTW_I_HI(a), QTW_I_MID(a), QTW_I_LO(a), \
            QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res)); \
} while (0)

// QTW × TW-scalar multiply: res = (sh,sm,sl) * a  (TW scalar × TW complex)
#define QTW_SCAL_TT(res, sh, sm, sl, a) do { \
    tw_mul((sh), (sm), (sl), QTW_R_HI(a), QTW_R_MID(a), QTW_R_LO(a), \
            QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res)); \
    tw_mul((sh), (sm), (sl), QTW_I_HI(a), QTW_I_MID(a), QTW_I_LO(a), \
            QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res)); \
} while (0)

// QTW complex multiply: complex-T scalar (T2 with .x, .y) × QTW complex.
// Used for gauge link × TW vector in the FP32-only DD-gauge path.
template<typename T2, typename TC>
__device__ __forceinline__ TC qtw_mult_uc(T2 u, TC v)
{
    typedef decltype(v.first.x) S;
    TC res;
    // real:  ur*vr - ui*vi
    S r1_h, r1_m, r1_l;
    tw_scal<S>(u.x, QTW_R_HI(v), QTW_R_MID(v), QTW_R_LO(v), r1_h, r1_m, r1_l);
    S r2_h, r2_m, r2_l;
    tw_scal<S>(u.y, QTW_I_HI(v), QTW_I_MID(v), QTW_I_LO(v), r2_h, r2_m, r2_l);
    tw_add<S>(r1_h, r1_m, r1_l, -r2_h, -r2_m, -r2_l,
              QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res));
    // imag: ur*vi + ui*vr
    tw_scal<S>(u.x, QTW_I_HI(v), QTW_I_MID(v), QTW_I_LO(v), r1_h, r1_m, r1_l);
    tw_scal<S>(u.y, QTW_R_HI(v), QTW_R_MID(v), QTW_R_LO(v), r2_h, r2_m, r2_l);
    tw_add<S>(r1_h, r1_m, r1_l, r2_h, r2_m, r2_l,
              QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res));
    return res;
}

// QTW × QTW complex multiply (TW gauge × TW vector).
template<typename TC>
__device__ __forceinline__ TC qtw_mult_cc(TC u, TC v)
{
    typedef decltype(v.first.x) S;
    TC res;
    S ar_h, ar_m, ar_l, ai_h, ai_m, ai_l;
    // real: ur*vr - ui*vi
    tw_mul<S>(QTW_R_HI(u), QTW_R_MID(u), QTW_R_LO(u),
              QTW_R_HI(v), QTW_R_MID(v), QTW_R_LO(v), ar_h, ar_m, ar_l);
    tw_mul<S>(QTW_I_HI(u), QTW_I_MID(u), QTW_I_LO(u),
              QTW_I_HI(v), QTW_I_MID(v), QTW_I_LO(v), ai_h, ai_m, ai_l);
    tw_add<S>(ar_h, ar_m, ar_l, -ai_h, -ai_m, -ai_l,
              QTW_R_HI(res), QTW_R_MID(res), QTW_R_LO(res));
    // imag: ur*vi + ui*vr
    tw_mul<S>(QTW_R_HI(u), QTW_R_MID(u), QTW_R_LO(u),
              QTW_I_HI(v), QTW_I_MID(v), QTW_I_LO(v), ar_h, ar_m, ar_l);
    tw_mul<S>(QTW_I_HI(u), QTW_I_MID(u), QTW_I_LO(u),
              QTW_R_HI(v), QTW_R_MID(v), QTW_R_LO(v), ai_h, ai_m, ai_l);
    tw_add<S>(ar_h, ar_m, ar_l, ai_h, ai_m, ai_l,
              QTW_I_HI(res), QTW_I_MID(res), QTW_I_LO(res));
    return res;
}

#endif // __CUDACC__


// real_t = double

void afield_init(double *data, const int size);
void afield_tidyup(double *data, const int size);
void afield_set(double *v,  double a, const int nin, const int nv2);

void set(double *v,  double a, int size);

void copy_to_device(double *v, int nv);
void copy_to_device(double *v, int nv1, int nv);
void copy_from_device(double *v, int nv);
void copy_from_device(double *v, int nv1, int nv);

void copy_to_qdw(double* qdw_v, const double* std_v, int nvol, int nin4);
void copy_from_qdw(double* std_v, const double* qdw_v, int nvol, int nin4);

void convert(double *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, double *w, int nin, int nvol, int nvol_pad);

void copy(double *v, double *w, int nin, int nvol);
void copy(double *v, int nv1, double *w, int nv2, int nin, int nvol);

void axpy(double* v, int nv1, double a,
          double* w, int nv2, int nin, int nvol, int mode = 0);
void axpy(double* v, int nv1, double ar, double ai,
          double* w, int nv2, int nin, int nvol, int mode = 0);

void aypx(double a, double* v, int nv1,
          double* w, int nv2, int nin, int nvol, int mode = 0);
void aypx(double ar, double ai, double* v, int nv1,
          double* w, int nv2, int nin, int nvol, int mode = 0);

void scal(double* v, int nv1, double a, int nin, int nvol);
void scal(double* v, int nv1, double ar, double ai, int nin, int nvol);

double norm2(double* v1, double* red, int nin, int nvol, int nex, int mode = 0);

double dot(double* v1, double* v2, double* red, int nin, int nvol, int nex, int mode = 0);

void dotc(double* ar, double* ai, double* v1, double* v2,
          double* red1, double* red2, int nin, int nvol, int nex, int mode = 0);

// DD-pair BLAS for QDW/QTW: scalars carried as (h, l) float pairs.
void axpy_pair(double* v, int nv1, double a_h, double a_l,
               double* w, int nv2, int nin, int nvol, int mode = 0);
void aypx_pair(double a_h, double a_l, double* v, int nv1,
               double* w, int nv2, int nin, int nvol, int mode = 0);
void norm2_pair(double* h, double* l, double* v1, double* red,
                int nin, int nvol, int nex, int mode = 0);
void dot_pair(double* h, double* l, double* v1, double* v2, double* red,
              int nin, int nvol, int nex, int mode = 0);
void dotc_pair(double* ar_h, double* ar_l, double* ai_h, double* ai_l,
               double* v1, double* v2, double* red1, double* red2,
               int nin, int nvol, int nex, int mode = 0);

// Triple-word (QTW) BLAS. Field layout: 6 reals/cplx {rh, ih, rm, im, rl, il}.
void axpy_tw3(double* y, int nv1, double a,
              double* x, int nv2, int nin, int nvol);
void axpy_tw3_pair(double* y, int nv1, double a_h, double a_m, double a_l,
                   double* x, int nv2, int nin, int nvol);
void aypx_tw3(double a, double* y, int nv1,
              double* x, int nv2, int nin, int nvol);
void aypx_tw3_pair(double a_h, double a_m, double a_l, double* y, int nv1,
                   double* x, int nv2, int nin, int nvol);
void norm2_tw3(double* h, double* m, double* l,
               double* v1, double* red, int nin, int nvol, int nex);
void dot_tw3(double* h, double* m, double* l,
             double* v1, double* v2, double* red, int nin, int nvol, int nex);
void dotc_tw3(double* ar_h, double* ar_m, double* ar_l,
              double* ai_h, double* ai_m, double* ai_l,
              double* v1, double* v2, double* red1, double* red2,
              int nin, int nvol, int nex);
void normalize_tw3(double* v, int nin, int nvol);

void xI(double* v1, int nin, int nvol);
void conjg(double* v1, int nin, int nvol);

void normalize(double* v, int nin, int nvol, int mode = 0);

// real_t = float

void afield_init(float *data, const int size);
void afield_tidyup(float *data, const int size);
void afield_set(float *v,  float a, const int nin, const int nv2);
void set(float *v,  float a, int size);

void copy_to_device(float *v, int nv);
void copy_to_device(float *v, int nv1, int nv);
void copy_from_device(float *v, int nv);
void copy_from_device(float *v, int nv1, int nv);

void copy_to_qdw(float* qdw_v, const float* std_v, int nvol, int nin4);
void copy_from_qdw(float* std_v, const float* qdw_v, int nvol, int nin4);

void convert(float *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, float *w, int nin, int nvol, int nvol_pad);

void copy(float *v, float *w, int nin, int nvol);
void copy(float *v, int nv1, float *w, int nv2, int nin, int nvol);

void axpy(float* v, int nv1, float a,
              float* w, int nv2, int nin, int nvol, int mode = 0);
void axpy(float* v, int nv1, float ar, float ai,
              float* w, int nv2, int nin, int nvol, int mode = 0);

void aypx(float a, float* v, int nv1,
          float* w, int nv2, int nin, int nvol, int mode = 0);
void aypx(float ar, float ai, float* v, int nv1,
          float* w, int nv2, int nin, int nvol, int mode = 0);

void scal(float* v, int nv1, float a, int nin, int nvol);
void scal(float* v, int nv1, float ar, float ai, int nin, int nvol);

float norm2(float* v1, float* red, int nin, int nvol, int nex, int mode = 0);

float dot(float* v1, float* v2, float* red, int nin, int nvol, int nex, int mode = 0);

void dotc(float* ar, float* ai, float* v1, float* v2,
          float* red1, float* red2, int nin, int nvol, int nex, int mode = 0);

void axpy_pair(float* v, int nv1, float a_h, float a_l,
               float* w, int nv2, int nin, int nvol, int mode = 0);
void aypx_pair(float a_h, float a_l, float* v, int nv1,
               float* w, int nv2, int nin, int nvol, int mode = 0);
void norm2_pair(float* h, float* l, float* v1, float* red,
                int nin, int nvol, int nex, int mode = 0);
void dot_pair(float* h, float* l, float* v1, float* v2, float* red,
              int nin, int nvol, int nex, int mode = 0);
void dotc_pair(float* ar_h, float* ar_l, float* ai_h, float* ai_l,
               float* v1, float* v2, float* red1, float* red2,
               int nin, int nvol, int nex, int mode = 0);

void axpy_tw3(float* y, int nv1, float a,
              float* x, int nv2, int nin, int nvol);
void axpy_tw3_pair(float* y, int nv1, float a_h, float a_m, float a_l,
                   float* x, int nv2, int nin, int nvol);
void aypx_tw3(float a, float* y, int nv1,
              float* x, int nv2, int nin, int nvol);
void aypx_tw3_pair(float a_h, float a_m, float a_l, float* y, int nv1,
                   float* x, int nv2, int nin, int nvol);
void norm2_tw3(float* h, float* m, float* l,
               float* v1, float* red, int nin, int nvol, int nex);
void dot_tw3(float* h, float* m, float* l,
             float* v1, float* v2, float* red, int nin, int nvol, int nex);
void dotc_tw3(float* ar_h, float* ar_m, float* ar_l,
              float* ai_h, float* ai_m, float* ai_l,
              float* v1, float* v2, float* red1, float* red2,
              int nin, int nvol, int nex);
void normalize_tw3(float* v, int nin, int nvol);

void xI(float* v1, int nin, int nvol);
void conjg(float* v1, int nin, int nvol);

void normalize(float* v, int nin, int nvol, int mode = 0);

// copy with double/float conversion

void copy(double *v, int nv1, float  *w, int nv2, int nin, int nvol);
void copy(float  *v, int nv1, double *w, int nv2, int nin, int nvol);

//====================================================================

}

#endif
