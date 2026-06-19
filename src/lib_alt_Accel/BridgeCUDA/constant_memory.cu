/*!
      @file    constant_memory.cu
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2025-02-26 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "constant_memory.h"
#include <cmath>   // std::fma for host double-double TW constant generation

namespace BridgeACC {

	__constant__ double const_e_double[NS_MAX]     = {0};
	__constant__ double const_f_double[NS_MAX]     = {0};
	__constant__ double const_dpinv_double[NS_MAX] = {0};
	__constant__ double const_dm_double[NS_MAX]    = {0};
	__constant__ double const_b_double[NS_MAX]     = {0};
	__constant__ double const_c_double[NS_MAX]     = {0};

	__constant__ float const_e_float[NS_MAX]     = {0};
	__constant__ float const_f_float[NS_MAX]     = {0};
	__constant__ float const_dpinv_float[NS_MAX] = {0};
	__constant__ float const_dm_float[NS_MAX]    = {0};
	__constant__ float const_b_float[NS_MAX]     = {0};
	__constant__ float const_c_float[NS_MAX]     = {0};

	// Float-float low words. Paired with const_*_float (which hold the
	// high words = (float)(double_coeff)). Allow a float run to carry
	// coefficients at ~48-bit (double-double-ish) precision using only
	// FP32 arithmetic via dw_add/dw_mul. Populated when extended_precision
	// is enabled in the fopr.
	__constant__ float const_e_ff_lo[NS_MAX]     = {0};
	__constant__ float const_f_ff_lo[NS_MAX]     = {0};
	__constant__ float const_dpinv_ff_lo[NS_MAX] = {0};
	__constant__ float const_dm_ff_lo[NS_MAX]    = {0};
	__constant__ float const_b_ff_lo[NS_MAX]     = {0};
	__constant__ float const_c_ff_lo[NS_MAX]     = {0};

	// Triple-word (QTW) coefficient arrays. Hi reuses const_*_float; mid/lo
	// are stored here. (hi + mid + lo) reproduces the original double to
	// ~72-bit precision using FP32 fma chains only.
	__constant__ float const_e_tw_mid[NS_MAX]     = {0};
	__constant__ float const_f_tw_mid[NS_MAX]     = {0};
	__constant__ float const_dpinv_tw_mid[NS_MAX] = {0};
	__constant__ float const_dm_tw_mid[NS_MAX]    = {0};
	__constant__ float const_b_tw_mid[NS_MAX]     = {0};
	__constant__ float const_c_tw_mid[NS_MAX]     = {0};

	__constant__ float const_e_tw_lo[NS_MAX]     = {0};
	__constant__ float const_f_tw_lo[NS_MAX]     = {0};
	__constant__ float const_dpinv_tw_lo[NS_MAX] = {0};
	__constant__ float const_dm_tw_lo[NS_MAX]    = {0};
	__constant__ float const_b_tw_lo[NS_MAX]     = {0};
	__constant__ float const_c_tw_lo[NS_MAX]     = {0};

	// Precomputed coupling a[j] = 0.5 * dm[j+1] * dpinv[j] (j=0..Ns-2), the
	// forward/backward-sweep coefficient that the LU(dag)inv kernels otherwise
	// recompute per-thread via a multiword multiply. Computed in host DD so the
	// stored words match the on-device product to full precision.
	__constant__ double const_a_double[NS_MAX] = {0};  // double-build hi word
	__constant__ float  const_a_float[NS_MAX]  = {0};  // hi word
	__constant__ float  const_a_ff_lo[NS_MAX]  = {0};  // float-float lo
	__constant__ float  const_a_tw_mid[NS_MAX] = {0};  // triple-word mid
	__constant__ float  const_a_tw_lo[NS_MAX]  = {0};  // triple-word lo

	// ---- double-build triple-word (double-triple) mid/lo words ----
	// For the double-base QDW/QTW path. hi reuses const_*_double. Phase 1:
	// sourced from the host double-double pipeline (mid = DD low word, lo = 0)
	// => genuine ~106-bit (double-double) reach, already ~1e16x past native
	// FP64. The lo slot is reserved for a future triple-double host upgrade
	// (full ~159-bit). Named "td" (triple-double) to stay distinct from the
	// float-build "tw" (triple-word) arrays above.
	__constant__ double const_e_td_mid[NS_MAX]     = {0};
	__constant__ double const_f_td_mid[NS_MAX]     = {0};
	__constant__ double const_dpinv_td_mid[NS_MAX] = {0};
	__constant__ double const_dm_td_mid[NS_MAX]    = {0};
	__constant__ double const_b_td_mid[NS_MAX]     = {0};
	__constant__ double const_c_td_mid[NS_MAX]     = {0};
	__constant__ double const_e_td_lo[NS_MAX]     = {0};
	__constant__ double const_f_td_lo[NS_MAX]     = {0};
	__constant__ double const_dpinv_td_lo[NS_MAX] = {0};
	__constant__ double const_dm_td_lo[NS_MAX]    = {0};
	__constant__ double const_b_td_lo[NS_MAX]     = {0};
	__constant__ double const_c_td_lo[NS_MAX]     = {0};
	__constant__ double const_a_td_mid[NS_MAX] = {0};  // double-build a mid (DD.lo)
	__constant__ double const_a_td_lo[NS_MAX]  = {0};  // double-build a lo (0; reserved)

	// ---- host double-double (DD) helpers for exact TW constant generation ----
	// The TW constants must carry >double precision, else they cap the QTW solve
	// at ~2^-53 (double) instead of the ~2^-72 the device triple-word arithmetic
	// can reach. So we recompute the preconditioner coefficients in double-double
	// and split the DD result into (hi,mid,lo) FP32 words whose lo now carries
	// genuine bits 54..72. Device stays FP32-only; DD is host constant prep only.
	namespace {
		struct DD { double hi, lo; };
		static inline DD quick_two_sum(double a, double b) {
			double s = a + b; double e = b - (s - a); return {s, e};
		}
		static inline DD two_sum(double a, double b) {
			double s = a + b; double v = s - a; double e = (a - (s - v)) + (b - v); return {s, e};
		}
		static inline DD two_prod(double a, double b) {
			double p = a * b; double e = std::fma(a, b, -p); return {p, e};
		}
		static inline DD dd_from(double a) { return DD{a, 0.0}; }
		static inline DD dd_add(DD a, DD b) {
			DD s = two_sum(a.hi, b.hi);
			DD t = two_sum(a.lo, b.lo);
			s.lo += t.hi; s = quick_two_sum(s.hi, s.lo);
			s.lo += t.lo; s = quick_two_sum(s.hi, s.lo);
			return s;
		}
		static inline DD dd_sub(DD a, DD b) { return dd_add(a, DD{-b.hi, -b.lo}); }
		static inline DD dd_mul(DD a, DD b) {
			DD p = two_prod(a.hi, b.hi);
			p.lo += a.hi * b.lo + a.lo * b.hi;
			p = quick_two_sum(p.hi, p.lo);
			return p;
		}
		static inline DD dd_div(DD a, DD b) {
			double q1 = a.hi / b.hi;
			DD r = dd_sub(a, dd_mul(b, dd_from(q1)));
			double q2 = r.hi / b.hi;
			r = dd_sub(r, dd_mul(b, dd_from(q2)));
			double q3 = r.hi / b.hi;
			DD q = quick_two_sum(q1, q2);
			return dd_add(q, dd_from(q3));
		}
		// Split a DD value into 3 FP32 words; lo carries genuine bits 54..72.
		static inline void split3_dd(DD d, float &h, float &m, float &l) {
			h = (float)d.hi;
			double r1 = d.hi - (double)h;          // exact (Sterbenz)
			m = (float)r1;
			double r2 = (r1 - (double)m) + d.lo;   // fold in the low word
			l = (float)r2;
		}
		// Split a DD value into 2 FP32 words (float-float); lo carries bits 25..48.
		static inline void split2_dd(DD d, float &h, float &l) {
			h = (float)d.hi;
			double r1 = (d.hi - (double)h) + d.lo; // exact part + low word
			l = (float)r1;
		}
		// Recompute the EO preconditioner coefficients in double-double from raw
		// inputs. Fills dm[Ns], dpinv[Ns], e[Ns-1], f[Ns-1] (the uploaded ones);
		// dp,g are intermediates. MUST mirror
		//   AFopr_Domainwall_5din_eo<>::set_precond_parameters()  (the *_d block).
		static inline void compute_precond_coeffs_dd(
				double M0, double mq, double alpha,
				const double* b, const double* c, int Ns,
				DD* dm, DD* dpinv, DD* e, DD* f) {
			DD dp[NS_MAX];
			const DD aDD = dd_from(alpha);
			const DD one = dd_from(1.0);
			const DD fourMinusM0 = dd_sub(dd_from(4.0), dd_from(M0));
			for (int is = 0; is < Ns; ++is) {
				dp[is] = dd_mul(aDD, dd_add(one, dd_mul(dd_from(b[is]), fourMinusM0)));
				dm[is] = dd_mul(aDD, dd_sub(one, dd_mul(dd_from(c[is]), fourMinusM0)));
			}
			const DD mqDD = dd_from(mq);
			e[0] = dd_div(dd_mul(mqDD, dm[Ns - 1]), dp[0]);
			f[0] = dd_div(dd_mul(mqDD, dm[0]),       aDD);
			for (int is = 1; is < Ns - 1; ++is) {
				e[is] = dd_div(dd_mul(e[is - 1], dm[is - 1]), dp[is]);
				f[is] = dd_div(dd_mul(f[is - 1], dm[is]),     dp[is - 1]);
			}
			const DD g = dd_mul(e[Ns - 2], dm[Ns - 2]);
			for (int is = 0; is < Ns - 1; ++is) dpinv[is] = dd_div(one, dp[is]);
			dpinv[Ns - 1] = dd_div(one, dd_add(dp[Ns - 1], g));
		}

		// ---- host triple-double (TD) for genuine double-TRIPLE constants ----
		// Phase 2: feed the double-base QTW kernels coefficients at ~159-bit so the
		// device triple-double arithmetic is NOT capped at double-double. The TD
		// top 106 bits are computed exactly via the dd_* primitives; the third limb
		// collects the leading sub-2^-106 cross/round terms. Inputs are O(1) and the
		// preconditioner formulas have no catastrophic cancellation, so this yields
		// ~150+ genuine bits — far past DD's 106. Correctness is checked empirically
		// by the true-residual floor dropping well below the double-double ~6e-57.
		struct TD { double m0, m1, m2; };
		static inline TD td_renorm(double a0, double a1, double a2) {
			DD r = quick_two_sum(a0, a1); double m0 = r.hi, c = r.lo;
			r = two_sum(c, a2);           double m1 = r.hi, m2 = r.lo;
			r = quick_two_sum(m0, m1);    m0 = r.hi; m1 = r.lo;
			r = quick_two_sum(m1, m2);    m1 = r.hi; m2 = r.lo;
			return { m0, m1, m2 };
		}
		static inline TD td_from(double a) { return { a, 0.0, 0.0 }; }
		static inline TD td_add(TD a, TD b) {
			DD sh = dd_add(DD{a.m0, a.m1}, DD{b.m0, b.m1});  // top ~106 bits
			double lo = a.m2 + b.m2;                          // third-limb sum
			return td_renorm(sh.hi, sh.lo, lo);
		}
		static inline TD td_neg(TD a) { return { -a.m0, -a.m1, -a.m2 }; }
		static inline TD td_sub(TD a, TD b) { return td_add(a, td_neg(b)); }
		static inline TD td_mul(TD a, TD b) {
			DD ph = dd_mul(DD{a.m0, a.m1}, DD{b.m0, b.m1});  // top ~106 bits
			// leading third-limb terms below the DD product (each ~2^-106):
			double cross = a.m0 * b.m2 + a.m2 * b.m0 + a.m1 * b.m1;
			return td_renorm(ph.hi, ph.lo, cross);
		}
		static inline TD td_div(TD a, TD b) {
			double q0 = a.m0 / b.m0;
			TD r = td_sub(a, td_mul(b, td_from(q0)));
			double q1 = r.m0 / b.m0;
			r = td_sub(r, td_mul(b, td_from(q1)));
			double q2 = r.m0 / b.m0;
			return td_renorm(q0, q1, q2);
		}
		static inline TD td_scale_half(TD a) { return { 0.5 * a.m0, 0.5 * a.m1, 0.5 * a.m2 }; }

		// TD mirror of compute_precond_coeffs_dd (same formula, ~159-bit).
		static inline void compute_precond_coeffs_td(
				double M0, double mq, double alpha,
				const double* b, const double* c, int Ns,
				TD* dm, TD* dpinv, TD* e, TD* f) {
			TD dp[NS_MAX];
			const TD aTD = td_from(alpha);
			const TD one = td_from(1.0);
			const TD fourMinusM0 = td_sub(td_from(4.0), td_from(M0));
			for (int is = 0; is < Ns; ++is) {
				dp[is] = td_mul(aTD, td_add(one, td_mul(td_from(b[is]), fourMinusM0)));
				dm[is] = td_mul(aTD, td_sub(one, td_mul(td_from(c[is]), fourMinusM0)));
			}
			const TD mqTD = td_from(mq);
			e[0] = td_div(td_mul(mqTD, dm[Ns - 1]), dp[0]);
			f[0] = td_div(td_mul(mqTD, dm[0]),       aTD);
			for (int is = 1; is < Ns - 1; ++is) {
				e[is] = td_div(td_mul(e[is - 1], dm[is - 1]), dp[is]);
				f[is] = td_div(td_mul(f[is - 1], dm[is]),     dp[is - 1]);
			}
			const TD g = td_mul(e[Ns - 2], dm[Ns - 2]);
			for (int is = 0; is < Ns - 1; ++is) dpinv[is] = td_div(one, dp[is]);
			dpinv[Ns - 1] = td_div(one, td_add(dp[Ns - 1], g));
		}
	} // anonymous namespace


	void initDomainwallConstantMemory(
		double* e, double* f,
		double* dpinv, double* dm,
		double* b, double* c,
		int Ns){

		CHECK(cudaMemcpyToSymbol(const_e_double,     e,     Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_f_double,     f,     Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_double, dpinv, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dm_double,    dm,    Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_b_double,     b,     Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_c_double,     c,     Ns * sizeof(double)));

		};

	void initDomainwallConstantMemory(
			float* e, float* f,
			float* dpinv, float* dm,
			float* b, float* c,
			int Ns){

			CHECK(cudaMemcpyToSymbol(const_e_float,     e,     Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_f_float,     f,     Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_dpinv_float, dpinv, Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_dm_float,    dm,    Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_b_float,     b,     Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_c_float,     c,     Ns * sizeof(float)));

		};

	// Upload ONLY the Moebius b/c coefficients (no e/f/dpinv/dm). Used by the
	// non-eo dd domainwall operator, whose 5dir/5dirdag kernels read const_b/
	// const_c directly. b/c are mass-independent, so multiple coexisting dd
	// operators (e.g. fineF and Pauli-Villars in a multigrid setup) all push the
	// same values -> const_b/const_c stay consistent. This deliberately does NOT
	// touch const_e/f/dpinv (which are mass-dependent and carry the QTW/FF
	// hi/mid/lo words for the eo preconditioner) so it cannot clobber an eo
	// operator's coefficients living in the same process.
	void setDomainwallConstantBC(
			double* b, double* c, int Ns){
			CHECK(cudaMemcpyToSymbol(const_b_double, b, Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_c_double, c, Ns * sizeof(double)));
		};

	void setDomainwallConstantBC(
			float* b, float* c, int Ns){
			CHECK(cudaMemcpyToSymbol(const_b_float, b, Ns * sizeof(float)));
			CHECK(cudaMemcpyToSymbol(const_c_float, c, Ns * sizeof(float)));
		};

	// Push float-float pairs derived from the double coefficients. const_*_float
	// gets the high words (= (float)d) and const_*_ff_lo gets the low words
	// ((float)(d - (float)d)). Together they encode each coefficient with
	// ~48-bit precision using only FP32 storage and arithmetic.
	void initDomainwallConstantMemoryFF(
			double M0, double mq, double alpha,
			const double* b, const double* c,
			int Ns){

		float e_hi[NS_MAX], f_hi[NS_MAX], dpinv_hi[NS_MAX], dm_hi[NS_MAX];
		float b_hi[NS_MAX], c_hi[NS_MAX];
		float e_lo[NS_MAX], f_lo[NS_MAX], dpinv_lo[NS_MAX], dm_lo[NS_MAX];
		float b_lo[NS_MAX], c_lo[NS_MAX];

		// DD coefficient pipeline (shared with the TW path), then 2-word split.
		// float-float caps at ~48-bit, so DD bits beyond 48 are dropped; this just
		// keeps FF consistent with TW and exact at the float-float level.
		DD dm[NS_MAX], dpinv[NS_MAX], e[NS_MAX], f[NS_MAX];
		compute_precond_coeffs_dd(M0, mq, alpha, b, c, Ns, dm, dpinv, e, f);

		auto split2 = [](double d, float &h, float &l) {
			h = (float)d;
			l = (float)(d - (double)h);
		};

		for (int i = 0; i < Ns; ++i) {
			split2_dd(dm[i],    dm_hi[i],    dm_lo[i]);
			split2_dd(dpinv[i], dpinv_hi[i], dpinv_lo[i]);
			split2(b[i], b_hi[i], b_lo[i]);   // exact input -> plain split
			split2(c[i], c_hi[i], c_lo[i]);
		}
		for (int i = 0; i < Ns - 1; ++i) {
			split2_dd(e[i], e_hi[i], e_lo[i]);
			split2_dd(f[i], f_hi[i], f_lo[i]);
		}
		e_hi[Ns - 1] = e_lo[Ns - 1] = 0.0f;
		f_hi[Ns - 1] = f_lo[Ns - 1] = 0.0f;

		// a[j] = 0.5 * dm[j+1] * dpinv[j], j=0..Ns-2 (float-float).
		float a_hi[NS_MAX], a_lo[NS_MAX];
		for (int j = 0; j < Ns - 1; ++j) {
			DD ap = dd_mul(dm[j + 1], dpinv[j]);
			DD aj = { 0.5 * ap.hi, 0.5 * ap.lo };
			split2_dd(aj, a_hi[j], a_lo[j]);
		}
		a_hi[Ns - 1] = a_lo[Ns - 1] = 0.0f;

		CHECK(cudaMemcpyToSymbol(const_e_float,     e_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_float,     f_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_float, dpinv_hi, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_float,    dm_hi,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_float,     b_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_float,     c_hi,     Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_a_float,  a_hi, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_a_ff_lo,  a_lo, Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_e_ff_lo,     e_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_ff_lo,     f_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_ff_lo, dpinv_lo, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_ff_lo,    dm_lo,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_ff_lo,     b_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_ff_lo,     c_lo,     Ns * sizeof(float)));
	}

	// Recompute the EO preconditioner coefficients in double-double, then push the
	// (hi,mid,lo) FP32 decomposition. Raw inputs (M0,mq,alpha,b,c) are taken so the
	// e/f recurrence and the reciprocals are evaluated at DD precision instead of
	// refining already-rounded doubles.  b,c are exact inputs -> plain split.
	// NOTE: this coefficient formula MUST mirror
	//   AFopr_Domainwall_5din_eo<>::set_precond_parameters()  (the *_d block).
	void initDomainwallConstantMemoryTW(
			double M0, double mq, double alpha,
			const double* b, const double* c,
			int Ns){

		float e_hi[NS_MAX],  f_hi[NS_MAX],  dpinv_hi[NS_MAX], dm_hi[NS_MAX];
		float b_hi[NS_MAX],  c_hi[NS_MAX];
		float e_md[NS_MAX],  f_md[NS_MAX],  dpinv_md[NS_MAX], dm_md[NS_MAX];
		float b_md[NS_MAX],  c_md[NS_MAX];
		float e_lo[NS_MAX],  f_lo[NS_MAX],  dpinv_lo[NS_MAX], dm_lo[NS_MAX];
		float b_lo[NS_MAX],  c_lo[NS_MAX];

		// --- DD coefficient pipeline (shared with the FF path) ---
		DD dm[NS_MAX], dpinv[NS_MAX], e[NS_MAX], f[NS_MAX];
		compute_precond_coeffs_dd(M0, mq, alpha, b, c, Ns, dm, dpinv, e, f);

		auto split3 = [](double d, float &h, float &m, float &l) {
			h = (float)d;
			double r1 = d - (double)h;
			m = (float)r1;
			double r2 = r1 - (double)m;
			l = (float)r2;
		};

		for (int i = 0; i < Ns; ++i) {
			split3_dd(dm[i],    dm_hi[i],    dm_md[i],    dm_lo[i]);
			split3_dd(dpinv[i], dpinv_hi[i], dpinv_md[i], dpinv_lo[i]);
			split3(b[i], b_hi[i], b_md[i], b_lo[i]);   // exact input -> plain split
			split3(c[i], c_hi[i], c_md[i], c_lo[i]);
		}
		for (int i = 0; i < Ns - 1; ++i) {
			split3_dd(e[i], e_hi[i], e_md[i], e_lo[i]);
			split3_dd(f[i], f_hi[i], f_md[i], f_lo[i]);
		}
		// e,f hold Ns-1 entries; zero the unused [Ns-1] slot (kernels never read it).
		e_hi[Ns - 1] = e_md[Ns - 1] = e_lo[Ns - 1] = 0.0f;
		f_hi[Ns - 1] = f_md[Ns - 1] = f_lo[Ns - 1] = 0.0f;

		// a[j] = 0.5 * dm[j+1] * dpinv[j], j=0..Ns-2 (triple-word). 0.5 is an exact
		// power-of-two scale (per-word), so no precision is lost vs the DD product.
		float a_hi[NS_MAX], a_md[NS_MAX], a_lo[NS_MAX];
		for (int j = 0; j < Ns - 1; ++j) {
			DD ap = dd_mul(dm[j + 1], dpinv[j]);
			DD aj = { 0.5 * ap.hi, 0.5 * ap.lo };
			split3_dd(aj, a_hi[j], a_md[j], a_lo[j]);
		}
		a_hi[Ns - 1] = a_md[Ns - 1] = a_lo[Ns - 1] = 0.0f;

		CHECK(cudaMemcpyToSymbol(const_e_float,     e_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_float,     f_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_float, dpinv_hi, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_float,    dm_hi,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_float,     b_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_float,     c_hi,     Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_e_tw_mid,     e_md,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_tw_mid,     f_md,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_tw_mid, dpinv_md, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_tw_mid,    dm_md,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_tw_mid,     b_md,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_tw_mid,     c_md,     Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_e_tw_lo,     e_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_tw_lo,     f_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_tw_lo, dpinv_lo, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_tw_lo,    dm_lo,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_tw_lo,     b_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_tw_lo,     c_lo,     Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_a_float,  a_hi, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_a_tw_mid, a_md, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_a_tw_lo,  a_lo, Ns * sizeof(float)));

		// ---- double-build (double-TRIPLE) constants, recomputed in host TD
		// (~159-bit). hi/mid/lo = TD.m0/m1/m2 -> genuine double-triple reach,
		// uncapped by double-double. We also OVERWRITE const_*_double with TD.m0
		// (the correctly-rounded hi) so hi+mid+lo reconstruct exactly; this is the
		// same value the FP native-double path wants, so it is not perturbed. The
		// float build is unaffected (its hi lives in the separate const_*_float).
		// const_a_double is populated here (was 0 for the double build -> the cause
		// of the earlier double+TW non-convergence). b,c are exact -> mid/lo = 0.
		TD dm_t[NS_MAX], dpinv_t[NS_MAX], e_t[NS_MAX], f_t[NS_MAX];
		compute_precond_coeffs_td(M0, mq, alpha, b, c, Ns, dm_t, dpinv_t, e_t, f_t);

		double dm_h[NS_MAX]={0}, dm_m[NS_MAX]={0}, dm_l[NS_MAX]={0};
		double dp_h[NS_MAX]={0}, dp_m[NS_MAX]={0}, dp_l[NS_MAX]={0};
		double e_h[NS_MAX]={0},  e_m[NS_MAX]={0},  e_l[NS_MAX]={0};
		double f_h[NS_MAX]={0},  f_m[NS_MAX]={0},  f_l[NS_MAX]={0};
		double b_h[NS_MAX]={0},  c_h[NS_MAX]={0},  zero_d[NS_MAX]={0};
		double a_h[NS_MAX]={0},  a_m[NS_MAX]={0},  a_l[NS_MAX]={0};
		for (int i = 0; i < Ns; ++i) {
			dm_h[i]=dm_t[i].m0; dm_m[i]=dm_t[i].m1; dm_l[i]=dm_t[i].m2;
			dp_h[i]=dpinv_t[i].m0; dp_m[i]=dpinv_t[i].m1; dp_l[i]=dpinv_t[i].m2;
			b_h[i]=b[i]; c_h[i]=c[i];   // exact inputs -> mid/lo = 0
		}
		for (int i = 0; i < Ns - 1; ++i) {
			e_h[i]=e_t[i].m0; e_m[i]=e_t[i].m1; e_l[i]=e_t[i].m2;
			f_h[i]=f_t[i].m0; f_m[i]=f_t[i].m1; f_l[i]=f_t[i].m2;
		}
		for (int j = 0; j < Ns - 1; ++j) {
			TD aj = td_scale_half(td_mul(dm_t[j + 1], dpinv_t[j]));
			a_h[j]=aj.m0; a_m[j]=aj.m1; a_l[j]=aj.m2;
		}
		// hi (overwrite const_*_double for consistency with mid/lo)
		CHECK(cudaMemcpyToSymbol(const_dm_double,    dm_h, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_double, dp_h, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_e_double,     e_h,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_f_double,     f_h,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_b_double,     b_h,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_c_double,     c_h,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_a_double,     a_h,  Ns * sizeof(double)));
		// mid
		CHECK(cudaMemcpyToSymbol(const_dm_td_mid,    dm_m, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_td_mid, dp_m, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_e_td_mid,     e_m,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_f_td_mid,     f_m,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_b_td_mid,     zero_d, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_c_td_mid,     zero_d, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_a_td_mid,     a_m,  Ns * sizeof(double)));
		// lo
		CHECK(cudaMemcpyToSymbol(const_dm_td_lo,     dm_l, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_td_lo,  dp_l, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_e_td_lo,      e_l,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_f_td_lo,      f_l,  Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_b_td_lo,      zero_d, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_c_td_lo,      zero_d, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_a_td_lo,      a_l,  Ns * sizeof(double)));
	}

};
