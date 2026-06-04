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
	__constant__ double const_a_double[NS_MAX] = {0};  // double-build stub (unused)
	__constant__ float  const_a_float[NS_MAX]  = {0};  // hi word
	__constant__ float  const_a_ff_lo[NS_MAX]  = {0};  // float-float lo
	__constant__ float  const_a_tw_mid[NS_MAX] = {0};  // triple-word mid
	__constant__ float  const_a_tw_lo[NS_MAX]  = {0};  // triple-word lo

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
	}

};
