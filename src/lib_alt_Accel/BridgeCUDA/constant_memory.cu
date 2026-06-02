/*!
      @file    constant_memory.cu
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2025-02-26 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "constant_memory.h"

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
			double* e, double* f,
			double* dpinv, double* dm,
			double* b, double* c,
			int Ns){

		float e_hi[NS_MAX], f_hi[NS_MAX], dpinv_hi[NS_MAX], dm_hi[NS_MAX];
		float b_hi[NS_MAX], c_hi[NS_MAX];
		float e_lo[NS_MAX], f_lo[NS_MAX], dpinv_lo[NS_MAX], dm_lo[NS_MAX];
		float b_lo[NS_MAX], c_lo[NS_MAX];

		auto split = [](double d, float &h, float &l) {
			h = (float)d;
			l = (float)(d - (double)h);
		};

		for (int i = 0; i < Ns; ++i) {
			split(e[i],     e_hi[i],     e_lo[i]);
			split(f[i],     f_hi[i],     f_lo[i]);
			split(dpinv[i], dpinv_hi[i], dpinv_lo[i]);
			split(dm[i],    dm_hi[i],    dm_lo[i]);
			split(b[i],     b_hi[i],     b_lo[i]);
			split(c[i],     c_hi[i],     c_lo[i]);
		}

		CHECK(cudaMemcpyToSymbol(const_e_float,     e_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_float,     f_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_float, dpinv_hi, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_float,    dm_hi,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_float,     b_hi,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_float,     c_hi,     Ns * sizeof(float)));

		CHECK(cudaMemcpyToSymbol(const_e_ff_lo,     e_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_f_ff_lo,     f_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dpinv_ff_lo, dpinv_lo, Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_dm_ff_lo,    dm_lo,    Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_b_ff_lo,     b_lo,     Ns * sizeof(float)));
		CHECK(cudaMemcpyToSymbol(const_c_ff_lo,     c_lo,     Ns * sizeof(float)));
	}

	// Push triple-word (h, m, l) decomposition of each double coefficient.
	// hi  = (float)d, mid = (float)(d - hi), lo = (float)(d - hi - mid).
	// Together they encode the coefficient at ~72-bit precision using only
	// FP32 storage. Used by the QTW (FP32-only triple-word) path.
	void initDomainwallConstantMemoryTW(
			double* e, double* f,
			double* dpinv, double* dm,
			double* b, double* c,
			int Ns){

		float e_hi[NS_MAX],  f_hi[NS_MAX],  dpinv_hi[NS_MAX], dm_hi[NS_MAX];
		float b_hi[NS_MAX],  c_hi[NS_MAX];
		float e_md[NS_MAX],  f_md[NS_MAX],  dpinv_md[NS_MAX], dm_md[NS_MAX];
		float b_md[NS_MAX],  c_md[NS_MAX];
		float e_lo[NS_MAX],  f_lo[NS_MAX],  dpinv_lo[NS_MAX], dm_lo[NS_MAX];
		float b_lo[NS_MAX],  c_lo[NS_MAX];

		auto split3 = [](double d, float &h, float &m, float &l) {
			h = (float)d;
			double r1 = d - (double)h;
			m = (float)r1;
			double r2 = r1 - (double)m;
			l = (float)r2;
		};

		for (int i = 0; i < Ns; ++i) {
			split3(e[i],     e_hi[i],     e_md[i],     e_lo[i]);
			split3(f[i],     f_hi[i],     f_md[i],     f_lo[i]);
			split3(dpinv[i], dpinv_hi[i], dpinv_md[i], dpinv_lo[i]);
			split3(dm[i],    dm_hi[i],    dm_md[i],    dm_lo[i]);
			split3(b[i],     b_hi[i],     b_md[i],     b_lo[i]);
			split3(c[i],     c_hi[i],     c_md[i],     c_lo[i]);
		}

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
	}

};
