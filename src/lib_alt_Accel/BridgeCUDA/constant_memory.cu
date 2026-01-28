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
	
};
