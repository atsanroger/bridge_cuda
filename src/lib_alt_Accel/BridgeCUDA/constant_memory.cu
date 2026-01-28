/*!
      @file    constant_memory.cu
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2025-02-26 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "constant_memory.h"

namespace BridgeACC {

	__constant__ double const_e[NS_MAX]     = {0};
	__constant__ double const_f[NS_MAX]     = {0};
	__constant__ double const_dpinv[NS_MAX] = {0};
	__constant__ double const_dm[NS_MAX]    = {0};
	__constant__ double const_b[NS_MAX]     = {0};
	__constant__ double const_c[NS_MAX]     = {0};

	void initDomainwallConstantMemory(
		double* e, double* f,
		double* dpinv, double* dm,
		double* b, double* c,
		int Ns){

		CHECK(cudaMemcpyToSymbol(const_e, e,         Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_f, f,         Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dpinv, dpinv, Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_dm, dm,       Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_b, b,         Ns * sizeof(double)));
		CHECK(cudaMemcpyToSymbol(const_c, c,         Ns * sizeof(double)));

		};

	void initDomainwallConstantMemory(
			float* e, float* f,
			float* dpinv, float* dm,
			float* b, float* c,
			int Ns){

			double* e_f     = new double[Ns];
			double* f_f     = new double[Ns];
			double* dpinv_f = new double[Ns];
			double* dm_f    = new double[Ns];
			double* b_f     = new double[Ns];
			double* c_f     = new double[Ns];
	
			for (int i = 0; i < Ns; i++) {
				e_f[i]     = static_cast<double>(e[i]);
				f_f[i]     = static_cast<double>(f[i]);
				dpinv_f[i] = static_cast<double>(dpinv[i]);
				dm_f[i]    = static_cast<double>(dm[i]);
				b_f[i]     = static_cast<double>(b[i]);
				c_f[i]     = static_cast<double>(c[i]);
			}
	
			CHECK(cudaMemcpyToSymbol(const_e, e_f,         Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_f, f_f,         Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_dpinv, dpinv_f, Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_dm, dm_f,       Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_b, b_f,         Ns * sizeof(double)));
			CHECK(cudaMemcpyToSymbol(const_c, c_f,         Ns * sizeof(double)));

			delete [] e_f;
			delete [] f_f;
			delete [] dpinv_f;
			delete [] dm_f;
			delete [] b_f;
			delete [] c_f;		

		};
	
};
