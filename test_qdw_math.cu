#include <iostream>
#include <cmath>
#include <cuda_runtime.h>

#define QDW_R_HI(val) (val.x)
#define QDW_I_HI(val) (val.y)
#define QDW_R_LO(val) (val.z)
#define QDW_I_LO(val) (val.w)

__device__ __forceinline__ void TwoSum(double a, double b, double &s, double &e) {
    s = a + b;
    double v = s - a;
    e = (a - (s - v)) + (b - v);
}

__device__ __forceinline__ void TwoProd(double a, double b, double &p, double &e) {
    p = a * b;
    e = __fma_rn(a, b, -p);
}

__device__ __forceinline__ void dw_add(double ah, double al, double bh, double bl, double &dh, double &dl) {
    double s1, e1, d_hi, d_lo;
    TwoSum(ah, bh, s1, e1);
    double s2 = al + bl;
    double s3 = e1 + s2;
    TwoSum(s1, s3, d_hi, d_lo);
    dh = d_hi;
    dl = d_lo;
}

__global__ void test_dw_add(double4 *res, double4 *a, double4 *b) {
    dw_add(QDW_R_HI(a[0]), QDW_R_LO(a[0]), QDW_R_HI(b[0]), QDW_R_LO(b[0]), QDW_R_HI(res[0]), QDW_R_LO(res[0]));
    dw_add(QDW_I_HI(a[0]), QDW_I_LO(a[0]), QDW_I_HI(b[0]), QDW_I_LO(b[0]), QDW_I_HI(res[0]), QDW_I_LO(res[0]));
}

int main() {
    double4 h_a = {1.0, 2.0, 1e-16, 2e-16};
    double4 h_b = {3.0, 4.0, 3e-16, 4e-16};
    double4 h_res;

    double4 *d_a, *d_b, *d_res;
    cudaMalloc(&d_a, sizeof(double4));
    cudaMalloc(&d_b, sizeof(double4));
    cudaMalloc(&d_res, sizeof(double4));

    cudaMemcpy(d_a, &h_a, sizeof(double4), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, &h_b, sizeof(double4), cudaMemcpyHostToDevice);

    test_dw_add<<<1, 1>>>(d_res, d_a, d_b);

    cudaMemcpy(&h_res, d_res, sizeof(double4), cudaMemcpyDeviceToHost);

    std::cout << "Result HI Re: " << h_res.x << " Im: " << h_res.y << std::endl;
    std::cout << "Result LO Re: " << h_res.z << " Im: " << h_res.w << std::endl;

    cudaFree(d_a); cudaFree(d_b); cudaFree(d_res);
    return 0;
}
