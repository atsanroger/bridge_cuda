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

namespace BridgeACC {

#ifdef __CUDACC__
// QDW (Quasi Double-Word) Helper Macros & Functions
// Data Layout Option C: double4 {r_hi, i_hi, r_lo, i_lo}

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

// QDW Add: res = a + b
#define QDW_ADD(res, a, b) \
    QDW_R_HI(res) = QDW_R_HI(a) + QDW_R_HI(b); \
    QDW_I_HI(res) = QDW_I_HI(a) + QDW_I_HI(b); \
    QDW_R_LO(res) = QDW_R_LO(a) + QDW_R_LO(b); \
    QDW_I_LO(res) = QDW_I_LO(a) + QDW_I_LO(b);

// QDW Subtract: res = a - b
#define QDW_SUB(res, a, b) \
    QDW_R_HI(res) = QDW_R_HI(a) - QDW_R_HI(b); \
    QDW_I_HI(res) = QDW_I_HI(a) - QDW_I_HI(b); \
    QDW_R_LO(res) = QDW_R_LO(a) - QDW_R_LO(b); \
    QDW_I_LO(res) = QDW_I_LO(a) - QDW_I_LO(b);

// QDW Negate: res = -a
#define QDW_NEG(res, a) \
    QDW_R_HI(res) = -QDW_R_HI(a); \
    QDW_I_HI(res) = -QDW_I_HI(a); \
    QDW_R_LO(res) = -QDW_R_LO(a); \
    QDW_I_LO(res) = -QDW_I_LO(a);


// QDW Multiply: res = u * v (Complex * Complex QDW)
// u is double2 (r, i), v is double4 (QDW)
__device__ __forceinline__ double4 qdw_mult_uc(double2 u, double4 v) {
    double4 res;
    double p_r_hi, e_r_hi;
    double p_i_hi, e_i_hi;
    
    // Real part: (ur*vr - ui*vi)
    double p1, e1, p2, e2;
    TwoProd(u.x, v.x, p1, e1); // ur * vr_hi
    TwoProd(u.y, v.y, p2, e2); // ui * vi_hi
    TwoSum(p1, -p2, p_r_hi, e_r_hi); // p1 - p2
    
    // Low part approximation
    double p_r_lo = u.x * v.z - u.y * v.w;
    
    res.x = p_r_hi;
    res.z = p_r_lo + e1 - e2 + e_r_hi;

    // Imaginary part: (ur*vi + ui*vr)
    TwoProd(u.x, v.y, p1, e1); // ur * vi_hi
    TwoProd(u.y, v.x, p2, e2); // ui * vr_hi
    TwoSum(p1, p2, p_i_hi, e_i_hi); // p1 + p2
    
    double p_i_lo = u.x * v.w + u.y * v.z;
    
    res.y = p_i_hi;
    res.w = p_i_lo + e1 + e2 + e_i_hi;
    
    return res;
}

#define QDW_LOAD(var, ptr, idx) var = ptr[idx]
#define QDW_STORE(ptr, idx, var) ptr[idx] = var

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

void convert(double *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, double *w, int nin, int nvol, int nvol_pad);

void copy(double *v, double *w, int nin, int nvol);
void copy(double *v, int nv1, double *w, int nv2, int nin, int nvol);

void axpy(double* v, int nv1, double a,
          double* w, int nv2, int nin, int nvol);
void axpy(double* v, int nv1, double ar, double ai,
          double* w, int nv2, int nin, int nvol);

void aypx(double a, double* v, int nv1,
          double* w, int nv2, int nin, int nvol);
void aypx(double ar, double ai, double* v, int nv1,
          double* w, int nv2, int nin, int nvol);

void scal(double* v, int nv1, double a, int nin, int nvol);
void scal(double* v, int nv1, double ar, double ai, int nin, int nvol);

double norm2(double* v1, double* red, int nin, int nvol, int nex);

double dot(double* v1, double* v2, double* red, int nin, int nvol, int nex);

void dotc(double* ar, double* ai, double* v1, double* v2,
          double* red1, double* red2, int nin, int nvol, int nex);

void xI(double* v1, int nin, int nvol);
void conjg(double* v1, int nin, int nvol);

// real_t = float

void afield_init(float *data, const int size);
void afield_tidyup(float *data, const int size);
void afield_set(float *v,  float a, const int nin, const int nv2);
void set(float *v,  float a, int size);

void copy_to_device(float *v, int nv);
void copy_to_device(float *v, int nv1, int nv);
void copy_from_device(float *v, int nv);
void copy_from_device(float *v, int nv1, int nv);

void convert(float *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, float *w, int nin, int nvol, int nvol_pad);

void copy(float *v, float *w, int nin, int nvol);
void copy(float *v, int nv1, float *w, int nv2, int nin, int nvol);

void axpy(float* v, int nv1, float a,
              float* w, int nv2, int nin, int nvol);
void axpy(float* v, int nv1, float ar, float ai,
              float* w, int nv2, int nin, int nvol);

void aypx(float a, float* v, int nv1,
          float* w, int nv2, int nin, int nvol);
void aypx(float ar, float ai, float* v, int nv1,
          float* w, int nv2, int nin, int nvol);

void scal(float* v, int nv1, float a, int nin, int nvol);
void scal(float* v, int nv1, float ar, float ai, int nin, int nvol);

float norm2(float* v1, float* red, int nin, int nvol, int nex);

float dot(float* v1, float* v2, float* red, int nin, int nvol, int nex);

void dotc(float* ar, float* ai, float* v1, float* v2,
          float* red1, float* red2, int nin, int nvol, int nex);

void xI(float* v1, int nin, int nvol);
void conjg(float* v1, int nin, int nvol);

// copy with double/float conversion

void copy(double *v, int nv1, float  *w, int nv2, int nin, int nvol);
void copy(float  *v, int nv1, double *w, int nv2, int nin, int nvol);

//====================================================================

}

#endif
