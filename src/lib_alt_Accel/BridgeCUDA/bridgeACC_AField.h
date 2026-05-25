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

// DW Addition: d = a + b 
// Algorithm by Knuth (requires |a.hi| >= |b.hi|, but usually TwoSum is safer without branching)
// However, since a and b are both DW, the accurate formula is:
// s1, e1 = TwoSum(a.hi, b.hi)
// s2 = a.lo + b.lo
// s3 = e1 + s2
// d.hi, d.lo = TwoSum(s1, s3)
__device__ __forceinline__ void dw_add(double ah, double al, double bh, double bl, double &dh, double &dl) {
    double s1, e1, d_hi, d_lo;
    TwoSum(ah, bh, s1, e1);
    double s2 = al + bl;
    double s3 = e1 + s2;
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


// QDW Scalar Multiply: res = u * a  where u is scalar double, a is QDW
// Note: u is double, a is QDW.
// s1, e1 = TwoProd(u, a.hi)
// s2 = u * a.lo
// s3 = e1 + s2
// d.hi, d.lo = TwoSum(s1, s3)
__device__ __forceinline__ void dw_scal(double u, double vh, double vl, double &dh, double &dl) {
    double s1, e1, d_hi, d_lo;
    TwoProd(u, vh, s1, e1);
    double s2 = u * vl;
    double s3 = e1 + s2;
    TwoSum(s1, s3, d_hi, d_lo);
    dh = d_hi;
    dl = d_lo;
}

#define QDW_SCAL(res, u, a) \
    dw_scal(u, QDW_R_HI(a), QDW_R_LO(a), QDW_R_HI(res), QDW_R_LO(res)); \
    dw_scal(u, QDW_I_HI(a), QDW_I_LO(a), QDW_I_HI(res), QDW_I_LO(res));

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

void xI(float* v1, int nin, int nvol);
void conjg(float* v1, int nin, int nvol);

void normalize(float* v, int nin, int nvol, int mode = 0);

// copy with double/float conversion

void copy(double *v, int nv1, float  *w, int nv2, int nin, int nvol);
void copy(float  *v, int nv1, double *w, int nv2, int nin, int nvol);

//====================================================================

}

#endif
