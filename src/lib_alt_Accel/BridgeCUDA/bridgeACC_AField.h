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
