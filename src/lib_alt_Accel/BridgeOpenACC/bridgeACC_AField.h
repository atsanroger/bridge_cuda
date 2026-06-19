/*!
      @file    bridgeACC_AField.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef BRIDGEACC_AFIELD_INCLUDED
#define BRIDGEACC_AFIELD_INCLUDED

namespace BridgeACC {

// real_t = double

void afield_init(double *data, const int size);
void afield_tidyup(double *data, const int size);
void afield_set(double *v,  double a, const int nin, const int nv2);
void copy_to_device(double *v, int nv);
void copy_to_device(double *v, int nv1, int nv);
void copy_from_device(double *v, int nv);
void copy_from_device(double *v, int nv1, int nv);

void convert(double *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, double *w, int nin, int nvol, int nvol_pad);

void copy(double *v, double *w, int nin, int nvol);
void copy(double *v, int nv1, double *w, int nv2, int nin, int nvol);

void axpy(double *restrict v, int nv1, double a,
          double *restrict w, int nv2, int nin, int nvol);
void axpy(double *restrict v, int nv1, double ar, double ai,
          double *restrict w, int nv2, int nin, int nvol);

void aypx(double a, double *restrict v, int nv1,
          double *restrict w, int nv2, int nin, int nvol);
void aypx(double ar, double ai, double *restrict v, int nv1,
          double *restrict w, int nv2, int nin, int nvol);

void scal(double *restrict v, int nv1, double a, int nin, int nvol);
void scal(double *restrict v, int nv1, double ar, double ai,
                                                     int nin, int nvol);

double norm2(double *restrict v1, int nin, int nvol);

double dot(double *restrict v1, double *restrict v2, int nin, int nvol);
void dotc(double* ar, double* ai,
          double *restrict v1, double *restrict v2, int nin, int nvv);

void xI(double *restrict v1, int nin, int nvol);
void conjg(double *restrict v1, int nin, int nvol);

// real_t = float

void afield_init(float *data, const int size);
void afield_tidyup(float *data, const int size);
void afield_set(float *v,  float a, const int nin, const int nv2);
void copy_to_device(float *v, int nv);
void copy_to_device(float *v, int nv1, int nv);
void copy_from_device(float *v, int nv);
void copy_from_device(float *v, int nv1, int nv);

void convert(float *v, double *w, int nin, int nvol, int nvol_pad);
void reverse(double *v, float *w, int nin, int nvol, int nvol_pad);

void copy(float *v, float *w, int nin, int nvol);
void copy(float *v, int nv1, float *w, int nv2, int nin, int nvol);

void axpy(float *restrict v, int nv1, float a,
              float *restrict w, int nv2, int nin, int nvol);
void axpy(float *restrict v, int nv1, float ar, float ai,
              float *restrict w, int nv2, int nin, int nvol);

void aypx(float a, float *restrict v, int nv1,
          float *restrict w, int nv2, int nin, int nvol);
void aypx(float ar, float ai, float *restrict v, int nv1,
          float *restrict w, int nv2, int nin, int nvol);

void scal(float *restrict v, int nv1, float a, int nin, int nvol);
void scal(float *restrict v, int nv1, float ar, float ai,
                                                   int nin, int nvol);
float norm2(float *restrict v1, int nin, int nvol);

float dot(float *restrict v1, float *restrict v2, int nin, int nvol);
void dotc(float* ar, float* ai,
              float *restrict v1, float *restrict v2, int nin, int nvol);

void xI(float *restrict v1, int nin, int nvol);
void conjg(float *restrict v1, int nin, int nvol);

// copy with double/float conversion

void copy(double *v, int nv1, float  *w, int nv2, int nin, int nvol);
void copy(float  *v, int nv1, double *w, int nv2, int nin, int nvol);

//====================================================================

}

#endif
