/*!
      @file    bridgeVec_afield.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_AFIELD_H
#define BRIDGEVEC_AFIELD_H

namespace BridgeVec{

void afield_set(double *v,  double a, const int nin, const int nv2);
void convert_dev(double *v, double *w, int nin, int nvol, int nvol_pad);
void reverse_dev(double *v, double *w, int nin, int nvol, int nvol_pad);
void copy_dev(double *v, double *w, int nin, int nvol);
void copy_dev(double *v, int nv1, double *w, int nv2, int nin, int nvol);
void axpy_dev(double *v, int nv1, double a,
              double *w, int nv2, int nin, int nvol);
void axpy_dev(double *v, int nv1, double ar, double ai,
              double *w, int nv2, int nin, int nvol);
void aypx_dev(double *v, int nv1, double a,
              double *w, int nv2, int nin, int nvol);
void aypx_dev(double *v, int nv1, double ar, double ai,
              double *w, int nv2, int nin, int nvol);
void scal_dev(double* v, int nv1, double a, int nin, int nvol);
double dot_dev(double *v1, double *v2, int nin, int nvol);
double norm2_dev(double *v1, int nin, int nvol);
void dotc_dev(double* ar, double* ai,
              double *v1, double *v2, int nin, int nvv);
void xI_dev(double *v1, int nin, int nvol);
void conjg_dev(double *v1, int nin, int nvol);

//====================================================================
void afield_set(float *v,  float a, const int nin, const int nv2);
void convert_dev(float *v, double *w, int nin, int nvol, int nvol_pad);
void reverse_dev(double *v, float *w, int nin, int nvol, int nvol_pad);
void copy_dev(float *v, float *w, int nin, int nvol);
void copy_dev(float *v, int nv1, float *w, int nv2, int nin, int nvol);
void axpy_dev(float *v, int nv1, float a,
              float *w, int nv2, int nin, int nvol);
void axpy_dev(float *v, int nv1, float ar, float ai,
              float *w, int nv2, int nin, int nvol);
void aypx_dev(float *v, int nv1, float a,
              float *w, int nv2, int nin, int nvol);
void aypx_dev(float *v, int nv1, float ar, float ai,
              float *w, int nv2, int nin, int nvol);

void scal_dev(float* v, int nv1, float a, int nin, int nvol);
float dot_dev(float *v1, float *v2, int nin, int nvol);
float norm2_dev(float *v1, int nin, int nvol);
void dotc_dev(float* ar, float* ai,
              float *v1, float *v2, int nin, int nvol);
void xI_dev(float *v1, int nin, int nvol);
void conjg_dev(float *v1, int nin, int nvol);

}

#endif
//====================================================================
