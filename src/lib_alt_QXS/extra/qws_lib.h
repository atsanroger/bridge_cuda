#ifndef QWS_LIB_H
#define QWS_LIB_H

#ifdef USE_QWSLIB

#include <qws.h>

namespace QWS_lib {
  extern int qws_is_setup;
  extern Bridge::VerboseLevel vl;
  void init_qws(const int *boundary);
  void tidyup_qws();
}

// declaration of functions in qws and its extension in qws_bridge
extern "C" {
extern int nt, nz, ny, nx;
extern int vols, vold;
// functions defined in qws_lib.cpp:
//  int std_xyzt2i_(int* j);
//  int e_o_(int* j);

// functions in libqws.a:
void qws_init_(int *lx, int *ly, int *lz, int *lt,
               int *npe_f, int *fbc_f, int *pce_f, int *pco_f,
               int *block_size);

void qws_finalize_();


void qws_loadg_bridge_(const double *u, const int *const volh_tot, const double *in_kappa);

//  void qws_loadclov_bridge_(const double *u, const int* const volh_tot, const double* in_kappa);


// set gauge field with domain decompsed indexing
void qws_loadg_dd_bridge_(const double *u, const double *in_kappa);
void qws_loadgs_dd_bridge_(const float *u, const double *in_kappa);

// inverse clover term
void qws_loadc_dd_bridge_(const double *c);
void qws_loadcs_dd_bridge_(const float *c);

// set fermion field with even-odd indexing
void qws_loadfd_bridge_(scd_t *out, const double *in);
void qws_loadfs_bridge_(scs_t *out, const float *in);

// get fermion gauge field domain decompsed indexing
void qws_storefd_bridge_(double *out, const scd_t *in);
void qws_storefs_bridge_(float *out, const scs_t *in);

// set fermion field with domain decompsed indexing
void qws_loadfd_dd_bridge_(scd_t *out, const double *in);
void qws_loadfs_dd_bridge_(scs_t *out, const float *in);

// get fermion gauge field domain decompsed indexing
void qws_storefd_dd_bridge_(double *out, const scd_t *in);
void qws_storefs_dd_bridge_(float *out, const scs_t *in);

void deo_vm_(int *pe, int *po, scd_t *out, scd_t *in);

void deo_dag_vm_(int *pe, int *po, scd_t *out, scd_t *in);

void deo_s_(int *pe, int *po, scs_t *out, scs_t *in);

// mult with domain deomposed data layout
void ddd_s_noprl_(scs_t * out, scs_t *in);
void ddd_d_(scd_t * out, scd_t *in);

void prec_s_noprl_(scs_t * out, scs_t *in, int *nsap, int *nm, scs_t *s, scs_t *q);

void clv_s_dd_(int *deo, scs_t *inout);
void clv_d_dd_(int *deo, scd_t *inout);

// not tested
void clv_matvec_s_dd_(int *deo, scs_t *out, clvs_t *clv, scs_t *in);
void clv_matvec_d_dd_(int *deo, scd_t *out, clvd_t *clv, scd_t *in);

void bicgstab_dd_mix_(scd_t *x, scd_t *b, double *tol, int *conviter, int *maxiter, double *tol_s, int *maxiter_s, int *nsap, int *nm);
}

#endif

#endif
