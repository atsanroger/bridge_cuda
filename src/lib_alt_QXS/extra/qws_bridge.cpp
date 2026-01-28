/*!
        @file    qws_bridge.cpp
        @brief   extension of QWS functionalities for Bridge++
        @author  Issaku Kanamori
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2#$
        @version $LastChangedRevision: 2499 $
*/

#ifdef USE_QWSLIB

#include <qws.h>
#include <clover_d.h>
#include <clover_s.h>
#include <wilson_d.h>
#include <wilson_s.h>
#include <mpi.h>


#include "lib_alt_QXS/inline/define_vlen.h"

// for thread parallelization a la Bridge
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"


#ifdef __cplusplus
extern "C" {
#endif

extern __attribute__((aligned(64))) pglud_t glud;
extern __attribute__((aligned(64))) pglus_t glus;
extern __attribute__((aligned(64))) pclvd_t clvd;
extern __attribute__((aligned(64))) pclvs_t clvs;
extern double kappa, kappa2, mkappa;
extern int    nt, nz, ny, nx, nxh, nxd, vold, nxs, vols;
extern int    px, py, pz, pt;
extern int    domain_e, domain_o;
extern int    npe[4];
//  extern int volse;

// cf. qws.cc
//   nt, nz, ny, nx:  local lattice size
//   nxh = nx/2
//   nxd = nx/2/VLEND;
//   nxs = nx/2/VLENS;
//   vold = nxd * ny * nz * nt;   ( = local_vol /2/VLEND )
//   vols = nxs * ny * nz * nt;   ( = local_vol /2/VLENS )

int std_xyzt2i_(int *j)
{
  return j[0] + nx * (j[1] + ny * (j[2] + nz * j[3]));
}


int std_xyzt2kd_(int *j)     // in-vector index (double)
{
  int kx = j[0] % VLENXD;
  int ky = j[1] % VLENYD;
  return kx + VLENXD * ky;
}


int std_xyzt2ks_(int *j)     // in-vector index (single)
{
  int kx = j[0] % VLENXS;
  int ky = j[1] % VLENYS;
  return kx + VLENXS * ky;
}


int std_xyzt2ivd_(int *j)     // vector index (double)
{
  int ivx  = j[0] / VLENXD;
  int ivy  = j[1] / VLENYD;
  int nvx  = nx / VLENXD;
  int nvxy = nx * ny / VLEND;
  return ivx + nvx * ivy + nvxy * (j[2] + nz * j[3]);
}


int std_xyzt2ivs_(int *j)     // vector index (single)
{
  int ivx  = j[0] / VLENXS;
  int ivy  = j[1] / VLENYS;
  int nvx  = nx / VLENXS;
  int nvxy = nx * ny / VLENS;
  return ivx + nvx * ivy + nvxy * (j[2] + nz * j[3]);
}


int std_xvyzt2ivd_(int x, int ivy, int zt)     // vector index (double)
{
  int ivx  = x / VLENXD;
  int nvx  = nx / VLENXD;
  int nvxy = nx * ny / VLEND;
  return ivx + nvx * ivy + nvxy * zt;
}


int std_xky2kd_(int x, int ky)     // in-vector index (double)
{
  int ikx = x % VLENXD;
  return ikx + VLENXD * ky;
}


int std_xvyzt2ivs_(int x, int ivy, int zt)     // vector index (single)
{
  int ivx  = x / VLENXS;
  int nvx  = nx / VLENXS;
  int nvxy = nx * ny / VLENS;
  return ivx + nvx * ivy + nvxy * zt;
}


int std_xky2ks_(int x, int ky)     // in-vector index (single)
{
  int ikx = x % VLENXS;
  return ikx + VLENXS * ky;
}


int e_o_(int *j, int eo_offset)
{
  return (j[0] + j[1] + j[2] + j[3] + eo_offset) % 2;
}


void qws_loadg_bridge_(const double *u, const int *volh_tot, const double *in_kappa)
{
  int x, y, z, t, j[4], mu, c1, c2, eo_offset;
  kappa     = *in_kappa;
  mkappa    = -kappa;
  kappa2    = kappa * kappa;
  eo_offset = (px * nx + py * ny + pz * nz + pt * nt) % 2;

  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          j[0] = x;
          j[1] = y;
          j[2] = z;
          j[3] = t;

          // index in bridge (lexical)
          //   Note that input AField is in double prec.
          int ikd = std_xyzt2kd_(j);
          int ivd = std_xyzt2ivd_(j);
          //int iks = std_xyzt2ks_(j);
          //int ivs = std_xyzt2ivs_(j);

          // index in qws (even-odd)
          int eo   = e_o_(j, eo_offset);
          int iwsd = x / 2 / VLEND + nxd * y + nxd * ny * z + nxd * ny * nz * t + NDIM * vold * eo;
          int ixxd = (x / 2) % VLEND;
          int iwss = x / 2 / VLENS + nxs * y + nxs * ny * z + nxs * ny * nz * t + NDIM * vols * eo;
          int ixxs = (x / 2) % VLENS;

          for (mu = 0; mu < NDIM; mu++) {
#if SU3_RECONSTRUCT_D == 18
            for (c1 = 0; c1 < NCOL; ++c1) {
              for (c2 = 0; c2 < NCOL; ++c2) {
                //                  double val=u[ikd + VLEND*(0 + 2*c1 + 6*c2 + 18*(ivd + mu*vold*2))];
                //                  if(val>0.1){
                //                    printf("hoge: val=%15.7e: (x,y,z,t)=(%d,%d,%d,%d), eo=%d, mu=%d, c1=%d, c2=%d, ikd=%d, ixxd=%d, ivd=%d, iwsd=%d\n",
                //                           val, x,y,z,t, eo, mu, c1, c2, ikd, ixxd, ivd, iwsd);
                //                  }

                glud[iwsd + vold * mu].c[c2][c1][0][ixxd]
                  = u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
                glud[iwsd + vold * mu].c[c2][c1][1][ixxd]
                  = u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              }
            }
#elif SU3_RECONSTRUCT_D == 12
            for (c1 = 0; c1 < NCOL; ++c1) {
              for (c2 = 0; c2 < NCOL - 1; ++c2) {
                glud[iwsd + vold * mu].c[c2][c1][0][ixxd]
                  = u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
                glud[iwsd + vold * mu].c[c2][c1][1][ixxd]
                  = u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              }
            }
#endif
#if SU3_RECONSTRUCT_S == 18
            for (c1 = 0; c1 < NCOL; ++c1) {
              for (c2 = 0; c2 < NCOL; ++c2) {
                glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                  = (float)u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
                glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                  = (float)u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              }
            }
#elif SU3_RECONSTRUCT_S == 12
            for (c1 = 0; c1 < NCOL; ++c1) {
              for (c2 = 0; c2 < NCOL - 1; ++c2) {
                glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                  = (float)u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
                glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                  = (float)u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              }
            }
#endif
          }
        }
      }
    }
  }
}


void qws_loadfd_bridge_(scd_t *out, const double *in)
{
  int x, y, z, t, j[4], mu, c, s, eo_offset;
  eo_offset = (px * nx + py * ny + pz * nz + pt * nt) % 2;
  printf("hoge: %s: eo_offset=%d: (px,py,pz,pt)=(%d,%d,%d,%d): (nx,ny,nz,nt)=((%d,%d,%d,%d), VLENXD=%d\n",
         __func__, eo_offset, px, py, pz, pt, nx, ny, nz, nt, VLENXD);
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          j[0] = x;
          j[1] = y;
          j[2] = z;
          j[3] = t;

          // index in bridge (lexical)
          int ikd = std_xyzt2kd_(j);
          int ivd = std_xyzt2ivd_(j);

          // index in qws (even-odd)
          int eo = e_o_(j, eo_offset);
          if (eo == 0) {
            int iwsd = x / 2 / VLEND + nxd * y + nxd * ny * z + nxd * ny * nz * t;
            int ixxd = (x / 2) % VLEND;
            //              printf("hoge: (x,y,z,t)=(%d,%d,%d,%d):  ivd=%d, iwsd=%d, ixxd=%d, ikd=%d\n",
            //                     x,y,z,t, ivd, iwsd, ixxd, ikd);

            for (c = 0; c < NCOL; c++) {
              for (s = 0; s < 4; s++) {
                out[iwsd].c[c][s][0][ixxd]
                  = in[VLEND * (0 + 2 * s + 8 * c + 24 * ivd) + ikd];
                out[iwsd].c[c][s][1][ixxd]
                  = in[VLEND * (1 + 2 * s + 8 * c + 24 * ivd) + ikd];
              }   //s
            }  // c
          }
        }
      }
    }
  }
}


void qws_loadfs_bridge_(scs_t *out, const float *in)
{
  int x, y, z, t, j[4], mu, c, s, eo_offset;
  eo_offset = (px * nx + py * ny + pz * nz + pt * nt) % 2;

  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          j[0] = x;
          j[1] = y;
          j[2] = z;
          j[3] = t;

          // index in bridge (lexical)
          int iks = std_xyzt2ks_(j);
          int ivs = std_xyzt2ivs_(j);

          // index in qws (even-odd)
          int eo = e_o_(j, eo_offset);
          if (eo == 0) {
            int iwss = x / 2 / VLENS + nxs * y + nxs * ny * z + nxs * ny * nz * t;
            int ixxs = (x / 2) % VLENS;

            for (c = 0; c < NCOL; c++) {
              for (s = 0; s < 4; s++) {
                out[iwss].c[c][s][0][ixxs] = in[VLENS * (0 + 2 * s + 8 * c + 24 * ivs) + iks];
                out[iwss].c[c][s][1][ixxs] = in[VLENS * (1 + 2 * s + 8 * c + 24 * ivs) + iks];
              } //s
            }  // c
          }     // eo==0
        }
      }
    }
  }
}


void qws_storefd_bridge_(double *out, const scd_t *in)
{
  int x, y, z, t, j[4], mu, c, s, eo_offset;
  eo_offset = (px * nx + py * ny + pz * nz + pt * nt) % 2;

  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          j[0] = x;
          j[1] = y;
          j[2] = z;
          j[3] = t;

          // index in bridge (lexical)
          int ikd = std_xyzt2kd_(j);
          int ivd = std_xyzt2ivd_(j);

          // index in qws (even-odd)
          int eo = e_o_(j, eo_offset);
          if (eo == 0) {
            int iwsd = x / 2 / VLEND + nxd * y + nxd * ny * z + nxd * ny * nz * t;
            int ixxd = (x / 2) % VLEND;

            for (c = 0; c < NCOL; c++) {
              for (s = 0; s < 4; s++) {
                out[VLEND * (0 + 2 * s + 8 * c + 24 * ivd) + ikd]
                  = in[iwsd].c[c][s][0][ixxd];
                out[VLEND * (1 + 2 * s + 8 * c + 24 * ivd) + ikd]
                  = in[iwsd].c[c][s][1][ixxd];
              } //s
            }  // c
          }     // eo==0
        }
      }
    }
  }
}


void qws_storefs_bridge_(float *out, const scs_t *in)
{
  int x, y, z, t, j[4], mu, c, s, eo_offset;
  eo_offset = (px * nx + py * ny + pz * nz + pt * nt) % 2;

  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          j[0] = x;
          j[1] = y;
          j[2] = z;
          j[3] = t;

          // index in bridge (lexical)
          int iks = std_xyzt2ks_(j);
          int ivs = std_xyzt2ivs_(j);

          // index in qws (even-odd)
          int eo = e_o_(j, eo_offset);
          if (eo == 0) {
            int iwss = x / 2 / VLENS + nxs * y + nxs * ny * z + nxs * ny * nz * t;
            int ixxs = (x / 2) % VLENS;

            for (c = 0; c < NCOL; c++) {
              for (s = 0; s < 4; s++) {
                out[VLENS * (0 + 2 * s + 8 * c + 24 * ivs) + iks]
                  = in[iwss].c[c][s][0][ixxs];
                out[VLENS * (1 + 2 * s + 8 * c + 24 * ivs) + iks]
                  = in[iwss].c[c][s][1][ixxs];
              } //s
            }  // c
          }     // eo==0
        }
      }
    }
  }
}


void qws_loadg_dd_bridge_(const double *u, const double *in_kappa)
{
  int x, y, z, t, mu, c1, c2;
  kappa  = *in_kappa;
  mkappa = -kappa;
  kappa2 = kappa * kappa;
  int nvy = ny / VLENYD;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYD; ky++) {
      y = ky + VLENYD * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int ikd = std_xky2kd_(x, ky);
        int ivd = std_xvyzt2ivd_(x, vy, zt);

        // index in qws (domain decomposed)
        int iwsd = (x % nxh) / VLEND + nxd * y + nxd * ny * zt + NDIM * vold * (x / nxh);
        int ixxd = (x % nxh) % VLEND;
        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + NDIM * vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;
        for (mu = 0; mu < NDIM; mu++) {
#if SU3_RECONSTRUCT_D == 18
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL; c2++) {
              double val = u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              glud[iwsd + vold * mu].c[c2][c1][0][ixxd]
                = u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              glud[iwsd + vold * mu].c[c2][c1][1][ixxd]
                = u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
            }
          }
#elif SU3_RECONSTRUCT_D == 12
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL - 1; c2++) {
              glud[iwsd + vold * mu].c[c2][c1][0][ixxd]
                = u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              glud[iwsd + vold * mu].c[c2][c1][1][ixxd]
                = u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
            }
          }
#endif
#if SU3_RECONSTRUCT_S == 18
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL; c2++) {
              glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                = (float)u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                = (float)u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
            }
          }
#elif SU3_RECONSTRUCT_S == 12
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL - 1; c2++) {
              glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                = (float)u[ikd + VLEND * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
              glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                = (float)u[ikd + VLEND * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vold * 2))];
            }
          }
#endif
        } //mu
      }   //x
    }     // ky
  }       // yvtv
#pragma omp barrier
}


void qws_loadgs_dd_bridge_(const float *u, const double *in_kappa)
{
  int x, y, z, t, mu, c1, c2;
  kappa  = *in_kappa;
  mkappa = -kappa;
  kappa2 = kappa * kappa;
  int nvy = ny / VLENYS;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYS; ky++) {
      y = ky + VLENYS * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int iks = std_xky2ks_(x, ky);
        int ivs = std_xvyzt2ivs_(x, vy, zt);

        // index in qws (domain decomposed)
        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + NDIM * vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;
        for (mu = 0; mu < NDIM; mu++) {
#if SU3_RECONSTRUCT_S == 18
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL; c2++) {
              glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                = u[iks + VLENS * (0 + 2 * c1 + 6 * c2 + 18 * (ivs + mu * vols * 2))];
              glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                = u[iks + VLENS * (1 + 2 * c1 + 6 * c2 + 18 * (ivs + mu * vols * 2))];
            }
          }
#elif SU3_RECONSTRUCT_S == 12
          for (c1 = 0; c1 < NCOL; c1++) {
            for (c2 = 0; c2 < NCOL - 1; c2++) {
              glus[iwss + vols * mu].c[c2][c1][0][ixxs]
                = u[iws + VLENS * (0 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vols * 2))];
              glus[iwss + vols * mu].c[c2][c1][1][ixxs]
                = u[iws + VLENS * (1 + 2 * c1 + 6 * c2 + 18 * (ivd + mu * vols * 2))];
            }
          }
#endif
        } // mu
      }   // x
    }     // ky
  }       // yztv
#pragma omp barrier
}


void qws_loadc_dd_bridge_(const double *c)
{
  // notation in qws
  //   c[2][36]   2 x (6x6 Hermitian matrix)
  //     c[ud][0], c[ud][1],..., c[ud][5]:  daigonal elements (real): C_ii
  //     (c[ud][6], c[ud][7]),..., (c[ud][14],c[ud][15])     C01, C02, C03, C04, C05
  //     (c[ud][16], c[ud][17]),..., (c[ud][22],c[ud][23])   C12, C13, C14, C15
  //     (c[ud][24], c[ud][25]),..., (c[ud][28],c[ud][29])   C23, C24, C25
  //     (c[ud][30], c[ud][31]),(c[ud][32],c[ud][33])        C34, C35
  //     (c[ud][34], c[ud][35] )                             C45
  // (spin, color) --> 3spin +color = 0,1,...,5
  // IMPORTANT:
  //   the stored colver term is in the Chiral representation while the fermion field
  //   is in the Dirac represenation.  The conversion of represenation of the fermion
  //   feild requires a factor 1/sqrt(2) for each of Dirac --> Chiral and Chiral --> Dirac,
  //   i.e., one needs to mulitply a factor 1/2.  This 1/2 multiplication is moved to
  //   the Clover term.
  //       M^dag Clov[Chiral] M psi[Dirac]
  //        = ( sqrt(2) M^dag ) Clov[Chial]/2  (sqrt(2) M) psi[Dirac]
  //   and we store Clov[Chiral]/2
  //
  //   0  6  8 10 12 14
  //      1 16 18 20 22
  //         2 24 26 28
  //            3 30 32
  //               4 34
  //                  5

  int x, y, ky, yztv, ud;
  int nvy = ny / VLENYD;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (ky = 0; ky < VLENYD; ky++) {
      y = ky + VLENYD * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int ikd = std_xky2kd_(x, ky);
        int ivd = std_xvyzt2ivd_(x, vy, zt);

        int iwsd = (x % nxh) / VLEND + nxd * y + nxd * ny * zt + vold * (x / nxh);
        int ixxd = (x % nxh) % VLEND;
        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;

        for (ud = 0; ud < 2; ud++) {
          int offset = 72 * ivd + 36 * ud;
          for (int ij = 0; ij < 36; ij++) {
            clvd[iwsd].c[ud][ij][ixxd] = c[VLEND * (ij + offset) + ikd];
          }

          // single
          for (int ij = 0; ij < 36; ij++) {
            clvs[iwss].c[ud][ij][ixxs] = (float)c[VLEND * (ij + offset) + ikd];
          }
        }   // ud
      }
    }
  }
}


void qws_loadcs_dd_bridge_(const float *c)
{
  // notation in qws
  //   c[2][36]   2 x (6x6 Hermitian matrix)
  //     c[ud][0], c[ud][1],..., c[ud][5]:  daigonal elements (real): C_ii
  //     (c[ud][6], c[ud][7]),..., (c[ud][14],c[ud][15])     C01, C02, C03, C04, C05
  //     (c[ud][16], c[ud][17]),..., (c[ud][22],c[ud][23])   C12, C13, C14, C15
  //     (c[ud][24], c[ud][25]),..., (c[ud][28],c[ud][29])   C23, C24, C25
  //     (c[ud][30], c[ud][31]),(c[ud][32],c[ud][33])        C34, C35
  //     (c[ud][34], c[ud][35] )                             C45

  // a naive 6x6 notation as real varibles
  //   0  2  4  6  8 10
  //  12 14 16 18 20 22
  //  24 26 28 30 32 34
  //  36 38 40 42 44 46
  //  48 50 52 52 56 58
  //  60 62 64 66 68 70

  int x, y, ky, yztv, ud;
  int nvy = ny / VLENYS;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (ky = 0; ky < VLENYS; ky++) {
      y = ky + VLENYS * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int iks = std_xky2ks_(x, ky);
        int ivs = std_xvyzt2ivs_(x, vy, zt);

        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;

        for (ud = 0; ud < 2; ud++) {
          int offset = 72 * ivs + 36 * ud;
          for (int ij = 0; ij < 36; ij++) {
            clvs[iwss].c[ud][ij][ixxs] = c[VLENS * (ij + offset) + iks];
          }
        }   // ud
      }
    }
  }
}


void qws_loadfd_dd_bridge_(scd_t *out, const double *in)
{
  int x, y, z, t, s, c;
  int nvy = ny / VLENYD;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYD; ky++) {
      y = ky + VLENYD * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int ikd = std_xky2kd_(x, ky);
        int ivd = std_xvyzt2ivd_(x, vy, zt);


        // index in qws (domain decomposed)
        int iwsd = (x % nxh) / VLEND + nxd * y + nxd * ny * zt + vold * (x / nxh);
        int ixxd = (x % nxh) % VLEND;
        for (c = 0; c < NCOL; c++) {
          for (s = 0; s < 4; s++) {
            out[iwsd].c[c][s][0][ixxd]
              = in[VLEND * (0 + 2 * s + 8 * c + 24 * ivd) + ikd];
            out[iwsd].c[c][s][1][ixxd]
              = in[VLEND * (1 + 2 * s + 8 * c + 24 * ivd) + ikd];
          } //s
        }  // c
      }     //x
    }       // ky
  }         // yztv
#pragma omp barrier
}


void qws_storefd_dd_bridge_(double *out, const scd_t *in)
{
  int x, y, z, t, s, c;
  int nvy = ny / VLENYD;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYD; ky++) {
      y = ky + VLENYD * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int ikd = std_xky2kd_(x, ky);
        int ivd = std_xvyzt2ivd_(x, vy, zt);

        // index in qws (domain decomposed)
        int iwsd = (x % nxh) / VLEND + nxd * y + nxd * ny * zt + vold * (x / nxh);
        int ixxd = (x % nxh) % VLEND;
        for (c = 0; c < NCOL; c++) {
          for (s = 0; s < 4; s++) {
            out[VLEND * (0 + 2 * s + 8 * c + 24 * ivd) + ikd] = in[iwsd].c[c][s][0][ixxd];
            out[VLEND * (1 + 2 * s + 8 * c + 24 * ivd) + ikd] = in[iwsd].c[c][s][1][ixxd];
          } //s
        }  // c
      }     // x
    }       // ky
  }         // yztv
#pragma omp barrier
}


void qws_loadfs_dd_bridge_(scs_t *out, const float *in)
{
  int x, y, z, t, s, c;
  int nvy = ny / VLENYS;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYS; ky++) {
      y = ky + VLENYS * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int iks = std_xky2ks_(x, ky);
        int ivs = std_xvyzt2ivs_(x, vy, zt);

        // index in qws (domain decomposed)
        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;
        for (c = 0; c < NCOL; c++) {
          for (s = 0; s < 4; s++) {
            out[iwss].c[c][s][0][ixxs]
              = in[VLENS * (0 + 2 * s + 8 * c + 24 * ivs) + iks];
            out[iwss].c[c][s][1][ixxs]
              = in[VLENS * (1 + 2 * s + 8 * c + 24 * ivs) + iks];
          }  //s
        }    // c
      }
    }
  }         // x,ky,yztv
#pragma omp barrier
}


void qws_storefs_dd_bridge_(float *out, const scs_t *in)
{
  int x, y, z, t, s, c;
  int nvy = ny / VLENYS;
  int nyztv = nt * nz * nvy;
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, nyztv);
#pragma omp barrier
  for (int yztv = is; yztv < ns; yztv++) {
    int vy = yztv % nvy;
    int zt = yztv / nvy;
    for (int ky = 0; ky < VLENYS; ky++) {
      y = ky + VLENYS * vy;
      for (x = 0; x < nx; x++) {
        // index in bridge (lexical)
        int iks = std_xky2ks_(x, ky);
        int ivs = std_xvyzt2ivs_(x, vy, zt);

        // index in qws (domain decomposed)
        int iwss = (x % nxh) / VLENS + nxs * y + nxs * ny * zt + vols * (x / nxh);
        int ixxs = (x % nxh) % VLENS;
        for (c = 0; c < NCOL; c++) {
          for (s = 0; s < 4; s++) {
            out[VLENS * (0 + 2 * s + 8 * c + 24 * ivs) + iks] = in[iwss].c[c][s][0][ixxs];
            out[VLENS * (1 + 2 * s + 8 * c + 24 * ivs) + iks] = in[iwss].c[c][s][1][ixxs];
          }   //s
        }  // c
      }
    }
  }         // x,ky,yztv
#pragma omp barrier
}


//---------------------------------------------------------------------------------------- for DD solver
void ddd_out_pre_d_(scd_t *in, int *domain);
void ddd_out_pos_d_(scd_t *out, scd_t *in, int *domain);
void ddd_in_d_(scd_t *out, scd_t *in, int *domain);

void ddd_out_pre_s_(scs_t *in, int *domain);
void ddd_out_pre_s_noprl_(scs_t *in, int *domain);
void ddd_out_pre_s_no_timer_(scs_t *in, int *domain);
void ddd_out_pre_s_noprl_no_timer_(scs_t *in, int *domain);
void ddd_out_pos_s_(scs_t *out, scs_t *in, int *domain, float factor);
void ddd_out_pos_s_no_timer_(scs_t *out, scs_t *in, int *domain, float factor);
void ddd_in_s_(scs_t *out, scs_t *in, int *domain);
void ddd_out_pos_s_noprl_(scs_t *out, scs_t *in, int *domain, float factor);
void ddd_out_pos_s_noprl_no_timer_(scs_t *out, scs_t *in, int *domain, float factor);
void ddd_in_s_noprl(scs_t *out, scs_t *in, int *domain);
void jinv_ddd_in_s_noprl_(scs_t *x, scs_t *b, int *domain, int *maxiter);

//----------------------------------------------------------------------------------------

/*
void ddd_d_(scd_t* out, scd_t* in){
  printf("hoge: %s: kappa=%f\n", __func__,kappa); fflush(stdout);
  ddd_out_pre_d_(in, &domain0);
  ddd_in_d_(     &out[vold*0], &in[vold*0], &domain0);
  ddd_out_pos_d_(&out[vold*0], &in[vold*0], &domain0);
  ddd_out_pre_d_(in, &domain1);
  ddd_in_d_(     &out[vold*1], &in[vold*1], &domain1);
  ddd_out_pos_d_(&out[vold*1], &in[vold*0], &domain1);
}
*/
//----------------------------------------------------------------------------------------
void ddd_eo_s_noprl_(scs_t *out, scs_t *in, int *domain)
{
  ddd_out_pre_s_noprl_no_timer_(in, domain);
  ddd_in_s_noprl(&out[vols * (*domain)], &in[vols * (*domain)], domain);
  ddd_out_pos_s_noprl_no_timer_(&out[vols * (*domain)], in, domain, (float)mkappa);
}


//----------------------------------------------------------------------------------------
void ddd_s_noprl_(scs_t *out, scs_t *in)
{
  ddd_eo_s_noprl_(out, in, &domain_e);
  ddd_eo_s_noprl_(out, in, &domain_o);
}


// inline functions defiend in qws.c
extern void copy_vols2_s_noprl_(scs_t *out, scs_t *in);
extern void zero_vols_s_noprl_(scs_t *out);
extern void accum_add_vols_s_noprl_(scs_t *out, scs_t *in);
extern void accum_sub_vols_s_noprl_(scs_t *out, scs_t *in);
extern void accum_addsub_vols_s_noprl_(scs_t *out, scs_t *in1, scs_t *in2);

//----------------------------------------------------------------------------------------
void prec_s_noprl_(scs_t *out, scs_t *in, int *nsap, int *nm, scs_t *s, scs_t *q)
{
#ifdef COMPILE_TIME_DIM_SIZE
  const int vols = VOLS;
#else
  const int vols = ::vols;
#endif
  float kappa = ::kappa;

  {
    copy_vols2_s_noprl_(s, in);
    if ((npe[1] == 1) || (npe[2] == 1) || (npe[3] == 1)) zero_vols_s_noprl_(&q[vols * domain_o]);
    for (int isap = 0; isap < *nsap; isap++) {
      jinv_ddd_in_s_noprl_(&q[vols * domain_e], &s[vols * domain_e], &domain_e, nm);
      ddd_out_pre_s_noprl_no_timer_(q, &domain_o);
      ddd_in_s_noprl(&q[vols * domain_o], &q[vols * domain_e], &domain_e);
      accum_addsub_vols_s_noprl_(&s[vols * domain_e], &in[vols * domain_e], &q[vols * domain_o]);
      ddd_out_pos_s_noprl_no_timer_(&s[vols * domain_o], q, &domain_o, (float)kappa);
      jinv_ddd_in_s_noprl_(&q[vols * domain_o], &s[vols * domain_o], &domain_o, nm);
      ddd_out_pre_s_noprl_no_timer_(q, &domain_e);
      ddd_in_s_noprl(&q[vols * domain_e], &q[vols * domain_o], &domain_o);
      accum_addsub_vols_s_noprl_(&s[vols * domain_o], &in[vols * domain_o], &q[vols * domain_e]);
      ddd_out_pos_s_noprl_no_timer_(&s[vols * domain_e], q, &domain_e, (float)kappa);
    }
    if ((npe[1] == 1) || (npe[2] == 1) || (npe[3] == 1)) zero_vols_s_noprl_(&out[vols * domain_o]);
    jinv_ddd_in_s_noprl_(&out[vols * domain_e], &s[vols * domain_e], &domain_e, nm);
    ddd_out_pre_s_noprl_no_timer_(out, &domain_o);
    ddd_out_pos_s_noprl_no_timer_(&s[vols * domain_o], out, &domain_o, (float)kappa);
    jinv_ddd_in_s_noprl_(&out[vols * domain_o], &s[vols * domain_o], &domain_o, nm);
  }
}


//---------------------------------------------------------------------
void clv_s_dd_(int *deo, scs_t *inout)
{
  __attribute__((aligned(64))) scs_t tmp;
  int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nxs; x++) {
          int i0 = x + nxs * y + nxs * ny * z + nxs * ny * nz * t + vols * (*deo);
          __load_sc_s(tmp.c, inout[i0].c);
          __mult_clvs(tmp.cv, clvs[i0].cv);
          __store_sc_s(inout[i0].c, tmp.c);
        }
      }
    }
  }
}


void clv_matvec_s_dd_(int *deo, scs_t *out, clvs_t *clv, scs_t *in)
{
  __attribute__((aligned(64))) scs_t tmp;
  int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nxs; x++) {
          int i0 = x + nxs * y + nxs * ny * z + nxs * ny * nz * t + vols * (*deo);
          __load_sc_s(tmp.c, in[i0].c);
          __mult_clvs(tmp.cv, clv[i0].cv);
          __store_sc_s(out[i0].c, tmp.c);
        }
      }
    }
  }
}


//--------------------------------------------------------------------
void clv_d_dd_(int *deo, scd_t *inout)
{
  __attribute__((aligned(64))) scd_t tmp;
  int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nxd; x++) {
          int i0 = x + nxd * y + nxd * ny * z + nxd * ny * nz * t + vold * (*deo);
          __load_sc(tmp.c, inout[i0].c);
          __mult_clvd(tmp.cv, clvd[i0].cv);
          __store_sc(inout[i0].c, tmp.c);
        }
      }
    }
  }
}


void clv_matvec_d_dd_(int *deo, scd_t *out, clvd_t *clv, scd_t *in)
{
  __attribute__((aligned(64))) scd_t tmp;
  int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nxd; x++) {
          int i0 = x + nxd * y + nxd * ny * z + nxd * ny * nz * t + vold * (*deo);
          __load_sc(tmp.c, in[i0].c);
          __mult_clvd(tmp.cv, clv[i0].cv);
          __store_sc(out[i0].c, tmp.c);
        }
      }
    }
  }
}


#ifdef __cplusplus
}  // extern "C"
#endif

#endif
//==============================================================END===
