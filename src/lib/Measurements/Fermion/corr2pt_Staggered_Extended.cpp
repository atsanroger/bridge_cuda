/*!
        @file    corr2pt_Staggered_Extended.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "corr2pt_Staggered_Extended.h"
#include "Communicator/communicator.h"
using Bridge::vout;

#if defined USE_GROUP_SU3
#define NC     3
#define NC2    6
#define C1     0
#define C2     2
#define C3     4
#elif defined USE_GROUP_SU2
#define NC     2
#define NC2    4
#endif

const std::string Corr2pt_Staggered_Extended::class_name
  = "Corr2pt_Staggered_Extended";

//====================================================================
void Corr2pt_Staggered_Extended::init()
{
  int Nc = CommonParameters::Nc();

  if (Nc == 3) {
    set_asymtensor_Nc3();
  } else {
    vout.crucial(m_vl, "%s: Nc = %d is not implemented.\n", Nc);
    exit(EXIT_FAILURE);
  }

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();

  // Note that local lattice sizes in x,y,z directions must be even.
  int k = (Nx % 2) + (Ny % 2) + (Nz % 2);
  if (k > 0) {
    vout.crucial("%s: local spatial lattice sizes must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_filename_output = "stdout";
}


//====================================================================
void Corr2pt_Staggered_Extended::meson_all(
  std::vector<double>& meson,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2)
{
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  int Lt       = CommonParameters::Lt();
  int Ndim     = CommonParameters::Ndim();
  int Ndim_spc = Ndim - 1;
  if (meson.size() != Lt) meson.resize(Lt);

  std::vector<dcomplex> corr(Lt);
  std::vector<dcomplex> corr1(Lt);
  std::vector<dcomplex> corr2(Lt);
  std::vector<double>   corr_r(Lt);

  std::vector<int> shft(Ndim_spc);

  vout.general(m_vl, "(0,0,0) <-- wall (qq+oo) correlator:\n");
  shft[0] = 0;
  shft[1] = 0;
  shft[2] = 0;
  meson_correlator(corr1, shft, sq1, sq2, 0, 0);
  meson_correlator(corr2, shft, sq1, sq2, 1, 1);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  // correlator returned to main, preserving previous definition.
  for (int t = 0; t < corr.size(); ++t) {
    meson[t] = real(corr1[t]);
  }

  vout.general(m_vl, "(0,1,1) <-- wall (qq+oo) correlator:\n");
  shft[0] = 0;
  shft[1] = 1;
  shft[2] = 1;
  meson_correlator(corr1, shft, sq1, sq2, 0, 0);
  meson_correlator(corr2, shft, sq1, sq2, 1, 1);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(1,0,1) <-- wall (qq+oo) correlator:\n");
  shft[0] = 1;
  shft[1] = 0;
  shft[2] = 1;
  meson_correlator(corr1, shft, sq1, sq2, 0, 0);
  meson_correlator(corr2, shft, sq1, sq2, 1, 1);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(1,1,0) <-- wall (qq+oo) correlator:\n");
  shft[0] = 1;
  shft[1] = 1;
  shft[2] = 0;
  meson_correlator(corr1, shft, sq1, sq2, 0, 0);
  meson_correlator(corr2, shft, sq1, sq2, 1, 1);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(0,0,0) <-- wall (qo+oq) correlator:\n");
  shft[0] = 0;
  shft[1] = 0;
  shft[2] = 0;
  meson_correlator(corr1, shft, sq1, sq2, 0, 1);
  meson_correlator(corr2, shft, sq1, sq2, 1, 0);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(0,1,1) <-- wall (qo+oq) correlator:\n");
  shft[0] = 0;
  shft[1] = 1;
  shft[2] = 1;
  meson_correlator(corr1, shft, sq1, sq2, 0, 1);
  meson_correlator(corr2, shft, sq1, sq2, 1, 0);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(1,0,1) <-- wall (qo+oq) correlator:\n");
  shft[0] = 1;
  shft[1] = 0;
  shft[2] = 1;
  meson_correlator(corr1, shft, sq1, sq2, 0, 1);
  meson_correlator(corr2, shft, sq1, sq2, 1, 0);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "(1,1,0) <-- wall (qo+oq) correlator:\n");
  shft[0] = 1;
  shft[1] = 1;
  shft[2] = 0;
  meson_correlator(corr1, shft, sq1, sq2, 0, 1);
  meson_correlator(corr2, shft, sq1, sq2, 1, 0);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void Corr2pt_Staggered_Extended::baryon_all(
  std::vector<double>& baryon,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2)
{
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  int Lt       = CommonParameters::Lt();
  int Ndim     = CommonParameters::Ndim();
  int Ndim_spc = Ndim - 1;

  if (baryon.size() != Lt) baryon.resize(Lt);

  std::vector<dcomplex> corr(Lt);
  std::vector<dcomplex> corr1(Lt);
  std::vector<dcomplex> corr2(Lt);

  std::vector<double> prty(Ndim_spc);

  vout.general(m_vl, "nucleon <-- wall (qqq+ooo) correlator:\n");

  nucleon_correlator(corr1, sq1, sq2, 0, 0);
  nucleon_correlator(corr2, sq1, sq2, 1, 1);
  for (int t = 0; t < corr.size(); ++t) {
    corr[t] = corr1[t] + corr2[t];
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  for (int t = 0; t < corr.size(); ++t) {
    baryon[t] = real(corr[t]);
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void Corr2pt_Staggered_Extended::meson_correlator(
  std::vector<dcomplex>& meson,
  const std::vector<int>& shft,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2,
  const int isrc1, const int isrc2)
{
  int Nc    = CommonParameters::Nc();
  int Nvol  = CommonParameters::Nvol();
  int Lt    = CommonParameters::Lt();
  int Nx    = CommonParameters::Nx();
  int Ny    = CommonParameters::Ny();
  int Nz    = CommonParameters::Nz();
  int Nt    = CommonParameters::Nt();
  int Nx2   = Nx / 2;
  int Ny2   = Ny / 2;
  int Nz2   = Nz / 2;
  int Nspc2 = Nx2 * Ny2 * Nz2;
  int Neta  = 8;

  Field corrF(2, Nvol, 1);
  int   ex = 0;

  for (int it = 0; it < Nt; ++it) {
    for (int site2 = 0; site2 < Nspc2; ++site2) {
      int y1 = site2 % Nx2;
      int y2 = (site2 / Nx2) % Ny2;
      int y3 = site2 / (Nx2 * Ny2);

      for (int eta = 0; eta < Neta; ++eta) {
        int eta1  = eta % 2;
        int eta2  = (eta / 2) % 2;
        int eta3  = eta / 4;
        int etap1 = (eta1 + shft[0]) % 2;
        int etap2 = (eta2 + shft[1]) % 2;
        int etap3 = (eta3 + shft[2]) % 2;

        int x1     = 2 * y1 + eta1;
        int x2     = 2 * y2 + eta2;
        int x3     = 2 * y3 + eta3;
        int site_x = m_index.site(x1, x2, x3, it);

        int z1     = 2 * y1 + etap1;
        int z2     = 2 * y2 + etap2;
        int z3     = 2 * y3 + etap3;
        int site_z = m_index.site(z1, z2, z3, it);

        double corr_r = 0.0;
        double corr_i = 0.0;
        for (int c0 = 0; c0 < Nc; ++c0) {
          for (int c1 = 0; c1 < Nc; ++c1) {
            corr_r += sq1[c0 + Nc * isrc1].cmp_r(c1, site_x, ex)
                      * sq2[c0 + Nc * isrc2].cmp_r(c1, site_z, ex)
                      + sq1[c0 + Nc * isrc1].cmp_i(c1, site_x, ex)
                      * sq2[c0 + Nc * isrc2].cmp_i(c1, site_z, ex);
            corr_i += sq1[c0 + Nc * isrc1].cmp_r(c1, site_x, ex)
                      * sq2[c0 + Nc * isrc2].cmp_i(c1, site_z, ex)
                      - sq1[c0 + Nc * isrc1].cmp_i(c1, site_x, ex)
                      * sq2[c0 + Nc * isrc2].cmp_r(c1, site_z, ex);
          }
        }
        corrF.set(0, site_x, ex, corr_r);
        corrF.set(1, site_x, ex, corr_i);
      }
    }
  }

  std::vector<dcomplex> corr_local(Nt);

  for (int t = 0; t < Nt; ++t) {
    corr_local[t] = cmplx(0.0, 0.0);

    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int    site = m_index.site(x, y, z, t);
          double cr   = corrF.cmp(0, site, ex);
          double ci   = corrF.cmp(1, site, ex);
          corr_local[t] += cmplx(cr, ci);
        }
      }
    }
  }

  int ipet = Communicator::ipe(3);
  std::vector<dcomplex> corr_tmp(Lt);

  for (int itime = 0; itime < Lt; ++itime) {
    corr_tmp[itime] = cmplx(0.0, 0.0);
  }

  for (int it = 0; it < Nt; ++it) {
    int itime = it + ipet * Nt;
    corr_tmp[itime] = corr_local[it];
  }

  for (int itime = 0; itime < Lt; ++itime) {
    double crr  = real(corr_tmp[itime]);
    double cri  = imag(corr_tmp[itime]);
    double crr2 = Communicator::reduce_sum(crr);
    double cri2 = Communicator::reduce_sum(cri);
    meson[itime] = cmplx(crr2, cri2);
  }
}


//====================================================================
void Corr2pt_Staggered_Extended::nucleon_correlator(
  std::vector<dcomplex>& baryon,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2,
  const int isrc1, const int isrc2)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();
  int Nt   = CommonParameters::Nt();
  int ipet = Communicator::ipe(3);

  if (baryon.size() != Lt) baryon.resize(Lt);

  std::vector<dcomplex> corr_local(Lt);

  int FactNc = 6;  // this is specific to Nc=3 case.

  for (int itime = 0; itime < Lt; ++itime) {
    corr_local[itime] = cmplx(0.0, 0.0);
  }

  for (int it = 0; it < Nt; ++it) {
    int itime = it + Nt * ipet;

    for (int icp = 0; icp < FactNc; ++icp) {
      int ic1 = epsilon_index(0, icp) + Nc * isrc1;
      int ic2 = epsilon_index(1, icp) + Nc * isrc1;
      int ic3 = epsilon_index(2, icp) + Nc * isrc2;

      dcomplex sum1;
      contract_at_t(sum1, sq1[ic1], sq1[ic2], sq2[ic3], it);

      corr_local[itime] += epsilon_value(icp) * sum1;
    }
  }

  for (int itime = 0; itime < Lt; ++itime) {
    double crr  = real(corr_local[itime]);
    double cri  = imag(corr_local[itime]);
    double crr2 = Communicator::reduce_sum(crr);
    double cri2 = Communicator::reduce_sum(cri);
    baryon[itime] = cmplx(crr2, cri2);
  }
}


//====================================================================
void Corr2pt_Staggered_Extended::contract_at_t(dcomplex& corr,
                                               const Field_F_1spinor& v1,
                                               const Field_F_1spinor& v2,
                                               const Field_F_1spinor& v3,
                                               const int it)
{   // for staggered nucleon
  corr = cmplx(0.0, 0.0);

#if defined USE_GROUP_SU3
  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();

  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);

  const double *w1 = v1.ptr(0);
  const double *w2 = v2.ptr(0);
  const double *w3 = v3.ptr(0);

  double cr = 0.0;
  double ci = 0.0;

  for (int iz = 0; iz < Nz; ++iz) {
    int z = iz + Nz * ipez;
    if (z % 2 == 1) continue;
    for (int iy = 0; iy < Ny; ++iy) {
      int y = iy + Ny * ipey;
      if (y % 2 == 1) continue;
      for (int ix = 0; ix < Nx; ++ix) {
        int x = ix + Nx * ipex;
        if (x % 2 == 1) continue;

        int site = NC2 * m_index.site(ix, iy, iz, it);

        int ic11r = C1 + site;
        int ic22r = C2 + site;
        int ic33r = C3 + site;

        int ic11i = C1 + 1 + site;
        int ic22i = C2 + 1 + site;
        int ic33i = C3 + 1 + site;

        int ic21r = C2 + site;
        int ic32r = C3 + site;
        int ic13r = C1 + site;

        int ic21i = C2 + 1 + site;
        int ic32i = C3 + 1 + site;
        int ic13i = C1 + 1 + site;

        int ic31r = C3 + site;
        int ic12r = C1 + site;
        int ic23r = C2 + site;

        int ic31i = C3 + 1 + site;
        int ic12i = C1 + 1 + site;
        int ic23i = C2 + 1 + site;

        cr +=
          (w1[ic11r] * w2[ic22r] - w1[ic11i] * w2[ic22i]) * w3[ic33r]
          - (w1[ic11r] * w2[ic22i] + w1[ic11i] * w2[ic22r]) * w3[ic33i];

        ci +=
          (w1[ic11r] * w2[ic22r] - w1[ic11i] * w2[ic22i]) * w3[ic33i]
          + (w1[ic11r] * w2[ic22i] + w1[ic11i] * w2[ic22r]) * w3[ic33r];

        cr +=
          (w1[ic21r] * w2[ic32r] - w1[ic21i] * w2[ic32i]) * w3[ic13r]
          - (w1[ic21r] * w2[ic32i] + w1[ic21i] * w2[ic32r]) * w3[ic13i];

        ci +=
          (w1[ic21r] * w2[ic32r] - w1[ic21i] * w2[ic32i]) * w3[ic13i]
          + (w1[ic21r] * w2[ic32i] + w1[ic21i] * w2[ic32r]) * w3[ic13r];

        cr +=
          (w1[ic31r] * w2[ic12r] - w1[ic31i] * w2[ic12i]) * w3[ic23r]
          - (w1[ic31r] * w2[ic12i] + w1[ic31i] * w2[ic12r]) * w3[ic23i];

        ci +=
          (w1[ic31r] * w2[ic12r] - w1[ic31i] * w2[ic12i]) * w3[ic23i]
          + (w1[ic31r] * w2[ic12i] + w1[ic31i] * w2[ic12r]) * w3[ic23r];

        cr -=
          (w1[ic31r] * w2[ic22r] - w1[ic31i] * w2[ic22i]) * w3[ic13r]
          - (w1[ic31r] * w2[ic22i] + w1[ic31i] * w2[ic22r]) * w3[ic13i];

        ci -=
          (w1[ic31r] * w2[ic22r] - w1[ic31i] * w2[ic22i]) * w3[ic13i]
          + (w1[ic31r] * w2[ic22i] + w1[ic31i] * w2[ic22r]) * w3[ic13r];

        cr -=
          (w1[ic21r] * w2[ic12r] - w1[ic21i] * w2[ic12i]) * w3[ic33r]
          - (w1[ic21r] * w2[ic12i] + w1[ic21i] * w2[ic12r]) * w3[ic33i];

        ci -=
          (w1[ic21r] * w2[ic12r] - w1[ic21i] * w2[ic12i]) * w3[ic33i]
          + (w1[ic21r] * w2[ic12i] + w1[ic21i] * w2[ic12r]) * w3[ic33r];

        cr -=
          (w1[ic11r] * w2[ic32r] - w1[ic11i] * w2[ic32i]) * w3[ic23r]
          - (w1[ic11r] * w2[ic32i] + w1[ic11i] * w2[ic32r]) * w3[ic23i];

        ci -=
          (w1[ic11r] * w2[ic32r] - w1[ic11i] * w2[ic32i]) * w3[ic23i]
          + (w1[ic11r] * w2[ic32i] + w1[ic11i] * w2[ic32r]) * w3[ic23r];
      }
    }
  }

  corr = cmplx(cr, ci);
#endif  // (USE_GROUP_SU3)
}


//====================================================================
void Corr2pt_Staggered_Extended::set_asymtensor_Nc3()
{
  assert(CommonParameters::Nc() == 3);

  int Nc = 3;
  int n;

  m_epsilon_index.resize(Nc * 6);

  n = 0;
  m_epsilon_index[Nc * n]     = 0;
  m_epsilon_index[1 + Nc * n] = 1;
  m_epsilon_index[2 + Nc * n] = 2;

  n = 1;
  m_epsilon_index[Nc * n]     = 1;
  m_epsilon_index[1 + Nc * n] = 2;
  m_epsilon_index[2 + Nc * n] = 0;

  n = 2;
  m_epsilon_index[Nc * n]     = 2;
  m_epsilon_index[1 + Nc * n] = 0;
  m_epsilon_index[2 + Nc * n] = 1;

  n = 3;
  m_epsilon_index[Nc * n]     = 2;
  m_epsilon_index[1 + Nc * n] = 1;
  m_epsilon_index[2 + Nc * n] = 0;

  n = 4;
  m_epsilon_index[Nc * n]     = 1;
  m_epsilon_index[1 + Nc * n] = 0;
  m_epsilon_index[2 + Nc * n] = 2;

  n = 5;
  m_epsilon_index[Nc * n]     = 0;
  m_epsilon_index[1 + Nc * n] = 2;
  m_epsilon_index[2 + Nc * n] = 1;
}


//============================================================END=====
