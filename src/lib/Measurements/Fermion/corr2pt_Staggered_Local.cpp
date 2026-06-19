/*!
        @file    corr2pt_Staggered_Local.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "lib/Measurements/Fermion/corr2pt_Staggered_Local.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Communicator/communicator.h"
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

const std::string Corr2pt_Staggered_Local::class_name
  = "Corr2pt_Staggered_Local";

//====================================================================
void Corr2pt_Staggered_Local::init()
{
  ThreadManager::assert_single_thread(class_name);

  int Nc = CommonParameters::Nc();

  if (Nc == 3) {
    set_asymtensor_Nc3();
  } else {
    vout.crucial(m_vl, "%s: Nc = %d is not implemented.\n", Nc);
    exit(EXIT_FAILURE);
  }

  // Note that global lattice sizes in x,y,z directions must be even.
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  int       k  = (Lx % 2) + (Ly % 2) + (Lz % 2);
  if (k > 0) {
    vout.crucial("%s: global spatial lattice sizes must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_filename_output = "stdout";
}


//====================================================================
void Corr2pt_Staggered_Local::meson_all(std::vector<double>& meson,
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

  std::vector<double> prty(Ndim_spc);

  vout.general(m_vl, "PS <-- wall_cube correlator:\n");
  prty[0] = 1.0;
  prty[1] = 1.0;
  prty[2] = 1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  for (int t = 0; t < corr.size(); ++t) {
    meson[t] = real(corr[t]);
  }

  vout.general(m_vl, "SC <-- wall_cube correlator:\n");
  prty[0] = -1.0;
  prty[1] = -1.0;
  prty[2] = -1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "VT1 <-- wall_cube correlator:\n");
  prty[0] = -1.0;
  prty[1] = 1.0;
  prty[2] = 1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "VT2 <-- wall_cube correlator:\n");
  prty[0] = 1.0;
  prty[1] = -1.0;
  prty[2] = 1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "VT3 <-- wall_cube correlator:\n");
  prty[0] = 1.0;
  prty[1] = 1.0;
  prty[2] = -1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "PV1 <-- wall_cube correlator:\n");
  prty[0] = 1.0;
  prty[1] = -1.0;
  prty[2] = -1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "PV2 <-- wall_cube correlator:\n");
  prty[0] = -1.0;
  prty[1] = 1.0;
  prty[2] = -1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  vout.general(m_vl, "PV3 <-- wall_cube correlator:\n");
  prty[0] = -1.0;
  prty[1] = -1.0;
  prty[2] = 1.0;
  meson_correlator(corr, prty, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void Corr2pt_Staggered_Local::baryon_all(std::vector<double>& baryon,
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

  std::vector<double> prty(Ndim_spc);

  vout.general(m_vl, "nucleon <-- wall_cube correlator:\n");

  nucleon_correlator(corr, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
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
void Corr2pt_Staggered_Local::meson_correlator(
  std::vector<dcomplex>& meson,
  const std::vector<double>& prty,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  int ipe_x = Communicator::ipe(0);
  int ipe_y = Communicator::ipe(1);
  int ipe_z = Communicator::ipe(2);
  int ipe_t = Communicator::ipe(3);

  int Nspc = Nx * Ny * Nz;

  Field corrF(2, Nvol, 1);

  int isrc = 0;
  int ex   = 0;

  for (int it = 0; it < Nt; ++it) {
    for (int ispc = 0; ispc < Nspc; ++ispc) {
      int ix   = ispc % Nx;
      int iy   = (ispc / Nx) % Ny;
      int iz   = ispc / (Nx * Ny);
      int lx   = ix + Nx * ipe_x;
      int ly   = iy + Ny * ipe_y;
      int lz   = iz + Nz * ipe_z;
      int eta1 = lx % 2;
      int eta2 = ly % 2;
      int eta3 = lz % 2;
      int site = ispc + Nspc * it;

      double sign1  = (1 - eta1) * 1.0 + eta1 * prty[0];
      double sign2  = (1 - eta2) * 1.0 + eta2 * prty[1];
      double sign3  = (1 - eta3) * 1.0 + eta3 * prty[2];
      double parity = sign1 * sign2 * sign3;

      double corr_r = 0.0;
      double corr_i = 0.0;
      for (int c0 = 0; c0 < Nc; ++c0) {   // source color
        for (int c1 = 0; c1 < Nc; ++c1) { // sink color
          corr_r += sq1[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    + sq1[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_i(c1, site, ex);

          corr_i += sq1[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    - sq1[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_r(c1, site, ex);
        }
      }
      corr_r *= parity;
      corr_i *= parity;
      corrF.set(0, site, ex, corr_r);
      corrF.set(1, site, ex, corr_i);
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

  std::vector<dcomplex> corr_tmp(Lt);

  int ipet = Communicator::ipe(3);
  for (int t = 0; t < Lt; ++t) {
    corr_tmp[t] = cmplx(0.0, 0.0);
  }

  for (int t = 0; t < Nt; ++t) {
    int t2 = t + ipet * Nt;
    corr_tmp[t2] = corr_local[t];
  }

  for (int t = 0; t < Lt; ++t) {
    double crr  = real(corr_tmp[t]);
    double cri  = imag(corr_tmp[t]);
    double crr2 = Communicator::reduce_sum(crr);
    double cri2 = Communicator::reduce_sum(cri);
    meson[t] = cmplx(crr2, cri2);
  }
}


//====================================================================
void Corr2pt_Staggered_Local::nucleon_correlator(
  std::vector<dcomplex>& baryon,
  const std::vector<Field_F_1spinor>& sq1,
  const std::vector<Field_F_1spinor>& sq2)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();
  int Nt   = CommonParameters::Nt();
  int ipet = Communicator::ipe(3);

  int isrc = 0;

  if (baryon.size() != Lt) baryon.resize(Lt);

  std::vector<dcomplex> corr_local(Lt);

  int FactNc = 6;  // this is specific to Nc=3 case.

  for (int itime = 0; itime < Lt; ++itime) {
    corr_local[itime] = cmplx(0.0, 0.0);
  }

  for (int it = 0; it < Nt; ++it) {
    int itime = it + Nt * ipet;

    for (int icp = 0; icp < FactNc; ++icp) {
      int ic1 = epsilon_index(0, icp) + Nc * isrc;
      int ic2 = epsilon_index(1, icp) + Nc * isrc;
      int ic3 = epsilon_index(2, icp) + Nc * isrc;

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
void Corr2pt_Staggered_Local::contract_at_t(dcomplex& corr,
                                            const Field_F_1spinor& v1,
                                            const Field_F_1spinor& v2,
                                            const Field_F_1spinor& v3,
                                            const int it)
{
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
void Corr2pt_Staggered_Local::set_asymtensor_Nc3()
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
