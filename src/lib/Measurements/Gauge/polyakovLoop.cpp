/*!
        @file    polyakovLoop.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "polyakovLoop.h"
#include "ResourceManager/threadManager.h"
#include "Field/field_thread-inc.h"

const std::string PolyakovLoop::class_name = "PolyakovLoop";

//====================================================================
void PolyakovLoop::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

#if 0
  //- fetch and check input parameters
  int Nspc_size, Ntype;

  int err = 0;
  err += params.fetch_int("spatial_correlator_size", Nspc_size);
  err += params.fetch_int("number_of_correlator_type", Ntype);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nspc_size, Ntype);
#endif
}


//====================================================================
void PolyakovLoop::get_parameters(Parameters& params) const
{
  params.set_int("spatial_correlator_size", m_Nspc_size);
  params.set_int("number_of_correlator_type", m_Ntype);

  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void PolyakovLoop::set_parameters(const int Nspc_size, const int Ntype)
{
#if 0
  //- print input parameters
  vout.general(m_vl, "Polyakov loop measurement:\n");
  vout.general(m_vl, "  Nspc_size = %d\n", Nspc_size);
  vout.general(m_vl, "  Ntype     = %d\n", Ntype);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nspc_size);
  err += ParameterCheck::non_negative(Ntype);

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  if (Ntype > 6) ++err;

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nspc_size = Nspc_size;
  m_Ntype     = Ntype;


  //- post-process
  m_Nx_ext   = CommonParameters::Nx() + m_Nspc_size + 1;
  m_Ny_ext   = CommonParameters::Ny() + m_Nspc_size + 1;
  m_Nz_ext   = CommonParameters::Nz() + m_Nspc_size + 1;
  m_Nt_ext   = 1;
  m_Nvol_ext = m_Nx_ext * m_Ny_ext * m_Nz_ext * m_Nt_ext;

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  m_Nmax[0] = Nspc_size;
  m_Nmax[1] = Nspc_size;
  m_Nmax[2] = Nspc_size / 2;
  m_Nmax[3] = Nspc_size;
  m_Nmax[4] = Nspc_size / 2;
  m_Nmax[5] = Nspc_size / 2;
#endif
}


//====================================================================
void PolyakovLoop::init()
{
  const int Ndim = CommonParameters::Ndim();

  assert(Ndim == 4);

  m_filename_output = "stdout";

  const int Nx       = CommonParameters::Nx();
  const int Ny       = CommonParameters::Ny();
  const int Nz       = CommonParameters::Nz();
  const int Nvol_spc = Nx * Ny * Nz;

  m_P.reset(Nvol_spc, 1);
  m_Pcp1.reset(Nvol_spc, 1);
  m_Pcp2.reset(Nvol_spc, 1);

#if 0
  m_Ntype_max = 6;
  const int Ndim_spc = Ndim - 1;

  m_Nunit.resize(m_Ntype_max);
  m_Nmax.resize(m_Ntype_max);

  for (int i = 0; i < m_Ntype_max; ++i) {
    m_Nunit[i].resize(Ndim_spc);
  }

  // The following setting explicitly depends on the definition
  // of unit vectors.
  assert(m_Ntype_max >= 6);

  m_Nunit[0][0] = 1;
  m_Nunit[0][1] = 0;
  m_Nunit[0][2] = 0;

  m_Nunit[1][0] = 1;
  m_Nunit[1][1] = 1;
  m_Nunit[1][2] = 0;

  m_Nunit[2][0] = 2;
  m_Nunit[2][1] = 1;
  m_Nunit[2][2] = 0;

  m_Nunit[3][0] = 1;
  m_Nunit[3][1] = 1;
  m_Nunit[3][2] = 1;

  m_Nunit[4][0] = 2;
  m_Nunit[4][1] = 1;
  m_Nunit[4][2] = 1;

  m_Nunit[5][0] = 2;
  m_Nunit[5][1] = 2;
  m_Nunit[5][2] = 1;
#endif
}


//====================================================================
dcomplex PolyakovLoop::measure_ploop(const Field_G& U)
{
  ThreadManager::assert_single_thread(class_name);

  dcomplex ploop;

#pragma omp parallel
  {
    dcomplex ploop1;
    calc_ploop(ploop1, U);

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) ploop = ploop1;
  }

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "PolyakovLoop(Re,Im) = %20.16e %20.16e\n",
               real(ploop), imag(ploop));

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return ploop;
}


//====================================================================
void PolyakovLoop::calc_ploop(dcomplex& ploop, const Field_G& U)
{
#pragma omp barrier

  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  const int Nvol = U.nvol();
  assert(Nvol == Nx * Ny * Nz * Nt);

  const int Nvol_spc = m_P.nvol();
  assert(Nvol_spc == Nx * Ny * Nz);

  copy(m_Ut, 0, U, Ndim - 1);
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol_spc);

  const Index_lex index;
  const Index_lex index_spc(Nx, Ny, Nz, 1);

  //- Definition of local Polyakov loops
  const int t = 0;

  for (int site2 = is; site2 < ns; ++site2) {
    int x = site2 % Nx;
    int y = (site2 % (Nx * Ny) - x) / Nx;
    int z = (site2 - Nx * y - x) / (Nx * Ny);

    int site = index.site(x, y, z, t);

    Mat_SU_N utmp1(Nc);
    m_Ut.mat(utmp1, site, 0);

    m_P.set_mat(site2, 0, utmp1);
  }
#pragma omp barrier

  for (int t = 1; t < Nt; ++t) {
    for (int site2 = is; site2 < ns; ++site2) {
      int x = site2 % Nx;
      int y = (site2 % (Nx * Ny) - x) / Nx;
      int z = (site2 - Nx * y - x) / (Nx * Ny);

      int site = index.site(x, y, z, t);

      Mat_SU_N utmp1(Nc);
      m_Ut.mat(utmp1, site, 0);

      Mat_SU_N utmp2(Nc);
      m_P.mat(utmp2, site2, 0);

      Mat_SU_N utmp3(Nc);
      utmp3.mult_nn(utmp2, utmp1);

      m_P.set_mat(site2, 0, utmp3);
    }
#pragma omp barrier
  }

  //- global Polyakov loops
  const int Npe_t = Communicator::npe(Ndim - 1);
  if (Npe_t > 1) {
    const int size_cp = m_P.nin() * Nvol_spc;

    for (int ipe_t = 1; ipe_t < Npe_t; ++ipe_t) {
      if (ipe_t == 1) {
        copy(m_Pcp1, 0, m_P, 0);
      } else {
        copy(m_Pcp1, 0, m_Pcp2, 0);
      }

#pragma omp barrier

      if (ith == 0) {
        Communicator::exchange(size_cp, m_Pcp2.ptr(0), m_Pcp1.ptr(0),
                               Ndim - 1, 1, 0);
      }
#pragma omp barrier

      mult_Field_Gnn(m_Pcp1, 0, m_P, 0, m_Pcp2, 0);
#pragma omp barrier

      copy(m_P, 0, m_Pcp1, 0);
#pragma omp barrier
    }
  }

  //- Take the trace
  dcomplex tr = 0.0;

  for (int site = is; site < ns; ++site) {
    Mat_SU_N utmp(Nc);
    m_P.mat(utmp, site, 0);

    tr += Tr(utmp);
  }
#pragma omp barrier

  ThreadManager::reduce_sum_global(tr, ith, nth);

  const int Lvol_spc = CommonParameters::Lvol() / CommonParameters::Lt();

  const double re_ploop = real(tr) / (Nc * Lvol_spc * Npe_t);
  const double im_ploop = imag(tr) / (Nc * Lvol_spc * Npe_t);

  ploop = cmplx(re_ploop, im_ploop);

#pragma omp barrier
}


//============================================================END=====
