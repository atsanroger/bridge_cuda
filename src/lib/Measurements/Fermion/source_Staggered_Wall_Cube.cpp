/*!
        @file    source_Staggered_Wall_Cube.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "lib/Measurements/Fermion/source_Staggered_Wall_Cube.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_thread-inc.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Source_Staggered_Wall_Cube::register_factory();
}
#endif

const std::string Source_Staggered_Wall_Cube::class_name
  = "Source_Staggered_Wall_Cube";

//====================================================================
void Source_Staggered_Wall_Cube::init()
{
  ThreadManager::assert_single_thread(class_name);

  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();

  // Note that global lattice sizes in x,y,z directions must be even.
  int k = (Lx % 2) + (Ly % 2) + (Lz % 2);
  if (k > 0) {
    vout.crucial("%s: global spatial lattice sizes must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Source_Staggered_Wall_Cube::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  int source_position;

  int err = 0;
  err += params.fetch_int("source_position", source_position);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(source_position);
}


//====================================================================
void Source_Staggered_Wall_Cube::set_parameters(const int source_position)
{
  const int Lt = CommonParameters::Lt();

  int err = 0;
  // 0 <= source_position < Lt must be satisfied.
  err += ParameterCheck::non_negative(source_position + 1);
  err += ParameterCheck::non_negative(Lt - source_position);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_t_src = source_position;

  vout.general(m_vl, "%s: Source for staggered fermion - Wall_Cube:\n",
               class_name.c_str());
  vout.general(m_vl, "  source_position[t] = %d\n", m_t_src);
}


//====================================================================
void Source_Staggered_Wall_Cube::get_parameters(Parameters& params) const
{
  params.set_int("source_position", m_t_src);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Source_Staggered_Wall_Cube::set(Field& src, const int idx)
{
  int Nc   = CommonParameters::Nc();
  int ic   = idx % Nc;
  int isrc = idx / Nc;

  set(src, ic, isrc);
}


//====================================================================
void Source_Staggered_Wall_Cube::set(Field_F_1spinor& src,
                                     const int ic, const int isrc)
{
  ThreadManager::assert_single_thread(class_name);

  int   Nin  = src.nin();
  int   Nvol = src.nvol();
  int   Nex  = src.nex();
  Field src_tmp(Nin, Nvol, Nex);

  set(src_tmp, ic, isrc);
  copy(src, src_tmp);
}


//====================================================================
void Source_Staggered_Wall_Cube::set(Field& src,
                                     const int ic, const int isrc)
{
#pragma omp barrier

  if ((isrc < 0) || (isrc > 7)) {
    vout.crucial(m_vl, "Error at %s: illegal source index\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  int ipe_x = Communicator::ipe(0);
  int ipe_y = Communicator::ipe(1);
  int ipe_z = Communicator::ipe(2);
  int ipe_t = Communicator::ipe(3);

  int ex = 0;

  // eta is position in each spatial 2x2x2 cube.
  int eta_x = isrc % 2;
  int eta_y = (isrc / 2) % 2;
  int eta_z = isrc / 4;

  int it_src   = m_t_src % Nt;
  int ipet_src = m_t_src / Nt;

  src.set(0.0);
#pragma omp barrier

  int Nspc = Nx * Ny * Nz;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nspc);

  if (ipet_src == ipe_t) {
    for (int ispc = is; ispc < ns; ++ispc) {
      int ix = ispc % Nx;
      int iy = (ispc / Nx) % Ny;
      int iz = ispc / (Nx * Ny);
      int lx = ix + Nx * ipe_x;
      int ly = iy + Ny * ipe_y;
      int lz = iz + Nz * ipe_z;
      int kx = lx % 2;
      int ky = ly % 2;
      int kz = lz % 2;
      if ((kx == eta_x) && (ky == eta_y) && (kz == eta_z)) {
        int site = m_index.site(ix, iy, iz, it_src);
        src.set(2 * ic, site, ex, 1.0);
      }
    }
  }

#pragma omp barrier
}


//============================================================END=====
