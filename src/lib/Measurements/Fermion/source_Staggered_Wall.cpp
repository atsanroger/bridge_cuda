/*!
        @file    source_Staggered_Wall.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "source_Staggered_Wall.h"

const std::string Source_Staggered_Wall::class_name = "Source_Staggered_Wall";

//====================================================================
void Source_Staggered_Wall::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int source_position;

  int err = 0;
  err += params.fetch_int("source_position", source_position);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(source_position);
}


//====================================================================
void Source_Staggered_Wall::get_parameters(Parameters& params) const
{
  params.set_int("source_position", m_t_src);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Source_Staggered_Wall::set_parameters(const int source_position)
{
  const int Lt = CommonParameters::Lt();

  //- print input parameters
  vout.general(m_vl, "Source for staggered spinor field - Wall:\n");
  vout.general(m_vl, "  source_position[t] = %d\n", source_position);

  //- range check
  int err = 0;
  // NB. Lt > abs(source_position)
  err += ParameterCheck::non_negative(Lt - abs(source_position));

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_t_src = source_position;

  //- post-process
  init();
}


//====================================================================
void Source_Staggered_Wall::init()
{
  const int Nvol = CommonParameters::Nvol();

  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  const int ipe_x = Communicator::ipe(0);
  const int ipe_y = Communicator::ipe(1);
  const int ipe_z = Communicator::ipe(2);
  const int ipe_t = Communicator::ipe(3);

  src_wall_1.reset(1, Nvol, 1);
  src_wall_2.reset(1, Nvol, 1);
  src_wall_1.set(0.0);
  src_wall_2.set(0.0);

  const int tl_src = m_t_src % Nt;
  const int ex     = 0;

  if (m_t_src / Nt == ipe_t) {
    for (int z = 0; z < Nz; ++z) {
      int z_global = z + ipe_z * Nz;

      for (int y = 0; y < Ny; ++y) {
        int y_global = y + ipe_y * Ny;

        for (int x = 0; x < Nx; ++x) {
          int x_global = x + ipe_x * Nx;

          int site = m_index.site(x, y, z, tl_src);

          if (((x_global + y_global + z_global) % 2) == 0) {
            src_wall_1.set(0, site, ex, 1.0);
            src_wall_2.set(0, site, ex, 1.0);
          } else {
            src_wall_1.set(0, site, ex, 1.0);
            src_wall_2.set(0, site, ex, -1.0);
          }
        }
      }
    }
  }
}


//====================================================================
void Source_Staggered_Wall::set(Field_F_1spinor& src,
                                const int ic, const int i_src)
{
  const int Nvol = CommonParameters::Nvol();

  src.set(0.0);

  const int ex = 0;

  if (i_src == 0) {
    for (int site = 0; site < Nvol; ++site) {
      src.set_r(ic, site, ex, src_wall_1.cmp(0, site, ex));
    }
  } else if (i_src == 1) {
    for (int site = 0; site < Nvol; ++site) {
      src.set_r(ic, site, ex, src_wall_2.cmp(0, site, ex));
    }
  } else {
    vout.crucial(m_vl, "Error at %s: illegal source index\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
//============================================================END=====
