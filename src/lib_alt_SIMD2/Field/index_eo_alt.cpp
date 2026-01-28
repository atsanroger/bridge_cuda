/*!
        @file    $Id:: index_eo_alt.cpp #$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$

        @version $LastChangedRevision: 2160 $
*/


#include "lib_alt_SIMD2/Field/index_eo_alt.h"

#include <assert.h>

#include "lib_alt_SIMD2/Field/afield.h"


//====================================================================
template<typename REALTYPE>
void Index_eo_alt<REALTYPE,SIMD2>::init()
{
  Nx = CommonParameters::Nx();
  Ny = CommonParameters::Ny();
  Nz = CommonParameters::Nz();
  Nt = CommonParameters::Nt();
  Nvol = CommonParameters::Nvol();
  Nx2 = Nx/2;
  Nvol2 = Nvol/2;
  m_vl = CommonParameters::Vlevel();


  Nc = CommonParameters::Nc();
  Nd = CommonParameters::Nd();
  Ndf  = 2 * Nc * Nc;
  Nvcd = 2 * Nc * Nd;

  if((Nx % 2) == 1){
    vout.crucial(m_vl, "Index_eo_alt: Nx is not even.\n");
    exit(EXIT_FAILURE);
  }

  Leo.resize(Ny * Nz * Nt);
  for(int t = 0; t < Nt; ++t) {
    int t2 = Communicator::ipe(3) * Nt + t;
    for(int z = 0; z < Nz; ++z) {
      int z2 = Communicator::ipe(2) * Nz + z;
      for(int y = 0; y < Ny; ++y) {
        int y2 = Communicator::ipe(1) * Ny + y;
        Leo[y + Ny * (z + Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

}

//====================================================================
// explicit instanciation.
template class Index_eo_alt<float,SIMD2>;
template class Index_eo_alt<double,SIMD2>;
//============================================================END=====
