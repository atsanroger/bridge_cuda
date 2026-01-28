/*!
        @file    aindex_eo-tmpl.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "aindex_eo.h"


template<typename REALTYPE>
const std::string AIndex_eo<REALTYPE, VECTOR>::class_name
                                = "AIndex_eo<REALTYPE, VECTOR>";

//====================================================================
template<typename REALTYPE>
void AIndex_eo<REALTYPE, VECTOR>::init()
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
    vout.crucial(m_vl, "AIndex_eo: Nx is not even.\n");
    exit(EXIT_FAILURE);
  }

  int ipe1 = Communicator::ipe(1);
  int ipe2 = Communicator::ipe(2);
  int ipe3 = Communicator::ipe(3);

  m_ieo_origin = (ipe1 * Ny + ipe2 * Nz + ipe3 * Nt) % 2;
  // note that Nx must be even.

  m_Nsize[0] = Nx;
  m_Nsize[1] = Ny;
  m_Nsize[2] = Nz;
  m_Nsize[3] = Nt;

  m_Nsize2[0] = Nx2;
  m_Nsize2[1] = Ny;
  m_Nsize2[2] = Nz;
  m_Nsize2[3] = Nt;

  Leo.resize(Ny * Nz * Nt);
  for(int t = 0; t < Nt; ++t) {
    int t2 = ipe3 * Nt + t;
    for(int z = 0; z < Nz; ++z) {
      int z2 = ipe2 * Nz + z;
      for(int y = 0; y < Ny; ++y) {
        int y2 = ipe1 * Ny + y;
        Leo[y + Ny * (z + Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

}

//============================================================END=====
