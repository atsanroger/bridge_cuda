/*!
      @file    aindex_eo-tmpl.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


//====================================================================
template<typename REALTYPE>
void AIndex_eo<REALTYPE, QXS>::init()
{
  Nx    = CommonParameters::Nx();
  Ny    = CommonParameters::Ny();
  Nz    = CommonParameters::Nz();
  Nt    = CommonParameters::Nt();
  Nvol  = CommonParameters::Nvol();
  Nx2   = Nx / 2;
  Nvol2 = Nvol / 2;
  m_vl  = CommonParameters::Vlevel();

  Nc   = CommonParameters::Nc();
  Nd   = CommonParameters::Nd();
  Ndf  = 2 * Nc * Nc;
  Nvcd = 2 * Nc * Nd;

  if ((Nx % 2) != 0) {
    vout.crucial(m_vl, "Index_eo_alt: Nx must be even.\n");
    exit(EXIT_FAILURE);
  }

  if ((Nx2 % VLENX) != 0) {
    vout.crucial(m_vl, "Index_eo_alt: Nx2 % VLENX must be 0.\n");
    exit(EXIT_FAILURE);
  }

  if ((Ny % VLENY) != 0) {
    vout.crucial(m_vl, "Index_eo_alt: Ny % VLENY must be 0.\n");
    exit(EXIT_FAILURE);
  }

  Leo.resize(Ny * Nz * Nt);
  for (int it = 0; it < Nt; ++it) {
    int it2 = Communicator::ipe(3) * Nt + it;
    for (int iz = 0; iz < Nz; ++iz) {
      int iz2 = Communicator::ipe(2) * Nz + iz;
      for (int iy = 0; iy < Ny; ++iy) {
        int iy2 = Communicator::ipe(1) * Ny + iy;
        Leo[iy + Ny * (iz + Nz * it)] = (iy2 + iz2 + it2) % 2;
      }
    }
  }
}


//============================================================END=====
