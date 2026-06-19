/*!
      @file    shiftAField_lex_vector-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

//====================================================================
void shift_lex_xp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  real_t bc2 = bc[0];

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;
    for(int in = 0; in < nin; ++in){
      buf[in + nin * iyzt] = bc2 * v1[IDXV(nin, Nst, in, ist, 0)];
    }
  }

}

//====================================================================
void shift_lex_xp2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;
    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = buf[in + nin * iyzt];
    }
  }

}

//====================================================================
void shift_lex_xpb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ((ix+1) % Nx) + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == Nx-1) bc2 = bc[0];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_xm1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;
    for(int in = 0; in < nin; ++in){
      buf[in + nin * iyzt] = v1[IDXV(nin, Nst, in, ist, 0)];
    }
  }

}

//====================================================================
void shift_lex_xm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  real_t bc2 = bc[0];

  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;
    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * buf[in + nin * iyzt];
    }
  }

}

//====================================================================
void shift_lex_xmb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  // int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ix-1 + Nx * iyzt;
    if(ix == 0) nei = Nx-1 + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == 0) bc2 = bc[0];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_yp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  real_t bc2 = bc[1];

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;
     for(int in = 0; in < nin; ++in){
       buf[in + nin * ixzt] = bc2 * v1[IDXV(nin, Nst, in, ist, 0)];
     }
   }
  }

}

//====================================================================
void shift_lex_yp2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;

     for(int in = 0; in < nin; ++in){
       v2[IDXV(nin, Nst, in, ist, 0)] = buf[in + nin * ixzt];
     }
   }
  }

}

//====================================================================
void shift_lex_ypb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy+1) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);
    real_t bc2 = 1.0;
    if(iy == Ny-1) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }
  }

}

//====================================================================
void shift_lex_ym1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  // int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny*izt);
     int ixzt = ix + Nx * izt;

     for(int in = 0; in < nin; ++in){
       buf[in + nin * ixzt] = v1[IDXV(nin, Nst, in, ist, 0)];
     }

   }
  }

}

//====================================================================
void shift_lex_ym2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  // int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  real_t bc2 = bc[1];

  for(int izt = 0; izt < Nzt; ++izt){
   for(int ix = 0; ix < Nx; ++ix){
     int iy = 0;
     int ist  = ix + Nx*(iy + Ny*izt);
     int ixzt = ix + Nx * izt;

     for(int in = 0; in < nin; ++in){
       v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * buf[in + nin*ixzt];
     }

   }
  }

}

//====================================================================
void shift_lex_ymb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst  = Nx * Ny * Nzt;

  for(int ist = 0; ist < Nst; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy-1+Ny) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);
    real_t bc2 = 1.0;
    if(iy == 0) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_zp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nxy * Nz * Nt;

  real_t bc2 = bc[idir];

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     for(int in = 0; in < nin; ++in){
       buf[in + nin * ixyt] = bc2 * v1[IDXV(nin, Nst, in, ist, 0)];
     }

   }
  }

}

//====================================================================
void shift_lex_zp2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy = Nsize[0] * Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     for(int in = 0; in < nin; ++in){
       v2[IDXV(nin, Nst, in, ist, 0)] = buf[in + nin * ixyt];
     }

   }
  }

}

//====================================================================
void shift_lex_zpb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy = Nsize[0] * Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz+1) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == Nz-1) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_zm1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy = Nsize[0] * Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     for(int in = 0; in < nin; ++in){
       buf[in + nin * ixyt] = v1[IDXV(nin, Nst, in, ist, 0)];
     }
   }
  }

}

//====================================================================
void shift_lex_zm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  real_t bc2 = bc[idir];

  for(int it = 0; it < Nt; ++it){
   for(int ixy = 0; ixy < Nxy; ++ixy){
     int iz = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     for(int in = 0; in < nin; ++in){
       v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * buf[in + nin*ixyt];
     }

   }
  }

}

//==================================================================== 
void shift_lex_zmb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 2;

  int Nxy = Nsize[0] * Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz-1+Nz) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == 0) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_tp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nxyz * Nt;

  real_t bc2 = bc[idir];

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;
    for(int in = 0; in < nin; ++in){
      buf[in + nin * ixyz] = bc2 * v1[IDXV(nin, Nst, in, ist, 0)];
    }
  }

}

//====================================================================
void shift_lex_tp2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;
    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = buf[in + nin * ixyz];
    }
  }

}

//====================================================================
void shift_lex_tpb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == Nt-1) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }

  }

}

//====================================================================
void shift_lex_tm1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;
    for(int in = 0; in < nin; ++in){
      buf[in + nin * ixyz] = v1[IDXV(nin, Nst, in, ist, 0)];
    }
  }

}

//====================================================================
void shift_lex_tm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst  = Nxyz * Nt;

  real_t bc2 = bc[idir];

  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;
    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * buf[in + nin * ixyz];
    }
  }

}

//==================================================================== 
void shift_lex_tmb(real_t *RESTRICT v2, real_t *RESTRICT v1,
                   int nin, int *Nsize, int *bc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  for(int ist = 0; ist < Nst; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == 0) bc2 = bc[idir];

    for(int in = 0; in < nin; ++in){
      v2[IDXV(nin, Nst, in, ist, 0)] = bc2 * v1[IDXV(nin, Nst, in, nei, 0)];
    }
  }

}

//====================================================================

