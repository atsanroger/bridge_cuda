/*!
        @file    qws_lib.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/inline/define_params.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"

#ifdef USE_QWSLIB
#include "qws_lib.h"

using Bridge::vout;

namespace QWS_lib {
  // static int qws_is_setup = 0;
  int qws_is_setup        = 0;
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  void init_qws(const int *boundary)
  {
    if (qws_is_setup > 0) {
      ++qws_is_setup;
      return;
    }

    int Nx = CommonParameters::Nx();
    int Ny = CommonParameters::Ny();
    int Nz = CommonParameters::Nz();
    int Nt = CommonParameters::Nt();

    int npe_f[4];
    npe_f[0] = CommonParameters::NPEx();
    npe_f[1] = CommonParameters::NPEy();
    npe_f[2] = CommonParameters::NPEz();
    npe_f[3] = CommonParameters::NPEt();
    int fbc_f[4];
    fbc_f[0] = boundary[0];
    fbc_f[1] = boundary[1];
    fbc_f[2] = boundary[2];
    fbc_f[3] = boundary[3];

    int pce_f = 0;
    int pco_f = 1 - pce_f;

    int block_size[4];  // now not used, dummy.
    block_size[0] = 2;
    block_size[1] = 1;
    block_size[2] = 1;
    block_size[3] = 1;

    vout.general(vl, "qws_is_setup = %d\n", qws_is_setup);
    vout.general(vl, "   boundary conditation: %d %d %d %d\n", fbc_f[0], fbc_f[1], fbc_f[2], fbc_f[3]);
    if (qws_is_setup == 0) {
      qws_init_(&Nx, &Ny, &Nz, &Nt, npe_f, fbc_f, &pce_f, &pco_f,
                block_size);
      ++qws_is_setup;
      vout.general(vl, "qws_init_ was called: new qws_is_setup=%d\n", qws_is_setup);
    }
  }


  void tidyup_qws()
  {
    if (qws_is_setup > 0) {
      qws_finalize_();
    }
    return;

    if (qws_is_setup == 0) {
      return;
    }
    --qws_is_setup;
    if (qws_is_setup == 0) {
      qws_finalize_();
    }
  }
}

extern "C" {
/*
  extern int nt, nz, ny, nx;

  int std_xyzt2i_(int* j){
  return j[0] + nx*(j[1] + ny*(j[2] + nz*j[3]));
  }

  int e_o_(int* j){
  return (j[0] + j[1] + j[2] + j[3]) % 2;
  }
*/
}

#endif
