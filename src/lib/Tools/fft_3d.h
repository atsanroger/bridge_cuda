/*!
        @file    fft_3d.h

        @brief

        @author  $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FFT_3D_INCLUDED
#define FFT_3D_INCLUDED

#ifdef USE_FFTWLIB

// implementations
#include "fft.h"
#include "fft_3d_local.h"
#ifdef USE_MPI
#include "fft_3d_parallel1d.h"
#include "fft_3d_parallel3d.h"
#endif

#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "factory.h"
#endif


//! Base class for FFT class family.

/*!
    This class implements a base class of Fast Fourier Transformation.
                               [7 Sep 2015 Y.Namekawa]
    Alternative implementations that allow parallelization in arbitrary
    directions.
                               [5 May 2017 T.Aoyama]
 */

class FFT_3d
{
#ifdef USE_FACTORY
 private:
  static FFT *create_object()
  {
    // auto-select

    int npe_x = Communicator::npe(0);
    int npe_y = Communicator::npe(1);
    int npe_z = Communicator::npe(2);
    int npe_t = Communicator::npe(3);

    if ((npe_x == 1) && (npe_y == 1) && (npe_z == 1)) {
      // no parallelization in xyz-directions
      return new FFT_3d_local;

#ifdef USE_MPI
    } else if ((npe_x == 1) && (npe_y == 1)) {
      // parallelization only in z-direction
      return new FFT_3d_parallel1d;
    } else {
      // most general case
      return new FFT_3d_parallel3d;
#endif
    }

    // default
    return NULL;
  }

  static FFT *create_object_with_params(const Parameters& params)
  {
    // auto-select

    int npe_x = Communicator::npe(0);
    int npe_y = Communicator::npe(1);
    int npe_z = Communicator::npe(2);
    int npe_t = Communicator::npe(3);

    if ((npe_x == 1) && (npe_y == 1) && (npe_z == 1)) {
      // no parallelization in xyz-directions
      return new FFT_3d_local(params);

#ifdef USE_MPI
    } else if ((npe_x == 1) && (npe_y == 1)) {
      // parallelization only in z-direction
      return new FFT_3d_parallel1d(params);
    } else {
      // most general case
      return new FFT_3d_parallel3d(params);
#endif
    }

    // default
    return NULL;
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= FFT::Factory::Register("auto", create_object);
    init &= FFT::Factory_params::Register("auto", create_object_with_params);
    return init;
  }
#endif
};

//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
