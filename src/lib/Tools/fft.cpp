/*!
        @file    fft.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "fft.h"

#ifdef USE_FFTWLIB


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "fft_xyz_1dim.h"
#include "fft_xyz_3dim.h"

#include "fft_3d_local.h"
#ifdef USE_MPI
#include "fft_3d_parallel1d.h"
#include "fft_3d_parallel3d.h"
#endif
#include "fft_3d.h"

bool FFT::init_factory()
{
  bool result = true;

  result &= FFT_xyz_1dim::register_factory();
  result &= FFT_xyz_3dim::register_factory();

  result &= FFT_3d_local::register_factory();
#ifdef USE_MPI
  result &= FFT_3d_parallel1d::register_factory();
  result &= FFT_3d_parallel3d::register_factory();
#endif
  result &= FFT_3d::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


#endif /* USE_FFTWLIB */
