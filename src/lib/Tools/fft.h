/*!
        @file    fft.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifndef FFT_INCLUDED
#define FFT_INCLUDED

#ifdef USE_FFTWLIB

#ifdef USE_MPI
#include <fftw3-mpi.h>
#include "Communicator/MPI/communicator_mpi.h"
#else
#include <fftw3.h>
#endif

#ifdef USE_OPENMP
#include "ResourceManager/threadManager.h"
#endif

#include "Field/field.h"
#include "Field/index_lex.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "factory.h"
#endif


//! Base class for FFT class family.

/*!
    This class implements a base class of Fast Fourier Transformation.
                               [7 Sep 2015 Y.Namekawa]
 */

class FFT
{
 public:
  enum Direction
  {
    FORWARD  = FFTW_FORWARD,
    BACKWARD = FFTW_BACKWARD,
    UNDEF
  };

 public:
  FFT() {}
  virtual ~FFT() {}

 private:
  // non-copyable
  FFT(const FFT&);
  FFT& operator=(const FFT&);

 public:
  virtual void fft(Field& field) = 0;  // field is overwritten
  virtual void fft(Field& field_out, const Field& field_in) = 0;
  virtual void fft(Field& field_out, const Field& field_in, const Direction dir) = 0;

  virtual void set_parameters(const std::string& str_fft_direction) {}


#ifdef USE_FACTORY
 public:
  typedef FFT *(*ProductCreator)();
  typedef FFT *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<FFT, ProductCreator>          Factory;
  typedef FactoryTemplate<FFT, ProductCreator_params>   Factory_params;

  static FFT *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static FFT *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
