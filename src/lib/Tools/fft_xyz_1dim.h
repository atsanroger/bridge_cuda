/*!
        @file    fft_xyz_1dim.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifndef FFT_XYZ_1DIM_INCLUDED
#define FFT_XYZ_1DIM_INCLUDED

#ifdef USE_FFTWLIB

#include "fft.h"

//! Fast Fourier Transformation in x,y,z directions.

/*!
   This class implements Fast Fourier Transformation
   in x,y,z directions 1 dim by 1 dim.
   NB. FFTW supports parallelization only in 1 direction.
                                   [06 Jun 2015 Y.Namekawa]
 */

class FFT_xyz_1dim : public FFT
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  bool m_is_forward;

  Index_lex m_index;

  fftw_complex *m_in;
  fftw_complex *m_out;

  fftw_plan m_plan;

#ifdef USE_MPI
  ptrdiff_t m_Nsize_in_p, m_start_in_p;
  ptrdiff_t m_Nsize_out_p, m_start_out_p;
#endif


 public:
  FFT_xyz_1dim()
    : m_vl(CommonParameters::Vlevel()) {}

  FFT_xyz_1dim(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  ~FFT_xyz_1dim() {}

  void fft(Field& field);  // field is overwritten
  void fft(Field& field_out, const Field& field_in);
  void fft(Field& field_out, const Field& field_in, const Direction dir);

  void set_parameters(const Parameters&);
  void set_parameters(const std::string& str_fft_direction);

  void get_parameters(Parameters&) const;

#ifdef USE_FACTORY
 private:
  static FFT *create_object()
  {
    return new FFT_xyz_1dim();
  }

  static FFT *create_object_with_params(const Parameters& params)
  {
    return new FFT_xyz_1dim(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= FFT::Factory::Register("FFT_xyz_1dim", create_object);
    init &= FFT::Factory_params::Register("FFT_xyz_1dim", create_object_with_params);
    return init;
  }
#endif
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
