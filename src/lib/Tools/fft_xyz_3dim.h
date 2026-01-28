/*!
        @file    fft_xyz_3dim.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifndef FFT_XYZ_3DIM_INCLUDED
#define FFT_XYZ_3DIM_INCLUDED

#ifdef USE_FFTWLIB

#include "fft.h"

//! Fast Fourier Transformation in x,y,z directions.

/*!
   This class implements Fast Fourier Transformation
   in x,y,z directions simultaneously.
   NB. FFTW supports parallelization only in z direction.
                                   [06 Jun 2015 Y.Namekawa]
 */

class FFT_xyz_3dim : public FFT
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
  ptrdiff_t m_Nz_p, m_z_start_p;
#endif


 public:
  FFT_xyz_3dim()
    : m_vl(CommonParameters::Vlevel())
  {
    init();
  }

  FFT_xyz_3dim(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    init();
    set_parameters(params);
  }

  ~FFT_xyz_3dim()
  {
    tidy_up();
  }

  void fft(Field& field);  // field is overwritten
  void fft(Field& field_out, const Field& field_in);
  void fft(Field& field_out, const Field& field_in, const Direction dir);

  void set_parameters(const Parameters&);
  void set_parameters(const std::string& str_fft_direction);

  void get_parameters(Parameters&) const;

 private:
  void init();
  void tidy_up();

#ifdef USE_FACTORY
 private:
  static FFT *create_object()
  {
    return new FFT_xyz_3dim();
  }

  static FFT *create_object_with_params(const Parameters& params)
  {
    return new FFT_xyz_3dim(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= FFT::Factory::Register("FFT_xyz_3dim", create_object);
    init &= FFT::Factory_params::Register("FFT_xyz_3dim", create_object_with_params);
    return init;
  }
#endif
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
