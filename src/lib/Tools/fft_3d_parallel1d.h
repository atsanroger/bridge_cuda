/*!
        @file    fft_3d_parallel1d.h

        @brief

        @author  $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifndef FFT_3D_PARALLEL1D_INCLUDED
#define FFT_3D_PARALLEL1D_INCLUDED

// requires FFTW library
#ifdef USE_FFTWLIB

// requires MPI
#ifdef USE_MPI

#include <fftw3-mpi.h>

#include "fft.h"
#include "Communicator/MPI/communicator_mpi.h"

class FFT_3d_parallel1d : public FFT
{
 public:
  static const std::string class_name;

 public:
  FFT_3d_parallel1d();
  FFT_3d_parallel1d(const Parameters& params);
  virtual ~FFT_3d_parallel1d();

  void fft(Field& dst, const Field& src, enum Direction dir);
  void fft(Field& dst, const Field& src);
  void fft(Field& field);

  void set_parameters(const Parameters& params);
  void set_parameters(const std::string& direction);

  void get_parameters(Parameters& params) const;

 private:

  void initialize();
  void finalize();

  bool check_ok();

  void initialize_plan(const Field& src);
  void clear_plan();

  Bridge::VerboseLevel m_vl;

  // size info
  int m_ndim;
  ptrdiff_t m_nsize[4]; // assume atmost 4dimensional
  ptrdiff_t m_vol;
  ptrdiff_t m_nv;       // number of complex elements

  fftw_complex *m_buf_in;
  fftw_complex *m_buf_out;
  fftw_plan m_plan_fw;
  fftw_plan m_plan_bw;

  MPI_Comm m_comm;

  Direction m_direction;


#ifdef USE_FACTORY
 private:
  static FFT *create_object()
  {
    return new FFT_3d_parallel1d();
  }

  static FFT *create_object_with_params(const Parameters& params)
  {
    return new FFT_3d_parallel1d(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= FFT::Factory::Register("FFT_3d_parallel_1dim", create_object);
    init &= FFT::Factory_params::Register("FFT_3d_parallel_1dim", create_object_with_params);
    return init;
  }
#endif
};

#endif /* USE_MPI */
#endif /* USE_FFTWLIB */

#endif /* FFT_3D_PARALLEL1D_INCLUDED */
