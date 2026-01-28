/*!
        @file    fft_3d_parallel3d.h

        @brief

        @author  $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifndef FFT_3D_PARALLEL3D_INCLUDED
#define FFT_3D_PARALLEL3D_INCLUDED

// requires FFTW library
#ifdef USE_FFTWLIB

// requires MPI
#ifdef USE_MPI

#include <fftw3.h>

#include "fft.h"
#include "Parameters/parameters.h"
#include "Communicator/MPI/communicator_mpi.h"

class FFT_3d_parallel3d : public FFT
{
 public:
  static const std::string class_name;

 public:
  FFT_3d_parallel3d();
  FFT_3d_parallel3d(const Parameters& params);
  virtual ~FFT_3d_parallel3d();

  void fft(Field& dst, const Field& src, enum Direction dir);
  void fft(Field& dst, const Field& src);
  void fft(Field& field);

  void set_parameters(const Parameters& params);
  void set_parameters(const std::string& direction);

  void get_parameters(Parameters& params) const;

 private:

  bool check_ok();

  void initialize();
  void finalize();

  bool need_create_plan(const Field&);

  void create_plan(int site_dof);
  void release_plan();

  void create_mpi_datatype(int site_dof);
  void release_mpi_datatype();

  void create_fft_plan(int site_dof);
  void release_fft_plan();

  // utility functions
  std::vector<int> grid_rank_to_coord(int r);
  int find_global_index(const std::vector<int>& coord);

  Bridge::VerboseLevel m_vl;

  int m_site_dof;

  // {x,y,z} subcommunicator
  MPI_Comm m_comm;

  int m_local_rank;
  int m_local_ipe_x;
  int m_local_ipe_y;
  int m_local_ipe_z;

  // geometry and grid info of xyz-plane
  int m_ndims;

  std::vector<int> m_grid_size;
  int m_grid_vol;

  std::vector<int> m_lattice_size;
  int m_lattice_vol;

  std::vector<int> m_local_size;
  int m_local_vol;

  // datatypes for gather/scatter
  MPI_Datatype m_site_vector_type;
  MPI_Datatype m_local_patch_type;
  MPI_Datatype m_subarray_type;

  std::vector<int> m_sendcounts;
  std::vector<int> m_subarray_displs;
  std::vector<int> m_local_patch_displs;

  // fft workspace. do fft on-the-fly.
  fftw_complex *m_buf;

  fftw_plan m_plan_fw;
  fftw_plan m_plan_bw;

  bool m_initialized;

  // for compatibility with FFT_xyz classes
  Direction m_direction;


#ifdef USE_FACTORY
 private:
  static FFT *create_object()
  {
    return new FFT_3d_parallel3d();
  }

  static FFT *create_object_with_params(const Parameters& params)
  {
    return new FFT_3d_parallel3d(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= FFT::Factory::Register("FFT_3d_parallel_3dim", create_object);
    init &= FFT::Factory_params::Register("FFT_3d_parallel_3dim", create_object_with_params);
    return init;
  }
#endif
};

#endif /* USE_MPI */
#endif /* USE_FFTWLIB */

#endif /* FFT_3D_PARALLEL3D_INCLUDED */
