/*!
        @file    fft_3d_parallel1d.cpp

        @brief

        @author  $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifdef USE_FFTWLIB
#ifdef USE_MPI

#include "fft_3d_parallel1d.h"
#include "Field/index_lex.h"
#include "Communicator/communicator.h"
#include "Communicator/MPI/communicator_mpi.h"

#ifdef USE_OPENMP
#include "ResourceManager/threadManager.h"
#endif


const std::string FFT_3d_parallel1d::class_name = "FFT_3d_parallel1d";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = FFT_3d_parallel1d::register_factory();
}
#endif

//====================================================================
FFT_3d_parallel1d::FFT_3d_parallel1d()
  : m_vl(CommonParameters::Vlevel())
  , m_ndim(0)
  , m_vol(0)
  , m_nv(0)
  , m_buf_in(NULL)
  , m_buf_out(NULL)
  , m_plan_fw(NULL)
  , m_plan_bw(NULL)
  , m_direction(UNDEF)
{
  if (check_ok()) {
    initialize();
  }
}


//====================================================================
FFT_3d_parallel1d::FFT_3d_parallel1d(const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
  , m_ndim(0)
  , m_vol(0)
  , m_nv(0)
  , m_buf_in(NULL)
  , m_buf_out(NULL)
  , m_plan_fw(NULL)
  , m_plan_bw(NULL)
  , m_direction(UNDEF)
{
  if (check_ok()) {
    initialize();
  }
  set_parameters(params);
}


//====================================================================
FFT_3d_parallel1d::~FFT_3d_parallel1d()
{
  finalize();
}


//====================================================================
void FFT_3d_parallel1d::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  std::string direction;
  if (!params.fetch_string("FFT_direction", direction)) {
    set_parameters(direction);
  }
}


//====================================================================
void FFT_3d_parallel1d::get_parameters(Parameters& params) const
{
  if (m_direction == FORWARD) {
    params.set_string("FFT_direction", "Forward");
  } else if (m_direction == BACKWARD) {
    params.set_string("FFT_direction", "Backward");
  } else {
    params.set_string("FFT_direction", "None");
  }

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void FFT_3d_parallel1d::set_parameters(const std::string& direction)
{
  if (direction == "Forward") {
    m_direction = FORWARD;
  } else if (direction == "Backward") {
    m_direction = BACKWARD;
  } else {
    m_direction = UNDEF;

    vout.crucial(m_vl, "Error at %s: unsupported FFT direction \"%s\"\n", class_name.c_str(), direction.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
bool FFT_3d_parallel1d::check_ok()
{
  int npe_xy = Communicator::npe(0) * Communicator::npe(1);

  if (npe_xy > 1) {
    vout.crucial(m_vl, "%s: incompatible with xy parallelization.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return true;
}


//====================================================================
void FFT_3d_parallel1d::initialize()
{
#ifdef USE_OPENMP
  int thread_ok = fftw_init_threads();

  if (thread_ok) {
    fftw_plan_with_nthreads(ThreadManager::get_num_threads());
  }
#endif

  fftw_mpi_init();

  // split communicator along t-axis to form 3dim subspaces
  int ipe_x = Communicator::ipe(0);
  int ipe_y = Communicator::ipe(1);
  int ipe_z = Communicator::ipe(2);
  int ipe_t = Communicator::ipe(3);

  int npe_x = Communicator::npe(0);
  int npe_y = Communicator::npe(1);
  int npe_z = Communicator::npe(2);
  int npe_t = Communicator::npe(3);

  int local_rank = ipe_x + npe_x * (ipe_y + npe_y * ipe_z);

  MPI_Comm_split(Communicator_impl::world(), ipe_t, local_rank, &m_comm);
}


//====================================================================
void FFT_3d_parallel1d::initialize_plan(const Field& src)
{
  int local_vol =
    CommonParameters::Lz() * CommonParameters::Ly() * CommonParameters::Lx() / Communicator::npe(2);

  if ((m_nv == src.nin() / 2) && (m_vol == local_vol)) {
    vout.general(m_vl, "%s: plan recycled.\n", class_name.c_str());
    return;
  } else {
    vout.general(m_vl, "%s: create plan.\n", class_name.c_str());
  }

  // first, clear pre-existing plan, if any.
  clear_plan();

  // assume 3dimensional
  m_ndim = 3;
  // row-major fortran order
  m_nsize[0] = CommonParameters::Lz();
  m_nsize[1] = CommonParameters::Ly();
  m_nsize[2] = CommonParameters::Lx();

  // local volume size
  m_vol = local_vol;

  // number of complex elements. run nin at a time.
  m_nv = src.nin() / 2;

  // local size
  ptrdiff_t local_n0, local_0_start;

  ptrdiff_t psize = fftw_mpi_local_size_many(m_ndim, m_nsize, m_nv,
                                             FFTW_MPI_DEFAULT_BLOCK, m_comm,
                                             &local_n0, &local_0_start);

  int nz    = CommonParameters::Nz();
  int ipe_z = Communicator::ipe(2);

  if ((local_n0 != nz) || (local_0_start != nz * ipe_z)) {
    vout.crucial(m_vl, "%s: data distribution plan not matched.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_buf_in  = fftw_alloc_complex(psize);
  m_buf_out = fftw_alloc_complex(psize);

  if ((!m_buf_in) || (!m_buf_out)) {
    vout.crucial(m_vl, "%s: memory allocation failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_plan_fw = fftw_mpi_plan_many_dft(m_ndim, m_nsize, m_nv,
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                     m_buf_in, m_buf_out,
                                     m_comm,
                                     FFTW_FORWARD,
                                     FFTW_ESTIMATE);

  m_plan_bw = fftw_mpi_plan_many_dft(m_ndim, m_nsize, m_nv,
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                     m_buf_in, m_buf_out,
                                     m_comm,
                                     FFTW_BACKWARD,
                                     FFTW_ESTIMATE);

  if ((!m_plan_fw) || (!m_plan_bw)) {
    vout.crucial(m_vl, "%s: create plan failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void FFT_3d_parallel1d::clear_plan()
{
  if (m_buf_in) fftw_free(m_buf_in);
  if (m_buf_out) fftw_free(m_buf_out);

  if (m_plan_fw) fftw_destroy_plan(m_plan_fw);
  if (m_plan_bw) fftw_destroy_plan(m_plan_bw);
}


//====================================================================
void FFT_3d_parallel1d::finalize()
{
  clear_plan();

  MPI_Comm_free(&m_comm);
}


//====================================================================
void FFT_3d_parallel1d::fft(Field& dst, const Field& src, enum Direction dir)
{
  if (not ((dir == FORWARD) || (dir == BACKWARD))) {
    vout.crucial(m_vl, "%s: unsupported direction. %d\n", class_name.c_str(), dir);
    exit(EXIT_FAILURE);
  }

  initialize_plan(src);

  int nex = src.nex();
  int nt  = CommonParameters::Nt();

  Index_lex index;

  size_t count = m_nv * m_vol;  // count in complex numbers

  for (int iex = 0; iex < nex; ++iex) {
    for (int it = 0; it < nt; ++it) {
      memcpy(m_buf_in, src.ptr(0, index.site(0, 0, 0, it), iex), sizeof(double) * count * 2);

      if (dir == FORWARD) {
        fftw_execute(m_plan_fw);
      } else if (dir == BACKWARD) {
        fftw_execute(m_plan_bw);
      } else {
        vout.crucial(m_vl, "%s: unsupported direction. %d\n", class_name.c_str(), dir);
        exit(EXIT_FAILURE);
      }

      memcpy(dst.ptr(0, index.site(0, 0, 0, it), iex), m_buf_out, sizeof(double) * count * 2);
    }
  }

  if (dir == BACKWARD) {
    size_t global_vol = CommonParameters::Lvol() / CommonParameters::Lt();
    scal(dst, 1.0 / global_vol);
  }
}


//====================================================================
void FFT_3d_parallel1d::fft(Field& dst, const Field& src)
{
  return fft(dst, src, m_direction);
}


//====================================================================
void FFT_3d_parallel1d::fft(Field& field)
{
  vout.crucial(m_vl, "Error at %s: fft on-the-fly unsupported.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
#endif /* USE_MPI */
#endif /* USE_FFTWLIB */
