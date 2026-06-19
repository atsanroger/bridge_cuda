/*!
        @file    fft_3d_local.cpp

        @brief

        @author  $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifdef USE_FFTWLIB

#include "fft_3d_local.h"
#include "Field/index_lex.h"
#include <cstring>

#ifdef USE_OPENMP
#include "ResourceManager/threadManager.h"
#endif

const std::string FFT_3d_local::class_name = "FFT_3d_local";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = FFT_3d_local::register_factory();
}
#endif

//====================================================================
FFT_3d_local::FFT_3d_local()
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
FFT_3d_local::FFT_3d_local(const Parameters& params)
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
FFT_3d_local::~FFT_3d_local()
{
  finalize();
}


//====================================================================
void FFT_3d_local::set_parameters(const Parameters& params)
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
void FFT_3d_local::get_parameters(Parameters& params) const
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
void FFT_3d_local::set_parameters(const std::string& direction)
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
bool FFT_3d_local::check_ok()
{
  int npe_xyz = Communicator::npe(0) * Communicator::npe(1) * Communicator::npe(2);

  if (npe_xyz > 1) {
    vout.crucial(m_vl, "%s: incompatible with xyz parallelization.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return true;
}


//====================================================================
void FFT_3d_local::initialize()
{
#ifdef USE_OPENMP
  int thread_ok = fftw_init_threads();

  if (thread_ok) {
    fftw_plan_with_nthreads(ThreadManager::get_num_threads());
  }
#endif
}


//====================================================================
void FFT_3d_local::initialize_plan(const Field& src)
{
  if ((m_nv == src.nin() / 2) && (m_vol == CommonParameters::Lvol() / CommonParameters::Lt())) {
    vout.detailed(m_vl, "%s: plan recycled.\n", class_name.c_str());
    return;
  } else {
    vout.detailed(m_vl, "%s: create plan.\n", class_name.c_str());
  }

  // first, clear pre-existing plan, if any.
  clear_plan();

  // assume 3dimensional
  m_ndim = 3;
  // row-major fortran order
  m_nsize[0] = CommonParameters::Lz();
  m_nsize[1] = CommonParameters::Ly();
  m_nsize[2] = CommonParameters::Lx();

  // local volume
  m_vol = 1;
  for (int i = 0; i < m_ndim; ++i) {
    m_vol *= m_nsize[i];
  }

  // number of complex elements. run nin at a time.
  m_nv = src.nin() / 2;

  m_buf_in  = fftw_alloc_complex(m_nv * m_vol);
  m_buf_out = fftw_alloc_complex(m_nv * m_vol);

  if ((!m_buf_in) || (!m_buf_out)) {
    vout.crucial(m_vl, "%s: memory allocation failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_plan_fw = fftw_plan_many_dft(m_ndim, m_nsize, m_nv,
                                 m_buf_in, m_nsize, m_nv, 1,
                                 m_buf_out, m_nsize, m_nv, 1,
                                 FFTW_FORWARD, FFTW_ESTIMATE);

  m_plan_bw = fftw_plan_many_dft(m_ndim, m_nsize, m_nv,
                                 m_buf_in, m_nsize, m_nv, 1,
                                 m_buf_out, m_nsize, m_nv, 1,
                                 FFTW_BACKWARD, FFTW_ESTIMATE);

  if ((!m_plan_fw) || (!m_plan_bw)) {
    vout.crucial(m_vl, "%s: create plan failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void FFT_3d_local::clear_plan()
{
  if (m_buf_in) fftw_free(m_buf_in);
  if (m_buf_out) fftw_free(m_buf_out);

  if (m_plan_fw) fftw_destroy_plan(m_plan_fw);
  if (m_plan_bw) fftw_destroy_plan(m_plan_bw);
}


//====================================================================
void FFT_3d_local::finalize()
{
  clear_plan();
}


//====================================================================
void FFT_3d_local::fft(Field& dst, const Field& src, enum Direction dir)
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
    scal(dst, 1.0 / m_vol);
  }
}


//====================================================================
void FFT_3d_local::fft(Field& dst, const Field& src)
{
  return fft(dst, src, m_direction);
}


//====================================================================
void FFT_3d_local::fft(Field& field)
{
  // return fft(field, field, m_direction);
  vout.crucial(m_vl, "Error at %s: fft on-the-fly unsupported.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
#endif /* USE_FFTWLIB */

//====================================================================
//============================================================END=====
