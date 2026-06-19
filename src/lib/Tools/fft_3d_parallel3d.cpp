/*!
        @file    fft_3d_parallel3d.cpp

        @brief

        @author  $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#ifdef USE_FFTWLIB
#ifdef USE_MPI

#include "fft_3d_parallel3d.h"
#include "Communicator/communicator.h"
#include "Communicator/MPI/communicator_mpi.h"

#ifdef USE_OPENMP
#include "ResourceManager/threadManager.h"
#endif


const std::string FFT_3d_parallel3d::class_name = "FFT_3d_parallel3d";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = FFT_3d_parallel3d::register_factory();
}
#endif

//====================================================================
void FFT_3d_parallel3d::set_parameters(const Parameters& params)
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
void FFT_3d_parallel3d::get_parameters(Parameters& params) const
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
void FFT_3d_parallel3d::set_parameters(const std::string& direction)
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
FFT_3d_parallel3d::FFT_3d_parallel3d()
  : m_vl(CommonParameters::Vlevel())
  , m_initialized(false)
  , m_direction(UNDEF)
{
  if (check_ok()) {
    initialize();
  }
}


//====================================================================
FFT_3d_parallel3d::FFT_3d_parallel3d(const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
  , m_initialized(false)
  , m_direction(UNDEF)
{
  if (check_ok()) {
    initialize();
  }
  set_parameters(params);
}


//====================================================================
FFT_3d_parallel3d::~FFT_3d_parallel3d()
{
  finalize();
}


//====================================================================
bool FFT_3d_parallel3d::check_ok()
{
  return true;
}


//====================================================================
void FFT_3d_parallel3d::initialize()
{
#ifdef USE_OPENMP
  int thread_ok = fftw_init_threads();

  if (thread_ok) {
    fftw_plan_with_nthreads(ThreadManager::get_num_threads());
  }
#endif

  // fft itself is not mpi-parallelized.

  int ipe_x = Communicator::ipe(0);
  int ipe_y = Communicator::ipe(1);
  int ipe_z = Communicator::ipe(2);
  int ipe_t = Communicator::ipe(3);

  int npe_x = Communicator::npe(0);
  int npe_y = Communicator::npe(1);
  int npe_z = Communicator::npe(2);
  int npe_t = Communicator::npe(3);

  // split communicator along t directions
  int ipe_xyz = ipe_x + npe_x * (ipe_y + npe_y * (ipe_z));

  MPI_Comm_split(Communicator_impl::world(), ipe_t, ipe_xyz, &m_comm);

  // find rank and coordinate in subcommunicator
  int local_rank;
  MPI_Comm_rank(m_comm, &local_rank);

  int local_ipe_x = local_rank % npe_x;
  int local_ipe_y = (local_rank / npe_x) % npe_y;
  int local_ipe_z = (local_rank / npe_x / npe_y) % npe_z;

  // just for check
  if ((local_ipe_x != ipe_x) || (local_ipe_y != ipe_y) || (local_ipe_z != ipe_z)) {
    vout.crucial(m_vl, "%s: split commnicator failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_local_rank = local_rank;

  m_local_ipe_x = local_ipe_x;
  m_local_ipe_y = local_ipe_y;
  m_local_ipe_z = local_ipe_z;

  // lattice and grid information
  m_ndims = 3;

  m_grid_size.resize(m_ndims);
  m_grid_size[0] = CommonParameters::NPEx();
  m_grid_size[1] = CommonParameters::NPEy();
  m_grid_size[2] = CommonParameters::NPEz();

  m_grid_vol = 1;
  for (int i = 0; i < m_ndims; ++i) {
    m_grid_vol *= m_grid_size[i];
  }

  m_lattice_size.resize(m_ndims);
  m_lattice_size[0] = CommonParameters::Lx();
  m_lattice_size[1] = CommonParameters::Ly();
  m_lattice_size[2] = CommonParameters::Lz();

  m_lattice_vol = 1;
  for (int i = 0; i < m_ndims; ++i) {
    m_lattice_vol *= m_lattice_size[i];
  }

  m_local_size.resize(m_ndims);
  m_local_size[0] = CommonParameters::Nx();
  m_local_size[1] = CommonParameters::Ny();
  m_local_size[2] = CommonParameters::Nz();

  m_local_vol = 1;
  for (int i = 0; i < m_ndims; ++i) {
    m_local_vol *= m_local_size[i];
  }
}


//====================================================================
void FFT_3d_parallel3d::finalize()
{
  if (m_initialized) {
    release_plan();
  }
}


//====================================================================
void FFT_3d_parallel3d::create_mpi_datatype(int site_dof)
{
  // MPI datatypes and gather/scatter parameters

  // assume site_dof in complex

  MPI_Type_contiguous(2 * site_dof,
                      MPI_DOUBLE,
                      &m_site_vector_type);

  MPI_Type_commit(&m_site_vector_type);

  MPI_Type_contiguous(m_local_vol,
                      m_site_vector_type,
                      &m_local_patch_type);

  MPI_Type_commit(&m_local_patch_type);

  int size_;
  MPI_Type_size(m_site_vector_type, &size_);

  MPI_Datatype     type_;
  std::vector<int> local_origin(m_ndims, 0);

  MPI_Type_create_subarray(m_ndims,
                           &m_lattice_size[0],
                           &m_local_size[0],
                           &local_origin[0],
                           MPI_ORDER_FORTRAN,
                           m_site_vector_type,
                           &type_);

  MPI_Type_create_resized(type_, 0, size_, &m_subarray_type);

  MPI_Type_commit(&m_subarray_type);


  m_sendcounts.resize(m_grid_vol);

  for (int r = 0; r < m_grid_vol; ++r) {
    m_sendcounts[r] = 1;
  }

  m_subarray_displs.resize(m_grid_vol);

  for (int r = 0; r < m_grid_vol; ++r) {
    std::vector<int> coord = grid_rank_to_coord(r);

    // find global coordinate of origin of each local patch
    for (int j = 0; j < m_ndims; ++j) {
      coord[j] *= m_local_size[j];
    }

    int idx = find_global_index(coord);

    m_subarray_displs[r] = idx;
  }

  m_local_patch_displs.resize(m_grid_vol);

  for (int r = 0; r < m_grid_vol; ++r) {
    m_local_patch_displs[r] = r;
  }
}


//====================================================================
void FFT_3d_parallel3d::release_mpi_datatype()
{
  int is_finalized = 0;

  MPI_Finalized(&is_finalized);

  if (is_finalized) {
    vout.crucial(m_vl, "%s: MPI has already gone...\n", class_name.c_str());
    return;
  }

  MPI_Type_free(&m_site_vector_type);
  MPI_Type_free(&m_subarray_type);
  MPI_Type_free(&m_local_patch_type);
}


//====================================================================
void FFT_3d_parallel3d::create_fft_plan(int site_dof)
{
  // allocate buffer (run on-the-fly)
  m_buf = fftw_alloc_complex(site_dof * m_lattice_vol);
  if (!m_buf) {
    vout.crucial(m_vl, "%s: buffer allocation failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // create plan

  m_plan_fw = fftw_plan_many_dft(m_ndims, &m_lattice_size[0], site_dof,
                                 m_buf, NULL, site_dof, 1,
                                 m_buf, NULL, site_dof, 1,
                                 FFTW_FORWARD, FFTW_ESTIMATE);

  m_plan_bw = fftw_plan_many_dft(m_ndims, &m_lattice_size[0], site_dof,
                                 m_buf, NULL, site_dof, 1,
                                 m_buf, NULL, site_dof, 1,
                                 FFTW_BACKWARD, FFTW_ESTIMATE);

  if (!m_plan_fw || !m_plan_bw) {
    vout.crucial(m_vl, "%s: create plan failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void FFT_3d_parallel3d::release_fft_plan()
{
  if (m_buf) fftw_free(m_buf);
  m_buf = NULL;
  if (m_plan_fw) fftw_destroy_plan(m_plan_fw);
  m_plan_fw = NULL;
  if (m_plan_bw) fftw_destroy_plan(m_plan_bw);
  m_plan_bw = NULL;
}


//====================================================================
void FFT_3d_parallel3d::create_plan(int site_dof)
{
  create_mpi_datatype(site_dof);
  create_fft_plan(site_dof);

  m_site_dof = site_dof;

  m_initialized = true;
}


//====================================================================
void FFT_3d_parallel3d::release_plan()
{
  release_fft_plan();
  release_mpi_datatype();

  m_initialized = false;
}


//====================================================================
bool FFT_3d_parallel3d::need_create_plan(const Field& field)
{
  if (field.nin() / 2 == m_site_dof) return false;

  return true;
}


//====================================================================
void FFT_3d_parallel3d::fft(Field& dst, const Field& src, enum Direction dir)
{
  if (not ((dir == FORWARD) || (dir == BACKWARD))) {
    vout.crucial(m_vl, "%s: unsupported direction. %d\n", class_name.c_str(), dir);
    exit(EXIT_FAILURE);
  }

  // check if mpi types and fft plans are recyclable.
  if (m_initialized == false) {
    vout.general(m_vl, "%s: create plan.\n", class_name.c_str());
    create_plan(src.nin() / 2);
  } else {
    if (need_create_plan(src)) {
      vout.general(m_vl, "%s: discard plan and create new.\n", class_name.c_str());
      release_plan();
      create_plan(src.nin() / 2);
    } else {
      vout.general(m_vl, "%s: plan recycled.\n", class_name.c_str());
    }
  }

  int nex = src.nex();
  int nt  = CommonParameters::Nt();

  int ndata = nt * nex;

  std::vector<dcomplex *> src_array(ndata, nullptr);
  std::vector<dcomplex *> dst_array(ndata, nullptr);

  int local_vol = m_local_vol;

  int k = 0;
  for (int iex = 0; iex < nex; ++iex) {
    for (int it = 0; it < nt; ++it) {
      src_array[k] = (dcomplex *)(src.ptr(0, local_vol * it, iex));
      dst_array[k] = (dcomplex *)(dst.ptr(0, local_vol * it, iex));
      ++k;
    }
  }

  int nblock = m_grid_vol;

  for (int k = 0; k < ndata; k += nblock) {
    bool do_full = (k + nblock <= ndata);
    int  nwork   = do_full ? nblock : (ndata % nblock);

    if (do_full) {
      MPI_Alltoallv(src_array[k], &m_sendcounts[0], &m_local_patch_displs[0], m_local_patch_type,
                    m_buf, &m_sendcounts[0], &m_subarray_displs[0], m_subarray_type,
                    m_comm);
    } else {
      for (int j = 0; j < nwork; ++j) {
        MPI_Gatherv(src_array[k + j], 1, m_local_patch_type,
                    m_buf, &m_sendcounts[0], &m_subarray_displs[0], m_subarray_type,
                    j, m_comm);
      }
    }

    if (m_local_rank < nwork) {
      fftw_execute(dir == FORWARD ? m_plan_fw : m_plan_bw);
    }

    if (do_full) {
      MPI_Alltoallv(m_buf, &m_sendcounts[0], &m_subarray_displs[0], m_subarray_type,
                    dst_array[k], &m_sendcounts[0], &m_local_patch_displs[0], m_local_patch_type,
                    m_comm);
    } else {
      for (int j = 0; j < nwork; ++j) {
        MPI_Scatterv(m_buf, &m_sendcounts[0], &m_subarray_displs[0], m_subarray_type,
                     dst_array[k + j], 1, m_local_patch_type,
                     j, m_comm);
      }
    }
  }

  if (dir == BACKWARD) {
    scal(dst, 1.0 / m_lattice_vol);
  }
}


//====================================================================
void FFT_3d_parallel3d::fft(Field& dst, const Field& src)
{
  return fft(dst, src, m_direction);
}


//====================================================================
void FFT_3d_parallel3d::fft(Field& field)
{
  // return fft(field, field, m_direction);
  vout.crucial(m_vl, "Error at %s: fft on-the-fly unsupported.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
std::vector<int> FFT_3d_parallel3d::grid_rank_to_coord(int r)
{
  std::vector<int> coord(m_ndims);

  for (int i = 0; i < m_ndims; ++i) {
    coord[i] = r % m_grid_size[i];
    r       /= m_grid_size[i];
  }

  return coord;
}


//====================================================================
int FFT_3d_parallel3d::find_global_index(const std::vector<int>& coord)
{
  assert(coord.size() == m_ndims);

  int idx = coord[m_ndims - 1];
  for (int i = m_ndims - 2; i >= 0; --i) {
    idx *= m_lattice_size[i];
    idx += coord[i];
  }

  return idx;
}


//====================================================================
#endif /* USE_MPI */
#endif /* USE_FFTWLIB */

//====================================================================
//============================================================END=====
