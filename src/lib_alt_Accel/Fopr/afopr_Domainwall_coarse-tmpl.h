/*!
      @file    afopr_Domainwall_coarse-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-03-26 21:23:45 #$
      @version $LastChangedRevision: 2591 $
*/

#include "lib/Tools/decompose_LU_Cmplx.h"

#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
//#include "lib_alt_Accel/Fopr/afopr_Domainwall_dd.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Domainwall_coarse.h"

template<typename AFIELD>
const std::string AFopr_Domainwall_coarse<AFIELD>::class_name
                                         = "AFopr_Domainwall_coarse";

//====================================================================
namespace{

inline void accum_mult_u(real_t *v2, real_t *u2, real_t *v1,
			 int site, int nei,
			 int Ndf, int Ncol, int Ncol2)
{
  for(int i = 0; i < Ncol; ++i){
    real_t vtr = 0.0;
    real_t vti = 0.0;
    for(int j = 0; j < Ncol; ++j){
      real_t wtr = v1[IDX2(Ncol2, 2*j,   nei)];
      real_t wti = v1[IDX2(Ncol2, 2*j+1, nei)];
      real_t utr = u2[IDX2(Ndf, Ncol2*j + 2*i,   site)];
      real_t uti = u2[IDX2(Ndf, Ncol2*j + 2*i+1, site)];
      vtr += utr * wtr - uti * wti;
      vti += utr * wti + uti * wtr;
    }
    v2[IDX2(Ncol2, 2*i,   site)] += vtr;
    v2[IDX2(Ncol2, 2*i+1, site)] += vti;
  }
}

inline int index_mat(int ic1, int ic2, int Ncol){
  return ic1 + Ncol * ic2;
}

} // namespace

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general("%s: initial setup:\n", class_name.c_str());

  m_repr = "Dirac";  // now only the Dirac repr is available.

  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 0 if communication called in necessary

  int Ndim = CommonParameters::Ndim();

  do_comm_any = 0;
  for(int mu = 0; mu < Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    m_bc2[mu] = 1;
    if(do_comm[mu] > 0) m_bc2[mu] = 0;
  }

  m_bdsize.resize(Ndim);

  /*
  int fine_nvol = CommonParameters::Nvol();
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int NinF = 2 * Nc * Nd;
  //  workvec1.reset(NinF, fine_nvol, 1);
  //  workvec2.reset(NinF, fine_nvol, 1);
  //  workvec3.reset(NinF, fine_nvol, 1);
  //  workvec4.reset(NinF, fine_nvol, 1);
  */

  vout.general("%s: setup finished.\n", class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  int Ndim = CommonParameters::Ndim();

  // openacc device memory clean up
  for(int mu = 0; mu < Ndim; ++mu){

#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_dn1, m_bdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn2, m_bdsize[mu]);
    BridgeACC::afield_tidyup(buf_up1, m_bdsize[mu]);
    BridgeACC::afield_tidyup(buf_up2, m_bdsize[mu]);
#else
    real_t* buf_up = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_up, m_bdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn, m_bdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_parameters(
                                            const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int num_testvectors;
  std::vector<int> coarse_lattice;
  int fine_fopr_nin, fine_fopr_nvol, fine_fopr_nex;

  int err = 0;
  err += params.fetch_int("number_of_testvectors", num_testvectors);
  err += params.fetch_int_vector("coarse_lattice_size", coarse_lattice);
  err += params.fetch_int("fine_operator_nin",  fine_fopr_nin);
  err += params.fetch_int("fine_operator_nvol", fine_fopr_nvol);
  err += params.fetch_int("fine_operator_nex",  fine_fopr_nex);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                                                class_name.c_str());
  exit(EXIT_FAILURE);
 }

  set_parameters(num_testvectors, coarse_lattice,
                 fine_fopr_nin, fine_fopr_nvol, fine_fopr_nex);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_parameters(
                             const int num_testvectors,
                             const std::vector<int> &coarse_lattice,
                             const int fine_fopr_nin,
                             const int fine_fopr_nvol,
                             const int fine_fopr_nex)
{
  ThreadManager::assert_single_thread(class_name);

  int Ndim = CommonParameters::Ndim();
  assert(coarse_lattice.size() == Ndim);

  m_num_testvectors = num_testvectors;
  m_ncol = 2 * num_testvectors;  // number of chirality is multiplied.
  m_Nc = m_ncol;
  m_Nc2 = m_ncol * m_ncol;
  m_Nvc = 2 * m_Nc;        // 2 for complex
  m_Ndf = 2 * m_Nc * m_Nc; // 2 for complex

  m_fine_fopr_nin  = fine_fopr_nin;
  m_fine_fopr_nvol = fine_fopr_nvol;
  m_fine_fopr_nex  = fine_fopr_nex;

  m_Nx = coarse_lattice[0];
  m_Ny = coarse_lattice[1];
  m_Nz = coarse_lattice[2];
  m_Nt = coarse_lattice[3];
  m_Nst = m_Nx * m_Ny * m_Nz * m_Nt;
  m_Nsize[0] = coarse_lattice[0];
  m_Nsize[1] = coarse_lattice[1];
  m_Nsize[2] = coarse_lattice[2];
  m_Nsize[3] = coarse_lattice[3];

  m_bdsize[0] = m_Nvc * m_Ny * m_Nz * m_Nt;
  m_bdsize[1] = m_Nvc * m_Nx * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * m_Nx * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * m_Nx * m_Ny * m_Nz;

  setup_channels();

  size_t coarse_nvol = m_Nst;
  //  m_coarse_lvol = coarse_nvol * CommonParameters::NPE();

  m_Clov.reset(m_Ndf, m_Nst, 1);   // on-site term
  m_Clov_inv.reset(m_Ndf, m_Nst, 1);   // on-site term

  m_U.reset(m_Ndf, m_Nst, 2 * Ndim);   // hopping term
  m_U_unprec.reset(m_Ndf, m_Nst, 2 * Ndim);   // hopping term

  tmp_buffer1.resize(coarse_nvol);
  tmp_buffer2.resize(coarse_nvol);
  tmp_buffer_real1.resize(coarse_nvol);

  workvec1.reset(m_fine_fopr_nin, m_fine_fopr_nvol, m_fine_fopr_nex);
  workvec2.reset(m_fine_fopr_nin, m_fine_fopr_nvol, m_fine_fopr_nex);
  workvec3.reset(m_fine_fopr_nin, m_fine_fopr_nvol, m_fine_fopr_nex);

  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  number of testvectors = %d\n",
               m_num_testvectors);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  coarse_lattice_size[%d] = %2d\n",
                 mu, coarse_lattice[mu]);
  }

  vout.general(m_vl, "Parameters of %s: setup finished\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::setup_channels()
{
  //check
  if( (m_Nst/m_Nx % NWP) != 0){
    vout.crucial(m_vl, "x-slice volume(%d) must be multiple of NWP(%d)\n",
                 m_Nst/m_Nx, NWP);
    exit(EXIT_FAILURE);
  }

  if( (m_Nst/m_Ny % NWP) != 0){
    vout.crucial(m_vl, "y-slice volume(%d) must be multiple of NWP(%d)\n",
                 m_Nst/m_Ny, NWP);
    exit(EXIT_FAILURE);
  }

  if( (m_Nst/m_Nz % NWP) != 0){
    vout.crucial(m_vl, "z-slice volume(%d) must be multiple of NWP(%d)\n",
                 m_Nst/m_Nz, NWP);
    exit(EXIT_FAILURE);
  }

  if( (m_Nst/m_Nt % NWP) != 0){
    vout.crucial(m_vl, "t-slice volume(%d) must be multiple of NWP(%d)\n",
                 m_Nst/m_Nt, NWP);
    exit(EXIT_FAILURE);
  }

  int Ndim = CommonParameters::Ndim();

  chsend_up.resize(Ndim);
  chrecv_up.resize(Ndim);
  chsend_dn.resize(Ndim);
  chrecv_dn.resize(Ndim);

  for(int mu = 0; mu < Ndim; ++mu){

    size_t Nvsize = m_bdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu,  1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if(do_comm[mu] == 1){
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_init(buf_dn1, m_bdsize[mu]);
    BridgeACC::afield_init(buf_dn2, m_bdsize[mu]);
    BridgeACC::afield_init(buf_up1, m_bdsize[mu]);
    BridgeACC::afield_init(buf_up2, m_bdsize[mu]);
#else
    BridgeACC::afield_init((real_t*)buf_up, m_bdsize[mu]);
    BridgeACC::afield_init((real_t*)buf_dn, m_bdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::generate_coarse_op(
                                AFopr_dd<AFIELD>* fine_afopr,
                                ASolver<AFIELD> *solver_PV,
                                const std::vector<AFIELD>& atestvec)
{
  Timer timer;
  timer.start();

  Timer timer2;

  const int coarse_nvol = m_U.nvol();

  int ith, nth, coarse_is, coarse_ns;
  set_threadtask_afopr(ith, nth, coarse_is, coarse_ns, coarse_nvol);

  real_t    *out_clov   = m_Clov.ptr(0);
  real_t    *out_gauge  = m_U_unprec.ptr(0);
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;
  assert(m_Nc2 == coarse_Nc2);

  m_Clov.set(0.0);
  m_Clov.update_host();

  m_U_unprec.set(0.0);
  m_U_unprec.update_host();

#pragma omp barrier

  std::vector<int> coarse_lattice(4);
  coarse_lattice[0] = m_Nx;
  coarse_lattice[1] = m_Ny;
  coarse_lattice[2] = m_Nz;
  coarse_lattice[3] = m_Nt;
  AIndex_block_lex<real_t, ACCEL> index_block(coarse_lattice);

  AIndex_coarse_lex<real_t, ACCEL> index_coarse(
                          m_Nx, m_Ny, m_Nz, m_Nt, num_vectors, 2);

#pragma omp barrier

  // coarse clover:
  //  I = 2 * i1 + chirality1
  //  J = 2 * i2 + chirality2
  //  JI = (2 * num_vectors) * J + I
  //  <J|D|I>

  for (int i1 = 0; i1 < num_vectors; ++i1) {
    for (int ch1 = -1; ch1 < 2; ch1 += 2) { // ch=-1,+1

      for(int ieo = 0; ieo < 2; ++ieo) {

        // diag block: "clover"
        //project_chiral(workvec2, atestvec[i1], ch1, fine_afopr);
        fine_afopr->mult_gm5(workvec2, atestvec[i1]);
        aypx(real_t(ch1), workvec2, atestvec[i1]);

        // calculate work1 := M_{PV})_{SAP} |i1, eo>
        int ieo_mask = 1-ieo;
        int Nconv;
        real_t diff;
        block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo_mask, index_block);
        //  vout.general("hoge: calling Solver_PV\n");
        //solver_PV->solve(workvec1, workvec2, Nconv, diff);
        solver_PV->solve(workvec1, workvec2, Nconv, diff, ieo);
        /*
        {
          double w1=workvec1.norm2();
          double w2=workvec2.norm2();
          vout.general("i1=%d, ch1=%d, ieo=%d, w1=%23.15e, w2=%23.15e\n",
                       i1, ch1, ieo, w1, w2);
        }
        */
        fine_afopr->mult(workvec2, workvec1);
        block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo_mask, index_block);

        real_t *out = out_clov;
        int    I    = 2 * i1 + (ch1 + 1) / 2;

        for (int i2 = 0; i2 < num_vectors; ++i2) {

          //project_chiral(workvec3, workvec2, -1, fine_afopr);
          fine_afopr->mult_gm5(workvec3, workvec2);
          aypx(real_t(-1.0), workvec3, workvec2);
#pragma omp barrier

          block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                        ieo, index_block);
#pragma omp barrier

          //project_chiral(workvec3, workvec2, 1, fine_afopr);
          fine_afopr->mult_gm5(workvec3, workvec2);
          aypx(real_t(1.0), workvec3, workvec2);
#pragma omp barrier

          block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                        ieo, index_block);

          int J = 2 * i2;
          for (int s = coarse_is; s < coarse_ns; s++) {
            int idx_r = index_coarse.idx_Gr(J, I, s, 0);
            int idx_i = index_coarse.idx_Gi(J, I, s, 0);
            out[idx_r] += real(tmp_buffer1[s]);
            out[idx_i] += imag(tmp_buffer1[s]);
          }

          ++J;
          for (int s = coarse_is; s < coarse_ns; s++) {
            int idx_r = index_coarse.idx_Gr(J, I, s, 0);
            int idx_i = index_coarse.idx_Gi(J, I, s, 0);
            out[idx_r] += real(tmp_buffer2[s]);
            out[idx_i] += imag(tmp_buffer2[s]);
          }

#pragma omp barrier
        } // i2


        // hopping block: "gauge "
        for (int mu = 0; mu < 8; mu++) {
          workvec2.set(0.0);
          if( mu<4 ) {
            fine_afopr->mult_dup(workvec2, workvec1, mu);
          } else {
            fine_afopr->mult_ddn(workvec2, workvec1, mu-4);
          }
          // block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo, index_block);

          real_t *out = out_gauge + mu * 2 * coarse_Nc2 * coarse_nvol;
          // mu comes last

          int I = 2 * i1 + (ch1 + 1) / 2;
          for (int i2 = 0; i2 < num_vectors; ++i2) {
            //project_chiral(workvec3, workvec2, -1, fine_afopr);
            fine_afopr->mult_gm5(workvec3, workvec2);
            aypx(real_t(-1.0), workvec3, workvec2);
#pragma omp barrier

	    block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                          1-ieo, index_block);
#pragma omp barrier

            //project_chiral(workvec3, workvec2, 1, fine_afopr);
            fine_afopr->mult_gm5(workvec3, workvec2);
            aypx(real_t(1.0), workvec3, workvec2);
#pragma omp barrier

            block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                          1-ieo, index_block);

            int J = 2 * i2;
            for (int s = coarse_is; s < coarse_ns; ++s) {
              int idx_r = index_coarse.idx_Gr(J, I, s, 0);
              int idx_i = index_coarse.idx_Gi(J, I, s, 0);
              out[idx_r] += real(tmp_buffer1[s]);
              out[idx_i] += imag(tmp_buffer1[s]);
            }

            ++J;
            for (int s = coarse_is; s < coarse_ns; ++s) {
              int idx_r = index_coarse.idx_Gr(J, I, s, 0);
              int idx_i = index_coarse.idx_Gi(J, I, s, 0);
              out[idx_r] += real(tmp_buffer2[s]);
              out[idx_i] += imag(tmp_buffer2[s]);
            }

#pragma omp barrier
          } // i2
        }   // mu
      } // ieo
    }
  }       // i1, ch1

  make_preconditioned();

  if(ith == 0){
    m_U_unprec.update_device();
    m_Clov.update_device();
    m_Clov_inv.update_device();
    m_U.update_device();
  }

  // rescale operator as project chiral does not have 1/2
  m_Clov.scal(0.25);
  m_U_unprec.scal(0.25);
  m_Clov_inv.scal(4.0);

  {
    double clv2  = m_Clov.norm2();
    double u2_un = m_U_unprec.norm2();
    vout.general("%s: |m_Clov|^2     = %23.15e\n", class_name.c_str(), clv2);
    vout.general("%s: |m_U_unprec|^2 = %23.15e\n", class_name.c_str(), u2_un);

    double clv2_inv = m_Clov_inv.norm2();
    double u2   = m_U.norm2();
    vout.general("%s: after preconditioned operator\n", class_name.c_str());
    vout.general("%s: |m_Clov_inv|^2 = %23.15e\n", class_name.c_str(), clv2_inv);
    vout.general("%s: |m_U|^2        = %23.15e\n", class_name.c_str(), u2);
  }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "%s: generate_coarse_op finished: %e sec\n",
               class_name.c_str(), elapsed_time);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::make_preconditioned()
{
  int ith, nth, is, ns;
  const int coarse_nvol = m_U.nvol();
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;

  AIndex_coarse_lex<real_t, ACCEL> index_coarse(
                           m_Nx, m_Ny, m_Nz, m_Nt, num_vectors, 2);

  std::vector<double> mat(2*coarse_Nc2);

  std::vector<double> mat2(2*coarse_Nc2);

  const real_t *clov  = m_Clov.ptr(0);
  const real_t *gauge = m_U_unprec.ptr(0);
  real_t *clov_inv  = m_Clov_inv.ptr(0);
  real_t *gauge_out = m_U.ptr(0);

  set_threadtask(ith, nth, is, ns, coarse_nvol);

  Decompose_LU_Cmplx lu(m_ncol);

  for(int site = is; site < ns; ++site){

    for(int i = 0; i < 2 * coarse_Nc2; ++i) {
      mat[i] = clov[index_coarse.idx_G(i, site, 0)];
    }

    lu.set_matrix(&mat[0]);

    // preconditioning matrix
    lu.get_inverse(&mat[0]);

    for(int i = 0; i < 2 * coarse_Nc2; ++i) {
      clov_inv[index_coarse.idx_G(i, site, 0)] = mat[i];
    }

    // preconditioned matrix
    for(int mu = 0; mu < 8; ++mu) {

      /*
      for(int i = 0; i < 2 * coarse_Nc2; ++i) {
        mat[i] = gauge[index_coarse.idx_G(i, site, mu)];
      }

      lu.mult_inverse(&mat[0]);

      for(int i = 0; i < 2 * coarse_Nc2; ++i) {
        gauge_out[index_coarse.idx_G(i, site, mu)] = mat[i];
      }
      */

      for(int i = 0; i < 2 * coarse_Nc2; ++i) {
        mat2[i] = gauge[index_coarse.idx_G(i, site, mu)];
      }

      for(int c1 = 0; c1 < m_ncol; ++c1){
        for(int c2 = 0; c2 < m_ncol; ++c2){
          double matr = 0.0;
          double mati = 0.0;
          for(int k = 0; k < m_ncol; ++k){
            int idx1 = index_mat(c1, k, m_ncol);
            int idx2 = index_mat(k, c2, m_ncol);
            double mat1r = mat[2*idx1];
            double mat1i = mat[2*idx1+1];
            double mat2r = mat2[2*idx2];
            double mat2i = mat2[2*idx2+1];
            matr += mat1r * mat2r - mat1i * mat2i;
            mati += mat1r * mat2i + mat1i * mat2r;
          }
          int idx = index_mat(c1, c2, m_ncol);
          gauge_out[index_coarse.idx_G(2*idx,   site, mu)] = matr;
          gauge_out[index_coarse.idx_G(2*idx+1, site, mu)] = mati;
        }	
      }	

    }
  }

  // debug
  /*
  int mu = 0;
  for(int site = 0; site < 4; ++site) {
    for(int i = 0; i < 4; ++i) {
      vout.general("coarse gauge: %d  %d  %d  %e\n", i, site, mu,
                   gauge_out[index_coarse.idx_G(i, site, mu)]);
    }
  }
  */
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_config(Field* u){
  vout.crucial(m_vl, "%s: set_config is called\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_up(int mu, AFIELD &v,
                                              const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  w.update_host();
  v.update_host();

  if(mu == 0){
    mult_xp(vp, wp);
  }else if(mu == 1){
    mult_yp(vp, wp);
  }else if(mu == 2){
    mult_zp(vp, wp);
  }else if(mu == 3){
    mult_tp(vp, wp);
  }else{
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

  v.update_device();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_dn(int mu, AFIELD &v,
                                              const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  w.update_host();
  v.update_host();

  if(mu == 0){
    mult_xm(vp, wp);
  }else if(mu == 1){
    mult_ym(vp, wp);
  }else if(mu == 2){
    mult_zm(vp, wp);
  }else if(mu == 3){
    mult_tm(vp, wp);
  }else{
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }

  v.update_device();

}

//====================================================================
template<typename AFIELD>
std::string AFopr_Domainwall_coarse<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                           const string mode)
{
  if (mode == "prec") {
    mult_prec(v, w);
    //mult_prec_alt(v, w);
  } else if (mode == "clov") {
    mult_csw(v, w);
    //mult_csw_alt(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
  //mult_D_alt(v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth;
  set_thread(ith, nth);

#pragma omp barrier

  if(ith == 0){
    BridgeACC::mult_domainwall_coarse_mult_gm5(vp, wp, m_ncol, m_Nsize);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_gm5_host(real_t *v, real_t *w)
{ // this is implementation on the host

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    for(int i = 0; i < m_ncol/2; ++i){
      int chm = 2 * i;      // negative chirality
      int chp = 2 * i + 1;  // positive chirality
      v[IDX2(m_Nvc, 2*chm,   site)] = -w[IDX2(m_Nvc, 2*chm,   site)];
      v[IDX2(m_Nvc, 2*chm+1, site)] = -w[IDX2(m_Nvc, 2*chm+1, site)];
      v[IDX2(m_Nvc, 2*chp,   site)] =  w[IDX2(m_Nvc, 2*chp,   site)];
      v[IDX2(m_Nvc, 2*chp+1, site)] =  w[IDX2(m_Nvc, 2*chp+1, site)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_D(AFIELD &v, const AFIELD &w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u  = m_U.ptr(0);
  real_t *ct = m_Clov.ptr(0);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){

    if(do_comm_any > 0){

      real_t *buf1xp = (real_t*)chsend_dn[0].ptr();
      real_t *buf1xm = (real_t*)chsend_up[0].ptr();

      real_t *buf1yp = (real_t*)chsend_dn[1].ptr();
      real_t *buf1ym = (real_t*)chsend_up[1].ptr();

      real_t *buf1zp = (real_t*)chsend_dn[2].ptr();
      real_t *buf1zm = (real_t*)chsend_up[2].ptr();

      real_t *buf1tp = (real_t*)chsend_dn[3].ptr();
      real_t *buf1tm = (real_t*)chsend_up[3].ptr();

      BridgeACC::mult_domainwall_coarse_1(
                              buf1xp, buf1xm, buf1yp, buf1ym,
                              buf1zp, buf1zm, buf1tp, buf1tm,
                              v1, m_ncol, m_Nsize, do_comm);
      chset_send.start();
      chset_recv.start();
    }

    BridgeACC::mult_domainwall_coarse_bulk(v2, u, v1,
                                       m_ncol, m_Nsize, m_bc2);

    if(do_comm_any > 0){

      chset_send.wait();
      chset_recv.wait();

      real_t *buf2xp = (real_t*)chrecv_up[0].ptr();
      real_t *buf2xm = (real_t*)chrecv_dn[0].ptr();

      real_t *buf2yp = (real_t*)chrecv_up[1].ptr();
      real_t *buf2ym = (real_t*)chrecv_dn[1].ptr();

      real_t *buf2zp = (real_t*)chrecv_up[2].ptr();
      real_t *buf2zm = (real_t*)chrecv_dn[2].ptr();

      real_t *buf2tp = (real_t*)chrecv_up[3].ptr();
      real_t *buf2tm = (real_t*)chrecv_dn[3].ptr();

      BridgeACC::mult_domainwall_coarse_2(v2, u,
                                      buf2xp, buf2xm, buf2yp, buf2ym,
                                      buf2zp, buf2zm, buf2tp, buf2tm,
                                      m_ncol, m_Nsize, do_comm);

    }

  } // if(ith == 0)

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_csw_alt(AFIELD& v, const AFIELD& w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  w.update_host();

  mult_csw_host(vp, wp);

  v.update_device();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_csw_host(real_t *v2, real_t *v1)
{
  real_t *u = m_Clov.ptr(0);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    accum_mult_u(v2, u, v1, site, site, m_Ndf, m_ncol, m_Nvc);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_csw(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *ct = m_Clov.ptr(0);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::mult_domainwall_coarse_prec(v2, ct, v1, m_ncol, m_Nsize);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_prec(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *ct = m_Clov_inv.ptr(0);

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    BridgeACC::mult_domainwall_coarse_prec(v2, ct, v1, m_ncol, m_Nsize);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_prec_alt(AFIELD& v,
                                                    const AFIELD& w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  w.update_host();

  mult_prec_host(vp, wp);

  v.update_device();
  
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_prec_host(real_t *v2, real_t *v1)
{
  real_t *u = m_Clov_inv.ptr(0);

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

#pragma omp barrier

  for(int site = is; site < ns; ++site){
    accum_mult_u(v2, u, v1, site, site, m_Ndf, m_ncol, m_Nvc);
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_D_alt(AFIELD &v, const AFIELD &w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  copy(v,w);

  w.update_host();

#pragma omp barrier

  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);

  v.update_device();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::H(AFIELD &v, const AFIELD &w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::clear_host(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  for(int site = is; site < ns; ++site){
    for(int ivc = 0; ivc < m_Nvc; ++ivc){
      v[IDX2(m_Nvc, ivc, site)] = 0.0;
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_dn[0].ptr();
  real_t *buf2 = (real_t*)chrecv_up[0].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if(do_comm[0] > 0){

    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nx;
      int iyzt = site/m_Nx;
      if(ix == 0){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, iyzt)] = v1[IDX2(m_Nvc, ivc, site)];
        }
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_dn[0].start();
    chrecv_up[0].start();
    chsend_dn[0].wait();
    chrecv_up[0].wait();
  }
#pragma omp barrier

  } // if(do_comm[0] == 1)

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx;
    int iyzt = site/m_Nx;

    if(ix < m_Nx-1 || do_comm[0] == 0){
      int nei = ix+1 + m_Nx * iyzt;
      if(ix == m_Nx-1) nei = 0 + m_Nx * iyzt;
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_ncol, m_Nvc);
    }else{
      accum_mult_u(v2, u, buf2, site, iyzt, m_Ndf, m_ncol, m_Nvc);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_up[0].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[0].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir + NDIM));

  int NcolH = m_ncol/2;

#pragma omp barrier

  if(do_comm[0] > 0){

    for(int site = is; site < ns; ++site){
      int ix   = site % m_Nx;
      int iyzt = site/m_Nx;
      if(ix == m_Nx-1){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, iyzt)] = 0.0;
	}
        accum_mult_u(buf1, u, v1, iyzt, site, m_Ndf, m_Nvc, NcolH);
      }
    }

#pragma omp barrier
#pragma omp master
 {
    chsend_up[0].start();
    chrecv_dn[0].start();
    chsend_up[0].wait();
    chrecv_dn[0].wait();
 }
#pragma omp barrier

  } // end of if(do_comm[0] > 0)

  for(int site = is; site < ns; ++site){
    int ix   = site % m_Nx;
    int iyzt = site/m_Nx;

    if(ix > 0 || do_comm[0] == 0){
      int nei = ix-1 + m_Nx * iyzt;
      if(ix == 0) nei = m_Nx-1 + m_Nx * iyzt;
        accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_Nvc, NcolH);
    }else{
      for(int ivc = 0; ivc < m_Nvc; ++ivc){
        v2[IDX2(m_Nvc, ivc, site)] += buf2[IDX2(m_Nvc, ivc, iyzt)];
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy = m_Nx * m_Ny;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_dn[1].ptr();
  real_t *buf2 = (real_t*)chrecv_up[1].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if(do_comm[1] > 0){

    for(int site = is; site < ns; ++site){
      int ix  = site % m_Nx;
      int iy  = (site/m_Nx) % m_Ny;
      int izt = site/Nxy;
      int ixzt = ix + m_Nx * izt;
      if(iy == 0){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixzt)] = v1[IDX2(m_Nvc, ivc, site)];
        }
      }
    }

#pragma omp barrier

#pragma omp master
   {
    chsend_dn[1].start();
    chrecv_up[1].start();
    chsend_dn[1].wait();
    chrecv_up[1].wait();
   }

#pragma omp barrier

  }  // end of if(do_comm[1] > 0)

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nx;
    int iy  = (site/m_Nx) % m_Ny;
    int izt = site/Nxy;

    if(iy < m_Ny-1 || do_comm[1] == 0){
      int iy2 = (iy+1) % m_Ny;
      int nei = ix + m_Nx * (iy2 + m_Ny * izt);
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_ncol, m_Nvc);
    }else{
      int ixzt = ix + m_Nx * izt;
      accum_mult_u(v2, u, buf2, site, ixzt, m_Ndf, m_ncol, m_Nvc);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxy = m_Nx * m_Ny;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_up[1].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[1].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir + NDIM));

  int NcolH = m_ncol/2;

#pragma omp barrier

  if(do_comm[1] > 0){


    for(int site = is; site < ns; ++site){
      int ix  = site % m_Nx;
      int iy  = (site/m_Nx) % m_Ny;
      int izt = site/Nxy;
      if(iy == m_Ny-1){
        int ixzt = ix + m_Nx * izt;
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixzt)] = 0.0;
	}
        accum_mult_u(buf1, u, v1, ixzt, site, m_Ndf, m_Nvc, NcolH);
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_up[1].start();
    chrecv_dn[1].start();
    chsend_up[1].wait();
    chrecv_dn[1].wait();
  }

#pragma omp barrier

  }

  for(int site = is; site < ns; ++site){
    int ix  = site % m_Nx;
    int iy  = (site/m_Nx) % m_Ny;
    int izt = site/Nxy;

    if(iy != 0 || do_comm[idir] == 0){
      int iy2 = (iy-1 + m_Ny) % m_Ny;
      int nei = ix + m_Nx * (iy2 + m_Ny * izt);
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_Nvc, NcolH);
    }else{
      int ixzt = ix + m_Nx * izt;
      for(int ivc = 0; ivc < m_Nvc; ++ivc){
        v2[IDX2(m_Nvc, ivc, site)] += buf2[IDX2(m_Nvc, ivc, ixzt)];
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy = m_Nx * m_Ny;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_dn[2].ptr();
  real_t *buf2 = (real_t*)chrecv_up[2].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if(do_comm[2] > 0){

    for(int site = is; site < ns; ++site){
      int ixy = site % Nxy;
      int iz  = (site/Nxy) % m_Nz;
      int it  = site/(Nxy * m_Nz);
      int ixyt = ixy + Nxy * it;
      if(iz == 0){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixyt)] = v1[IDX2(m_Nvc, ivc, site)];
        }
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_dn[2].start();
    chrecv_up[2].start();
    chsend_dn[2].wait();
    chrecv_up[2].wait();
  }

#pragma omp barrier

  }

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);

    if(iz != m_Nz-1 || do_comm[2] == 0){
      int iz2 = (iz+1) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz*it);
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_ncol, m_Nvc);
    }else{
      int ixyt = ixy + Nxy * it;
      accum_mult_u(v2, u, buf2, site, ixyt, m_Ndf, m_ncol, m_Nvc);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxy = m_Nx * m_Ny;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_up[2].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[2].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir + NDIM));

  int NcolH = m_ncol/2;

#pragma omp barrier

  if(do_comm[2] > 0){

    for(int site = is; site < ns; ++site){
      int ixy = site % Nxy;
      int iz  = (site/Nxy) % m_Nz;
      int it  = site/(Nxy * m_Nz);
      if(iz == m_Nz-1){
        int ixyt = ixy + Nxy * it;
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixyt)] = 0.0;
	}
        accum_mult_u(buf1, u, v1, ixyt, site, m_Ndf, m_Nvc, NcolH);
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_up[2].start();
    chrecv_dn[2].start();
    chsend_up[2].wait();
    chrecv_dn[2].wait();
  }

#pragma omp barrier

  }

  for(int site = is; site < ns; ++site){
    int ixy = site % Nxy;
    int iz  = (site/Nxy) % m_Nz;
    int it  = site/(Nxy * m_Nz);

    if(iz > 0 || do_comm[2] == 0){
      int iz2 = (iz-1 + m_Nz) % m_Nz;
      int nei = ixy + Nxy * (iz2 + m_Nz*it);
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_Nvc, NcolH);
    }else{
      int ixyt = ixy + Nxy * it;
      for(int ivc = 0; ivc < m_Nvc; ++ivc){
        v2[IDX2(m_Nvc, ivc, site)] += buf2[IDX2(m_Nvc, ivc, ixyt)];
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nx * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_dn[3].ptr();
  real_t *buf2 = (real_t*)chrecv_up[3].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if(do_comm[3] > 0){

    for(int site = is; site < ns; ++site){
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if(it == 0){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixyz)] = v1[IDX2(m_Nvc, ivc, site)];
        }
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_dn[3].start();
    chrecv_up[3].start();
    chsend_dn[3].wait();
    chrecv_up[3].wait();
  }

#pragma omp barrier

  }

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    if(it < m_Nt-1 || do_comm[3] == 0){
      int it2 = (it+1) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_ncol, m_Nvc);
    }else{
      accum_mult_u(v2, u, buf2, site, ixyz, m_Ndf, m_ncol, m_Nvc);
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir = 3;
  int Nxyz = m_Nx * m_Ny * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_afopr(ith, nth, is, ns, m_Nst);

  real_t *buf1 = (real_t*)chsend_up[3].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[3].ptr();
  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir + NDIM));

  int NcolH = m_ncol/2;

#pragma omp barrier

  if(do_comm[3] > 0){

    for(int site = is; site < ns; ++site){
      int ixyz = site % Nxyz;
      int it   = site / Nxyz;
      if(it == m_Nt-1){
        for(int ivc = 0; ivc < m_Nvc; ++ivc){
          buf1[IDX2(m_Nvc, ivc, ixyz)] = 0.0;
	}
        accum_mult_u(buf1, u, v1, ixyz, site, m_Ndf, m_Nvc, NcolH);
      }
    }

#pragma omp barrier

#pragma omp master
  {
    chsend_up[3].start();
    chrecv_dn[3].start();
    chsend_up[3].wait();
    chrecv_dn[3].wait();
  }
#pragma omp barrier

  }

  for(int site = is; site < ns; ++site){
    int ixyz = site % Nxyz;
    int it   = site / Nxyz;

    if(it > 0 || do_comm[3] == 0){
      int it2 = (it-1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyz * it2;
      accum_mult_u(v2, u, v1, site, nei, m_Ndf, m_Nvc, NcolH);
    }else{
      for(int ivc = 0; ivc < m_Nvc; ++ivc){
        v2[IDX2(m_Nvc, ivc, site)] += buf2[IDX2(m_Nvc, ivc, ixyz)];
      }
    }

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_coarse<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int NPE  = CommonParameters::NPE();
  int Lvol = m_Nx * m_Ny * m_Nz * m_Nt * NPE;
  int Ndim = CommonParameters::Ndim();

  double flop_site, flop;

  // Note that m_ncol = 2 * m_num_testvectors;
  flop_site = static_cast<double>(8 * m_ncol * m_ncol * (2 * Ndim + 1));
  flop = flop_site * static_cast<double>(Lvol);

  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}

//============================================================END=====
