/*!
      @file    MultiGrid_Clover-tmpl.h
      @brief   MultiGrid operation for Clover fermion
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2#$
      @version $LastChangedRevision: 2499 $
 */

#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/randomNumbers_Mseries.h"
#include "lib/Tools/timer.h"

//====================================================================
namespace {
  // memory saving implementation: no work array of fields are needed
  template<typename INDEX, typename AFIELD>
  void gramschmidt_in_block_impl(std::vector<AFIELD>& v,
                                 AFopr<AFIELD> *afopr,
                                 const INDEX& block_index,
                                 AFIELD& tmp1, AFIELD& tmp2,
                                 typename AFIELD::real_t *block_norm,
                                 typename AFIELD::complex_t *block_prod)
  {
    typedef typename AFIELD::real_t      real_t;
    typedef typename AFIELD::complex_t   complex_t;

    vout.detailed("gramschimdt_in_block, called (general version with smaller memory): size of(real_t)=%d\n", sizeof(real_t));

#ifdef DEBUG_GS
    {
      double norm2 = v[0].norm2();
      vout.general("initial: n2[0]=%e\n", norm2);
    }
#endif

    int nvector = v.size();

    int coarse_nvol = block_index.coarse_nvol();
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, coarse_nvol);


    // make <v[j]|v[i]> = 0 for j < i
    for (int i = 0; i < nvector; i++) {
      // vout.general("i = %d  norm2 = %f\n", i, v[i].norm2());

#pragma omp barrier

      // linear algebras on the fine lattice
      copy(tmp1, v[i]);
#pragma omp barrier

      afopr->mult_gm5(tmp2, v[i]);

      axpy(tmp1, real_t(1.0), tmp2);  // tmp1  = (1 + gm5) v[i]
      axpy(v[i], real_t(-1.0), tmp2); // v[i] := (1 - gm5) v[i]
#pragma omp barrier                   // work on the fine lattice, done

      // linear algebras with block
      // ch=-1
      for (int j = 0; j < i; j++) {
        // P|v[i]> := P|v[i]> -  |v[j]> <v[j]|P|v[i]>
        block_dotc(block_prod, v[j], v[i], block_index);
        block_axpy(v[i], block_prod, v[j], real_t(-1.0), block_index);
      }
#pragma omp barrier  // work on the blocks, done

      // normalize
      // linear algebras on the fine lattice
      afopr->mult_gm5(tmp2, v[i]);
      axpy(v[i], real_t(-1.0), tmp2); // apply P(-)
#pragma omp barrier                   // work on the fine lattice, done

      // linear algebras with block
      block_norm2(block_norm, v[i], block_index);
#pragma omp barrier  // work on the blocks, done

      for (int k = is; k < ns; k++) {
        double n2 = block_norm[k];
        block_norm[k] = 1.0 / sqrt(n2);
      }
#pragma omp barrier

      block_scal(v[i], block_norm, block_index);


      // ch=+1
      for (int j = 0; j < i; ++j) {
        // |Pv[i]> := |Pv[i]> - |v[j]> <v[j]|Pv[i]>
        block_dotc(block_prod, v[j], tmp1, block_index);
        block_axpy(tmp1, block_prod, v[j], real_t(-1.0), block_index);
      }
#pragma omp barrier // work on blocks, done


      // normalize

      // linear algebras on the fine lattice
      afopr->mult_gm5(tmp2, tmp1);
      axpy(tmp1, real_t(1.0), tmp2); // apply P
#pragma omp barrier                  // work on the fine lattice, done

      // linear algebra with blocks
      block_norm2(block_norm, tmp1, block_index);
#pragma omp barrier  // work on the blocks, done

      for (int k = is; k < ns; ++k) {
        double n2 = block_norm[k];
        block_norm[k] = 1.0 / sqrt(n2);
      }
#pragma omp barrier

      // merge both chiralities
      block_axpy(v[i], block_norm, tmp1, real_t(1.0), block_index);
    } // i

#pragma omp barrier

#ifdef DEBUG
    {
      double norm2 = v[0].norm2();
      vout.general("  after Gramschmidt n2[0]=%e\n", norm2);
    }
#endif
  }


  //====================================================================
  template<typename INDEX, class AFIELD1, class AFIELD2>
  void make_coarse_vector_impl(AFIELD1& coarse_vector, const AFIELD2& fine_vector,
                               const std::vector<AFIELD2>& testvectors,
                               AFopr<AFIELD2> *afopr, const INDEX& block_index,
                               AFIELD2& tmp1, AFIELD2& tmp2, typename AFIELD1::complex_t *inner_prod)
  {
#ifdef DEBUG
    vout.detailed("%s is called\n", __func__);
#endif
    typedef typename AFIELD1::complex_t   complex_t;
    typedef typename AFIELD1::real_t      real_t;

    const int num_vectors = testvectors.size();
    assert(coarse_vector.nin() == 4 * num_vectors);  // 4: 2 for complex, 2 for chirality

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, block_index.coarse_nvol());
    int coarse_Nx = block_index.coarse_lattice_size(0);
    int coarse_Ny = block_index.coarse_lattice_size(1);
    int coarse_Nz = block_index.coarse_lattice_size(2);
    int coarse_Nt = block_index.coarse_lattice_size(3);
    AIndex_coarse_lex<real_t, AFIELD1::IMPL> index(coarse_Nx, coarse_Ny, coarse_Nz, coarse_Nt, num_vectors, 2);



    for (int i = 0; i < num_vectors; i++) {
#pragma omp barrier

      //linear algebra on fine lattice
      copy(tmp2, testvectors[i]);
      afopr->mult_gm5(tmp1, tmp2);
      axpy(tmp2, -1.0, tmp1);            // tmp2=(1-gm5)v[i]
      axpy(tmp1, +1.0, testvectors[i]);  // tmp1=(1+gm5)v[i]

#pragma omp barrier

      // linear algerbar with blocks

      // ch=-1
      block_dotc(inner_prod, tmp2, fine_vector, block_index);
      for (int block_id = is; block_id < ns; block_id++) {
        int index_re = index.idx_SPr(i, 0, block_id, 0);
        int index_im = index.idx_SPi(i, 0, block_id, 0);
        coarse_vector.set(index_re, 0.5 * real(inner_prod[block_id]));
        coarse_vector.set(index_im, 0.5 * imag(inner_prod[block_id]));
      }

      // ch=+1
      block_dotc(inner_prod, tmp1, fine_vector, block_index);
      for (int block_id = is; block_id < ns; block_id++) {
        int index_re = index.idx_SPr(i, 1, block_id, 0);
        int index_im = index.idx_SPi(i, 1, block_id, 0);
        coarse_vector.set(index_re, 0.5 * real(inner_prod[block_id]));
        coarse_vector.set(index_im, 0.5 * imag(inner_prod[block_id]));
      }
    } // i
#pragma omp barrier
  }


  //====================================================================
  template<typename INDEX, class AFIELD1, class AFIELD2>
  void make_fine_vector_impl(AFIELD2& fine_vector, const AFIELD1& coarse_vector,
                             const std::vector<AFIELD2>& testvectors,
                             AFopr<AFIELD2> *afopr, const INDEX& block_index,
                             AFIELD2& tmp1, AFIELD2& tmp2, typename AFIELD1::complex_t *complex_array)
  {
#ifdef DEBUG
    vout.detailed("%s is called\n", __func__);
#endif

    typedef typename AFIELD1::complex_t   complex_t;
    typedef typename AFIELD1::real_t      real_t;

    const int num_vectors = testvectors.size();
    assert(coarse_vector.nin() == 4 * num_vectors);  // 4: 2 for complex, 2 for chirality

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, block_index.coarse_nvol());
    int coarse_Nx = block_index.coarse_lattice_size(0);
    int coarse_Ny = block_index.coarse_lattice_size(1);
    int coarse_Nz = block_index.coarse_lattice_size(2);
    int coarse_Nt = block_index.coarse_lattice_size(3);
    AIndex_coarse_lex<real_t, AFIELD1::IMPL> index(coarse_Nx, coarse_Ny, coarse_Nz, coarse_Nt, num_vectors, 2);

#pragma omp barrier

    // linear algebra on fine lattice

    fine_vector.set(0.0);

    for (int i = 0; i < num_vectors; i++) {
      copy(tmp2, testvectors[i]);
      afopr->mult_gm5(tmp1, tmp2);
      axpy(tmp2, -1.0, tmp1);            // tmp2=(1-gm5)v[i]
      axpy(tmp1, +1.0, testvectors[i]);  // tmp1=(1+gm5)v[i]

#pragma omp barrier

      // linear algebra with bolcks
      for (int block_id = is; block_id < ns; block_id++) { // ch=-1
        int    index_re = index.idx_SPr(i, 0, block_id, 0);
        int    index_im = index.idx_SPi(i, 0, block_id, 0);
        real_t re       = coarse_vector.cmp(index_re);
        real_t im       = coarse_vector.cmp(index_im);
        complex_array[block_id] = complex_t(re, im);
      }
      block_axpy(fine_vector, &(complex_array[0]), tmp2, real_t(0.5),
                 block_index);

      for (int block_id = is; block_id < ns; block_id++) { // ch=+1
        int    index_re = index.idx_SPr(i, 1, block_id, 0);
        int    index_im = index.idx_SPi(i, 1, block_id, 0);
        real_t re       = coarse_vector.cmp(index_re);
        real_t im       = coarse_vector.cmp(index_im);
        complex_array[block_id] = complex_t(re, im);
      }
      block_axpy(fine_vector, &(complex_array[0]), tmp1, real_t(0.5),
                 block_index);
    }//i
  }


  //====================================================================
  template<class AFIELD>
  void set_testvectors_impl_global(std::vector<AFIELD>& testvectors, AFopr<AFIELD> *afopr, Field& vec_tmp)
  {
    vout.general("using global random numbers for the initial testvectors\n");

    typedef typename AFIELD::real_t real_t;
    int nvec = testvectors.size();

    RandomNumbers *random = RandomNumberManager::getInstance();
    for (int i = 0; i < nvec; i++) {
      random->uniform_lex_global(vec_tmp);
      //#pragma omp parallel
      {
        if (afopr->needs_convert()) {
          if (i == 0) {
            vout.detailed("convert in rep. required.\n");
          }
          afopr->convert(testvectors[i], vec_tmp);
        } else {
          if (i == 0) {
            vout.detailed("convert in rep. not required.\n");
          }
          AIndex_lex<real_t, AFIELD::IMPL> index_f;
          convert(index_f, testvectors[i], vec_tmp);
        }
      } // omp parallel
#pragma omp barrier
    }   // i
  }


  //====================================================================
  template<class AFIELD>
  void set_testvectors_impl_global(std::vector<AFIELD>& testvectors, AFopr<AFIELD> *afopr)
  {
    std::string class_name = __func__;
    ThreadManager::assert_single_thread(class_name);

    int    nin  = testvectors[0].nin();
    size_t nvol = testvectors[0].nvol();
    Field  vec_tmp(nin, nvol, 1);
#pragma omp parallel
    {
      set_testvectors_impl_global(testvectors, afopr, vec_tmp);
    }
  }


  //====================================================================
  template<class AFIELD>
  void set_testvectors_impl_local(std::vector<AFIELD>& testvectors, AFopr<AFIELD> *afopr, Field& vec_tmp)
  {
    vout.general("using local random numbers for the initial testvectors\n");

    typedef typename AFIELD::real_t real_t;

    int    nvec = testvectors.size();
    int    nin  = testvectors[0].nin();
    size_t nvol = testvectors[0].nvol();
    assert(nin == vec_tmp.nin() && nvol == vec_tmp.nvol());
    {
      //      const int seed_base=13957039; // charged pion: 139.57039 MeV
      const int seed_base = 1395704; // charged pion: 139.57039 MeV
      //      const int seed_base=493677; // charged Kaon: 493.677 MeV
      //      const int seed_base=1234567;
      int myrank = Communicator::nodeid();
      int is0, ns0, ith, nth;
      set_threadtask(ith, nth, is0, ns0, nvol);
      int seed = seed_base + nth * myrank + ith;
      RandomNumbers_Mseries random(seed);
      size_t                is = is0 * nin;
      size_t                ns = ns0 * nin;
      for (int i = 0; i < nvec; i++) {
        double *pv = vec_tmp.ptr(0);
        for (size_t s = is; s < ns; s++) {
          pv[s] = random.get();
        }
#pragma omp barrier
        if (afopr->needs_convert()) {
          if (i == 0) {
            vout.detailed("convert in rep. required.\n");
          }
          afopr->convert(testvectors[i], vec_tmp);
        } else {
          if (i == 0) {
            vout.detailed("convert in rep. not required.\n");
          }
          AIndex_lex<real_t, AFIELD::IMPL> index_f;
          convert(index_f, testvectors[i], vec_tmp);
        }
#pragma omp barrier
      } // i
    }
  }


  //====================================================================
  template<class AFIELD>
  void set_testvectors_impl_local(std::vector<AFIELD>& testvectors, AFopr<AFIELD> *afopr)
  {
    std::string class_name = __func__;
    ThreadManager::assert_single_thread(class_name);

    int    nin  = testvectors[0].nin();
    size_t nvol = testvectors[0].nvol();
    Field  vec_tmp(nin, nvol, 1);

#pragma omp parallel
    {
      set_testvectors_impl_local(testvectors, afopr, vec_tmp);
    }
  }
} // end of namespace

//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::init_resources()
{
  m_tmp1.reset(m_nin, m_fine_nvol, 1);
  m_tmp2.reset(m_nin, m_fine_nvol, 1);
  m_field_tmp.reset(m_nin, m_fine_nvol, 1);

  int coarse_nvol = m_block_index.coarse_nvol();
  m_real_array.resize(coarse_nvol);
  m_complex_array.resize(coarse_nvol);
  m_coarse_array.resize(coarse_nvol * 2 * m_nvec);
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::set_testvectors(const std::vector<AFIELD2>& w)
{
  assert(w.size() == m_nvec);
  for (int ivec = 0; ivec < m_nvec; ++ivec) {
    copy(m_testvectors[ivec], w[ivec]);
  }
}


template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::set_testvectors()
{
  Timer timer;
  timer.start();

  //set_testvectors_impl_global(m_testvectors, m_afopr_fine);
  set_testvectors_impl_local(m_testvectors, m_afopr_fine, m_field_tmp);

  timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.detailed(m_vl, "%s: time for set_testvector = %e sec\n",
                class_name.c_str(), elapsed_time);
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::set_coarse_array(
  const AFIELD1& coarse_vector) const

{
  typedef typename AFIELD1::real_t    real_t;
  typedef typename AFIELD1::complex_t complex_t;

  const int coarse_nvol = m_block_index.coarse_nvol();
  const int num_vectors = m_testvectors.size();

  const int Nxc = m_block_index.coarse_lattice_size(0);
  const int Nyc = m_block_index.coarse_lattice_size(1);
  const int Nzc = m_block_index.coarse_lattice_size(2);
  const int Ntc = m_block_index.coarse_lattice_size(3);

  AIndex_coarse_lex<real_t, AFIELD1::IMPL> index_c(
    Nxc, Nyc, Nzc, Ntc, num_vectors, 2);

  int ith, nth, is_coarse, ns_coarse;
  set_threadtask(ith, nth, is_coarse, ns_coarse, coarse_nvol);

#pragma omp barrier

  for (int ivec = 0; ivec < num_vectors; ++ivec) {
    for (int block = is_coarse; block < ns_coarse; ++block) {
      int    index_re1 = index_c.idx_SPr(ivec, 0, block, 0);
      int    index_im1 = index_c.idx_SPi(ivec, 0, block, 0);
      real_t re1       = coarse_vector.cmp(index_re1);
      real_t im1       = coarse_vector.cmp(index_im1);

      int    index_re2 = index_c.idx_SPr(ivec, 1, block, 0);
      int    index_im2 = index_c.idx_SPi(ivec, 1, block, 0);
      real_t re2       = coarse_vector.cmp(index_re2);
      real_t im2       = coarse_vector.cmp(index_im2);

      // impl-2
      //  |f> += cm (1-gm5)/2 |i> + cp (1+gm5)/2 |i>
      //       = (cm+cp) |i>/2  + (-cm+cp) gm5|i>/2
      //   cm=(re1,im1), cp=(re2,im2)
      m_coarse_array[block + coarse_nvol * (2 * ivec)]
        = complex_t(re2 + re1, im2 + im1);
      m_coarse_array[block + coarse_nvol * (2 * ivec + 1)]
        = complex_t(re2 - re1, im2 - im1);

      /*
      // impl-1
      m_coarse_array[block + coarse_nvol * (2*ivec)]
                            = complex_t(re1, im1);
      m_coarse_array[block + coarse_nvol * (2*ivec+1)]
                            = complex_t(re2, im2);
      */
    }
  }

#pragma omp barrier
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::set_coarse_vector(
  AFIELD1& coarse_vector) const
{
  typedef typename AFIELD1::real_t    real_t;
  typedef typename AFIELD1::complex_t complex_t;

  const int coarse_nvol = m_block_index.coarse_nvol();
  const int num_vectors = m_testvectors.size();

  const int Nxc = m_block_index.coarse_lattice_size(0);
  const int Nyc = m_block_index.coarse_lattice_size(1);
  const int Nzc = m_block_index.coarse_lattice_size(2);
  const int Ntc = m_block_index.coarse_lattice_size(3);

  AIndex_coarse_lex<real_t, AFIELD1::IMPL> index_c(
    Nxc, Nyc, Nzc, Ntc, num_vectors, 2);

  complex_t *array1 = &m_coarse_array[0];

  int ith, nth, is_coarse, ns_coarse;
  set_threadtask(ith, nth, is_coarse, ns_coarse, coarse_nvol);

#pragma omp barrier

  for (int ivec = 0; ivec < num_vectors; ++ivec) {
    for (int block = is_coarse; block < ns_coarse; ++block) {
      real_t re1       = real(array1[block + coarse_nvol * (2 * ivec)]);
      real_t im1       = imag(array1[block + coarse_nvol * (2 * ivec)]);
      int    index_re1 = index_c.idx_SPr(ivec, 0, block, 0);
      int    index_im1 = index_c.idx_SPi(ivec, 0, block, 0);
      coarse_vector.set(index_re1, re1);
      coarse_vector.set(index_im1, im1);

      real_t re2       = real(array1[block + coarse_nvol * (2 * ivec + 1)]);
      real_t im2       = imag(array1[block + coarse_nvol * (2 * ivec + 1)]);
      int    index_re2 = index_c.idx_SPr(ivec, 1, block, 0);
      int    index_im2 = index_c.idx_SPi(ivec, 1, block, 0);
      coarse_vector.set(index_re2, re2);
      coarse_vector.set(index_im2, im2);
    }
  }
#pragma omp barrier
  coarse_vector.scal(real_t(0.5));
#pragma omp barrier
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::make_fine_vector(AFIELD2& fine_vector, const AFIELD1& coarse_vector) const
{
  vout.general("%s: make_fine_vector is called\n", class_name.c_str());

  // fine_vec
  // = coarse[ic, ch=-1, block] * P_-  testvecr[ic, block]
  //    + coarse[ic, ch=+1, block]*  P_+ testvec[ic, block]
  // = (coarse[ic, ch=-1, block] + coarse[ic, ch=+1, block])
  //            *testvec[ic, block] / 2
  //    + (- coarse[ic, ch=-1, block] + coarse[ic, ch=+1, block])
  //            * gamma5 * testvec[ic, block] /2

  //  unique_ptr<Timer> timer(new Timer);
  //  timer->start();

  //  make_fine_vector_impl(fine_vector, coarse_vector,
  //                        m_testvectors, m_afopr_fine, m_block_index,
  //                        m_tmp1, m_tmp2,  &(m_complex_array[0]));

  typedef typename AFIELD1::complex_t complex_t;
  typedef typename AFIELD1::real_t    real_t;

  const int num_vectors = m_testvectors.size();
  assert(coarse_vector.nin() == 4 * num_vectors);
  // 4: 2 for complex * 2 for chirality

  const int coarse_nvol = m_block_index.coarse_nvol();

  set_coarse_array(coarse_vector);

#pragma omp barrier

  fine_vector.set(0.0);

  for (int ivec = 0; ivec < num_vectors; ++ivec) {
    // imple-2
    m_afopr_fine->mult_gm5(m_tmp1, m_testvectors[ivec]);

    // impl-1
    //copy(m_tmp2, m_testvectors[ivec]);
    //m_afopr_fine->mult_gm5(m_tmp1, m_tmp2);
    //axpy(m_tmp2, -1.0, m_tmp1);               // tmp2 = (1 - gm5) v[i]
    //axpy(m_tmp1, +1.0, m_testvectors[ivec]);  // tmp1 = (1 + gm5) v[i]

#pragma omp barrier

    complex_t *array1 = &m_coarse_array[coarse_nvol * (2 * ivec)];

    // impl-2
    block_axpy(fine_vector, array1, m_testvectors[ivec], real_t(0.5),
               m_block_index);
    // impl-1
    //block_axpy(fine_vector, array1, m_tmp2, real_t(0.5),
    //                                                  m_block_index);

    complex_t *array2 = &m_coarse_array[coarse_nvol * (2 * ivec + 1)];
    block_axpy(fine_vector, array2, m_tmp1, real_t(0.5),
               m_block_index);
  }

  //  timer->stop();
  //  double elapsed_time = timer->elapsed_sec();
  //  vout.general("  time for make_fine_vector = %e sec\n", elapsed_time);
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::make_coarse_vector(AFIELD1& coarse_vector, const AFIELD2& fine_vector) const
{
  vout.general("%s: make_coarse_vector is called\n", class_name.c_str());

  typedef typename AFIELD1::complex_t complex_t;
  typedef typename AFIELD1::real_t    real_t;

  //  unique_ptr<Timer> timer(new Timer);
  //  timer->start();

  //  make_coarse_vector_impl(coarse_vector, fine_vector,
  //                          m_testvectors, m_afopr_fine, m_block_index,
  //                          m_tmp1, m_tmp2,  &(m_complex_array[0]));


  const int num_vectors = m_testvectors.size();
  assert(coarse_vector.nin() == 4 * num_vectors);
  // 4: 2 for complex, 2 for chirality

  const int coarse_nvol = m_block_index.coarse_nvol();

  for (int ivec = 0; ivec < num_vectors; ++ivec) {
    copy(m_tmp1, m_testvectors[ivec]);
    m_afopr_fine->mult_gm5(m_tmp2, m_testvectors[ivec]);
#pragma omp barrier
    axpy(m_tmp1, -1.0, m_tmp2);              // tmp1 = (1 - gm5) v[i]
    axpy(m_tmp2, +1.0, m_testvectors[ivec]); // tmp2 = (1 + gm5) v[i]
#pragma omp barrier

    // ch=-1
    complex_t *array1 = &m_coarse_array[coarse_nvol * (2 * ivec)];
    block_dotc(array1, m_tmp1, fine_vector, m_block_index);

    // ch=+1
    complex_t *array2 = &m_coarse_array[coarse_nvol * (2 * ivec + 1)];
    block_dotc(array2, m_tmp2, fine_vector, m_block_index);

#pragma omp barrier
  } // ivec

  set_coarse_vector(coarse_vector);

  //  timer->stop();
  //  double elapsed_time = timer->elapsed_sec();
  //  vout.general("  time for make_coarse_vector = %e sec\n", elapsed_time);
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::gramschmidt()
{
  gramschmidt(m_testvectors);
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::gramschmidt(
  std::vector<AFIELD2>& vectors) const
{
  vout.general("%s: ::gramschmidt is called\n", class_name.c_str());

  gramschmidt_in_block_impl(vectors, m_afopr_fine, m_block_index,
                            m_tmp1, m_tmp2,
                            &(m_real_array[0]), &(m_complex_array[0]));

  gramschmidt_in_block_impl(vectors, m_afopr_fine, m_block_index,
                            m_tmp1, m_tmp2,
                            &(m_real_array[0]), &(m_complex_array[0]));
}


//====================================================================
template<class AFIELD1, class AFIELD2>
void MultiGrid_Clover<AFIELD1, AFIELD2>::gramschmidt_double(
  std::vector<AFIELD2>& vectors) const
{
  vout.crucial("%s: gramschmed_double is called\n", class_name.c_str());
  abort();
}


//============================================================END=====
