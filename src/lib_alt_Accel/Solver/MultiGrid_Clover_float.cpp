/*!
        @file    MultiGrid_Clover_float.cpp
        @brief   MultiGrid operation for Clover fermion (Accel version)
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2258 $
 */

//====================================================================

#include "lib_alt_Accel/inline/define_params.h"

typedef float real_t;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_dd-inc.h"
#include "lib_alt_Accel/Fopr/afopr_Clover.h"
#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"
#include "lib_alt_Accel/Field/aindex_coarse_lex.h"
#include "lib_alt_Accel/Field/aindex_block_lex.h"

// template for MultiGrid_Clover
#include "lib_alt/Solver/MultiGrid_Clover.h"
#include "lib_alt/Solver/MultiGrid_Clover-tmpl.h"

typedef AField<float, ACCEL> AField_f;


template<>
const std::string MultiGrid_Clover< AField_f, AField_f >::class_name
 = "MultiGrid_Clover< AField<float,ACCEL>,  AField<float,ACCEL> >";

//====================================================================

template<>
void MultiGrid_Clover< AField_f, AField_f >::set_coarse_array(
                                 const AField_f& coarse_vector) const
{
  typedef typename AField_f::real_t real_t;
  typedef typename AField_f::complex_t complex_t;

  real_t* wp = const_cast<AField_f*>(&coarse_vector)->ptr(0);

  const int coarse_nvol = m_block_index.coarse_nvol();
  const int num_vectors = m_testvectors.size();

  coarse_vector.update_host();

  int ith, nth, is_coarse, ns_coarse;
  set_threadtask(ith, nth, is_coarse, ns_coarse, coarse_nvol);

#pragma omp barrier

  int Nin = 2 * 2 * num_vectors;

  for(int ivec = 0; ivec < num_vectors; ++ivec){
    for(int block = is_coarse; block < ns_coarse; ++block){

      real_t re1 = wp[IDX2(Nin, 2*(2*ivec)  , block)];
      real_t im1 = wp[IDX2(Nin, 2*(2*ivec)+1, block)];

      real_t re2 = wp[IDX2(Nin, 2*(2*ivec+1)  , block)];
      real_t im2 = wp[IDX2(Nin, 2*(2*ivec+1)+1, block)];

      // impl-2
      m_coarse_array[block + coarse_nvol * (2*ivec)]
                          = complex_t(re2 + re1, im2 + im1);
      m_coarse_array[block + coarse_nvol * (2*ivec+1)]
                          = complex_t(re2 - re1, im2 - im1);
      // impl-1
      //m_coarse_array[block + coarse_nvol * (2*ivec)]
      //                        = complex_t(re1, im1);
      //m_coarse_array[block + coarse_nvol * (2*ivec+1)]
      //                        = complex_t(re2, im2);
    }
  }

#pragma omp barrier

}

//====================================================================

template<>
void MultiGrid_Clover< AField_f, AField_f >::set_coarse_vector(
                                       AField_f& coarse_vector) const
{
  typedef typename AField_f::real_t real_t;
  typedef typename AField_f::complex_t complex_t;

  real_t* wp = coarse_vector.ptr(0);

  const int coarse_nvol = m_block_index.coarse_nvol();
  const int num_vectors = m_testvectors.size();

  int ith, nth, is_coarse, ns_coarse;
  set_threadtask(ith, nth, is_coarse, ns_coarse, coarse_nvol);

  complex_t *array1 = &m_coarse_array[0];

#pragma omp barrier

  int Nin = 2 * 2 * num_vectors;

  for(int ivec = 0; ivec < num_vectors; ++ivec){
    for(int block = is_coarse; block < ns_coarse; ++block){

      real_t re1 = real(array1[block + coarse_nvol * (2*ivec)]);
      real_t im1 = imag(array1[block + coarse_nvol * (2*ivec)]);
      wp[IDX2(Nin, 2*(2*ivec)  , block)] = 0.5 * re1;
      wp[IDX2(Nin, 2*(2*ivec)+1, block)] = 0.5 * im1;

      real_t re2 = real(array1[block + coarse_nvol * (2*ivec+1)]);
      real_t im2 = imag(array1[block + coarse_nvol * (2*ivec+1)]);
      wp[IDX2(Nin, 2*(2*ivec+1)  , block)] = 0.5 * re2;
      wp[IDX2(Nin, 2*(2*ivec+1)+1, block)] = 0.5 * im2;

    }
  }

#pragma omp barrier

  coarse_vector.update_device();

}

//====================================================================
// class instantiation

template class MultiGrid_Clover< AField_f, AField_f >;

//============================================================END=====
