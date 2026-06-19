/*!

        @file    $Id: MultiGrid_Clover.h #$

        @brief   MultiGrid operation for Clover fermion (template version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================
#ifndef MULTIGRID_CLOVER_H_INCLUDED
#define MULTIGRID_CLOVER_H_INCLUDED

#include "MultiGrid.h"


//====================================================================
template<class AFIELD1, class AFIELD2>
class MultiGrid_Clover : public MultiGrid<AFIELD1, AFIELD2>
{
 public:
  typedef AFIELD1                                                     AField_coarse_t;
  typedef AFIELD2                                                     AField_fine_t;
  typedef AIndex_block_lex<typename AFIELD2::real_t, AFIELD2::IMPL>   Index_t;
  //  typedef AFopr_Clover<AFIELD2> Afopr_fine_t;
  typedef AFopr_Clover_dd<AFIELD2>                                    Afopr_fine_dd_t;
  typedef void                                                        Afopr_coarse_t;
  static const std::string class_name;

 protected:
  using MultiGrid<AFIELD1, AFIELD2>::m_vl;
  //  AFopr<AFIELD1> *m_afopr_coarse; // only the fine op. is needed
  AFopr<AFIELD2> *m_afopr_fine = nullptr;
  Index_t m_block_index;
  std::vector<AFIELD2> m_testvectors;
  mutable AFIELD2 m_tmp1, m_tmp2;
  Field m_field_tmp;

  int m_nin;
  int m_nvec;
  size_t m_fine_nvol;
  mutable std::vector<typename AFIELD2::real_t> m_real_array;
  mutable std::vector<typename AFIELD2::complex_t> m_complex_array;
  mutable std::vector<typename AFIELD2::complex_t> m_complex_array2;

  mutable std::vector<typename AFIELD2::complex_t> m_coarse_array;

 public:
  MultiGrid_Clover() {}
  MultiGrid_Clover(const std::vector<int>& coarse_lattice,
                   const std::vector<int>& fine_lattice,
                   const int nin, const int nvec)
  {
    init(coarse_lattice, fine_lattice, nin, nvec);
  }

  ~MultiGrid_Clover()
  {
    tidyup();
  }

  void init(const std::vector<int>& coarse_lattice, const std::vector<int>& fine_lattice,
            const int nin, const int nvec)
  {
    m_block_index.init(coarse_lattice, fine_lattice);
    size_t fine_nvol = 1;
    for (int i = 0; i < fine_lattice.size(); i++) {
      fine_nvol *= fine_lattice[i];
    }
    m_testvectors.resize(nvec);
    for (int i = 0; i < nvec; i++) {
      m_testvectors[i].reset(nin, fine_nvol, 1);
    }

    m_fine_nvol = fine_nvol;
    m_nin       = nin;
    m_nvec      = nvec;
    init_resources();
  }

  void init_resources();

  void tidyup() {}

  const Index_t *get_block_index() const
  {
    return &m_block_index;
  }

  std::vector<AFIELD2> *get_testvectors()
  {
    return &m_testvectors;
  }

  const std::vector<AFIELD2> *get_testvectors() const
  {
    return &m_testvectors;
  }

  void set_testvectors();

  void set_testvectors(const std::vector<AFIELD2>&);

  //! array <- vector (function for optimization)
  void set_coarse_array(const AFIELD1& coarse_vector) const;

  //! vector <- array (function for optimization)
  void set_coarse_vector(AFIELD1& coarse_vector) const;

  void set_afopr_fine(AFopr_dd<AFIELD2> *afopr)
  {
    // check if downcasting is safe
    { // Clover_dd
      Afopr_fine_dd_t *afopr_tmp = dynamic_cast<Afopr_fine_dd_t *>(afopr);
      if (afopr_tmp != nullptr) {
        m_afopr_fine = afopr_tmp;
        return;
      }
    }
    vout.crucial("MultiGrid_Clover: bad afopr: only AFopr_Clover_dd is vaild\n");
    exit(EXIT_FAILURE);
  }

  void set_afopr_coarse(AFopr<AFIELD1> *afopr)
  {
    // do nothing
  }

  void gramschmidt();
  void gramschmidt(std::vector<AFIELD2>& fine_vectors) const;
  void gramschmidt_double(std::vector<AFIELD2>& fine_vectors) const;

  void make_fine_vector(AFIELD2& fine_vector, const AFIELD1& coarse_vector) const;
  void make_coarse_vector(AFIELD1& coarse_vector, const AFIELD2& fine_vector) const;
};

#endif // MULTIGRID_CLOVER_H_INCLUDED
//============================================================END=====
