/*!

        @file    $Id: MultiGrid.h #$

        @brief   base class for MultiGrid

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================
#pragma once

#include "lib_alt/Field/aindex_block_lex_base.h"
#include "lib/Fopr/afopr.h"


template<class AFIELD1, class AFIELD2>
class MultiGrid
{
 public:
  typedef AFIELD1                                                     Afield_coarse_t;
  typedef AFIELD2                                                     Afield_fine_t;
  typedef AIndex_block_lex<typename AFIELD2::real_t, AFIELD2::IMPL>   Index_t;

 protected:
  static const std::string class_name;
  Bridge::VerboseLevel m_vl;

 public:
  MultiGrid()
    : m_vl(CommonParameters::Vlevel()) {}

  MultiGrid(const std::vector<int>& coarse_lattice,
            const std::vector<int>& fine_lattice,
            const int nin  = 0,
            const int nvec = 0)
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~MultiGrid() {}

 private:
  MultiGrid(const MultiGrid&);
  MultiGrid& operator=(const MultiGrid&);

 public:
  virtual void init(const std::vector<int>& coarse_lattice, const std::vector<int>& fine_lattice, const int nin, const int nvec) = 0;

  virtual void set_afopr_coarse(AFopr<AFIELD1> *afopr) {}
  virtual void set_afopr_fine(AFopr<AFIELD2> *afopr) {}

  virtual const Index_t *get_block_index() const = 0;

  virtual std::vector<AFIELD2> *get_testvectors()             = 0;
  virtual const std::vector<AFIELD2> *get_testvectors() const = 0;

  virtual void set_testvectors() = 0;
  virtual void set_testvectors(const std::vector<AFIELD2>&) = 0;

  virtual void gramschmidt() = 0;
  virtual void gramschmidt(std::vector<AFIELD2>& fine_vectors) const = 0;

  virtual void make_fine_vector(AFIELD2&, const AFIELD1&) const   = 0;
  virtual void make_coarse_vector(AFIELD1&, const AFIELD2&) const = 0;
};


//============================================================END=====
