/*!
      @file    index_block_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AINDEX_BLOCK_LEX_INCLUDED
#define QXS_AINDEX_BLOCK_LEX_INCLUDED

#include <string>
#include <vector>

#include "lib/Parameters/commonParameters.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"
#include "lib_alt/Field/aindex_block_lex_base.h"


//! Blocked site index.

/*!
  This class defines lexicographycal site index for alternative
  code set: QXS version.
                                       [10 Mar 2021 H.Matsufuru]
*/
template<typename REALTYPE>
class AIndex_block_lex<REALTYPE, QXS>
{
 private:
  int m_coarse_lattice[4];
  int m_fine_lattice[4];
  int m_block_size[4];
  std::vector<int> m_block_eo;
  size_t m_coarse_nvol;
  size_t m_fine_nvol;
  int m_block_nvol;

 public:
  AIndex_block_lex() { /* do nothing */ }

  AIndex_block_lex(const std::vector<int>& coarse_lattice)
  {
    std::vector<int> fine_lattice(4);
    fine_lattice[0] = CommonParameters::Nx();
    fine_lattice[1] = CommonParameters::Ny();
    fine_lattice[2] = CommonParameters::Nz();
    fine_lattice[3] = CommonParameters::Nt();
    init(coarse_lattice, fine_lattice);
  }

  AIndex_block_lex(const std::vector<int>& coarse_lattice,
                   const std::vector<int>& fine_lattice)
  { init(coarse_lattice, fine_lattice); }

 public:

  int coarse_lattice_size(const int mu) const
  { return m_coarse_lattice[mu]; }

  int fine_lattice_size(const int mu) const
  { return m_fine_lattice[mu]; }

  int block_size(const int mu) const
  { return m_block_size[mu]; }

  size_t fine_nvol() const { return m_fine_nvol; }

  size_t coarse_nvol() const { return m_coarse_nvol; }

  size_t block_nvol() const { return m_block_nvol; }

  int block_eo(const int block_idx) const
  { return m_block_eo[block_idx]; }

  void get_coarse_coord(int& ibx, int& iby, int& ibz, int& ibt,
                        const int block_idx) const
  {
    int idx_tmp = block_idx;
    ibx     = idx_tmp % m_coarse_lattice[0];
    idx_tmp = idx_tmp / m_coarse_lattice[0];
    iby     = idx_tmp % m_coarse_lattice[1];
    idx_tmp = idx_tmp / m_coarse_lattice[1];
    ibz     = idx_tmp % m_coarse_lattice[2];
    ibt     = idx_tmp / m_coarse_lattice[2];
  }

  // private:

  void init(const std::vector<int>& coarse_lattice,
            const std::vector<int>& fine_lattice)
  {
    m_coarse_nvol = 1;
    m_fine_nvol   = 1;
    m_block_nvol  = 1;
    for (int mu = 0; mu < 4; ++mu) {
      m_coarse_lattice[mu] = coarse_lattice[mu];
      m_fine_lattice[mu]   = fine_lattice[mu];
      if (m_fine_lattice[mu] % m_coarse_lattice[mu] != 0) {
        vout.crucial("AIndex_block_lex: fine_lattice is");
        vout.crucial(" not multiple of coarse lattice\n");
        exit(EXIT_FAILURE);
      }
      m_block_size[mu] = fine_lattice[mu] / coarse_lattice[mu];
      m_coarse_nvol   *= m_coarse_lattice[mu];
      m_fine_nvol     *= m_fine_lattice[mu];
      m_block_nvol    *= m_block_size[mu];
    }
    // sanity check
    {
      int vlen = 1;
      if (sizeof(REALTYPE) == 4) {
        vlen = VLENS;
      } else if (sizeof(REALTYPE) == 8) {
        vlen = VLEND;
      }  else {
        vout.crucial("index_block_lex_alt: unkown REALTYPE\n");
        exit(EXIT_FAILURE);
      }
      if (m_coarse_nvol % vlen != 0) {
        vout.crucial("index_block_lex_alt: bad coarse_novl  (too small?), must be multiple of VLEN\n");
        vout.crucial("    coarse_nvol=%d, coarse_lattice = %d %d %d %d\n", m_coarse_nvol, m_coarse_lattice[0], m_coarse_lattice[1], m_coarse_lattice[2], m_coarse_lattice[3]);
        exit(EXIT_FAILURE);
      }
    }
    m_block_eo.resize(m_coarse_nvol);
    int NBx  = m_coarse_lattice[0];
    int NBy  = m_coarse_lattice[1];
    int NBz  = m_coarse_lattice[2];
    int NBt  = m_coarse_lattice[3];
    int ipex = Communicator::ipe(0);
    int ipey = Communicator::ipe(1);
    int ipez = Communicator::ipe(2);
    int ipet = Communicator::ipe(3);
    int Ieo  = (NBx * ipex + NBy * ipey + NBz * ipez + NBt * ipet) % 2;
    for (int block_idx = 0; block_idx < m_coarse_nvol; ++block_idx) {
      int kx, ky, kz, kt;
      get_coarse_coord(kx, ky, kz, kt, block_idx);
      m_block_eo[block_idx] = (Ieo + kx + ky + kz + kt) % 2;
    }
  }
};

#endif
