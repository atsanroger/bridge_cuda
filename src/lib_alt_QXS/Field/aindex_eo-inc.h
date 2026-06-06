/*!
      @file    aindex_eo-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
      @version $LastChangedRevision: 2529 $
*/

#include <assert.h>

#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

//====================================================================
template<typename REALTYPE>
template<typename AFIELD>
void AIndex_eo<REALTYPE, QXS>::split(AFIELD& field_e,
                                     AFIELD& field_o,
                                     const AFIELD& field_lex,
                                     bool /*qdw*/)
{
#pragma omp barrier

  int Nin   = field_lex.nin();
  int Nex   = field_lex.nex();
  int Nvol  = field_lex.nvol();
  int Nvol2 = Nvol / 2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  AIndex_lex<REALTYPE, QXS> index_lex;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int ist = is; ist < ns; ++ist) {
      int ix   = ist % Nx;
      int iyzt = ist / Nx;
      int ist2 = ist / 2;
      int ieo  = (ix + Leo[iyzt]) % 2;
      if (ieo == 0) {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_e.set(index2, field_lex.cmp(index1));
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_o.set(index2, field_lex.cmp(index1));
        }
      }
    }
  }
#pragma omp barrier

}


//====================================================================
template<typename REALTYPE>
template<typename AFIELD>
void AIndex_eo<REALTYPE, QXS>::merge(AFIELD& field_lex,
                                     const AFIELD& field_e,
                                     const AFIELD& field_o,
                                     bool /*qdw*/)
{
#pragma omp barrier

  int Nin   = field_lex.nin();
  int Nex   = field_lex.nex();
  int Nvol  = field_lex.nvol();
  int Nvol2 = Nvol / 2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  AIndex_lex<REALTYPE, QXS> index_lex;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int ist = is; ist < ns; ++ist) {
      int ix   = ist % Nx;
      int iyzt = ist / Nx;
      int ist2 = ist / 2;
      int ieo  = (ix + Leo[iyzt]) % 2;
      if (ieo == 0) {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set(index1, field_e.cmp(index2));
        }
      } else {
        for (int in = 0; in < Nin; ++in) {
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set(index1, field_o.cmp(index2));
        }
      }
    }
  }
#pragma omp barrier

}


//============================================================END=====
