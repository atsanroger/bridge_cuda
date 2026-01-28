/*!
        @file    afield_copy.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-04-04 12:51:46 #$
        @version $LastChangedRevision: 2597 $
*/

#include "lib_alt_QXS/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"

//====================================================================
void copy(AField<double, QXS>& v, const AField<float, QXS>& w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  AIndex_lex<double, QXS> index_d;
  AIndex_lex<float,  QXS> index_f;
  convert(index_d, v, index_f, w);
}

//====================================================================
void copy(AField<float, QXS>& v, const AField<double, QXS>& w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  AIndex_lex<float,  QXS> index_f;
  AIndex_lex<double, QXS> index_d;
  convert(index_f, v, index_d, w);
}

//============================================================END=====
