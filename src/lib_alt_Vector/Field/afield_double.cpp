/*!
      @file    afield_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/aindex_lex.h"

typedef double real_t;

#include "lib_alt_Vector/inline/afield_th-inc.h"


#include "lib_alt_Vector/BridgeVec/bridgeVec_afield.h"

#include "lib_alt_Vector/Field/afield-tmpl.h"


template<>
const std::string AField<double, VECTOR>::class_name = "AField<double, VECTOR>";


// explicit instanciation.
template class AField<double, VECTOR>;

// copy between fields with different precisions.
//====================================================================
void copy(AField<double, VECTOR>& v, const AField<float, VECTOR> &w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  float *wp = const_cast<AField<float, VECTOR>*>(&w)->ptr(0);
  //  copy_from_device(wp, w.size_pad());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, v.size_pad());

  for (int i = is; i < ns; ++i) {
    v.set_host(i, double(w.cmp(i)) );
  }

  //  copy_to_device(v.ptr(0), v.size_pad());

}

//============================================================END=====
