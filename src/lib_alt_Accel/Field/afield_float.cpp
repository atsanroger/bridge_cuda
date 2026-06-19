/*!
      @file    afield_float.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#include "lib_alt_Accel/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/aindex_lex.h"

typedef float real_t;

#include "lib_alt_Accel/inline/afield_th-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_Gauge.h"


// class template file
#include "lib_alt_Accel/Field/afield-tmpl.h"

template<>
const std::string AField<float, ACCEL>::class_name = "AField<float, ACCEL>";

// explicit instanciation.
template class AField<float, ACCEL>;

// copy between fields with different precisions.
//====================================================================
void copy(AField<float, ACCEL>& v, const AField<double, ACCEL> &w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  double *wp = const_cast<AField<double, ACCEL>*>(&w)->ptr(0);
  BridgeACC::copy_from_device(wp, w.size_pad());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, v.size_pad());

  for (int i = is; i < ns; ++i) {
    v.set_host(i, float(w.cmp_host(i)) );   // w already bulk-synced above;
  }                                         // cmp_host avoids a per-element D2H

  BridgeACC::copy_to_device(v.ptr(0), v.size_pad());

}

//============================================================END=====
