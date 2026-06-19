/*!
      @file    afopr_Domainwall_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2533 $
*/

#include "lib/Fopr/afopr_Domainwall_eo.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_eo.h"
#include "lib_alt_Accel/Field/aindex_eo-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,ACCEL> >::Factory_params::Register(
                      "Domainwall_eo", create_object_with_params);

  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                      "Domainwall_eo", create_object_with_params);
}
#endif


template<>
class Index_eo_Domainwall<AField<double, ACCEL> > {
  typedef AField<double, ACCEL> AFIELD;
public:
  void split(AFIELD &xe, AFIELD &xo, const AFIELD &x){
    m_idx.split(xe, xo, x);
  }
  void merge(AFIELD &x, const AFIELD &xe, const AFIELD &xo){
    m_idx.merge(x, xe, xo);
  }
private:
  AIndex_eo<AFIELD::real_t, AFIELD::IMPL>  m_idx;
};

template<>
class Index_eo_Domainwall<AField<float, ACCEL> > {
  typedef AField<float, ACCEL> AFIELD;
public:
  void split(AFIELD &xe, AFIELD &xo, const AFIELD &x){
    m_idx.split(xe, xo, x);
  }
  void merge(AFIELD &x, const AFIELD &xe, const AFIELD &xo){
    m_idx.merge(x, xe, xo);
  }
private:
  AIndex_eo<AFIELD::real_t, AFIELD::IMPL>  m_idx;
};


template<>
const std::string AFopr_Domainwall_eo<AField<double,ACCEL> >::
class_name = "AFopr_Domainwall_eo<AField<double> >";

template<>
const std::string AFopr_Domainwall_eo<AField<float,ACCEL> >::
class_name = "AFopr_Domainwall_eo<AField<float> >";


#include "lib/Fopr/afopr_Domainwall_eo-tmpl.h"


// class instanciation.
template class AFopr_Domainwall_eo<AField<double,ACCEL> >;
template class AFopr_Domainwall_eo<AField<float,ACCEL> >;

//============================================================END=====
