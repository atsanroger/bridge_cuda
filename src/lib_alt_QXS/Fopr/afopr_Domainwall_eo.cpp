/*!
      @file    afopr_Domainwall_eo.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2514 $
*/

#include "lib/Fopr/afopr_Domainwall_eo.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/aindex_eo-inc.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, QXS> >::Factory_params::Register(
    "Domainwall_eo", create_object_with_params);

  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Domainwall_eo", create_object_with_params);
}
#endif


template<>
class Index_eo_Domainwall<AField<double, QXS> > {
  typedef AField<double, QXS> AFIELD;
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
class Index_eo_Domainwall<AField<float, QXS> > {
  typedef AField<float, QXS> AFIELD;
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
const std::string AFopr_Domainwall_eo<AField<double, QXS> >::
class_name = "AFopr_Domainwall_eo<AField<double,QXS> >";

template<>
const std::string AFopr_Domainwall_eo<AField<float, QXS> >::
class_name = "AFopr_Domainwall_eo<AField<float,QXS> >";


#include "lib/Fopr/afopr_Domainwall_eo-tmpl.h"

// class instanciation.
template class AFopr_Domainwall_eo<AField<float, QXS> >;
template class AFopr_Domainwall_eo<AField<double, QXS> >;

//============================================================END=====
