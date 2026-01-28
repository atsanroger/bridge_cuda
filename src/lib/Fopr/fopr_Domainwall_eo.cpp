/*!
        @file    fopr_Domainwall_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#include "Fopr/fopr_Domainwall_eo.h"
#include "Field/index_eo.h"

template<>
class Index_eo_Domainwall<Field> {
public:
  void split(Field &xe, Field &xo, const Field &x){
    m_idx.convertField(xe, x, 0);
    m_idx.convertField(xo, x, 1);
  }
  void merge(Field &x, const Field &xe, const Field &xo){
    m_idx.reverseField(x, xe, 0);
    m_idx.reverseField(x, xo, 1);
  }
private:
  Index_eo m_idx;
};


#include "Fopr/afopr_Domainwall_eo-tmpl.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Domainwall::register_factory();
}
#endif
template<>
const std::string AFopr_Domainwall_eo<Field>::class_name
  = "Fopr_Domainwall_eo";

template class AFopr_Domainwall_eo<Field>;


//============================================================END=====
