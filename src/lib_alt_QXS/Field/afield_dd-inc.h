/*!
        @file    afield_dd-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-05-07 13:24:28 #$
        @version $LastChangedRevision: 2602 $
*/

#ifndef QXS_AFIELD_DD_INC_INCLUDED
#define QXS_AFIELD_DD_INC_INCLUDED

#include <cstdlib>

#include "complexTraits.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/afield_th-inc.h"


//====================================================================

// declaration of implementation functions

template<typename INDEX, typename AFIELD>
void block_dotc_eo_impl(typename AFIELD::complex_t *out,
                        const AFIELD& v, const AFIELD& w,
                        const int ieo, const INDEX& block_index);

template<typename INDEX, typename AFIELD>
void block_norm2_eo_impl(typename AFIELD::real_t *out,
                         const AFIELD& v, const int ieo,
                         const INDEX& block_index);

template<typename INDEX, typename AFIELD>
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
                        const int ieo, const INDEX& block_index);

template<typename INDEX, typename AFIELD>
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
                        const int ieo, const INDEX& block_index);

template<typename INDEX, typename AFIELD>
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
                        const INDEX& block_index);

template<typename INDEX, typename AFIELD>
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
                        const INDEX& block_index);

//====================================================================
template<typename INDEX, typename AFIELD>
void block_dotc(typename AFIELD::complex_t *out, const AFIELD& v,
                const AFIELD& w, const INDEX& block_index)
{
  block_dotc_eo_impl(out, v, w, -1, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_dotc_eo(typename AFIELD::complex_t *out,
                   const AFIELD& v, const AFIELD& w,
                   const int ieo, const INDEX& block_index)
{
  block_dotc_eo_impl(out, v, w, ieo, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_norm2(typename AFIELD::real_t *out, const AFIELD& v,
                 const INDEX& block_index)
{
  block_norm2_eo_impl(out, v, -1, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_norm2_eo(typename AFIELD::real_t *out,
                    const AFIELD& v, const int ieo,
                    const INDEX& block_index)
{
  block_norm2_eo_impl(out, v, ieo, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal(AFIELD& v, const typename AFIELD::real_t *a,
                const INDEX& block_index)
{
  block_scal_eo_impl(v, a, -1, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal_eo(AFIELD& v, const typename AFIELD::real_t *a,
                   const int ieo, const INDEX& block_index)
{
  block_scal_eo_impl(v, a, ieo, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal(AFIELD& v, const typename AFIELD::complex_t *a,
                const INDEX& block_index)
{
  block_scal_eo_impl(v, a, -1, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_scal_eo(AFIELD& v, const typename AFIELD::complex_t *a,
                   const int ieo, const INDEX& block_index)
{
  block_scal_eo_impl(v, a, ieo, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy(AFIELD& v, const typename AFIELD::real_t *a,
                const AFIELD& w, const typename AFIELD::real_t fac,
                const INDEX& block_index)
{
  block_axpy_eo_impl(v, a, w, -1, fac, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy_eo(AFIELD& v, const typename AFIELD::real_t *a,
                   const AFIELD& w, const int ieo,
                   const typename AFIELD::real_t fac,
                   const INDEX& block_index)
{
  block_axpy_eo_impl(v, a, w, ieo, fac, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy(AFIELD& v, const typename AFIELD::complex_t *a,
                const AFIELD& w, const typename AFIELD::real_t fac,
                const INDEX& block_index)
{
  block_axpy_eo_impl(v, a, w, -1, fac, block_index);
}


//====================================================================
template<typename INDEX, typename AFIELD>
void block_axpy_eo(AFIELD& v, const typename AFIELD::complex_t *a,
                   const AFIELD& w, const int ieo,
                   const typename AFIELD::real_t fac,
                   const INDEX& block_index)
{
  block_axpy_eo_impl(v, a, w, ieo, fac, block_index);
}


//============================================================END=====
#endif
