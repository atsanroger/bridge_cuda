/*!
        @file    afield_dd_double.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-05-07 13:24:28 #$
        @version $LastChangedRevision: 2602 $
*/

#include "lib_alt_QXS/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"

// template definition
#include "lib_alt_QXS/Field/afield_dd-tmpl.h"


// explicit instanciation.
//template class AField<double, QXS>;

typedef AField<double, QXS> AFIELD;
typedef AIndex_block_lex<double, QXS> INDEX;

template
void block_dotc_eo_impl(typename AFIELD::complex_t *out,
                        const AFIELD& v, const AFIELD& w,
                        const int ieo, const INDEX& block_index);
  
template
void block_norm2_eo_impl(typename AFIELD::real_t *out,
                         const AFIELD& v, const int ieo,
                         const INDEX& block_index);

template
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
			const int ieo, const INDEX& block_index);

template
void block_scal_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
			const int ieo, const INDEX& block_index);

template
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::real_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
			const INDEX& block_index);

template
void block_axpy_eo_impl(AFIELD& v, const typename AFIELD::complex_t *a,
                        const AFIELD& w, const int ieo,
                        const typename AFIELD::real_t fac,
                        const INDEX& block_index);

//============================================================END=====
