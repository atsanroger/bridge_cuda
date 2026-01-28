/*!
        @file    afield_Gauge_index-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-17 13:51:32 #$
        @version $LastChangedRevision: 2550 $
*/

#ifndef QXS_AFIELD_GAUGE_INDEX_INC_INCLUDED
#define QXS_AFIELD_GAUGE_INDEX_INC_INCLUDED

namespace {
  inline int idxGr(int ic1, int ic2){ return     2*(ic1 + NC * ic2); }
  inline int idxGi(int ic1, int ic2){ return 1 + 2*(ic1 + NC * ic2); }
}

#endif
