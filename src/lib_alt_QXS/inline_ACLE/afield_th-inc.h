/*!
        @file    afield_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

#ifndef QXS_AFIELD_TH_INCLUDED
#define QXS_AFIELD_TH_INCLUDED

namespace {

//====================================================================
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
    nth = ThreadManager::get_num_threads();
    ith = ThreadManager::get_thread_id();

    size_t is2 = size_t(size) * size_t(ith) / nth;
    size_t ns2 = size_t(size) * size_t(ith + 1) / nth;
    is = int(is2);
    ns = int(ns2);

}

//====================================================================
} // namespace

#endif
