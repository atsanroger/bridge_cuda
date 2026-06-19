/*!
        @file    threadManager.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef THREADMANAGER_INCLUDED
#define THREADMANAGER_INCLUDED

#include <string>

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"


//! Thread manager with OpenMP.

/*!
   This class wraps OpenMP API functions.
   The number of available threads is specified as an argument
   of init() function, which is to be called at the begining of
   main, soon after Communicator and CommonParameter are
   initialized.
   If the specified number of threads (nthread) is not the same
   as that of OMP_NUM_THREADS environment variable, the number
   of threads is set to the specified one.
   If Nthread = 0 is specified, the number of threads is set
   to OMP_NUM_THREADS or maximum number of platform.
                                      [29 Aug 2013 H.Matsufuru]

   Added reduce_sum_global for array  [30 Sep 2019 I.Kanamori]
   Add complex args                   [08 Aug 2020 Y.Namekawa]
   Renamed: ThreadManager_OpenMP -> ThreadManager
                                      [27 Aug 2022 H.Matsufuru]
 */

class ThreadManager {
 private:
  static int m_Nthread;              //!< number of threads.
  static Bridge::VerboseLevel m_vl;  //!< verbose level.
  static std::vector<dcomplex> m_darray_reductionDC;
  static std::vector<double> m_darray_reduction;
  static std::vector<float> m_darray_reductionF;

  // buffer size for global reduction
#ifdef REDUCTION_EACH_BUF_SIZE_DC
  static const int each_buf_size = REDUCTION_EACH_BUF_SIZE_DC; //!<  reduction buffer size for each thread (dcomplex)
#else
  static const int each_buf_sizeDC = 16;                       //!<  reduction buffer size for each thread (dcomplex)
#endif

#ifdef REDUCTION_EACH_BUF_SIZE
  static const int each_buf_size = REDUCTION_EACH_BUF_SIZE; //!<  reduction buffer size for each thread (double)
#else
  static const int each_buf_size = 8;                       //!<  reduction buffer size for each thread (double)
#endif

#ifdef REDUCTION_EACH_BUF_SIZE_F
  static const int each_buf_sizeF = REDUCTION_EACH_BUF_SIZE_F; //!<  reduction buffer size for each thread (float)
#else
  static const int each_buf_sizeF = 2 * each_buf_size;         //!<  reduction buffer size for each thread (float)
#endif

 private:
  // non-copyable
  ThreadManager(const ThreadManager&);
  ThreadManager& operator=(const ThreadManager&);

 public:
  static const std::string class_name;

  //! setup: called in main only once.
  static void init(int Nthread);

  //! finalization.
  static void finalize();

  //! returns number of threads (works outside of parallel region).
  static int get_num_threads_available() { return m_Nthread; }

  //! returns available number of threads.
  static int get_num_threads();

  //! returns thread id.
  static int get_thread_id();

  //! barrier among threads inside a node.
  static void wait();

  //! barrier among threads inside a node.
  static void barrier(const int Nthread);

  //! barrier among all the threads and nodes.
  static void sync_barrier_all();

  //! global reduction with summation: dcomplex values are assumed thread local.
  static void reduce_sum_global(dcomplex& value,
                                const int i_thread, const int Nthread);

  //! global reduction with summation for an array: dcomplex values are assumed thread local.
  static void reduce_sum_global(dcomplex *value, const int num,
                                const int i_thread, const int Nthread);

  //! global reduction with summation: double values are assumed thread local.
  static void reduce_sum_global(double& value,
                                const int i_thread, const int Nthread);

  //! global reduction with summation for an array: double values are assumed thread local.
  static void reduce_sum_global(double *value, const int num,
                                const int i_thread, const int Nthread);

  //! global reduction with summation: float values are assumed thread local.
  static void reduce_sum_global(float& value,
                                const int i_thread, const int Nthread);

  //! global reduction with summation for an array: float values are assumed thread local.
  static void reduce_sum_global(float *value, const int num,
                                const int i_thread, const int Nthread);

  //! global reduction with max for an array: double values are assumed thread local.
  static void reduce_max_global(double *value, const int num,
                                const int i_thread, const int Nthread);

  //! global reduction with max: double value is assumed thread local.
  static void reduce_max_global(double& value,
                                const int i_thread, const int Nthread);

  //! assert currently running on single thread.
  static void assert_single_thread(const std::string& class_name);
};


//typedef ThreadManager_OpenMP ThreadManager;
typedef ThreadManager ThreadManager_OpenMP;


#endif //THREADMANAGER_INCLUDED
