/*!
        @file    threadManager.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#include "ResourceManager/threadManager.h"

#include <omp.h>

#include "Communicator/communicator.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

// fix bug in global reduction/maximization for large size of data
//   [2023.04.14 I.Kanamori]

//====================================================================
namespace ThreadManager_Reduce {
  template<typename REALTYPE>
  void sum_global(REALTYPE *a,
                  const int num,
                  std::vector<REALTYPE>& array_reduction,
                  const int each_buf_size,
                  const int ith, const int nth)
  {
    typedef REALTYPE real_t;
    int                 remaining = num;
    std::vector<real_t> sum;
    real_t              *psum = nullptr;

#pragma omp master
    {
      sum.resize(num);
      for (int i = 0; i < num; i++) {
        sum[i] = 0;
      }
      psum = &sum[0];
    }

    real_t *pa = a;
    while (remaining > 0) // sum over threads; shared buffer size is each_buf_Size
    {
      const int n = (remaining < each_buf_size) ? remaining : each_buf_size;
      for (int i = 0; i < n; ++i) {
        array_reduction[ith * each_buf_size + i] = pa[i];
      }
#pragma omp barrier
#pragma omp master
      {
        for (int i = 0; i < nth; ++i) {
          for (int j = 0; j < n; ++j) {
            psum[j] += array_reduction[i * each_buf_size + j];
          }
        }
        psum += each_buf_size;
      } // master
      pa        += each_buf_size;
      remaining -= each_buf_size;
#pragma omp barrier
    } // sum over threads, done

#pragma omp master
    {
      Communicator::reduce_sum(num, a, &sum[0], 0);
    } // a in the master threads knows the global sum

    remaining = num;
    pa        = a;
    const int total_buf_size = each_buf_size * nth;
    while (remaining > 0) // distributes the sum to each thread
    {
      const int n = (remaining < total_buf_size) ? remaining : total_buf_size;
#pragma omp master
      {
        for (int i = 0; i < n; ++i) { // copy to the common buffer
          array_reduction[i] = pa[i];
        }
      } // master

      // ensures to read updated m_darray_reduction
#pragma omp barrier
      //#ifdef NECSX
      //#pragma omp flush
      //#else
      //if(sizeof(real_t)==4){
      //#pragma omp flush (ThreadManager::m_darray_reductionF)
      //      } else {
      //#pragma omp flush (ThreadManager::m_darray_reduction)
      // }
      //#endif
      for (int i = 0; i < n; ++i) { // copy from the common buffer
        pa[i] = array_reduction[i];
      }
      pa        += total_buf_size;
      remaining -= total_buf_size;
#pragma omp barrier
    } // distributes the sum to each thread, done
  }


  // global maximization is added [2021.12.25 H.Matsufuru]
  template<typename REALTYPE>
  void max_global(REALTYPE *a,
                  const int num,
                  std::vector<REALTYPE>& array_reduction,
                  const int each_buf_size,
                  const int ith, const int nth)
  {
    typedef REALTYPE real_t;
    int                 remaining = num;
    std::vector<real_t> vmax;
    real_t              *pmax = nullptr;

#pragma omp master
    {
      vmax.resize(num);
      for (int i = 0; i < num; ++i) {
        vmax[i] = 0.0;
      }
      pmax = &vmax[0];
    }

    real_t *pa = a;
    while (remaining > 0) // max over threads; shared buffer size is each_buf_size
    {
      const int n = (remaining < each_buf_size) ? remaining : each_buf_size;
      for (int i = 0; i < n; ++i) {
        array_reduction[ith * each_buf_size + i] = pa[i];
      }
#pragma omp barrier
#pragma omp master
      {
        for (int i = 0; i < nth; ++i) {
          for (int j = 0; j < n; ++j) {
#ifdef NECSX
            if (array_reduction[i * each_buf_size + j] > pmax[j]) {
              pmax[j] = array_reduction[i * each_buf_size + j];
            }
#else
            pmax[j] = std::max(pmax[j], array_reduction[i * each_buf_size + j]);
#endif
          }
        }
        pmax += each_buf_size;
      } // master
      pa        += each_buf_size;
      remaining -= each_buf_size;
#pragma omp barrier
    } // maximize over threads, done

#pragma omp master
    {
      Communicator::reduce_max(num, a, &vmax[0], 0);
    } // a in the master threads knows the global max

    remaining = num;
    pa        = a;
    const int total_buf_size = each_buf_size * nth;
    while (remaining > 0) // distributes the max to each thread
    {
      const int n = (remaining < total_buf_size) ? remaining : total_buf_size;
#pragma omp master
      {
        for (int i = 0; i < n; ++i) { // copy to the common buffer
          array_reduction[i] = pa[i];
        }
      } // master

      // ensures to read updated m_darray_reduction
#pragma omp barrier

      for (int i = 0; i < n; ++i) { // copy from the common buffer
        pa[i] = array_reduction[i];
      }
      pa        += total_buf_size;
      remaining -= total_buf_size;
#pragma omp barrier
    } // distributes the sum to each thread, done
  }
}     // namespace ThreadManager_Reduce

//====================================================================
// initialization of static member variables.

int ThreadManager::m_Nthread = 0;
Bridge::VerboseLevel  ThreadManager::m_vl = Bridge::CRUCIAL;
std::vector<dcomplex> ThreadManager::m_darray_reductionDC(0);
std::vector<double>   ThreadManager::m_darray_reduction(0);
std::vector<float>    ThreadManager::m_darray_reductionF(0);

const std::string ThreadManager::class_name = "ThreadManager";

//====================================================================
void ThreadManager::init(int Nthread)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: initialization\n", class_name.c_str());

  int Nthread_env = 0;

#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) {
      Nthread_env = omp_get_num_threads();
    }
  }


  if ((Nthread == Nthread_env) || (Nthread == 0)) {
    m_Nthread = Nthread_env;
  } else {
    vout.general(m_vl, "Warning at %s: Nthread(env) != Nthread(input)\n", class_name.c_str());
    vout.general(m_vl, "  Number of threads(env)   = %d\n", Nthread_env);
    vout.general(m_vl, "  Number of threads(input) = %d\n", Nthread);
    vout.general(m_vl, "  reset Nthread = Nthread(input).\n");

    omp_set_num_threads(Nthread);
    m_Nthread = Nthread;
  }

  vout.general(m_vl, "  Number of threads = %d\n", m_Nthread);

  m_darray_reductionDC.resize(m_Nthread * each_buf_sizeDC);
  m_darray_reduction.resize(m_Nthread * each_buf_size);
  m_darray_reductionF.resize(m_Nthread * each_buf_sizeF);
}


//====================================================================
void ThreadManager::finalize()
{
  vout.paranoiac(m_vl, "%s: finalize.\n", class_name.c_str());
}


//====================================================================
int ThreadManager::get_num_threads()
{
  return omp_get_num_threads();
}


//====================================================================
int ThreadManager::get_thread_id()
{
  return omp_get_thread_num();
}


//====================================================================
void ThreadManager::wait()
{
  int nth = get_num_threads();

  barrier(nth);
}


//====================================================================
void ThreadManager::barrier(int nth)
{
#pragma omp barrier
}


//====================================================================
void ThreadManager::sync_barrier_all()
{
#pragma omp barrier
#pragma omp master
  {
    Communicator::sync();
  }
#pragma omp barrier
}


//====================================================================
void ThreadManager::reduce_sum_global(dcomplex& a,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(&a, 1,
                                   m_darray_reductionDC, each_buf_sizeDC,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(dcomplex *a,
                                      const int num,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(a, num,
                                   m_darray_reductionDC, each_buf_sizeDC,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(double& a,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(&a, 1,
                                   m_darray_reduction, each_buf_size,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(double *a,
                                      const int num,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(a, num,
                                   m_darray_reduction, each_buf_size,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(float& a,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(&a, 1,
                                   m_darray_reductionF, each_buf_sizeF,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(float *a,
                                      const int num,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::sum_global(a, num,
                                   m_darray_reductionF, each_buf_sizeF,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_max_global(double& a,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::max_global(&a, 1,
                                   m_darray_reduction, each_buf_size,
                                   ith, nth);
}


//====================================================================
void ThreadManager::reduce_max_global(double *a,
                                      const int num,
                                      const int ith, const int nth)
{
  ThreadManager_Reduce::max_global(a, num,
                                   m_darray_reduction, each_buf_size,
                                   ith, nth);
}


//====================================================================
void ThreadManager::assert_single_thread(const std::string& name)
{
  int nth = get_num_threads();

  if (nth != 1) {
    vout.crucial(m_vl, "\n");
    vout.crucial(m_vl, "##### Caution #####\n");
    vout.crucial(m_vl, "Single-thread %s is called in parallel region.\n",
                 name.c_str());
    vout.crucial(m_vl, "Current number of thread = %d.\n", nth);

    exit(EXIT_FAILURE);
  }
}


//============================================================END=====
