/*!
        @file    threadManager.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "ResourceManager/threadManager.h"

#include "Communicator/communicator.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

// These are stub routine when OpenMP is not used.
//                                          [11 Mar 2014 H.Matsufuru]
//====================================================================
// initialization of static member variables.

int ThreadManager::m_Nthread             = 0;
Bridge::VerboseLevel ThreadManager::m_vl = Bridge::CRUCIAL;
std::vector<double>  ThreadManager::m_darray_reduction(0);
std::vector<float>   ThreadManager::m_darray_reductionF(0);

//====================================================================
void ThreadManager::init(int Nthread)
{
  m_vl = CommonParameters::Vlevel();

  m_Nthread = 1;

  vout.general(m_vl, "ThreadManager being setup.\n");
  vout.general(m_vl, "  Number of thread = %d\n", m_Nthread);
  vout.general(m_vl, "  OpenMP is not used: stub implementation.\n");
}


//====================================================================
void ThreadManager::finalize()
{
  vout.paranoiac(m_vl, "Thread manager says good-bye.\n");
}


//====================================================================
int ThreadManager::get_num_threads()
{
  return 1;
}


//====================================================================
int ThreadManager::get_thread_id()
{
  return 0;
}


//====================================================================
void ThreadManager::wait()
{
}


//====================================================================
void ThreadManager::barrier(int Nthread)
{
}


//====================================================================
void ThreadManager::sync_barrier_all()
{
  Communicator::sync();
}


//====================================================================
void ThreadManager::reduce_sum_global(dcomplex& a,
                                      const int ith, const int nth)
{
  reduce_sum_global(&a, 1, ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(dcomplex *a,
                                      const int num,
                                      const int ith, const int nth)
{
  std::vector<dcomplex> sum(num);
  Communicator::reduce_sum(num, &sum[0], a, 0);
  for (int i = 0; i < num; i++) {
    a[i] = sum[i];
  }
}


//====================================================================
void ThreadManager::reduce_sum_global(double& a,
                                      const int ith, const int nth)
{
  reduce_sum_global(&a, 1, ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(double *a,
                                      const int num,
                                      const int ith, const int nth)
{
  std::vector<double> sum(num);
  Communicator::reduce_sum(num, &sum[0], a, 0);
  for (int i = 0; i < num; i++) {
    a[i] = sum[i];
  }
}


//====================================================================
void ThreadManager::reduce_sum_global(float& a,
                                      const int ith, const int nth)
{
  reduce_sum_global(&a, 1, ith, nth);
}


//====================================================================
void ThreadManager::reduce_sum_global(float *a,
                                      const int num,
                                      const int ith, const int nth)
{
  std::vector<float> sum(num);
  Communicator::reduce_sum(num, &sum[0], a, 0);
  for (int i = 0; i < num; i++) {
    a[i] = sum[i];
  }
}


//====================================================================
void ThreadManager::reduce_max_global(double *a, const int num,
                                      const int ith, const int nth)
{
  std::vector<double> vmax(num);
  Communicator::reduce_max(num, &vmax[0], a, 0);
  for (int i = 0; i < num; ++i) {
    a[i] = vmax[i];
  }
}


//====================================================================
void ThreadManager::reduce_max_global(double& a,
                                      const int ith, const int nth)
{
  reduce_max_global(&a, 1, ith, nth);
}


//====================================================================
void ThreadManager::assert_single_thread(const std::string& class_name)
{
  // do nothing.
}


//============================================================END=====
