/*!
        @file    channel.h
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: wlchen $
        @date    $LastChangedDate: 2025-04-04 16:42:21 #$
        @version $LastChangedRevision: 2612 $
*/

/*
   std::vector version with the allocator templated
                            [5 May 2022 I.Kanamori]
 */
#ifndef CHANNEL_H_INCLUDED
#define CHANNEL_H_INCLUDED

#include <mpi.h>
#include <cassert>
#include <vector>

#include "lib/Communicator/communicator.h"
#include "aligned_allocator_impl.h"

#ifdef USE_ACCEL_CUDA
#include <cuda_runtime.h>
#endif

#ifdef USE_ACCEL_OPENACC
#include <openacc.h>
#endif


// forward declaration
class ChannelSet;

//! Channel class for asynchronous communication

/**
   Channel_impl class defines communication channel between sender and receiver
   for asynchronous data transfer. It actually provides a half-channel of
   sender part or receiver part.
   This class wraps MPI persistent communication, as well as send/receive
   buffer.  The template free communicator part is defined in
   Channel_communicator class.

   To access the buffer, operator[] is provided as if it is an ordinary
   container or array.

   start() method to start async data transfer, and wait() method to wait
   for the completion of operation. the channel instance will be re-used
   repeatedly for the same communication channel once it is created.
*/

class Channel_communicator {

 public:
  typedef char element_type;
  static constexpr int max_dimension = 8;

  Channel_communicator();          //!< constructor.
  //  Channel_communicator(const int count); //!< constructor with buffer size (count bytes)
  virtual ~Channel_communicator(); //!< destructor

  int start();                     //!< start asynchronous communication
  int wait();                      //!< wait for completion

  int send_init(int count, int idir, int ipm, void *buf);
  int recv_init(int count, int idir, int ipm, void *buf);

  //! accessor to MPI_Request
  MPI_Request& get_request() { return m_request; }

  int dest(){ return m_dest; }

 private:

  MPI_Request m_request;  //!< handler to MPI persistent communication
  MPI_Status m_status;    //!< handler to MPI status information

  int m_dest;             //!< destination rank

  friend class Communicator_impl;
  //  friend class ChannelSet;
};



template<typename ALLOCATOR>
class Channel_impl {

 public:
  typedef typename Channel_communicator::element_type element_type;
  typedef ALLOCATOR                                   allocator_t;
  typedef element_type*                               container_type;

 private:
  container_type m_buf;         //!< buffer
  element_type *m_ptr;
  Channel_communicator m_comm;  //!< template independent implementation
  
 public:
  Channel_impl() : m_buf(nullptr), m_ptr(nullptr) { }                     //!< constructor.
  Channel_impl(const int count) : m_buf(nullptr), m_ptr(nullptr) {
    allocate(count);
  }
  ~Channel_impl() {
    freeBuffer();
  };       //!< destructor

  void allocate(int count) {
    int size = count * sizeof(element_type);
    if (m_buf != nullptr) {
        freeBuffer();
    }
#ifdef USE_ACCEL_CUDA
    #ifdef USE_DIRECT_COMM
        cudaError_t err = cudaMalloc((void**)&m_buf, size);
        assert(err == cudaSuccess && "cudaMalloc failed");
    #else
        cudaError_t err = cudaMallocHost((void**)&m_buf, size);
        assert(err == cudaSuccess && "cudaMallocHost failed");
    #endif
#elif defined(USE_ACCEL_OPENACC)
        m_buf = (element_type*)malloc(size);
#endif
    m_ptr = m_buf;
  };

void freeBuffer() {
    if (m_buf != nullptr) {
#ifdef USE_ACCEL_CUDA
    #ifdef USE_DIRECT_COMM
        cudaError_t err = cudaFree(m_buf);
        assert(err == cudaSuccess && "cudaFree failed");
    #else
        cudaError_t err = cudaFreeHost(m_buf);
        assert(err == cudaSuccess && "cudaFreeHost failed");
    #endif
#elif defined(USE_ACCEL_OPENACC)
        free(m_buf);
#endif
        m_buf = nullptr;
        m_ptr = nullptr;
    }
}

  //! accessor to buffer
  inline element_type& operator[](unsigned int idx) { return m_buf[idx]; }
  //! accessor to buffer
  inline element_type operator[](unsigned int idx) const { return m_buf[idx]; }
  //! accessor to buffer; returns pointer to the first element.
  inline element_type *ptr() const { return m_ptr; }

  // communication
  int start() { return m_comm.start(); }  //!< start asynchronous communication
  int wait()  { return m_comm.wait(); }   //!< wait for completion

  int send_init(int count, int idir, int ipm){
    allocate(count);
    return m_comm.send_init(count, idir, ipm, (void*)m_ptr);
  }

  int recv_init(int count, int idir, int ipm){
    allocate(count);
    return m_comm.recv_init(count, idir, ipm, (void*)m_ptr);
  }

  MPI_Request& get_request() { return m_comm.get_request(); }

  //! return destination rank
  int destination(){ return m_comm.dest(); }

};



// alias
template<typename T>
constexpr int alignment_size();

template<int ALIGNMENT>
using Channel_aligned = Channel_impl<aligned_allocator_impl<char, ALIGNMENT> >;

// instance with default allocator
using Cannel = Channel_impl<std::allocator<char> >;

//! ChannelSet class for a collection of channels

/**
   ChannelSet defines a collection of channel class instances
   to invoke start and wait methods collectively for a set of channels.
*/
class ChannelSet {
 public:
  ChannelSet(int nchannel = 8); //!< constructor. default number of channels is 8 for upward and downward in 4 dimensions.

  template<typename T>
  int append(Channel_impl<T>& c)   //!< append channel to the set. there is no way to remove a channel.
  {
    return append(c.get_request());
  }

  int start();                  //!< collective start
  int wait();                   //!< collective wait

 private:
  int append(const MPI_Request& r); //!< implementation of adding the channel (template independent)

  std::vector<MPI_Request> m_array; //!< a collection of MPI request held in channels.
  std::vector<MPI_Status> m_status; //!< a collection of MPI status.
  unsigned int m_nreq;              //!< number of channels to hold.
};
#endif /* _CHANNEL_H_ */
