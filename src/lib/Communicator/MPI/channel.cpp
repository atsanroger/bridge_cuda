/*!
        @file    channel.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 2612 $
*/

#include "channel.h"
#include "layout.h"


/*
       std::vector version with the allocator templated
       this class is a template independent part
                               [5 May 2022 I.Kanamori]
 */

//====================================================================
// class Channel_communicator

//====================================================================

/**
   send_init() method
   creates a channel instance that wraps persistent communication object
   to send data to neighbour node located at upward or downward (specified
   by ipm) in the direction idir.
   count represents size of buffer by number of elements (bytes).
 */
int Channel_communicator::send_init(int count, int idir, int ipm, void *buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Communicator_impl::Layout::m_ndim <= max_dimension);

  //  m_buf.resize(count);

  int dest = (ipm == Forward) ? Communicator_impl::Layout::m_ipe_up[idir] : Communicator_impl::Layout::m_ipe_dn[idir];

  m_dest = dest;
  
  int tag  = idir + max_dimension * (((ipm == Forward) ? 0 : 1) + 2 * Communicator::nodeid());

  int retv = MPI_Send_init(buf, sizeof(element_type) * count, MPI_BYTE, dest, tag, Communicator_impl::world(), &m_request);

  return retv;
}


//====================================================================

/**
   recv_init() method
   creates a channel instance that wraps persistent communication object
   to receive data from neighbour node located at upward or downward (specified
   by ipm) in the direction idir.
   count represents size of buffer by number of elements (bytes).
 */
int Channel_communicator::recv_init(int count, int idir, int ipm, void *buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Communicator_impl::Layout::m_ndim <= max_dimension);

  //  m_buf.resize(count);

  int src = (ipm == Forward) ? Communicator_impl::Layout::m_ipe_up[idir] : Communicator_impl::Layout::m_ipe_dn[idir];
  int tag = idir + max_dimension * (((ipm == Forward) ? 1 : 0) + 2 * src);

  int retv = MPI_Recv_init(buf, sizeof(element_type) * count, MPI_BYTE, src, tag, Communicator_impl::world(), &m_request);

  return retv;
}


//====================================================================
Channel_communicator::Channel_communicator()
{
  LOG;
}


//====================================================================
Channel_communicator::~Channel_communicator()
{
  LOG;
}


//====================================================================
int Channel_communicator::start()
{
  LOG;
  return MPI_Start(&m_request);
}


//====================================================================
int Channel_communicator::wait()
{
  LOG;
  return MPI_Wait(&m_request, &m_status);
}


//====================================================================
// class ChannelSet
ChannelSet::ChannelSet(int count)
  : m_array(count), m_status(count), m_nreq(0)
{
  LOG;
}


//====================================================================
int ChannelSet::append(const MPI_Request& r)
{
  LOG;

  if (m_nreq >= m_array.size()) {
    return MPI_ERR_BUFFER;
  }
  m_array[m_nreq++] = r;

  return MPI_SUCCESS;
}


//====================================================================
int ChannelSet::start()
{
  LOG;
  return MPI_Startall(m_nreq, &m_array[0]);
}


//====================================================================
int ChannelSet::wait()
{
  LOG;
  //  return MPI_Waitall(m_nreq, &m_array[0], (MPI_Status *)0);
  return MPI_Waitall(m_nreq, &m_array[0], &m_status[0]);
}


//====================================================================
// Channel with defaul allocator
template class Channel_impl<std::allocator<char> >;

//====================================================================
//============================================================END=====
