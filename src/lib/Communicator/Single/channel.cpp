/*!
        @file    channel.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 2499 $
*/

#include "channel.h"

/*
       std::vector version with the allocator templated
       this class is a template independent part
                               [5 May 2022 I.Kanamori]
 */

//====================================================================
// class Channel_communicator

//====================================================================
int Channel_communicator::send_init(int count, int idir, int ipm, void *buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  return 0;
}


//====================================================================
int Channel_communicator::recv_init(int count, int idir, int ipm, void *buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  return 0;
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
  return 0;
}


//====================================================================
int Channel_communicator::wait()
{
  LOG;
  return 0;
}


//====================================================================
// class ChannelSet
ChannelSet::ChannelSet(int count)
{
  LOG;
}


//====================================================================
int ChannelSet::append()
{
  LOG;
  return 0;
}


//====================================================================
int ChannelSet::start()
{
  LOG;
  return 0;
}


//====================================================================
int ChannelSet::wait()
{
  LOG;
  return 0;
}


//====================================================================
// Channel with defaul allocator
template class Channel_impl<std::allocator<char> >;

//====================================================================
//============================================================END=====
