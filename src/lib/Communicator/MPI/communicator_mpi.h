/*!
        @file    communicator_mpi.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef COMMUNICATOR_MPI_INCLUDED
#define COMMUNICATOR_MPI_INCLUDED

#include <cstdarg>
#include <cstring>
#include <exception>

#include "Communicator/communicator.h"
#include "channel.h"

//! MPI-realisation of communicator class implementation

/**
   Communicator_impl class provides implementation of communication
   between parallel processes using MPI.
   Interface is basically the same as the abstract communicator class.
 */

class Communicator_impl {
 public:
  static int init(int *pargc, char ***pargv);
  static int finalize();
  static void abort();

  static int setup(int Ninstance = 1);

  static int setup(const std::vector<int>& lattice_size,
                   std::vector<int>& grid_size,
                   int Ninstance = 1);

  // info about rank
  static bool is_primary();

#ifdef ENABLE_MULTI_INSTANCE
  static bool is_primary_master();
#endif

  static int self();   //< rank within small world.
  static int size();   //< size of small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int self_global();
  static int world_id();
#endif

  // world
  static MPI_Comm& world() { return m_comm; }  //!< retrieves current communicator.

  // synchronize
  static int sync();   //< synchronize within small world.

  // synchronize (uses usleep to avoid busy wait)
  static int sync_usleep();   //< synchronize within small world. (w/o busy wait)

#ifdef ENABLE_MULTI_INSTANCE
  static int sync_global();   //< synchronize all processes.
#endif

  // info
  static double get_time();

  // debug
  static int status();

  // base case
  class Base {
   public:
    static int reduce(int count, void *recv_buf, void *send_buf, MPI_Datatype type, MPI_Op op, int pattern);
    static int broadcast(size_t size, void *data, int sender);
    static int exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag);

    static int send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag);
  };

  // for specific datatypes
  static int broadcast_string(int count, string& data, int sender);

  // logical and physical layout
  class Layout;

 private:
  Communicator_impl() {}
  Communicator_impl(const Communicator_impl&) {}
  Communicator_impl& operator=(const Communicator_impl&);

  ~Communicator_impl() {}

#ifdef ENABLE_MULTI_INSTANCE
  static int m_Ninstance;   // number of instances
  static int m_instance_id; // id of present instance

  static int m_global_rank;
  static int m_global_size;
#endif

  static int m_grid_rank;
  static int m_grid_size;

  static MPI_Comm m_comm;

  static char default_grid_map[16];
};
#endif /* COMMUNICATOR_MPI_INCLUDED */
