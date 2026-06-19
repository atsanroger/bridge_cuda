#ifndef COMMUNICATOR_IMPL_INCLUDED
#define COMMUNICATOR_IMPL_INCLUDED


#ifdef USE_MPI
#include <mpi.h>

#ifdef USE_FJMPI
#include "FJMPI/communicator_mpi.h"
#else
#include "MPI/communicator_mpi.h"
#endif

#else  // USE_MPI not defined
#include "Single/channel.h"
#endif

#endif
