/*!
        @file    fieldIO_NERSC.cpp

        @brief

        @author  Issaku Kanamori
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2014-06-06 18:09:02 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fieldIO_NERSC.h"

#include <fstream>
#include <strings.h>

#include "bridgeIO.h"
using Bridge::vout;

//typedef unsigned short int  n_uint16_t;
//typedef unsigned int        n_uint32_t;
//typedef unsigned long int   n_uint64_t;

const std::string FieldIO_NERSC::class_name = "FieldIO_NERSC";

//====================================================================
void FieldIO_NERSC::read_file(Field& v, const std::string& filename)
{
  const int nin_field = v.nin();
  const int nex_field = v.nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  // temporary buffer: allocated only at I/O node.
  Field vtmp;

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading field data from %s\n", filename.c_str());

    const long_t Lvol = CommonParameters::Lvol();
    vtmp.reset(nin_field, Lvol, nex_field);

    const bool do_swap = (is_bigendian() == false);
    if (do_swap) {
      vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
    }


    // read header and get header size
    std::ifstream infile_header(filename.c_str());
    if (!infile_header) {
      vout.crucial(m_vl, "Error at %s: file open failed, %s may not exist.\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);
    }
    // assumes the header size is smaller than this value
    //   N.B. w/o this kind of limit, whole configuration might be dumped
    //   to stdout
    const size_t max_header_size = 10000; // 25 lines x 400 bytes

    std::string line;
    getline(infile_header, line);
    if (line != "BEGIN_HEADER") {
      vout.crucial(m_vl, "Error at %s: invalid header\n", class_name.c_str());
      vout.crucial(m_vl, "  a bad configuration %s\n", filename.c_str());
      exit(EXIT_FAILURE);
    }
    vout.general(m_vl, "%s: %s\n", class_name.c_str(), line.c_str());
    while (line.find("END_HEADER") == std::string::npos)
    {
      if (infile_header.eof() ||
          (max_header_size < infile_header.tellg())) {
        vout.crucial(m_vl,
                     "%s: the header is too large,  maybe END_HEADER is missing\n",
                     class_name.c_str());
        vout.crucial(m_vl, "  filename: %s\n", filename.c_str());
        vout.crucial(m_vl, "  max_header_size = %zu\n", max_header_size);
        exit(EXIT_FAILURE);
      }
      getline(infile_header, line);
      vout.general(m_vl, "%s: %s\n", class_name.c_str(), line.c_str());
    }
    size_t offset = infile_header.tellg();
    vout.general("%s: offset=%zu\n", class_name.c_str(), offset);
    infile_header.close();

    const int block_size = nin_file;
    char      buf[sizeof(double) * block_size];

    std::ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
    if (!infile) {
      vout.crucial(m_vl, "Error at %s: file open failed, %s may not exist.\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);
    }
    infile.seekg(offset);

    for (int j = 0; j < nex_file; ++j) {
      for (long_t isite = 0; isite < Lvol; ++isite) {
        // read 1 block
        infile.read(buf, sizeof(double) * block_size);

        if (!infile) {
          if (infile.eof()) {  // short file
            vout.crucial(m_vl, "Error at %s: file size of %s is too small.\n", class_name.c_str(), __func__);
          } else {
            vout.crucial(m_vl, "Error at %s: io error of %s.\n", class_name.c_str(), __func__);
          }

          exit(EXIT_FAILURE);
        }

        if (do_swap) {
          byte_swap(buf, sizeof(double), block_size);
        }

        double *ptr = (double *)buf;

        for (int i = 0; i < nin_file; ++i) {
          int s, t;
          m_format->file_to_field(s, t, i, j);

          vtmp.set(s, isite, t, ptr[i]);
        }
      }
    }

    infile.close();
  }

  FieldIO::deliver(&v, &vtmp);

  vout.detailed(m_vl, "read successful\n");
}


//====================================================================
void FieldIO_NERSC::write_file(Field& v, const std::string& filename)
{
  vout.crucial(m_vl, "%s: write_file is not ready\n", class_name.c_str());
  exit(EXIT_FAILURE);

  /*
    Todo:
    see https://arxiv.org/pdf/1202.4813.pdf for the NERSC format
(here is the 3x3 format from 1202.4813)
BEGIN_HEADER
HDR_VERSION = 1.0
DATATYPE = 4D_SU3_GAUGE_3x3
DIMENSION_1 = %(NX)i
DIMENSION_2 = %(NY)i
DIMENSION_3 = %(NZ)i
DIMENSION_4 = %(NT)i
CHECKSUM = %( checksum)s
LINK_TRACE = %( linktrace)f
PLAQUETTE = %( plaquette)f
CREATOR = %( creator)s
ARCHIVE_DATE = %( archive_date)s
ENSEMBLE_LABEL = %( label)s
FLOATING_POINT = %( precision)s
ENSEMBLE_ID = %( ensemble_id)s
SEQUENCE_NUMBER = %( sequence_number)i
BETA = %(beta)f
MASS = %(mass)f
END_HEADER

The definition of CHECKSUM and LINK_TRACE are not explicitly given.

Here is the implementation in GRID for reference (see Grid/parallelIO/* )
 * CHECKSUM: global sum as uint32
 * LINK_TRACE: sum of the trace of all link variable, divided by (V4 * 4dim * 3color)
 * PLAQUETTE: average of all possible plaquette, normalized to [-1:1]
In reading,
 * always checks:  CHECKSUM
 * checks for gauge field:  LINK_TRACE, PLAQUETTE
   */

  /*
  const int nin_field = v->nin();
  const int nex_field = v->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  const long_t Lvol = CommonParameters::Lvol();

  Field vtmp;
  if (Communicator::is_primary()) {
    vtmp.reset(nin_field, Lvol, nex_field);
  }

  FieldIO::gather(&vtmp, v);

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing field data to %s\n", filename.c_str());

    const bool do_swap = (is_bigendian() == false);
    if (do_swap) {
      vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
    }

    const int block_size = nin_file;
    char      buf[sizeof(double) * block_size];

    std::ofstream outfile(filename.c_str(), std::ios::out | std::ios::binary);
    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: file open of %s failed\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);
    }

    for (int j = 0; j < nex_file; ++j) {
      for (long_t isite = 0; isite < Lvol; ++isite) {
        double *ptr = (double *)buf;

        for (int i = 0; i < nin_file; ++i) {
          int s, t;
          m_format->file_to_field(s, t, i, j);

          ptr[i] = vtmp.cmp(s, isite, t);
        }

        if (do_swap) {
          byte_swap(buf, sizeof(double), block_size);
        }

        outfile.write(buf, sizeof(double) * block_size);

        if (!outfile) {  // error
          vout.crucial(m_vl, "Error at %s: io error of %s.\n", class_name.c_str(), __func__);
          exit(EXIT_FAILURE);
        }
      }
    }

    outfile.close();
  }

  vout.detailed(m_vl, "write succeeded.\n");
  */
}


//====================================================================
//============================================================END=====
