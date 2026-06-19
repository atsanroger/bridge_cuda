#include "job_Utils.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
string filename_config(const string filename_base, const int iconf)
{
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  vout.paranoiac(vl, "\n");
  vout.paranoiac(vl, "filename_config: setting filename with iconf.\n");
  vout.paranoiac(vl, "  file name base = %s.\n", filename_base.c_str());
  vout.paranoiac(vl, "  iconf = %d.\n", iconf);

  string::size_type pos1 = filename_base.find_first_of("$");
  string::size_type pos2 = filename_base.find_last_of("$");
  vout.paranoiac(vl, "  pos1 = %d, pos2 =%d.\n", pos1, pos2);
  size_t length = pos2 - pos1 + 1;

  std::ostringstream ost;
  ost << std::setw(length) << std::setfill('0') << iconf;
  string num2 = ost.str();

  string filename = filename_base;
  filename.replace(pos1, length, num2);
  vout.paranoiac(vl, "  file name = %s.\n", filename.c_str());

  return filename;
}


//====================================================================
void read_input(int& iconf, int& nconf, std::string filename)
{
  Bridge::VerboseLevel vl = CommonParameters::Vlevel();

  if (Communicator::nodeid() == 0) {
    std::fstream field_in;
    field_in.open(filename.c_str(), std::ios::in);
    if (!field_in.is_open()) {
      vout.crucial(vl, "Failed to open the input file.\n");
      abort();
    }
    field_in >> iconf >> nconf;
    field_in.close();
  }

  Communicator::broadcast(1, &iconf, 0);
  Communicator::broadcast(1, &nconf, 0);

  //vout.general(vl,"iconf = %d  nconf = %d\n", iconf, nconf);
}


//====================================================================
void write_input(int iconf_next, int nconf, std::string filename)
{
  Bridge::VerboseLevel vl = CommonParameters::Vlevel();

  if (Communicator::nodeid() == 0) {
    std::fstream field_out;
    field_out.open(filename.c_str(), std::ios::out);
    if (!field_out.is_open()) {
      vout.crucial(vl, "Failed to open the text file.\n");
    } else {
      field_out << iconf_next << "  " << nconf << std::endl;
      field_out.close();
    }
  }
}


//============================================================END=====
