#ifndef JOB_UTILS_H
#define JOB_UTILS_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>


//! read configuration number from file.
void read_input(int& iconf, int& nconf, std::string input_file);

//! write next configuration number to file.
void write_input(int iconf, int nconf, std::string input_file);

//! generate file name for gauge configuration.
std::string filename_config(const std::string filename_base,
                            const int iconf);

#endif  // TEST_SPECTRUM_UTILS_H
