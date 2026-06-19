/*!
        @file    parameterManager.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "parameterManager_YAML.h"

#ifdef USE_XML
#include "parameterManager_XML.h"
#endif

const std::string ParameterManager::class_name = "ParameterManager";

//====================================================================
Parameters ParameterManager::read(const std::string& params_file)
{
  Parameters params;

  read(params_file, params);
  return params;
}


//====================================================================
void ParameterManager::read(const std::string& params_file, Parameters& params)
{
  if (params_file.size() == 0) return;

  vout.general("file = %s\n", params_file.c_str());

  std::string ext = params_file.substr(params_file.find_last_of('.'));

  vout.paranoiac("ext = %s\n", ext.c_str());

  if (ext == ".yaml") {
    ParameterManager_YAML().read_params(params_file, params);
#ifdef USE_XML
  } else if (ext == ".xml") {
    ParameterManager_XML().read_params(params_file, params);
#endif
  } else {
    vout.crucial("Error at %s: unrecognized file type: %s\n", class_name.c_str(), params_file.c_str());
    exit(EXIT_FAILURE);
  }
}


//============================================================END=====
