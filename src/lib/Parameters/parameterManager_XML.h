/*!
        @file    parameterManager_XML.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2020-03-17 18:03:39 #$

        @version $LastChangedRevision: 2061 $
 */

#ifndef PARAMETERMANAGER_XML_INCLUDED
#define PARAMETERMANAGER_XML_INCLUDED

#include "parameterManager.h"
#include <string>

//! Parameter manager with YAML parser.

/*!
   This is a simple parser to read parameters from a file
   prepared with YAML format.
   Only simple cases were checked.
                                      [17 Jul 2012 H.Matsufuru]

   read and set parameters from XML file, using tinyxml-2 parser.
   [16 Mar 2015 T.Aoyama]

 */

#ifdef USE_XML

class ParameterManager_XML : public ParameterManager
{
 public:
  static const std::string class_name;

  ParameterManager_XML() {}

  //! read parameters from file.
  void read_params(const std::string& params_file, Parameters& params);
};

#endif /* USE_XML */

#endif
