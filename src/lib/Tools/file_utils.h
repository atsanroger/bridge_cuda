/*!
        @file    file_utils.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef FILE_UTILS_INCLUDED
#define FILE_UTILS_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include "Parameters/commonParameters.h"

//! File utility.

/*!
    FileUtils provides file utilities, made by Aoyama-san.
                                   [28 May 2013 Y.Namekawa]
 */

namespace FileUtils
{
  std::string generate_filename(const char *fmt, ...);
}
#endif
