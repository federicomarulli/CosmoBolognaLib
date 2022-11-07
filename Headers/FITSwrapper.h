/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/FITSwrapper.h
 *
 *  @brief class FITSwrapper that wrap CCfits routines to manage FITS
 *  files
 *
 *  This file contains the wrappers of CCfits routines to manage FITS
 *  files
 *
 *  @author Alfonso Veropalumbo 
 *
 *  @author alfonso.veropalumbo@unibo.it
 */

#ifndef __FITSwrap__
#define __FITSwrap__

#include "Kernel.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


namespace cbl {
  
  namespace wrapper {
  
    /**
     *  @brief The namespace of the <B> FITS wrappers </B>
     *  
     *  The \e FITS namespace contains all the wrapper functions of the
     *  FITS routines.
     *
     *  The CCfits headers are expected to be installed in a
     *  subdirectory of the include path. The \<CCfits\> header file
     *  contains all that is necessary to use both the CCfits library
     *  and the cfitsio library (for example, it includes fitsio.h),
     *  thus making all of cfitsio's macro definitions available.
     *
     *  \<CCfits/CCfits\> includes 12 of the CCfits headers and will
     *  support all CCfits operations. The installed location of the
     *  library headers is $(ROOT)/include/CCfits to use the library
     *  either add -I$(ROOT)/include/CCfits or \# include
     *  \<CCfits/CCfits\> in the compilation target.
     */
    namespace ccfits {

      /**
       *  @brief function to read a table from a fits file 
       *
       *  This function searches for columns indicated in column_names
       *  and reads them. If one or more columns are not found in the
       *  selected extension, the output is filled with user-defined
       *  value(s). If no column is found, an error is raised.
       *
       *  @param input_fits name of the input fits file
       *
       *  @param column_names vector containing the columns to read
       *
       *  @param next the number of the extension
       *
       *  @param fill_value if a column is not found, the column
       *  vector is filled with the value specified in fill_value; if
       *  no column is found, an error is raised
       *
       *  @return vector containing the columns
       */
      std::vector<std::vector<double>> read_table_fits (const std::string input_fits, const std::vector<std::string> column_names, const int next=1, const double fill_value=cbl::par::defaultDouble);

      /**
       *  @brief function that writes a table on a fits file
       *
       *  This function writes columns in a binary table on a FITS file
       *
       *  @param output_dir output directory
       *  @param file_fits name of the file fits
       *  @param column_names vector containing the column names to write
       *  @param table vector containing the columns
       *  @param column_units vector containing the column units
       */
      void write_table_fits (const std::string output_dir, const std::string file_fits, const std::vector<std::string> column_names, const std::vector<std::vector<double>> table, const std::vector<std::string> column_units={});

    }
  }
}

#endif
