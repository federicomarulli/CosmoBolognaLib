/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file Headers/Data.h
 *
 *  @brief The class Data
 *
 *  This file defines the interface of the class Table
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TABLE__
#define __TABLE__

#include "Func.h"

namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> data </B> of any kind
   *  
   *  The \e data namespace contains all the main functions and
   *  classes to handle data of any kind
   */
  namespace data {

    /// Map type used in Table class
    typedef std::map<std::string, std::vector<double>> table_map;

    /**
     *  @class Table Table.h
     *  "Headers/Table.h"
     *
     *  @brief The class Table
     *
     *  This is the base class used generalize datasets
     */
    class Table
    {
    protected:

      /// the table in form of std::map
      table_map m_table;
      
      /**
       * @brief function to set internal members
       *
       * @param names the column names
       *
       * @param nrows the columns lenght
       *
       * 
       */
      void m_set (const std::vector<std::string> names, const size_t nrows);

      /**
       * @brief function to set internal members
       *
       * @param names the column names
       *
       * @param values the table values
       *
       * 
       */
      void m_set (const std::vector<std::string> names, const std::vector<std::vector<double>> values);

      /**
       * @brief get the map keys
       *
       * @return vector of string
       */
      std::vector<std::string> m_keys ();

    public:

      /**
       * @brief default constructor
       *
       * 
       */
      Table () = default;

      /**
       * @brief Constructor of object Table
       *
       * This constructor initializes the table, taking in input the
       * name of the columns and the values. The size of the two
       * arrays must be tha same.
       *
       * @param names the column names
       *
       * @param values the table values
       *
       * 
       */
      Table (const std::vector<std::string> names, const std::vector<std::vector<double>> values);

      /**
       * @brief Constructor of object Table
       *
       * This constructor initializes the table, creating empty
       * columns of size nrows
       *
       * @param names the column names
       *
       * @param nrows the columns lenght
       *
       * 
       */
      Table (const std::vector<std::string> names, const size_t nrows);

      /**
       * @brief  Constructor of object Table
       *
       * This constructor reads the table from an ASCII file
       *
       * @param input_dir directory of input file
       *
       * @param input_file input file name
       *
       * @param names the column names
       *
       * @param use_cols the index of the columns to be read 
       *
       * @param header_lines_to_skip number of header lines
       *
       * 
       */
      Table (const std::string input_dir, const std::string input_file, const std::vector<std::string> names, const std::vector<size_t> use_cols={}, const size_t header_lines_to_skip=1);

      /**
       * @brief default destructor
       *
       * 
       */
      virtual ~Table() = default;

      /**
       * @brief function to get a column
       *
       * @param name the column name
       *
       * @return vector containing the column values
       */
      std::vector<double> & operator[] (const std::string name);

      /**
       * @brief function to get some columns
       *
       * @param names the column names
       *
       * @return vector containing the columns
       */
      std::vector<std::vector<double>> operator[] (const std::vector<std::string> names);

      /**
       * @brief function to insert a column
       *
       * @param name the column name
       *
       * @param values the column values
       *
       * @param replace 
       *
       * 
       */
      void insert (const std::string name, const std::vector<double> values, const bool replace=false);

      /**
       * @brief function to insert some columns
       *
       * @param names the column names
       *
       * @param values the columns
       *
       * @param replace 
       *
       * 
       */
      void insert (const std::vector<std::string> names, const std::vector<std::vector<double>> values, const bool replace=false);

      /**
       * @brief function to read from an ascii file
       *
       * @param input_dir directory of input file
       *
       * @param input_file input file name
       *
       * @param names the column names
       *
       * @param use_cols the index of the columns to be read 
       *
       * @param header_lines_to_skip number of header lines
       *
       * 
       */
      void read (const std::string input_dir, const std::string input_file, const std::vector<std::string> names, const std::vector<size_t> use_cols = {}, const size_t header_lines_to_skip=1);

      /**
       * @brief function to write columns
       *
       * @param output_dir directory of output file
       *
       * @param output_file output file name
       *
       * @param names names of the colums to be wrote
       *
       * 
       */
      void write (const std::string output_dir, const std::string output_file, const std::vector<std::string> names={});
    };
  }
}

#endif
