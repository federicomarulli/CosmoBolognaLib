/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/EnumCast.h
 *
 *  @brief Classes used to cast integers and std::string into the enums
 *  used in the CosmoBolognaLib
 *
 *  This file contains the function used to cast integers and
 *  std::string into the enums used in the CosmoBolognaLib
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __ENUMCAST__
#define __ENUMCAST__ 

namespace cbl {

  /**
   * @brief cast an object of type enum class
   * from its index
   * @param i the enum index
   * @return object of class T
   */
  template<typename T, typename std::enable_if<std::is_enum<T>::value>::type* = nullptr>
    T castFromValue (const int i) { return static_cast<T>(i); }

  /**
   * @brief cast an object of type enum class
   * from its name
   * @param name the enum name
   * @param list std::vector containing the names of the 
   * enum
   * @return object of class T
   */
  template<typename T, typename std::enable_if<std::is_enum<T>::value>::type* = nullptr>
    T castFromName (const std::string name, const std::vector<std::string> list) { 
    std::vector<std::string>::const_iterator it = std::find(list.begin(), list.end(), name);
      
    if (it!=list.end()) 
      return castFromValue<T>(int(it-list.begin()));
      
    else {
      std::cerr << par::col_red << std::endl <<"*** Error in castFromName of EnumCast! The name doensn't correspond to any element in list ***" << std::endl << par::col_default << std::endl;
      exit(1);
    }

  }

  /**
   * @brief cast objects of type enum class
   * from indeces
   * @param ii the enum indeces
   * @return object of class T
   */
  template<typename T, typename std::enable_if<std::is_enum<T>::value>::type* = nullptr>
    std::vector<T> castFromValues (const std::vector<int> ii) 
    {
      std::vector<T> en(ii.size());
      for (size_t i=0; i<ii.size(); i++)
	en[i] = castFromValue<T>(ii[i]);
      return en;
    }

  /**
   * @brief cast an object of type enum class
   * from names
   * @param names the enum names
   * @param list std::vector containing the names of the 
   * enum
   * @return object of class T
   */
  template<typename T, typename std::enable_if<std::is_enum<T>::value>::type* = nullptr>
    std::vector<T> castFromNames (const std::vector<std::string> names, const std::vector<std::string> list) 
    {
      std::vector<T> en(names.size());

      for (size_t i=0; i<names.size(); i++)
	en[i] = castFromName<T>(names[i], list);
      return en;
    }

}

#endif
