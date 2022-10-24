/********************************************************************
 *  Copyright (C) 2022 by Sofia Contarini                           *
 *  sofia.contarini3@unibo.it                                       *
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
 *  @file Path/Path.cpp

 *  @brief The class Path used to handle the Cosmobolognalib paths
 *
 *  This file defines the interface of the class Path, used to handle
 *  the Cosmobolognalib paths
 *
 *  @authors Sofia Contarini
 *
 *  @authors sofia.contarini3@unibo.it
 */

#include "Path.h"

// ============================================================================

cbl::Path::Path ()
{
  SetDirs(DIRCOSMO, "./");
}

// ============================================================================

std::string cbl::Path::fullpath (std::string path, const bool isDir)
{ 
  const std::string find = "~";
  const std::string replace = getenv("HOME");
  char buff[PATH_MAX];

  size_t pos = 0;
  while ((pos=path.find(find, pos))!=std::string::npos) {
    path.replace(pos, find.length(), replace);
    pos += replace.length();
  }

  return std::string(realpath(path.c_str(),buff))+((isDir) ? "/" : "");
}
