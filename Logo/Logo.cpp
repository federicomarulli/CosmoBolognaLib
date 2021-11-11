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
 *  @file Logo/Logo.cpp
 *
 *  @brief Generic methods of the class Cosmology  
 *
 *  This file contains the logo of the CosmoBolognaLib
 *
 *  @author Sofia Contarini
 */

#include <iostream>
#include "Kernel.h"

using namespace std;

/**
 *  @brief main function to create the logo of the CosmoBolognaLib
 *
 *  @return none
 */
int main ()
{
  cout << " \033[0;33m \n ------------------------------------------------------------------------------ \n\n" << endl;
  cout << " CosmoBolognaLib                      \033[0m*******   " << endl;
  cout << " \033[0;33mcompiled with success!           \033[0m***\033[0;34m/////////\033[0m***   " << endl;
  cout << "                   \033[1;31m       _______ \033[0;34m/////////////////\033[0m*** " << endl;
  cout << "                   \033[1;31m      /       \\ \033[0;34m///////////////////\033[0m**" << endl;
  cout << "                   \033[1;31m     /  ------' \033[0;34m//////////////////////\033[0m*" << endl;
  cout << "                   \033[1;31m    /  / \033[0;34m////////" << " \033[1;31m________  \033[0;34m////////////\033[0m*" << endl; 
  cout << "                   \033[1;31m   (  ( \033[0;34m////////" << " \033[1;31m|   __   \\  \033[0;34m/////////////" << endl;
  cout << "                   \033[1;31m    \\  \\ \033[0;34m///////" << " \033[1;31m|  |__)   | \033[0;34m/////////////\033[0m* " << endl;
  cout << "                      \033[0;34m//\033[1;31m\\  ------- " << " |       _/ \033[0;34m//" << "  \033[1;31m__  \033[0;34m///////" << endl;
  cout << "                      \033[0;34m///\033[1;31m\\_______/ " << " |   ___   \\ \033[0;34m/" << " \033[1;31m|  | \033[0;34m///////\033[0m* " << endl;
  cout << "                      \033[0;34m/////////////" << " \033[1;31m|  |   \\   \\ " << " |  | \033[0;34m///////" << endl;
  cout << "                    \033[0m *\033[0;34m/////////////" << " \033[1;31m|  |___/   / " << " |  | \033[0;34m///////" << endl;
  cout << "                       \033[0;34m////////////" << " \033[1;31m|_________/ \033[0;34m/" <<  " \033[1;31m|  | \033[0;34m//////" << endl;
  cout << "                       \033[0m*\033[0;34m/////////////////////////\033[1;31m |  | \033[0;34m/////" << endl;
  cout << "                        \033[0m *\033[0;34m///////////////////////\033[1;31m |  --------" << endl;
  cout << "                           \033[0m**\033[0;34m////////////////////\033[1;31m |_________/ " << endl;
  cout << "                              \033[0m***\033[0;34m///////////////////      " << endl;
  cout << "                                 \033[0m ***\033[0;34m//////////\033[0m***   " << endl;
  cout << "                                     \033[0m *******                \033[0;33m...Have fun!" << endl;
  cout << " \033[0;33m \n ------------------------------------------------------------------------------ \n\n \033[0m" << endl;

  cbl::Beep("Congratulations! The CosmoBolognaLib compiled with success!");
}
 
