/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli, Tommaso Ronconi         *
 *  federico.marulli3@unibo.it tommaso.ronconi@studio.unibo.it      *
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
 *  @file Catalogue/GadgetCatalogue.cpp
 *
 *  @brief Methods of the class Catalogue to construct catalogues
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue used to create catalogues directly from GADGET-2.0 files
 *
 *  @author Tommaso Ronconi
 *
 *  @author tommaso.ronconi@studio.unibo.it
 */

#include "Func.h"
#include "Catalogue.h"
#include "Object.h"

using namespace cosmobl;


cosmobl::catalogue::Gadget_Header cosmobl::catalogue::Catalogue::m_swap_header (cosmobl::catalogue::Gadget_Header header)
{
  cosmobl::catalogue::Gadget_Header temp;
  for (int i=0; i<6; i++) temp.npart[i] = IntSwap(header.npart[i]);
  for (int i=0; i<6; i++) temp.massarr[i] = DoubleSwap(header.massarr[i]);
  temp.time = DoubleSwap(header.time);
  temp.redshift = DoubleSwap(header.redshift);
  temp.flag_sfr = IntSwap(header.flag_sfr);
  temp.flag_feedback = IntSwap(header.flag_feedback);
  for (int i=0; i<6; i++) temp.npartTotal[i] = IntSwap(header.npartTotal[i]);
  temp.flag_cool = IntSwap(header.flag_cool);
  temp.nfiles = IntSwap(header.nfiles);
  temp.boxsize = DoubleSwap(header.boxsize);
  temp.omega0 = DoubleSwap(header.omega0);
  temp.omegaLambda = DoubleSwap(header.omegaLambda);
  temp.hubblePar = DoubleSwap(header.hubblePar);
  temp.flag_stAge = IntSwap(header.flag_stAge);
  temp.flag_Metals = IntSwap(header.flag_Metals);
  temp.npart_totHW = IntSwap(header.npart_totHW);
  temp.flag_entr_ics = IntSwap(header.flag_entr_ics);
  for (int i=0; i<40; i++) temp.la[i] = ShortSwap(header.la[i]);
  
  return temp;
}


//==============================================================================================


void cosmobl::catalogue::Catalogue::m_check_it_in (ifstream& finr, bool swap)
{
  finr.read((char *)&m_blockheader, sizeof(m_blockheader));
  if (swap) m_blockheader = IntSwap(m_blockheader);
}

void cosmobl::catalogue::Catalogue::m_check_it_out (ifstream& finr, bool swap)
{
  int check;
  finr.read((char *)&check, sizeof(check));
  if (swap) check = IntSwap(check);
  if (check != m_blockheader) ErrorCBL("Error in reading gadget snapshot, block-headers do not match!");
}

//==============================================================================================

cosmobl::catalogue::Catalogue::Catalogue (const ObjType objType, const string file_cn, const bool swap, const double fact, const bool read_catalogue, const double nSub, const int seed)
{
  Gadget_Header header;
  string gdgt_head = file_cn+".0";
  ifstream finhead(gdgt_head.c_str()); checkIO(finhead, gdgt_head);
  m_check_it_in(finhead, swap);
  finhead.read((char *)&header, sizeof(header));
  if (swap) header = m_swap_header(header);
  m_check_it_out(finhead, swap);
  finhead.clear(); finhead.close();
  
  vector<string> components_name = {"Gas", "Halo", "Disk", "Bulge", "Stars", "Boundary"};

  cout << endl;
  coutCBL << "------- Total Particles -------" << endl << endl;
  for (int i = 0; i<6; i++) if (header.npartTotal[i] != 0) coutCBL << components_name[i]+": " << header.npartTotal[i] << endl;
  cout << endl;
  coutCBL << "------- Mass Resolution -------" << endl << endl;
  for (int i = 0; i<6; i++) if (header.massarr[i] != 0) coutCBL << components_name[i]+": " << header.massarr[i] << "" << endl;
  cout << endl;
  coutCBL << "Age (normalized): " << header.time << endl;
  coutCBL << "Redshift: " << header.redshift << endl;
  coutCBL << "Box Size: " << header.boxsize << " kpc/h" << endl;
  coutCBL << "Omega_M,0 = " << header.omega0 << endl;
  coutCBL << "Omega_Lambda = " << header.omegaLambda << endl;
  coutCBL << "h_0 = " << header.hubblePar << endl;
  cout << endl;
  coutCBL << "Snapshot divided in " << header.nfiles << " files." << endl;
  cout << endl;

  if (read_catalogue) {
    
    // parameters for random numbers used in case nSub!=1
    
    random::UniformRandomNumbers ran(0., 1., seed);
    
    float num_float1, num_float2, num_float3;
    for (int i = 0; i<header.nfiles; i++) {
      string gdgt_snap = file_cn+"."+conv(i, par::fINT);
      ifstream finsnap(gdgt_snap.c_str()); checkIO(finsnap, gdgt_snap);
      coutCBL << "Reading file " << i << " of " << header.nfiles << " ..." << endl;
      
      m_check_it_in(finsnap,swap);
      Gadget_Header data;
      finsnap.read((char *)&data, sizeof(data));
      if (swap) data = m_swap_header(data);
      m_check_it_out(finsnap,swap);
      int dimsnap = 0;
      for (int j = 0; j < 6; j++) dimsnap += data.npart[j];

      //reading particle positions
      m_check_it_in(finsnap,swap);
      
      for (int h = 0; h<dimsnap; h++) {
	
	comovingCoordinates coords;

	finsnap.read((char *)&num_float1, sizeof(num_float1));
	if (swap) num_float1 = FloatSwap(num_float1);
	coords.xx=(num_float1)*fact;

	finsnap.read((char *)&num_float2, sizeof(num_float2));
	if (swap) num_float2 = FloatSwap(num_float2);
	coords.yy=(num_float2)*fact;

	finsnap.read((char *)&num_float3, sizeof(num_float3));
	if (swap) num_float3 = FloatSwap(num_float3);
	coords.zz=(num_float3)*fact;

	if (ran()<nSub) m_object.push_back(move(Object::Create(objType, coords)));
      }
      m_check_it_out(finsnap,swap);
      finsnap.clear(); finsnap.close();
    }

  } //read_catalogue=true
  
  else WarningMsg("Warning: The catalogue is empty!"); //read_catalogue=false
    
}


