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
#include "Halo.h"
#include "HostHalo.h"
using namespace cbl;

cbl::catalogue::SubFindTab_Header cbl::catalogue::Catalogue::m_read_header (std::ifstream& finh, const bool swap)
{
  cbl::catalogue::SubFindTab_Header header;
  //int intvar;
  uint32_t intvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.Ngroups = (swap) ? IntSwap(intvar) : intvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.totNgroups = (swap) ? IntSwap(intvar) : intvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.Nids = (swap) ? IntSwap(intvar) : intvar;
  //long longvar;
  uint64_t longvar;
  finh.read((char *)&longvar, sizeof(longvar));
  header.totNids = (swap) ? LongSwap(longvar) : longvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.Ntask = (swap) ? IntSwap(intvar) : intvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.Nsubs = (swap) ? IntSwap(intvar) : intvar;
  finh.read((char *)&intvar, sizeof(intvar));
  header.totNsubs = (swap) ? IntSwap(intvar) : intvar;
  
  return header;
}

//==============================================================================================

cbl::catalogue::Gadget_Header cbl::catalogue::Catalogue::m_swap_header (cbl::catalogue::Gadget_Header header)
{
  cbl::catalogue::Gadget_Header temp;
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

cbl::catalogue::SubFindTab_Header cbl::catalogue::Catalogue::m_swap_header (cbl::catalogue::SubFindTab_Header header)
{
  cbl::catalogue::SubFindTab_Header temp;
  temp.Ngroups = IntSwap(header.Ngroups);
  temp.totNgroups = IntSwap(header.totNgroups);
  temp.Nids = IntSwap(header.Nids);
  temp.totNids = LongSwap(header.totNids);
  temp.Ntask = IntSwap(header.Ntask);
  temp.Nsubs = IntSwap(header.Nsubs);
  temp.totNsubs = IntSwap(header.totNsubs);

  return temp;
}

//==============================================================================================

void cbl::catalogue::Catalogue::m_check_it_in (std::ifstream& finr, const bool swap)
{
  finr.read((char *)&m_blockheader, sizeof(m_blockheader));
  if (swap) m_blockheader = IntSwap(m_blockheader);
}

void cbl::catalogue::Catalogue::m_check_it_out (std::ifstream& finr, const bool swap)
{
  int check;
  finr.read((char *)&check, sizeof(check));
  if (swap) check = IntSwap(check);
  if (check != m_blockheader) ErrorCBL("block-headers of gadget snapshot do not match!", "m_check_it_in", "GadgetCatalogue.cpp");
}

//==============================================================================================

cbl::catalogue::Catalogue::Catalogue (const ObjectType objectType, const std::string file_cn, const bool snapformat, const bool swap, const double fact, const bool read_catalogue, const double nSub, const std::string component_to_read, const std::vector<std::vector<double>> edges)
{
  std::string gdgt_head = file_cn+".0";
  std::ifstream finhead(gdgt_head.c_str(), std::ios::binary|std::ios::in); checkIO(finhead, gdgt_head);
  finhead.seekg(finhead.beg);
	bool cut=false;
	if (edges[0][0]!=par::defaultDouble && edges[0][1]!=par::defaultDouble && edges[1][0]!=par::defaultDouble && edges[1][1]!=par::defaultDouble
							&& edges[2][0]!=par::defaultDouble && edges[2][1]!=par::defaultDouble) cut=true;

  // for snapformat = 2:
  if (snapformat) {
    m_check_it_in(finhead, swap);
    char charvar;
    std::string stringvar;
    for (int i=0; i<4; i++) {
      finhead.read((char *)&charvar, sizeof(charvar));
      stringvar.push_back(charvar);
    }
    int intvar;
    finhead.read((char *)&intvar, sizeof(intvar));
    m_check_it_out(finhead, swap);
  }

  m_check_it_in(finhead,swap);
  Gadget_Header header;
  finhead.read((char *)&header, sizeof(header));
  if (swap) header = m_swap_header(header);
  m_check_it_out(finhead, swap);
  finhead.clear(); finhead.close();
  
  std::vector<std::string> components_name = {"Gas", "Halo", "Disk", "Bulge", "Stars", "Boundary"};

  coutCBL<< std::endl;
  coutCBL << "------- Total Particles -------" << std::endl << std::endl;
  for (int i = 0; i<6; i++) if (header.npartTotal[i] != 0) coutCBL << components_name[i]+": " << header.npartTotal[i] << std::endl;
  coutCBL<< std::endl;
  coutCBL << "------- Mass Resolution -------" << std::endl << std::endl;
  for (int i = 0; i<6; i++) if (header.massarr[i] != 0) coutCBL << components_name[i]+": " << header.massarr[i] << "" << std::endl;
  coutCBL<< std::endl;
  coutCBL << "Age (normalized): " << header.time << std::endl;
  coutCBL << "Redshift: " << header.redshift << std::endl;
  coutCBL << "Box Size: " << header.boxsize << " kpc/h" << std::endl;
  coutCBL << "Omega_M,0 = " << header.omega0 << std::endl;
  coutCBL << "Omega_Lambda = " << header.omegaLambda << std::endl;
  coutCBL << "h_0 = " << header.hubblePar << std::endl;
  coutCBL<< std::endl;
  coutCBL << "Snapshot divided in " << header.nfiles << " files." << std::endl;
  coutCBL<< std::endl;

  if (read_catalogue) {
    
    // parameters for random numbers used in case nSub!=1
    std::default_random_engine gen;
    std::uniform_real_distribution<float> ran(0., 1.);

    float num_float1, num_float2, num_float3;
    for (int i = 0; i<header.nfiles; i++) {
      std::string gdgt_snap = file_cn+"."+conv(i, par::fINT);
      std::ifstream finsnap(gdgt_snap.c_str(), std::ios::binary); checkIO(finsnap, gdgt_snap);
      coutCBL << "Reading file " << i+1 << " of " << header.nfiles << " ..." << std::endl;

      // for snapformat = 2:
      if (snapformat) {
				m_check_it_in(finsnap, swap);
				char charvar;
				std::string stringvar;
				for (int i=0; i<4; i++) {
					finsnap.read((char *)&charvar, sizeof(charvar));
					stringvar.push_back(charvar);
				}
				int intvar;
				finsnap.read((char *)&intvar, sizeof(intvar));
				m_check_it_out(finsnap, swap);
      }
      
      Gadget_Header data;
      m_check_it_in(finsnap,swap);
      finsnap.read((char *)&data, sizeof(data));
      if (swap) data = m_swap_header(data);
      m_check_it_out(finsnap,swap);
      
      int dimsnap = 0;
      for (int j = 0; j < 6; j++) dimsnap += data.npart[j];
      std::vector<bool> read (dimsnap, false);
      if (component_to_read == "ALL") std::replace(read.begin(), read.end(), false, true);
      else {
				int offset = 0;
				bool wrong_component_name = true;
				for (size_t jj = 0; jj<components_name.size(); jj++) {
	  			if (component_to_read == components_name[jj]) std::replace(read.begin()+offset,
								     read.begin()+offset+data.npart[jj],
								     false, true);
	  			offset += data.npart[jj];
	  			if (component_to_read == components_name[jj]) wrong_component_name = false;
				}
				if (wrong_component_name) WarningMsgCBL("selected component is not available, available components are ALL, Gas, Halo, Disk, Bulge, Stars, Boundary", "Catalogue", "GadgetCatalogue.cpp");
				if (offset != dimsnap) ErrorCBL("something horrible happened...", "Catalogue", "GadgetCatalogue.cpp");
      }
	
      // reading particle positions
      
      // for snapformat = 2:
      if (snapformat) {
				m_check_it_in(finsnap, swap);
				char charvar;
				std::string stringvar;
				for (int i=0; i<4; i++) {
	  			finsnap.read((char *)&charvar, sizeof(charvar));
	  			stringvar.push_back(charvar);
				}
				int intvar;
				finsnap.read((char *)&intvar, sizeof(intvar));
				m_check_it_out(finsnap, swap);
      }
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

				if (read[h]) {
					if (cut) {
						if (ran(gen)<nSub && coords.xx>edges[0][0] && coords.xx<edges[0][1] &&
									coords.yy>edges[1][0] && coords.yy<edges[1][1] &&
									coords.zz>edges[2][0] && coords.zz<edges[2][1] ) m_object.push_back(move(Object::Create(objectType, coords)));
					}
					else 
						if (ran(gen)<nSub) m_object.push_back(move(Object::Create(objectType, coords)));
				}
      }
      m_check_it_out(finsnap,swap);
      finsnap.clear(); finsnap.close();
    }
  }//read_catalogue=true
  
  else WarningMsgCBL("the catalogue is empty!", "Catalogue", "GadgetCatalogue.cpp"); //read_catalogue=false
    
}

//==============================================================================================

cbl::catalogue::Catalogue::Catalogue (const int snap, const std::string basedir, const bool swap, const bool long_ids, const double scaleFact, const double massFact, const EstimateCriterion estimate_crit, const bool veldisp, const bool masstab, const bool add_satellites, const bool verbose) 
{
  
  // defines the name of the catalogue files
  std::string snap_str = conv(snap, par::fINT);
  while (snap_str.size()!= 3) snap_str = "0"+snap_str;
  std::string file_base = basedir+"/groups_"+snap_str+"/subhalo_tab_"+snap_str+".";
  
  //typedef conditional<long_ids,
  //		      uint64_t,
  //		      uint32_t>::type ids_type;
  
  // begin a cicle to read and store data from the group files
  size_t NG=0, totNG=0, Nfiles=0, NS=0, totNS=0;
  std::vector<std::shared_ptr<Object>> groups, subgroups;
  unsigned int count_G = 0, count_S = 0;
  
  unsigned int filenum = 0;
  bool doneflag = false;
  while (!doneflag) {

    // try open file and check it works
    std::string current_file = file_base + conv(filenum, par::fINT);
    std::ifstream fincur (current_file.c_str(), std::ifstream::binary); checkIO(fincur, current_file);

    // read file header
    SubFindTab_Header header = m_read_header(fincur, swap);
    
    // store infos from header
    NG = header.Ngroups;
    NS = header.Nsubs;
    if (filenum == 0) {
      totNG = header.totNgroups;
      totNS = header.totNsubs;
      Nfiles = header.Ntask;
    }

    coutCBL << "Building catalogue " << std::setprecision(2) << std::setiosflags(std::ios::fixed) << std::setw(8) << 100.*(filenum+1)/Nfiles << " % done \r"; std::cout.flush();

    if (NG > 0) {
      // read data blocks, store those necessary
      fincur.seekg(sizeof(uint32_t)*NG, fincur.cur); // jumps the 'lenght' block
      fincur.seekg(sizeof(uint32_t)*NG, fincur.cur); // jumps the 'offset' block

      std::vector<float> vec_tot_mass;
      vectorReadFromBinary(fincur, vec_tot_mass, NG);
    
      std::vector<std::vector<float>> vec_pos;
      for (size_t ii = 0; ii < (size_t) NG; ii++) { 
	std::vector<float> pos;
	vectorReadFromBinary(fincur, pos, 3);
	vec_pos.push_back(pos);
      }

      std::vector<float> vec_mass_estimate, vec_radius_estimate, vec_veldisp_estimate;
      
      switch (estimate_crit) {

      case EstimateCriterion::_m200_:
	vectorReadFromBinary(fincur, vec_mass_estimate, NG);
	vectorReadFromBinary(fincur, vec_radius_estimate, NG);
	fincur.seekg(sizeof(float)*4*NG, fincur.cur);
	if (veldisp) {
	  vectorReadFromBinary(fincur, vec_veldisp_estimate, NG);
	  fincur.seekg(sizeof(float)*2*NG, fincur.cur);
	}
	break;

      case EstimateCriterion::_c200_:
	fincur.seekg(sizeof(float)*2*NG, fincur.cur);
	vectorReadFromBinary(fincur, vec_mass_estimate, NG);
	vectorReadFromBinary(fincur, vec_radius_estimate, NG);
	fincur.seekg(sizeof(float)*2*NG, fincur.cur);
	if (veldisp) {
	  fincur.seekg(sizeof(float)*NG, fincur.cur);
	  vectorReadFromBinary(fincur, vec_veldisp_estimate, NG);
	  fincur.seekg(sizeof(float)*NG, fincur.cur);
	}
	break;

      case EstimateCriterion::_t200_:
	fincur.seekg(sizeof(float)*4*NG, fincur.cur);
	vectorReadFromBinary(fincur, vec_mass_estimate, NG);
	vectorReadFromBinary(fincur, vec_radius_estimate, NG);
	if (veldisp) {
	  fincur.seekg(sizeof(float)*2*NG, fincur.cur);
	  vectorReadFromBinary(fincur, vec_veldisp_estimate, NG);
	}
	break;
	
      default:
	ErrorCBL("wrong estimate criterion selected; available criteria are: _m200_, _c200_ and _t200_", "Catalogue", "GadgetCatalogue.cpp");
      } // endswitch (estimate_crit)

      fincur.seekg(sizeof(uint32_t)*NG, fincur.cur); // jumps the 'contamination count' block
      fincur.seekg(sizeof(float)*NG, fincur.cur); // jumps the 'contamination mass' block

      std::vector<uint32_t> vec_nsubs;
      vectorReadFromBinary(fincur, vec_nsubs, NG);

      fincur.seekg(sizeof(uint32_t)*NG, fincur.cur); // end of Group data block

      // add data to group catalogue:
      size_t offset = groups.size();
      if (verbose) coutCBL << "Current dimension of catalogue: " << offset << std::endl;
      for (size_t ii = offset; ii < offset+(size_t)NG; ii++) {
	
	// add ii-th object to group catalogue
	comovingCoordinates coords = {scaleFact*vec_pos[ii-offset][0],
				      scaleFact*vec_pos[ii-offset][1], 
				      scaleFact*vec_pos[ii-offset][2]};
	groups.push_back(move(Object::Create(ObjectType::_HostHalo_, coords)));
	
	// set _Halo_ variables of the ii-th object
	groups[ii]->set_tot_mass(massFact*vec_tot_mass[ii-offset]);
	
	// set _HostHalo_ variables of the ii-th object
	groups[ii]->set_mass(massFact*vec_tot_mass[ii-offset]); // temporarily = to total mass of the halo
	groups[ii]->set_mass_estimate(massFact*vec_mass_estimate[ii-offset]);
	groups[ii]->set_radius_estimate(scaleFact*vec_radius_estimate[ii-offset]);
	if (veldisp) groups[ii]->set_veldisp_estimate(vec_veldisp_estimate[ii-offset]);
	groups[ii]->set_ID((int) ii);
	groups[ii]->set_nsub((int) vec_nsubs[ii-offset]);
	groups[ii]->set_parent(-1);
	count_G++;
      } // endfor (int ii = offset; ii < offset+NG; ii++)    
    } // endif (NG > 0)
    
    if (NS > 0) {

      // read data blocks, store those necessary
      fincur.seekg(sizeof(uint32_t)*NS, fincur.cur); // jumps the 'lenght' block
      fincur.seekg(sizeof(uint32_t)*NS, fincur.cur); // jumps the 'offset' block
      fincur.seekg(sizeof(uint32_t)*NS, fincur.cur); // jumps the 'parent' block
      
      std::vector<float> vec_sub_mass;
      vectorReadFromBinary(fincur, vec_sub_mass, NS); // reads the 'mass' block

      std::vector<std::vector<float>> vec_sub_pos;
      // reads the 'pos' block:
      for (size_t ii = 0; ii<NS; ii++) {
	std::vector<float> vec;
	vectorReadFromBinary(fincur, vec, 3);
	vec_sub_pos.push_back(vec);
      }

      fincur.seekg(3*sizeof(float)*NS, fincur.cur); // jumps the 'vel' block
      
      std::vector<std::vector<float>> vec_sub_cm;
      // reads the 'cm' block:
      for (size_t ii = 0; ii<NS; ii++) {
	std::vector<float> vec;
	vectorReadFromBinary(fincur, vec, 3);
	vec_sub_cm.push_back(vec);
      }
      std::vector<std::vector<float>> vec_sub_spin;
      // reads the 'spin' block:
      for (size_t ii = 0; ii<NS; ii++) {
	std::vector<float> vec;
	vectorReadFromBinary(fincur, vec, 3);
	vec_sub_spin.push_back(vec);
      }
      std::vector<float> vec_sub_veldisp;
      vectorReadFromBinary(fincur, vec_sub_veldisp, NS); // reads the 'veldisp' block
      std::vector<float> vec_sub_vmax;
      vectorReadFromBinary(fincur, vec_sub_vmax, NS); // reads the 'vmax' block
      std::vector<float> vec_sub_vmax_rad;
      vectorReadFromBinary(fincur, vec_sub_vmax_rad, NS); // reads the 'vmax_rad' block
      std::vector<float> vec_halfmass_rad;
      vectorReadFromBinary(fincur, vec_halfmass_rad, NS); //reads the 'halfmass_rad' block
      //fincur.seekg(sizeof(float)*NS, fincur.cur); // jumps the 'halfmass_rad' block

      // jumps the 'ID_most_bound' block:
      if (!long_ids) fincur.seekg(sizeof(uint32_t)*NS, fincur.cur); 
      else fincur.seekg(sizeof(uint64_t)*NS, fincur.cur);
      
      std::vector<uint32_t> vec_grnr;
      vectorReadFromBinary(fincur, vec_grnr, NS); // reads the 'GrNr' block
      
      if (masstab) {
	fincur.seekg(6*sizeof(float)*NS, fincur.cur); // jumps the 'masstab' block, if present
      } // endif (masstab)

      // add data to catalogue:
      size_t offset_sub = subgroups.size();
      for (unsigned int jj = offset_sub; jj < offset_sub+NS; jj++) {
	
	// add jj-th object to sub-group catalogue
	comovingCoordinates coords = {scaleFact*vec_sub_pos[jj-offset_sub][0],
				      scaleFact*vec_sub_pos[jj-offset_sub][1],
				      scaleFact*vec_sub_pos[jj-offset_sub][2]};
	subgroups.push_back(move(Object::Create(ObjectType::_HostHalo_, coords)));
	
	// set _Halo_ variables of the jj-th object
	subgroups[jj]->set_mass(massFact*vec_sub_mass[jj-offset_sub]);

	// set _HostHalo_ variables of the jj-th object
	subgroups[jj]->set_xcm(scaleFact*vec_sub_cm[jj-offset_sub][0]);
	subgroups[jj]->set_ycm(scaleFact*vec_sub_cm[jj-offset_sub][1]);
	subgroups[jj]->set_zcm(scaleFact*vec_sub_cm[jj-offset_sub][2]);
	subgroups[jj]->set_spin_x(vec_sub_spin[jj-offset_sub][0]);
	subgroups[jj]->set_spin_y(vec_sub_spin[jj-offset_sub][1]);
	subgroups[jj]->set_spin_z(vec_sub_spin[jj-offset_sub][2]);
	subgroups[jj]->set_veldisp(vec_sub_veldisp[jj-offset_sub]);
	subgroups[jj]->set_vmax(vec_sub_vmax[jj-offset_sub]);
	subgroups[jj]->set_vmax_rad(vec_sub_vmax_rad[jj-offset_sub]);
	subgroups[jj]->set_radius(vec_halfmass_rad[jj-offset_sub]);
	subgroups[jj]->set_ID((int) jj);
	subgroups[jj]->set_nsub(-1);
	subgroups[jj]->set_parent((int) vec_grnr[jj-offset_sub]);
	count_S++;
      } // endfor (int jj = offset_sub; jj < offset_sub+NS; jj++)   
    } // endif (NS > 0)

    // final check
    int curpos = fincur.tellg();
    fincur.seekg(0, fincur.end);
    if (curpos != fincur.tellg()) WarningMsgCBL("finished reading before EOF for file "+current_file, "Catalogue", "GadgetCatalogue.cpp");
    fincur.clear(); fincur.close();
    filenum ++;
    if (filenum >= Nfiles) doneflag = true;

  } // endwhile (!doneflag)
  if (count_G != totNG) ErrorCBL("groups read do not match expected!", "Catalogue", "GadgetCatalogue.cpp");
  if (count_S != totNS) ErrorCBL("subGroups read do not match expected!", "Catalogue", "GadgetCatalogue.cpp");

  // sorting the sub-groups according to their parent ID
  coutCBL << "Now sorting sub-groups ..." << std::endl;
  std::vector<int> parentIDs;
  for (size_t ii = 0; ii<subgroups.size(); ii++) parentIDs.push_back(subgroups[ii]->parent());
  bool bubbleswap = true;
  while (bubbleswap) {
    bubbleswap = false;
    for (size_t jj = 0; jj<subgroups.size()-1; ++jj) {
      if (parentIDs[jj] > parentIDs[jj+1]) {
	int tempVar = parentIDs[jj];
	parentIDs[jj] = parentIDs[jj+1];
	parentIDs[jj+1] = tempVar;
	std::shared_ptr<Object> tempObj = subgroups[jj];
	subgroups[jj] = subgroups[jj+1];
	subgroups[jj+1] = tempObj;
	bubbleswap = true;
      }
    }
  }
  coutCBL << "... done!" << std::endl;

  // building the catalogue with dependencies central/satellite
  for (size_t ii = 0; ii<groups.size(); ii++) {
    
    coutCBL << "Identifying centrals and satellites " << std::setprecision(2) << std::setiosflags(std::ios::fixed) << std::setw(8) << 100.*(ii+1)/groups.size() << " % done \r"; std::cout.flush(); 

    // groups with sub-groups within are linked to their satellites:
    if (groups[ii]->nsub() > 0) {
      int currentID = groups[ii]->ID();
      std::pair<std::vector<int>::iterator, std::vector<int>::iterator> bounds = equal_range(parentIDs.begin(), 
									      parentIDs.end(), 
									      currentID);
      int index = bounds.first - parentIDs.begin();
      if (verbose) {
	coutCBL << "check coordinates: " << std::endl;
	coutCBL << "\t Central: " 
		<< groups[ii]->xx() << "  "
		<< groups[ii]->yy() << "  "
		<< groups[ii]->zz() << std::endl;
	coutCBL << "\t Satellite: " 
		<< subgroups[index]->xx() << "  " 
		<< subgroups[index]->yy() << "  " 
		<< subgroups[index]->zz() << std::endl;
      }
      groups[ii]->set_mass(subgroups[index]->mass());
      groups[ii]->set_xcm(subgroups[index]->xcm());
      groups[ii]->set_ycm(subgroups[index]->ycm());
      groups[ii]->set_zcm(subgroups[index]->zcm());
      groups[ii]->set_spin_x(subgroups[index]->spin_x());
      groups[ii]->set_spin_y(subgroups[index]->spin_y());
      groups[ii]->set_spin_z(subgroups[index]->spin_z());
      groups[ii]->set_veldisp(subgroups[index]->veldisp());
      groups[ii]->set_vmax(subgroups[index]->vmax());
      groups[ii]->set_vmax_rad(subgroups[index]->vmax_rad());
      groups[ii]->set_radius(subgroups[index]->radius());
      m_object.push_back(groups[ii]); // adds a central halo to the catalogue.
      int last_central = m_object.size() -1; // remembers the central index
      for (int jj = 1; jj < bounds.second - bounds.first; jj++) {
	
	index = bounds.first - parentIDs.begin() + jj;
	subgroups[index]->set_tot_mass(groups[ii]->tot_mass());
	subgroups[index]->set_mass_estimate(groups[ii]->mass_estimate());
	subgroups[index]->set_radius_estimate(groups[ii]->radius_estimate());
	if (veldisp) subgroups[index]->set_veldisp_estimate(groups[ii]->veldisp_estimate());
	m_object[last_central]->set_satellite(subgroups[index]);
	// - if add_satellites = true adds the satellite to the catalogue as an HostHalo object, 
	// - if add_satellites = false, satellites can be only retrieved by the 'satellites' 
	//   member of the central HostHalo object:
	if (add_satellites) m_object.push_back(subgroups[index]); 
      } // endfor jj      
      if ((int) m_object[last_central]->satellites().size() != m_object[last_central]->nsub()-1) 
	ErrorCBL("number of satellites do not match!", "Catalogue", "GadgetCatalogue.cpp");
    } // endif (m_object[ii]->nsub() > 0)

    // groups without sub-groups added as-they-are to the catalogue
    // [will lack of sub-group properties (e.g. CM coords, Spin, ...)]
    else if (groups[ii]->nsub() == 0) m_object.push_back(groups[ii]);

  } //endfor ii
  
}
