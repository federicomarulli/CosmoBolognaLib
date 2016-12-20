/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/Catalogue.h
 *
 *  @brief The class Catalogue  
 *
 *  This file defines the interface of the class Catalogue, used
 *  handle catalogues of astronomical sources
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#ifndef __CATALOGUE__
#define __CATALOGUE__ 

#include "Field3D.h"
#include "ChainMesh.h"
#include "Object.h"
#include "RandomObject.h"
#include "Halo.h"
#include "Mock.h"
#include "Galaxy.h"
#include "Cluster.h"
#include "Void.h"


// ============================================================================================


namespace cosmobl {

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> catalogues of astronomical sources </B>
   *  
   *  The \e catalogue namespace contains all the functions and
   *  classes used to handle catalogues of astronomical sources
   */
  namespace catalogue {

    /**
     *  @enum Var
     *  @brief the catalogue variables
     */
    enum Var {
    
      /// coordinate x
      _X_,
    
      /// coordinate y
      _Y_, 

      /// coordinate z
      _Z_,

      /// Right Ascension
      _RA_, 

      /// Declination
      _Dec_, 

      /// redshift
      _Redshift_, 

      /// comoving distance
      _Dc_, 

      /// weight
      _Weight_,

      /// mass
      _Mass_, 

      /// magnitude
      _Magnitude_,

      /// star formation rate
      _SFR_,

      /// specific star formation rate
      _sSFR_, 

      /// richness
      _Richness_,
      
      /// velocity along the x direction
      _Vx_, 

      /// velocity along the y direction
      _Vy_, 

      /// velocity along the z direction
      _Vz_, 

      /// region
      _Region_,
      
      /// radius 
      _Radius_,

      /// densityContrast
      _densityContrast_,

      /// centralDensity
      _centralDensity_,

      /// ID
      _ID_,
      
      /// xx displacement
      _X_displacement_,

      /// yy displacement
      _Y_displacement_,   
   
      /// zz displacement
      _Z_displacement_,

      /// generic property
      _Generic_
      
    };
    
    /**
     *  @enum RandomType
     *  @brief the type of random catalogue
     */
    enum RandomType {

      /// random catalogue with cubic geometry (or parallelepiped) in comoving coordinates
      _createRandom_box_,

      /// random catalogue with square geometry in observed coordinates (R.A., Dec)
      _createRandom_square_,   

      /// random catalogue obtained with shuffling in observed coordinates (R.A., Dec)
      _createRandom_shuffle_,
      
      /// random catalogue with conic geometry
      _createRandom_cone_,

      /// random catalogue using mangle
      _createRandom_MANGLE_,

      /// random catalogue for VIPERS
      _createRandom_VIPERS_
      
    };
    
    /**
     *  @enum CharEncode
     *  @brief character encoding of input file
     */
    enum CharEncode {
    
      /// Format ASCII file
      _ascii_,
      
      /// Format binary file
      _binary_
      
    };
    
    /**
     *  @class Catalogue Catalogue.h "Headers/Lib/Catalogue.h"
     *
     *  @brief The class Catalogue
     *
     *  This class is used to handle objects of type <EM> Catalogue
     *  </EM>
     */
    class Catalogue {
      
    private :
      
      /// vector containing the objects of the catalogue
      vector<shared_ptr<Object> > m_object;
      
      /// vector containing the object indexes
      vector<int> m_index;      

      
    public :
      
      /**
       *  @name Constructors/destructors
       */
      ///@{
    
      /**
       *  @brief default constructor
       *  @return object of class Catalogue
       */
      Catalogue () = default;

      /**
       *  @brief default copy constructor
       *  @param cat object of class Catalogue
       *  @return object of class Catalogue
       */
      Catalogue (const Catalogue &cat)
	: m_object(cat.sample()), m_index(cat.index()) {}
      
      /**
       *  @brief constructor
       *
       *  @param objType the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration 
       *
       *  @param coordType the coordinate type, specified in the
       *  cosmobl::CoordType enumeration
       *
       *  @param coord1 vector containing the first coordinates, that
       *  can be either the x comoving coordinates, or the Right
       *  Ascensions (depending on coordtype)
       *
       *  @param coord2 vector containing the second coordinates, that
       *  can be either the y comoving coordinates, or the
       *  Declinations (depending on coordtype)
       *
       *  @param coord3 vector containing the third coordinates, that
       *  can be either the z comoving coordinates, or the redshits
       *  (depending on coordtype)
       *
       *  @param weight vector containing the weights
       *
       *  @param cosm object of class Cosmology
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @return object of type catalogue
       */
      Catalogue (const ObjType objType, const CoordType coordType, const vector<double> coord1, const vector<double> coord2, const vector<double> coord3, const vector<double> weight={}, const cosmology::Cosmology &cosm={}, const CoordUnits inputUnits=_radians_);

      /**
       *  @brief constructor
       *
       *  @param objType the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration 
       *
       *  @param coordType the coordinate type, specified in the
       *  cosmobl::CoordType enumeration
       *
       *  @param coord1 vector containing the first coordinates, that
       *  can be either the x comoving coordinates, or the Right
       *  Ascensions (depending on coordtype)
       *
       *  @param coord2 vector containing the second coordinates, that
       *  can be either the y comoving coordinates, or the
       *  Declinations (depending on coordtype)
       *
       *  @param coord3 vector containing the third coordinates, that
       *  can be either the z comoving coordinates, or the redshits
       *  (depending on coordtype)
       *
       *  @param cosm object of class Cosmology
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @return object of type catalogue
       */
      Catalogue (const ObjType objType, const CoordType coordType, const vector<double> coord1, const vector<double> coord2, const vector<double> coord3, const cosmology::Cosmology &cosm, const CoordUnits inputUnits=_radians_)
	: Catalogue(objType, coordType, coord1, coord2, coord3, {}, cosm, inputUnits) {}
      
      /**
       *  @brief constructor, reading a file with coordinates
       *
       *  @param objType the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *
       *  @param coordType the coordinate type, specified in the
       *  cosmobl::CoordType enumeration
       *
       *  @param file vector containing the files where the input
       *  catalogues are stored
       *
       *  @param col1 column of the input file containing the first
       *  coordinates, that can be either the x comoving coordinates,
       *  or the Right Ascensions (depending on coordtype)
       *
       *  @param col2 column of the input file containing the second
       *  coordinates, that can be either the y comoving coordinates,
       *  or the Declinations (depending on coordtype)
       *
       *  @param col3 column of the input file containing the third
       *  coordinates, that can be either the z comoving coordinates,
       *  or the redshits (depending on coordtype)
       *
       *  @param colWeight column of the input file containing the
       *  weights
       *
       *  @param colRegion column of the input file containing the
       *  regions (used for jackknife or bootstrap)
       *
       *  @param nSub the fracton of objects that will be randomly
       *  selected (nSub=1 &rArr; all objects are selected)
       *
       *  @param fact a factor used to multiply the coordinates,
       *  i.e. coordinate_i=coordinate_i*fact
       *
       *  @param cosm object of class Cosmology 
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param charEncode character encoding of input file,
       *  ascii or binary
       * 
       *  @return an object of class Catalogue
       */
      Catalogue (const ObjType objType, const CoordType coordType, const vector<string> file, const int col1=1, const int col2=2, const int col3=3, const int colWeight=-1, const int colRegion=-1, const double nSub=1.1, const double fact=1., const cosmology::Cosmology &cosm={}, const CoordUnits inputUnits=_radians_, const CharEncode charEncode=_ascii_);

      /**
       *  @brief constructor, reading a file with coordinates
       *
       *  @param objType the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *
       *  @param coordType the coordinate type, specified in the
       *  cosmobl::CoordType enumeration
       *
       *  @param file vector containing the files where the input
       *  catalogues are stored
       *
       *  @param cosm object of class Cosmology 
       *
       *  @param inputUnits the units of the input coordinates
       * 
       *  @return an object of class Catalogue
       */
      Catalogue (const ObjType objType, const CoordType coordType, const vector<string> file, const cosmology::Cosmology &cosm, const CoordUnits inputUnits=_radians_)
	: Catalogue(objType, coordType, file, 1, 2, 3, -1, -1, 1.1, 1., cosm, inputUnits, _ascii_) {}

      /**
       *  @brief constructor, using vectors of generic objects
       *  @param object objects of class T, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *  @return objects of type Catalogue
       */ 
      template<typename T> Catalogue (vector<T> object) {
	for (size_t i=0; i<object.size(); i++)
	  m_object.push_back(move(make_shared<T>(T(object[i]))));
      }

      /**
       *  @brief constructor, using vectors of pointers to generic
       *  objects
       *  @param sample vector of objects of type \e Object, specified
       *  in the cosmobl::catalogue::ObjType enumeration
       *  @return object of class Catalogue
       */
      Catalogue (vector<shared_ptr<Object> > sample) {
	for (auto &&i : sample)
	  m_object.push_back(move(i));
      }

      /**
       *  @brief constructor, creating a catalogue by matching the
       *  distribution of one quantity from a target catalogue
       *
       *  @param input_catalogue the input catalogue
       *
       *  @param target_catalogue the target catalogue
       *
       *  @param var_name the type of variable, specified
       *  cosmobl::catalogue::Var enumeration
       *
       *  @param nbin the binning for the variable
       *
       *  @param seed the seed for random number generation
       *
       *  @return an object of class Catalogue
       */
      Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const Var var_name, const int nbin, const int seed=3213);

      /**
       *  @brief constructor, creating a catalogue by matching the
       *  distributions of two quantities from a target catalogue
       *
       *  @param input_catalogue the input catalogue
       *
       *  @param target_catalogue the target catalogue
       *
       *  @param var_name1 the type of variable, specified
       *  cosmobl::catalogue::Var enumeration
       *
       *  @param nbin1 the binning for the variable
       *
       *  @param var_name2 the type of variable, specified
       *  cosmobl::catalogue::Var enumeration
       *
       *  @param nbin2 the binning for the variable
       *
       *  @param seed the seed for random number generation
       *
       *  @return an object of class Catalogue
       */
      Catalogue (const Catalogue input_catalogue, const Catalogue target_catalogue, const cosmobl::catalogue::Var var_name1, const int nbin1, const cosmobl::catalogue::Var var_name2, const int nbin2, const int seed=3213);

      /**
       * @brief default destructor
       * @return none
       */
      ~Catalogue () = default;

      ///@}

      
      /**
       *  @name Constructors of random catalogues 
       */
      ///@{
      
      /**
       *  @brief constructor that creates a random catalogue in a
       *  cubic box, warped by geometric distortions
       *
       *  this function reads a cubic random catalogue from a file,
       *  generated in a given cosmology, and trasforms it into a new
       *  one in a different cosmology
       *
       *  @param type the type of random catalogue, that must be set
       *  to \_createRandom_box\_
       *
       *  @param real_cosm object of class Cosmology representing the \e
       *  real (or \e assumed) cosmology
       *
       *  @param test_cosm object of class Cosmology representing the \e
       *  test cosmology
       *
       *  @param dir_in the input directory where the original random
       *  catalogue is stored
       *     
       *  @param Zguess_min minimum redshift used to search the redshift
       *
       *  @param Zguess_max maximum redshift used to search the redshift
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const cosmology::Cosmology &real_cosm, const cosmology::Cosmology &test_cosm, const string dir_in, const double Zguess_min, const double Zguess_max);

      /**
       *  @brief constructor that creates a random catalogue with
       *  either the square geometry or with 'the shuffle' method
       *
       *  constructor that creates a random catalogue with
       *  either the square geometry in observed coordinates (R.A.,
       *  Dec), or with the 'shuffle' method, i.e. using the R.A. and
       *  Dec coordinates of the input catalogue
       *
       *  @param type the type of random catalogue, that must be set
       *  to either \_createRandom_box\_, \_createRandom_square\_ or
       *  \_createRandom_shuffle\_
       *
       *  @param catalogue object of class Catalogue
       *
       *  @param N_R fraction of random objects, i.e.
       *  N<SUB>R</SUB>=N<SUB>random</SUB>/N<SUB>objects</SUB>
       *
       *  @param nbin number of redshift bins used to compute the
       *  redshift distribution
       *
       *  @param cosm object of class Cosmology
       *
       *  @param conv 1 &rarr; compute the Gaussian convolvolution of
       *  the distribution; 0 &rarr; do not convolve
       *
       *  @param sigma the standard deviation, &sigma;, of the
       *  Gaussian kernel
       *
       *  @param seed the seed for random number generation
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin=10, const cosmology::Cosmology &cosm={}, const bool conv=false, const double sigma=0., const int seed=3213);
      
      /**
       *  @brief constructor that creates a random catalogue in a cone
       *
       *  @param type the type of random catalogue, that must be
       *  set to \_createRandom_cone\_
       *
       *  @param catalogue object of class Catalogue
       *
       *  @param N_R fraction of random objects, i.e.
       *  N<SUB>R</SUB>=N<SUB>random</SUB>/N<SUB>objects</SUB>
       *
       *  @param nbin number of redshift bins used to compute the
       *  redshift distribution
       *
       *  @param Angle angle of the cone 
       *
       *  @param redshift vector containing the redshift of the
       *  objects in the catalogue
       *
       *  @param cosm object of class Cosmology 
       *
       *  @param conv 1 &rarr; compute the Gaussian convolvolution of
       *  the distribution; 0 &rarr; do not convolve
       *
       *  @param sigma the standard deviation, &sigma;, of the
       *  Gaussian kernel
       *
       *  @param seed the seed for random number generation
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const double N_R, const int nbin, const double Angle, const vector<double> redshift, const cosmology::Cosmology &cosm={}, const bool conv=false, const double sigma=0., const int seed=3213);

      /**
       *  @brief constructor that creates a random catalogue using the 
       *  a mask in the MANGLE format for the angular distribution
       *  and taking the redshift distribution from an input catalogue
       * 
       *  @param type the type of random catalogue, that must be
       *  set to \_createRandom_MANGLE\_
       *
       *  @param mangle_mask vector containing the input masks in MANGLE format
       *
       *  @param catalogue object of class Catalogue
       *
       *  @param N_R fraction of random objects, i.e.
       *  N<SUB>R</SUB>=N<SUB>random</SUB>/N<SUB>objects</SUB>
       *
       *  @param nbin number of redshift bins used to compute the
       *  redshift distribution
       *
       *  @param cosm object of class Cosmology
       *
       *  @param conv 1 &rarr; compute the Gaussian convolvolution of
       *  the distribution; 0 &rarr; do not convolve
       *
       *  @param sigma the standard deviation, &sigma;, of the
       *  Gaussian kernel
       *
       *  @param seed the seed for random number generation
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const vector<string> mangle_mask, const Catalogue catalogue, const double N_R, const int nbin, const cosmology::Cosmology cosm, const bool conv=false, const double sigma=0., const int seed=3213);   

      /// @cond extrandom
      
      Catalogue (const RandomType type, const string WField, const bool isSpectroscopic, const Catalogue catalogue, const Catalogue catalogue_for_nz, const double N_R, const cosmology::Cosmology &cosm, const int step_redshift, const vector<double> lim, const double redshift_min, const double redshift_max, const bool do_convol, const double sigma, const bool use_venice, const bool do_zdistr_with_venice, const string file_random, const string mask, const string pointing_file, const string dir_venice, const int seed); 
      
      /// @endcond
  
      ///@}

    
      /**
       *  @name Member functions used to get the private members and thier properties
       */
      ///@{

      /**
       *  @brief get the private member Catalogue::m_object
       *  @return the vector containing the objects of the catalogue
       */
      vector<shared_ptr<Object> > sample () const { return m_object; };
      
      /**
       *  @brief get the private member Catalogue::m_index
       *  @return the vector containing the object indexes
       */
      vector<int> index () const { return m_index; };
      
      /**
       * @brief get the private member Catalogue::m_object[i]->m_xx
       * @param i the object index
       * @return the coordinate x of the i-th object 
       */
      double xx (const int i) const { return m_object[i]->xx(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_yy
       * @param i the object index
       * @return the coordinate y of the i-th object 
       */
      double yy (const int i) const { return m_object[i]->yy(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_zz
       * @param i the object index
       * @return the coordinate z of the i-th object 
       */
      double zz (const int i) const { return m_object[i]->zz(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_vx
       * @param i the object index
       * @return the velocity along the x direction of the i-th object
       */
      double vx (const int i) const { return m_object[i]->vx(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_vy
       * @param i the object index
       * @return the velocity along the y direction of the i-th object
       */
      double vy (const int i) const { return m_object[i]->vy(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_vz
       * @param i the object index
       * @return the velocity along the z direction of the i-th object
       */
      double vz (const int i) const { return m_object[i]->vz(); }; 

      /**
       * @brief get the private member Catalogue::m_object[i]->m_dc
       * @param i the object index
       * @return the comoving distance of the i-th object
       */
      double dc (const int i) const { return m_object[i]->dc(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_ra
       * @param i the object index
       * @return the Right Ascension of the i-th object
       */
      double ra (const int i) const { return m_object[i]->ra(); };
    
      /**
       * @brief get the private member Catalogue::m_object[i]->m_dec
       * @param i the object index
       * @return the Declination of the i-th object
       */
      double dec (const int i) const { return m_object[i]->dec(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_redshift
       * @param i the object index
       * @return the redshift of the i-th object
       */
      double redshift (const int i) const { return m_object[i]->redshift(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_weight
       * @param i the object index
       * @return the weight of the i-th object
       */
      double weight (const int i) const { return m_object[i]->weight(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_region
       * @param i the object index
       * @return the index of the region of the i-th object
       */
      long region (const int i) const { return m_object[i]->region(); };

      /**
       * @brief get the private member Catalogue::m_object[i]->m_field
       * @param i the object index
       * @return the field where the i-th object has been observed
       */
      string field (const int i) const { return m_object[i]->field(); };

      /**
       * @brief get the private member
       * Catalogue::m_object[i]->m_x_displacement
       * @param i the object index
       * @return the displacement of the i-th object along the x-axis
       */
      double x_displacement (const int i) const { return m_object[i]->x_displacement(); };

      /**
       * @brief get the private member
       * Catalogue::m_object[i]->m_y_displacement
       * @param i the object index
       * @return the displacement of the i-th object along the x-axis
       */
      double y_displacement (const int i) const { return m_object[i]->y_displacement(); };

      /**
       * @brief get the private member
       * Catalogue::m_object[i]->m_z_displacement
       * @param i the object index
       * @return the displacement of the i-th object along the x-axis
       */
      double z_displacement (const int i) const { return m_object[i]->z_displacement(); };

      /**
       *  @brief get the list of regions in which the catalogue is
       *  divided
       *
       *  @return the list of regions 
       */
      vector<long> region_list () const { return different_elements(region()); }

      /**
       *  @brief get the total number of regions by which the
       *  Catalogues is divided
       *
       *  @return the total number of regions
       */
      size_t nRegions () const { return N_different_elements(region()); }
      
      /**
       *  @brief get the list of fields where the objects have been
       *  observed
       *
       *  @return the list of fields 
       */
      vector<string> field_list () const { return different_elements(field()); }
      
      /**
       *  @brief get the total number of fields where the objects have
       *  been observed
       *
       *  @return the total number of fields
       */
      size_t nFields () const { return N_different_elements(field()); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_mass
       * @param i the object index
       * @return the mass of the i-th object
       */
      double mass (const int i) const { return m_object[i]->mass(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_magnitude
       * @param i the object index
       * @return the magnitude of the i-th object
       */
      double magnitude (const int i) const { return m_object[i]->magnitude(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_radius
       * @param i the object index
       * @return radius of the i-th object
       */
      double radius (const int i) const { return m_object[i]->radius(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_densityContrast
       * @param i the object index
       * @return density contrast of the i-th object
       */
      double densityContrast (const int i) const { return m_object[i]->densityContrast(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_centralDensity
       * @param i the object index
       * @return central density of the i-th object
       */
      double centralDensity (const int i) const { return m_object[i]->centralDensity(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_ID
       * @param i the object index
       * @return ID of the i-th object
       */
      int ID (const int i) const { return m_object[i]->ID(); }
      
      /**
       * @brief get the private member
       * Catalogue::m_object[i]->m_richness
       * @param i the object index
       * @return the richness of the i-th object
       */
      double richness (const int i) const { return m_object[i]->richness(); }

      /**
       * @brief get the private member Catalogue::m_object[i]->m_generic
       * @param i the object index
       * @return generic properties of the i-th object
       */
      double generic (const int i) const { return m_object[i]->generic(); }

      /**
       * @brief get the values of the object regions  
       * @return the object regions
       */
      vector<long> region () const;

       /**
       * @brief get the values of the object fields  
       * @return the object fields
       */
      vector<string> field () const;
      
      /**
       * @brief get the values of the object variables  
       * @param var_name the variable name
       * @return the vector of the variable Var
       */
      vector<double> var (const Var) const;

      /**
       * @brief get the i-th object of the catalogue
       * @param i the object index
       * @return pointer to an object of the catalogue
       */
      shared_ptr<Object> catalogue_object (const int i) const { return m_object[i]; }

      /**
       * @brief get the object vector
       * @return vector of pointers to objects of the catalogue
       */
      vector<shared_ptr<Object>> catalogue_object () const { return m_object; }   

      /**
       * @brief get the X, Y, Z coordinates of the i-th object of the
       * catalogue
       *
       * @param i the object index
       * @return vector containing the three coordinates
       */
      vector<double> coordinate (const int i) const { return m_object[i]->coords(); }
    
      /**
       * @brief get the number of objects of the catalogue
       * @return the number of objects
       */
      size_t nObjects () const { return m_object.size(); }
  
      /**
       * @brief get the minimum value of a variable of the catalogue
       * objects
       * @param var_name the variable name
       * @return the minimum value of the variable
       */
      double Min (const Var var_name) const { return cosmobl::Min(var(var_name)); }

      /**
       * @brief get the maximum value of a variable of the catalogue
       * objects
       * @param var_name the variable name
       * @return the maximum value of the variable
       */
      double Max (const Var var_name) const { return cosmobl::Max(var(var_name)); }

      /**
       * @brief get the mean, the median, the standard deviation, and
       * the difference between the third and first quartiles of a
       * variable
       * @param [in] var_name the variable name
       * @param [out] stats 4 dimensional vector containing the mean,
       * the median, the standard deviation, and the difference between
       * the third and first quartiles of the variable
       * @return none
       */
      void stats_var (const Var var_name, vector<double> & stats) const;

      /**
       * @brief get the mean, the median, the standard deviation, and
       * the difference between the third and first quartiles of a
       * variable
       * @param [in] var_name vector of variable names
       * @param [out] stats vector of 4 dimensional vector containing
       * the mean, the median, the standard deviation, and the
       * difference between the third and first quartiles of the
       * variable 
       * @return none
       */
      void stats_var (const vector<Var> var_name, vector<vector<double>> &stats) const;
  
      /**
       * @brief get the distribution of a variable
       * @param [in] var_name the variable name
       * @param [out] _var vector of variables
       * @param [out] dist vector of values of f(varibles)
       * @param [out] err vector of Poissonian errors of f(varibles)
       * @param [in] nbin number of bins
       * @param [in] linear 1 &rarr; linear binning; 0 &rarr;
       * logarithmic binning
       * @param [in] file_out the output file where the distribution is
       * stored
       * @param [in] Volume the volume of the catalogue
       * @param [in] norm 1 &rarr; normalize to the number of objects;
       * 0 &rarr; do not normalize
       * @param [in] V1 the minimum limit of the distribution
       * @param [in] V2 the maximum limit of the distribution
       * @param [in] bin_type 1 &rarr; dn/dvar; 0 &rarr; dn/dlogvar; 
       * @param [in] convolution 0 &rarr; don't convolve the
       * distribution; 1 &rarr; convolve the distribution with a
       * gaussian function
       * @param [in] sigma &sigma;: the standard deviation of the
       * gaussian function used to convolve the distribution
       * @return none
       */
      void var_distr (const Var var_name, vector<double> &_var, vector<double> &dist, vector<double> &err, const int nbin, const bool linear=1, const string file_out=par::defaultString, const double Volume=1., const bool norm=0, const double V1=par::defaultDouble, const double V2=par::defaultDouble, const bool bin_type=1, const bool convolution=0, const double sigma=0.) const;
    
      /**
       * @brief get the total weight of the objects of the catalogue
       * @return the total weight
       */
      double weightedN () const;
    
      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{

      /**
       * @brief set a private variable
       * @param region vector containing the object regions
       * @return none
       */
      void set_region (const vector<long> region);

      /**
       * @brief set a private variable
       * @param field vector containing the object fields
       * @return none
       */
      void set_field (const vector<string> field);
      
      /**
       * @brief set a private variable
       * @param var_name name of the variable
       * @param var vector of variables
       * @return none
       */
      void set_var (const Var var_name, const vector<double> var);

      ///@}

    
      /**
       *  @name Member functions used to add or remove ojects to the catalogue
       */
      ///@{
    
      /**
       * @brief add one single object to the catalogue
       * @param object pointer to an object of type \e Object
       * @return none
       */
      void add_object (shared_ptr<Object> object) { m_object.push_back(move(object)); }

      /**
       * @brief add one single object to the catalogue
       * @param object object of type \e T
       * @return none
       */
      template<typename T>
	void add_object (T object) { m_object.push_back(move(make_shared<T>(T(object)))); }

      /**
       * @brief add some objects to the catalogue
       * @param sample vector of pointers to objects of type \e Object
       * @return none
       */
      void add_objects (vector<shared_ptr<Object> > sample) { 
	for (auto &&i : sample)
	  m_object.push_back(move(i));
      }

      /**
       * @brief add some objects to the catalogue
       * @param sample vector of objects of type \e T
       * @return none
       */
      template<typename T>
	void add_objects (vector<T> sample) { 
	for (auto &&i : sample)
	  add_object(i);
      }
    
      /**
       * @brief replace existing objects with new ones 
       * @param sample vector of objects of type \e T
       * @return none
       */
      template<typename T>
	void replace_objects(vector<T> sample) {
	m_object.erase(m_object.begin(), m_object.end());
	add_objects(sample);
      }

      /**
       * @brief replace existing objects with new ones 
       * @param sample vector of pointers to objects of type \e Object
       * @return none
       */
      void replace_objects (vector<shared_ptr<Object> > sample) {
	m_object.erase(m_object.begin(), m_object.end());
	for (auto &&i : sample)
	  m_object.push_back(move(i));
      }

      /**
       * @brief remove all objects 
       * @return none
       */
      void remove_objects () { m_object.erase(m_object.begin(), m_object.end()); }
      
      /**
       * @brief remove an existing object
       * @param index the index of the object to be removed
       * @return none
       */
      void remove_object (const int index) { m_object.erase(m_object.begin()+index); }
      
      ///@}

    
      /**
       *  @name Member functions used to operate on the ojects to the catalogue
       */
      ///@{
    
      /**
       *  @brief compute the comoving coordinates (x, y, z) from the
       *  observed coordinates (R.A., Dec, redshift)
       *
       *  @param cosm object of class Cosmology
       *  @param inputUnits the units of the input coordinates
       *  @return none
       */
      void computeComovingCoordinates (const cosmology::Cosmology &cosm, const CoordUnits inputUnits=_radians_); 

      /**
       *  @brief compute the polar coordinates (R.A., Dec,
       *  d<SUB>c</SUB>) from the comoving coordinates (x, y, z)
       *
       *  @param outputUnits the units of the output coordinates
       *  @return none
       */
      void computePolarCoordinates (const CoordUnits outputUnits=_radians_); 

      /**
       *  @brief compute the polar coordinates (R.A., Dec,
       *  d<SUB>c</SUB>, redshift) from the comoving (x, y, z), and
       *  assuming a cosmological model
       *
       *  @param cosm object of class Cosmology
       *  @param z1 the minimum redshift used in the computation
       *  @param z2 the maximum redshift used in the computation
       *  @param outputUnits the units of the output coordinates
       *  @return none
       */
      void computePolarCoordinates (const cosmology::Cosmology &cosm, const double z1=0., const double z2=10., const CoordUnits outputUnits=_radians_); 

      /**
       * @brief normalize (x, y, z) (i.e. &rarr; (x/dc, y/dc, z/dc))
       * @return none
       */
      void normalizeComovingCoordinates (); 
  
      /**
       * @brief back to (x, y, z) (i.e. the inverse of norm_xyv())
       * @return none
       */
      void restoreComovingCoordinates ();  

      /**
       * @brief order the catalogue according to the input vector
       * @param vv vector used to order the catalogue
       * @return none
       */
      void Order (const vector<int> vv); 

      /**
       * @brief restore the original vector (i.e. the opposite of
       * Order(vector<int>))
       * @return none
       */
      void Order ();  
    
      /**
       * @brief write the comoving coordinates of the catalogue to an
       * output file
       * @param outputFile the name of the output file
       * @return none
       */
      void write_comoving_coordinates (const string outputFile) const;

      /**
       * @brief write the polar coordinates of the catalogue to an
       * output file
       * @param outputFile the name of the output file
       * @return none
       */
      void write_obs_coordinates (const string outputFile) const;

      /**
       *  @brief write both the comoving and polar coordinates, and the
       *  regions (if present) of the catalogue to an output file
       *
       *  @param outputFile the name of the output file
       *
       *  @param var_name vector containing the variable names to be
       *  written
       *
       *  @return none
       */
      void write_data (const string outputFile, const vector<Var> var_name={}) const;
      
      /**
       * @brief get the distrance between the i-th object of the
       * catalogue and another object
       * @param i the object index
       * @param obj pointer to an object
       * @return distance between the i-th object of the catalogue and
       * the object obj
       */
      double distance (const int i, shared_ptr<Object> obj) const;
    
      /**
       * @brief get the angular distrance between the i-th object of the
       * catalogue and another object
       * @param i the object index
       * @param obj pointer to an object
       * @return distance between the i-th object of the catalogue and
       * the object obj
       */
      double angsep_xyz (const int i, shared_ptr<Object> obj) const;
    
      /**
       * @brief overloading of the += operator, to sum two catalogues
       * @param cc object of class Catalogue 
       * @return object of class catalogue
       */
      Catalogue operator += (shared_ptr<Catalogue> cc)
      {
	for (auto &&ss : cc->m_object)
	  m_object.push_back(shared_ptr<Object>(new Object(*ss)));
	return *this;
      }

      /**
       * @brief overloading of the += operator, to sum two catalogues
       * @param cc object of class Catalogue 
       * @return object of class catalogue
       */
      Catalogue operator += (const Catalogue cc)
      {
	for (auto &&ss : cc.m_object)
	  m_object.push_back(shared_ptr<Object>(new Object(*ss)));
	return *this;
      }

      /**
       * @brief create a sub-catalogue
       * @param var_name the variable name
       * @param down minimum variable used to cut the catalogue
       * @param up maximum variable used to cut the catalogue
       * @param excl 0 &rarr; creates a subcatalogue inside down-up; 1
       * &rarr; creates a subcatalogue outside down-up;
       * @return object of class catalogue
       */
      Catalogue cut (const Var, const double, const double, const bool excl=0);

      /**
       * @brief create a smoothed version of the catalogue
       * averaging quantities on a X, Y, Z grid 
       *
       * defalut averaged quantities are X, Y, Z, RA, DEC, REDSHIFT,
       * WEIGHT; others quantities must be passed trough a vector
       *
       * @param gridsize the cell size 
       * @param vars the vector of variable to average on
       * @param SUB the number of sub-catalogue used to create the
       * chain-mesh (use SUB>1 when there could be memory problems)
       * @return object of class catalogue
       */
      shared_ptr<Catalogue> smooth (const double, const vector<Var> vars={}, const int SUB=1);

      /**
       * @brief return the number of objectes following a condition
       * on the variable VAR
       * @param var_name the variable name
       * @param down minimum variable used to cut the catalogue
       * @param up maximum variable used to cut the catalogue
       * @param excl 0 &rarr; count objects inside down-up; 1
       * &rarr; count objects outside down-up;
       * @return number of objects following the condition
       */
      int nObjects_condition (const Var, const double, const double, const bool excl=0);

      /**
       * @brief return the weighted number of objectes following a condition
       * on the variable VAR
       * @param var_name the variable name
       * @param down minimum variable used to cut the catalogue
       * @param up maximum variable used to cut the catalogue
       * @param excl 0 &rarr; count objects inside down-up; 1
       * &rarr; count objects outside down-up;
       * @return weighted number of objects following the condition
       */
      double weightedN_condition (const Var, const double, const double, const bool excl=0);

      /**
       * @brief return the density field from object position
       * @param cell_size the minimum size of the density field
       * @param interpolation_type the type of interpolation 0 &rarr; nearest-grid-point;
       * 1 &rarr; cloud-in-cell
       * @param kernel_radius size of the kernel for the gaussian smoothing
       * @param useMass generate the density field using the mass information
       * @return the density field
       */
      data::ScalarField3D density_field (const double cell_size, const int interpolation_type=0, const double kernel_radius=0., const bool useMass=0) const;

      /**
       * @brief return the density field from object position
       * @param cell_size the minimum size of the density field
       * @param mask_catalogue catalogue containing points sampling the selecion function
       * of the catalogue
       * @param interpolation_type the type of interpolation 0 &rarr; nearest-grid-point;
       * 1 &rarr; cloud-in-cell
       * @param kernel_radius size of the kernel for the gaussian smoothing
       * @param useMass generate the density field using the mass information
       * @return the density field
       */
      data::ScalarField3D density_field (const double cell_size, const Catalogue mask_catalogue, const int interpolation_type = 0, const double kernel_radius=0., const bool useMass = 0) const;

      ///@}
      
    };
    
  }
}

#endif

