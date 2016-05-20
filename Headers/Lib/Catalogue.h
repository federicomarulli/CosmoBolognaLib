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
 * @file Headers/Lib/Catalogue.h
 *
 * @brief The class Catalogue 
 *
 * This file defines the interface of the class Catalogue, used
 * handle catalogues of astronomical sources
 *
 * @authors Federico Marulli, Alfonso Veropalumbo 
 *
 * @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
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
   *  @brief The namespace of the catalogue 
   *  
   *  The \e catalogue namespace contains all the functions and
   *  classes used to handle catalogues of astronomical sources
   */
  namespace catalogue {

    /**
     * @enum Var
     * @brief the catalogue variables
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

      /// richness
      _Richness_, 

      /// magnitude
      _Magnitude_, 

      /// velocity along the x direction
      _Vx_, 

      /// velocity along the y direction
      _Vy_, 

      /// velocity along the z direction
      _Vz_, 

      /// region
      _Region_,

      /// radius properties
      _Radius_,
      
      /// generic properties
      _Generic_,
    };

    
    /**
     * @enum RandomType
     * @brief the type of random catalogue
     */
    enum RandomType {

      /// random catalogue with cubic geometry (or parallelepiped)
      _createRandom_box_,

      /// random catalogue with square geometry (or parallelepiped)
      _createRandom_box_RADec_,   

      /// random catalogue obtained with shuffling
      _createRandom_shuffle_,
      
      /// random catalogue with conic geometry
      _createRandom_cone_,

      /// random catalogue for mocks
      _createRandom_mock_,

      /// random catalogue for VIPERS
      _createRandom_VIPERS_
      
    };
    
    /**
     * @class Catalogue Catalogue.h "Headers/Lib/Catalogue.h"
     *
     * @brief The class Catalogue
     *
     * This class is used to handle objects of type <EM> Catalogue
     * </EM>
     */
    class Catalogue {

    private :
    
      /// vector containing the objects of the catalogue
      vector<shared_ptr<Object> > m_sample;
    
      /// vector containing the object indexes
      vector<int> m_index;      

    public :

      /**
       *  @name Constructors/destructors
       */
      ///@{
    
      /**
       * @brief default constructor
       * @return object of class Catalogue
       */
      Catalogue () = default;


      /**
       * @brief default copy constructor
       * @param cat object of class Catalogue
       * @return object of class Catalogue
       */
      Catalogue (const Catalogue& cat);
      
      /**
       * @brief constructor, using vectors with Cartesian coordinates
       * @param type the object type, specified in the
       * cosmobl::catalogue::ObjType enumeration
       * @param xx vector containing the x coordinates
       * @param yy vector containing the y coordinates
       * @param zz vector containing the z coordinates
       * @param weight vector containing the weights
       * @return object of type catalogue
       */
      Catalogue (const ObjType type, const vector<double> xx, const vector<double> yy, const vector<double> zz, vector<double> weight={});

      /**
       *  @brief constructor, using vectors with polar coordinates
       *  @param type the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *  @param ra vector containing the Right Ascensions
       *  @param dec vector containing the Declinations
       *  @param redshift vector containing the redshifts
       *  @param cosm object of class Cosmology, used to estimate comoving distances 
       *  @param inputUnits the units of the input coordinates
       *  @param weight vector containing the weights
       *  @return object of type catalogue
       */ 
      Catalogue (const ObjType type, const vector<double> ra, const vector<double> dec, const vector<double> redshift, const Cosmology &cosm, const CoordUnits inputUnits=_radians_, vector<double> weight={}); 

      /**
       * @brief constructor, using vectors of generic objects
       * @param object objects of class T, specified in the
       * cosmobl::catalogue::ObjType enumeration
       * @return objects of type Catalogue
       */ 
      template<typename T> Catalogue (vector<T> object) {
	for (size_t i=0; i<object.size(); i++)
	  m_sample.push_back(move(make_shared<T>(T(object[i]))));
      }

      /**
       * @brief constructor, using vectors of pointers to generic
       * objects
       * @param sample vector of objects of type \e Object, specified
       * in the cosmobl::catalogue::ObjType enumeration
       * @return object of class Catalogue
       */
      Catalogue (vector<shared_ptr<Object> > sample) {
	for (auto &&i : sample)
	  m_sample.push_back(move(i));
      }

      /**
       *  @brief constructor, reading a file with Cartesian coordinates
       *  @param type the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *  @param file vector containing the files where the input
       *  catalogues are stored
       *  @param col_X column of the input file containing the X
       *  coordinates
       *  @param col_Y column of the input file containing the Y
       *  coordinates
       *  @param col_Z column of the input file containing the Z
       *  coordinates
       *  @param col_Weight column of the input file containing the
       *  weights
       *  @param nSub the fracton of objects that will be randomly
       *  selected (nSub=1 &rArr; all objects are selected)
       *  @return an object of class Catalogue
       */
      Catalogue (const ObjType type, const vector<string> file, const int col_X=0, const int col_Y=1, const int col_Z=2, const int col_Weight=-1, const double nSub=1.1);

      /**
       *  @brief constructor, reading a file with polar coordinates [R.A., Dec,
       *  redshift]
       *  @param type the object type, specified in the
       *  cosmobl::catalogue::ObjType enumeration
       *  @param file vector containing the files where the input
       *  catalogues are stored
       *  @param cosm object of class Cosmology
       *  @param col_RA column of the input file containing the R.A.
       *  coordinates
       *  @param col_Dec column of the input file containing the Dec
       *  coordinates
       *  @param col_redshift column of the input file containing the
       *  redshift coordinates
       *  @param col_Weight column of the input file containing the
       *  weights
       *  @param inputUnits the units of the input coordinates
       *  @param nSub the fracton of objects that will be randomly
       *  selected (nSub=1 &rArr; all objects are selected)
       *  @param fact the factor used to convert R.A. and Dec
       *  coordinates, i.e. &rArr; R.A.=R.A.*fact, Dec=Dec*fact
       *  @return an object of class Catalogue
       */
      Catalogue (const ObjType type, const vector<string> file, const Cosmology &cosm, const int col_RA=0, const int col_Dec=1, const int col_redshift=2, const int col_Weight=-1, const CoordUnits inputUnits=_radians_, const double nSub=1.1, const double fact=1.);

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
       *  cubic box
       *  @param type the type of random catalogue, that must be set
       *  to \_createRandom_box\_
       *  @param catalogue object of class Catalogue
       *  @param N_R fraction of random objects, i.e.
       *  N<SUB>R</SUB>=N<SUB>random</SUB>/N<SUB>objects</SUB>
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const double N_R);

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
       *  @param dir_out the output directory where the new random
       *  catalogue will be stored
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
      Catalogue (const RandomType type, const Cosmology &real_cosm, const Cosmology &test_cosm, const string dir_in, const string dir_out, const double Zguess_min, const double Zguess_max);

      /**
       *  @brief constructor that creates a random catalogue with the
       *  'shuffle' method
       *  @param type the type of random catalogue, that must be set
       *  to \_createRandom_shuffle\_
       *  @param catalogue object of class Catalogue
       *  @param cosm object of class Cosmology
       *  @param N_R fraction of random objects, i.e.
       *  N<SUB>R</SUB>=N<SUB>random</SUB>/N<SUB>objects</SUB>
       *  @param nbin numbers of bin to obtain the redshift distribution
       *  @param conv 1 &rarr; compute the Gaussian convolvolution of
       *  the distribution; 0 &rarr; do not convolve
       *  @param sigma the standard deviation, &sigma;, of the
       *  Gaussian kernel
       *  @param seed the seed for random number generation
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const Cosmology &cosm, const int N_R, const int nbin, const bool conv=0, const double sigma=0., const int seed = 3213);
      
      /**
       *  @brief constructor that creates a random catalogue in a cone
       *
       *  @param [in] type the type of random catalogue, that must be
       *  set to \_createRandom_cone\_
       *
       *  @param [in] catalogue object of class Catalogue
       *
       *  @param [in] nRandom the number of random objects
       *
       *  @param [in] cosm object of class Cosmology 
       *
       *  @param [in] Angle angle of the cone 
       *
       *  @param [in] step_redshift the number of steps in redshift
       *  used to redshift distribution of the random object; if
       *  step_redshift=0 the redshift distribution is estimated from
       *  the convolvolution of N(D<SUB>C</SUB>)
       *
       *  @param [in] redshift vector containing the redshift of the
       *  object in the real catalogue
       *
       *  @param [out] dc vector containing the central values of the
       *  binned comoving distances of the random objects
       *
       *  @param [out] convol vector containing the central values of
       *  the binned smoothed distribution of comoving distances of
       *  the random objects
       *
       *  @param [in] idum the random seed
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const int nRandom, const Cosmology &cosm, const double Angle, const int step_redshift, const vector<double> redshift, vector<double> &dc, vector<double> &convol, const int idum=13);
      
      /**
       *  @brief constructor that creates a random catalogue for a
       *  mock sample (with polar coordinates [R.A., Dec, redshift])
       *
       *  @param [in] type the type of random catalogue, that must be
       *  set to \_createRandom_mock\_
       *
       *  @param [in] catalogue object of class Catalogue
       *
       *  @param [in] nRandom the number of random objects
       *
       *  @param [in] cosm object of class Cosmology 
       *
       *  @param [in] dir the directory where the random catalogue is
       *  stored
       *
       *  @param [in] step_redshift the number of steps in redshift
       *  used to redshift distribution of the random object; if
       *  step_redshift=0 the redshift distribution is estimated from
       *  the convolvolution of N(D<SUB>C</SUB>)
       *
       *  @param [in] redshift vector containing the redshift of the
       *  object in the real catalogue
       *
       *  @param [out] dc vector containing the central values of the
       *  binned comoving distances of the random objects
       *
       *  @param [out] convol vector containing the central values of
       *  the binned smoothed distribution of comoving distances of
       *  the random objects
       *
       *  @param [in] idum the random seed
       *
       *  @return an object of class Catalogue
       *
       *  @warning the input parameter \e type is used only to make
       *  the constructor type explicit
       */
      Catalogue (const RandomType type, const Catalogue catalogue, const int nRandom, const Cosmology &cosm, const string dir, const int step_redshift, const vector<double> redshift, vector<double> &dc, vector<double> &convol, const int idum=13);

      /// @cond extrandom
      
      Catalogue (const RandomType type, const int nRandom, const Cosmology &cosm, const string dir_out, const int step_redshift, const vector<double> dc, const vector<double> convol, const vector<double> lim, const vector<double> redshift, const bool venice, string file_random, const string mask, const string dir_venice, const int idum=13); 
      
      /// @endcond
  
      ///@}

    
      /**
       *  @name Member functions used to get the protected members and thier properties
       */
      ///@{
    
      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_xx
       * @param i the object index
       * @return the coordinate x of the i-th object 
       */
      double xx (const int i) const { return m_sample[i]->xx(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_yy
       * @param i the object index
       * @return the coordinate y of the i-th object 
       */
      double yy (const int i) const { return m_sample[i]->yy(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_zz
       * @param i the object index
       * @return the coordinate z of the i-th object 
       */
      double zz (const int i) const { return m_sample[i]->zz(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_vx
       * @param i the object index
       * @return the velocity along the x direction of the i-th object
       */
      double vx (const int i) const { return m_sample[i]->vx(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_vy
       * @param i the object index
       * @return the velocity along the y direction of the i-th object
       */
      double vy (const int i) const { return m_sample[i]->vy(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_vz
       * @param i the object index
       * @return the velocity along the z direction of the i-th object
       */
      double vz (const int i) const { return m_sample[i]->vz(); }; 

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_dc
       * @param i the object index
       * @return the comoving distance of the i-th object
       */
      double dc (const int i) const { return m_sample[i]->dc(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_ra
       * @param i the object index
       * @return the Right Ascension of the i-th object
       */
      double ra (const int i) const { return m_sample[i]->ra(); };
    
      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_dec
       * @param i the object index
       * @return the Declination of the i-th object
       */
      double dec (const int i) const { return m_sample[i]->dec(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_redshift
       * @param i the object index
       * @return the redshift of the i-th object
       */
      double redshift (const int i) const { return m_sample[i]->redshift(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_weight
       * @param i the object index
       * @return the weight of the i-th object
       */
      double weight (const int i) const { return m_sample[i]->weight(); };

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_region
       * @param i the object index
       * @return the index of the region of the i-th object
       */
      long region (const int i) const { return m_sample[i]->region(); };

      /**
       * @brief get the total number of region the Catalogues is divided
       * @return the total number of regions
       */
      int Nregion () const;

      /**
       * @brief get the list of regions in which the Catalogue is
       divided     
       * @return the list of regions of regions
       */
      vector<long> get_region_list () const;

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_mass
       * @param i the object index
       * @return the mass of the i-th object
       */
      double mass (const int i) const { return m_sample[i]->mass(); }

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_magnitude
       * @param i the object index
       * @return the magnitude of the i-th object
       */
      double magnitude (const int i) const { return m_sample[i]->magnitude(); }

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_richness
       * @param i the object index
       * @return the richness of the i-th object
       */
      double richness (const int i) const { return m_sample[i]->richness(); }

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_generic
       * @param i the object index
       * @return generic properties of the i-th object
       */
      double generic (const int i) const { return m_sample[i]->generic(); }

      /**
       * @brief get the protected member Catalogue::m_sample[i]->m_radius
       * @param i the object index
       * @return radius of the i-th object
       */
      double radius (const int i) const { return m_sample[i]->radius(); }
  
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
      shared_ptr<Object> catalogue_object (const int i) const { return m_sample[i]; }

      /**
       * @brief get the object vector
       * @return vector of pointers to objects of the catalogue
       */
      vector<shared_ptr<Object>> catalogue_objects () const { return m_sample; }   

      /**
       * @brief get the X, Y, Z coordinates of the i-th object of the
       * catalogue
       *
       * @param i the object index
       * @return vector containing the three coordinates
       */
      vector<double> coordinates (const int i) const { return m_sample[i]->coords(); }
    
      /**
       * @brief get the number of objects of the catalogue
       * @return the number of objects
       */
      int nObjects () const { return m_sample.size(); }
  
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
      void stats_var (const Var, vector<double> &) const;

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
      void stats_var (const vector<Var>, vector<vector<double> > &) const;
  
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
       * @param var_name name of the variable
       * @param _var vector of variables
       * @return none
       */
      void set_var (const Var, const vector<double>); 

      /**
       * @brief change the number of objects of the catalogue
       * @param newN the new number of objects
       * @return none
       */
      void resize (const int newN) { m_sample.resize(newN); }

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
      void add_object (shared_ptr<Object> object) { m_sample.push_back(move(object)); }

      /**
       * @brief add one single object to the catalogue
       * @param object object of type \e T
       * @return none
       */
      template<typename T>
	void add_object (T object) { m_sample.push_back(move(make_shared<T>(T(object)))); }

   
      /**
       * @brief add some objects to the catalogue
       * @param sample vector of pointers to objects of type \e Object
       * @return none
       */
      void add_objects (vector<shared_ptr<Object> > sample) { 
	for (auto &&i : sample)
	  m_sample.push_back(move(i));
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
       * @brief remove all objects 
       * @return none
       */
      void remove_objects () {
	m_sample.erase(m_sample.begin(), m_sample.end());
      }
    
      /**
       * @brief replace existing objects with new ones 
       * @param sample vector of objects of type \e T
       * @return none
       */
      template<typename T>
	void remove_objects (vector<T > sample) {
	m_sample.erase(m_sample.begin(), m_sample.end());
	add_objects(sample);
      }

      /**
       * @brief replace existing objects with new ones 
       * @param sample vector of pointers to objects of type \e Object
       * @return none
       */
      void remove_objects (vector<shared_ptr<Object> > sample) {
	m_sample.erase(m_sample.begin(), m_sample.end());
	for (auto &&i : sample)
	  m_sample.push_back(move(i));
      }

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
      void computeComovingCoordinates (const Cosmology &cosm, const CoordUnits inputUnits=_radians_); 

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
      void computePolarCoordinates (const Cosmology &, const double z1=0., const double z2=10., const CoordUnits outputUnits=_radians_); 

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
      void Order (const vector<int>); 

      /**
       * @brief restore the original vector (i.e. the opposite of
       * Order(vector<int>))
       * @return none
       */
      void Order ();  
    
      /**
       * @brief write the comoving coordinates of the catalogue to an
       * output file
       * @param file_output the name of the output file
       * @return none
       */
      void write_comoving_coordinates (const string) const;

      /**
       * @brief write the polar coordinates of the catalogue to an
       * output file
       * @param file_output the name of the output file
       * @return none
       */
      void write_obs_coordinates (const string) const;

      /**
       * @brief write both the comoving and polar coordinates of the
       * catalogue to an output file
       * @param file_output the name of the output file
       * @return none
       */
      void write_coordinates (const string) const;

      /**
       * @brief get the distrance between the i-th object of the
       * catalogue and another object
       * @param i the object index
       * @param obj pointer to an object
       * @return distance between the i-th object of the catalogue and
       * the object obj
       */
      double distance (const int, shared_ptr<Object>) const;
    
      /**
       * @brief get the angular distrance between the i-th object of the
       * catalogue and another object
       * @param i the object index
       * @param obj pointer to an object
       * @return distance between the i-th object of the catalogue and
       * the object obj
       */
      double angsep_xyz (const int, shared_ptr<Object>) const;
    
      /**
       * @brief overloading of the += operator, to sum two catalogues
       * @param cc object of class Catalogue 
       * @return object of class catalogue
       */
      Catalogue operator += (shared_ptr<Catalogue> cc)
      {
	for (auto &&ss : cc->m_sample)
	  m_sample.push_back(shared_ptr<Object>(new Object(*ss)));
	return *this;
      }

      /**
       * @brief overloading of the += operator, to sum two catalogues
       * @param cc object of class Catalogue 
       * @return object of class catalogue
       */
      Catalogue operator += (const Catalogue cc)
      {
	for (auto &&ss : cc.m_sample)
	  m_sample.push_back(shared_ptr<Object>(new Object(*ss)));
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
      shared_ptr<Catalogue> cut (const Var, const double, const double, const bool excl=0);

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
      ScalarField3D density_field (const double cell_size, const int interpolation_type=0, const double kernel_radius=0., const bool useMass=0) const;

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
      ScalarField3D density_field (const double cell_size, const Catalogue mask_catalogue, const int interpolation_type = 0, const double kernel_radius=0., const bool useMass = 0) const;
      ///@}
    
    };
    
  }
}

#endif

