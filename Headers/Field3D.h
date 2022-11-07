/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *******************************************************************/

/**
 *  @file Headers/Field3D.h
 *
 *  @brief The class field3D
 *
 *  This file defines the interface of the class Field3D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __Field3D__
#define __Field3D__

#include "Data.h"


namespace cbl {

  namespace data {
  
    /**
     *  @class Field3D Field3D.h "Headers/Field3D.h"
     *
     *  @brief The class Field3D
     *
     *  This class is used to handle objects of type <EM> Field3D
     *  </EM>
     */
    class Field3D {
    
    protected:

      /// number of cells along the x-axis
      int m_nX;

      /// number of cells along the y-axis
      int m_nY;

      /// number of cells along the z-axis
      int m_nZ;

      /// number of cells along the z-axis, Fourier space
      int m_nZF;
    
      /// number of cells
      int m_nCells;

      /// number of cells, Fourier space
      int m_nCells_Fourier;

      /// X cell size 
      double m_deltaX;

      /// Y cell size
      double m_deltaY;

      /// Z cell size
      double m_deltaZ;

      /// lower x bound
      double m_MinX;

      /// lower y bound
      double m_MinY;

      /// lower z bound
      double m_MinZ;

      /// upper x bound
      double m_MaxX;

      /// upper y bound
      double m_MaxY;

      /// upper z bound
      double m_MaxZ;

      /// box volume
      double m_Volume;

      /// coordinates of the cells along the x-axis
      std::vector<double>  m_X;

      /// coordinates of the cells along the y-axis
      std::vector<double>  m_Y;

      /// coordinates of the cells along the z-axis
      std::vector<double>  m_Z;

      /// coordinates of the cells along the kx-axis
      std::vector<double>  m_kX;

      /// coordinates of the cells along the ky-axis
      std::vector<double>  m_kY;

      /// coordinates of the cells along the kz-axis
      std::vector<double>  m_kZ;

      /**
       * @brief contract 3 indeces into one
       *
       * @param i index of the i-th x-axis cell
       * @param j index of the j-th y-axis cell
       * @param k index of the k-th z-axis cell
       *
       * @return k+nZ*(j+nY*i)
       */
      long int inds_to_index (const int i, const int j, const int k) const { return k+m_nZ*(j+m_nY*i); } 

      /**
       * @brief contract 3 indeces into one, Fourier space
       *
       * @param i index of the i-th x-axis cell
       * @param j index of the j-th y-axis cell
       * @param k index of the k-th z-axis cell
       *
       * @return k+nZF*(j+nY*i)
       */
      long int inds_to_index_Fourier (const int i, const int j, const int k) const { return k+m_nZF*(j+m_nY*i); }

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      Field3D () = default;

      /**
       *  @brief constructor
       *
       *  @param deltaR size of the cubic cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */
      Field3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief constructor
       *
       *  @param nx number of x-axis cells
       *  @param ny number of y-axis cells
       *  @param nz number of z-axis cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */
      Field3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief default destructor
       *  
       */
      virtual ~Field3D () = default;

      ///@}


      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{
      
      /**
       *  @brief set the parameters
       *
       *  @param deltaR size of the cubic cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */    
      void set_parameters (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief set the parameters
       *
       *  @param nx number of x-axis cells
       *  @param ny number of y-axis cells
       *  @param nz number of z-axis cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */
      void set_parameters (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);
      
      ///@}


      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       * @brief get the private member m_nX
       * @return the number of cells along the x-axis 
       */
      int nx () const { return m_nX; }

      /**
       * @brief get the private member m_nY
       * @return the number of cells along the y-axis 
       */   
      int ny() const {return m_nY;}

      /**
       * @brief get the private member m_nZ
       * @return the number of cells along the Z-axis 
       */
      int nz () const { return m_nZ; }

      /**
       * @brief get the private member m_nZF
       * @return the number of cells along the z-axis, Fourier space
       */
      int nzFourier () const { return m_nZF; }
    
      /**
       * @brief get the private member m_nCells
       * @return the number of cells
       */
      int nCells () const { return m_nCells; }

      /**
       * @brief get the private member m_nCells_Fourier
       * @return the number of cells, Fourier space
       */
      int nCellsFourier () const { return m_nCells_Fourier; }

      /**
       * @brief get the private member m_MinX
       * @return the lower x bound 
       */
      double MinX () const { return m_MinX; }

      /**
       * @brief get the private member m_MinY
       * @return the lower y bound 
       */
      double MinY () const { return m_MinY; }

      /**
       * @brief get the private member m_MinZ
       * @return the lower z bound 
       */   
      double MinZ () const { return m_MinZ; }

      /**
       * @brief get the private member m_MaxX
       * @return the upper x bound 
       */                                    
      double MaxX () const { return m_MaxX; }
   
      /**
       * @brief get the private member m_MaxY
       * @return the upper y bound 
       */ 
      double MaxY() const { return m_MaxY; }
   
      /**
       * @brief get the private member m_MaxZ
       * @return the upper z bound 
       */
      double MaxZ () const { return m_MaxZ; }
   
      /**
       * @brief get the private member m_deltaX
       * @return the X cell size 
       */
      double deltaX () const { return m_deltaX; }
    
      /**
       * @brief get the private member m_deltaY
       * @return the Y cell size 
       */   
      double deltaY () const { return m_deltaY; }
    
      /**
       * @brief get the private member m_deltaZ
       * @return the Z cell size
       */
      double deltaZ () const { return m_deltaZ; }

      /**
       * @brief get the private member m_Volume
       * @return the box volume 
       */
      double Volume () const { return m_Volume; }

      /**
       * @brief get the value of the X coordinates at the i-th cell
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the x-axis 
       */
      double XX(const int i) const { return m_X[i];}

      /**
       * @brief get the value of the Y coordinates at the i-th cell
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the y-axis 
       */
      double YY(const int i) const { return m_Y[i];}

      /**
       * @brief get the value of the Z coordinates at the i-th cell
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the z-axis 
       */  
      double ZZ(const int i) const { return m_Z[i];}

      /**
       * @brief get the value of the X coordinates at the i-th cell,
       * Fourier space
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the x-axis 
       */   
      double kX(const int i) const { return m_kX[i];}

      /**
       * @brief get the value of the Y coordinates at the i-th cell,
       * Fourier space
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the y-axis 
       */   
      double kY(const int i) const { return m_kY[i];}
    
      /**
       * @brief get the value of the Z coordinates at the i-th cell,
       * Fourier space
       *
       * @param i the index of the cell
       *
       * @return the value of the center of the cell along the z-axis 
       */    
      double kZ(const int i) const { return m_kZ[i];}
      
      /**
       * @brief get the value of the scalar field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field
       */
      virtual double ScalarField (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "Scalarield", "Field3D.h"); return 0; }

        
      /**
       * @brief get the value of the scalar field
       * 
       * @param pos vector containing the point coordinates
       *
       * @return the value of the vector field
       */
      virtual double ScalarField (const std::vector<double> pos) const
      { (void)pos; ErrorCBL("", "Scalarield", "Field3D.h"); return 0; }
    
      /**
       * @brief get the value of the scalar field
       *
       * @return the values of the scalar field
       */
      virtual std::vector<double> ScalarField () const
      { ErrorCBL("", "Scalarield", "Field3D.h"); std::vector<double> vv; return vv; }

      /**
       * @brief get the value of the vector field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the values of the vector field
       */
      virtual std::vector<double> VectorField (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "VectorField", "Field3D.h"); return {0}; }

      /**
       * @brief get the value of the vector field
       * 
       * @param pos vector containing the point coordinates
       *
       * @return vector containing the values of the vector field
       */
      virtual std::vector<double> VectorField (const std::vector<double> pos) const
      { (void)pos; ErrorCBL("", "VectorField", "Field3D.h"); double vv = 0.; return {vv}; }

      /**
       * @brief get the value of the scalar field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the values of the scalar field, Fourier space, real part
       */
      virtual double ScalarField_FourierSpace_real (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "ScalarField_FourierSpace_real", "Field3D.h"); return 0; }

      /**
       * @brief get the value of the scalar field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the scalar field, Fourier space, complex part
       */
      virtual double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "ScalarField_FourierSpace_complex", "Field3D.h"); return 0; }

      /**
       * @brief get the value of the scalar field, Fourier space, real part
       *
       * @return the values of the scalar field
       */
      virtual std::vector<double> ScalarField_FourierSpace_real () const
      { ErrorCBL("", "Scalarield_FourierSpace_real", "Field3D.h"); std::vector<double> vv; return vv; }

      /**
       * @brief get the value of the scalar field, Fourier space, complex part
       *
       * @return the value of the vector field
       */
      virtual std::vector<double> ScalarField_FourierSpace_complex () const
      { ErrorCBL("", "Scalarield_FourierSpace_complex", "Field3D.h"); std::vector<double> vv; return vv; }

      /**
       * @brief get the value of the vector field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, real part
       */
      virtual std::vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "VectorField_FourierSpace_real", "Field3D.h"); return {0}; }

      /**
       * @brief get the value of the vector field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, complex part
       */
      virtual std::vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorCBL("", "VectorField_FourierSpace_complex", "Field3D.h"); return {0}; }

      ///@}


      /**
       *  @name Member functions to compute data properties
       */
      ///@{
      
      /**
       * @brief perform the Fourier transform on the field
       */
      virtual void FourierTransformField ()
      { ErrorCBL("", "FourierTransformField", "Field3D.h"); }

      /**
       * @brief perform the anti-Fourier transform on the field
       */
      virtual void FourierAntiTransformField ()
      { ErrorCBL("", "FourierAntiTransformField", "Field3D.h"); }

      /**
       * @brief perform a smoothing of the field with a gaussian kernel
       * @param kernel_size size of the gaussian kernel
       */
      virtual void GaussianConvolutionField (const double kernel_size)
      { (void)kernel_size; ErrorCBL("", "GaussianConvolutionField", "Field3D.h"); }

      /**
       * @brief set the value of the scalar field
       * 
       * @param value value of the scalar field
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add 1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       */
      virtual void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_ScalarField", "Field3D.h"); }
    
      /**
       * @brief set the value of the vectorr field
       * 
       * @param value vector containing values of the vector field
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add 1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       */
      virtual void set_VectorField (const std::vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_Vectorield", "Field3D.h"); }

      /**
       * @brief set the value of the scalar field in Fourier space, real part
       * 
       * @param value value of the scalar field in Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @param add 1 &rarr; add to the current value; 0 &rarr;
       * overwrite the value
       */
      virtual void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_ScalarField_FourierSpace_real", "Field3D.h"); }
    
      /**
       * @brief set the value of the scalar field in Fourier space, complex part
       * 
       * @param value value of the scalar field in Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @param add 1 &rarr; add to the current value; 0 &rarr;
       * overwrite the value
       */
      virtual void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_ScalarField_FourierSpace_complex", "Field3D.h"); }

      /**
       *  @brief set the value of the vector field, Fourier space,
       *  real part
       * 
       *  @param value vector containing values of the vector field,
       *  Fourier space, real part
       *  @param i the i-th cell along the x-axis
       *  @param j the j-th cell along the y-axis
       *  @param k the k-th cell along the z-axis
       *  @param add 1 &rarr; add to the current value; 0 &rarr;
       *  overwrite the value
       */
      virtual void set_VectorField_FourierSpace_real (const std::vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_VectorField_FourierSpace_real", "Field3D.h"); }

      /**
       * @brief set the value of the vector field, Fourier space,
       * complex part
       * 
       * @param value vector containing values of the vector field,
       * Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add 1 &rarr; add to the current value; 0 &rarr;
       * overwrite the value
       */
      virtual void set_VectorField_FourierSpace_complex (const std::vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorCBL("", "set_VectorField_FourierSpace_complex", "Field3D.h"); }

      /**
       * @brief set to 0 the fields
       */
      virtual void reset()
      { ErrorCBL("", "reset", "Field3D.h"); }

      ///@}
      
    };

    
    // ==================================================================================================================
    // ==================================================================================================================
    // ==================================================================================================================

    
    /**
     *  @class ScalarField3D Field3D.h "Headers/Field3D.h"
     *
     *  @brief The class ScalarField3D
     *
     *  This class is used to handle objects of type <EM> ScalarField3D
     *  </EM>
     */
    class ScalarField3D : public Field3D {
      
    protected:

      /// scalar field
      double *m_field;

      /// fourier transform of the scalar field
      fftw_complex *m_field_FourierSpace;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      ScalarField3D () = default;

      /**
       *  @brief constructor
       *
       *  @param deltaR size of the cubic cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */
      ScalarField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief constructor
       *
       *  @param nx number of x-axis cells
       *  @param ny number of y-axis cells
       *  @param nz number of z-axis cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       *  
       */
      ScalarField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);
      /**
       *  @brief default destructor
       *  
       */
      virtual ~ScalarField3D () = default;

      ///@}


      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{
      
      /**
       * @brief set the value of the scalar field
       * 
       * @param value value of the scalar field
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * 
       */
      void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0);

      /**
       *  @brief set the value of the scalar field in Fourier space,
       *  real part
       * 
       *  @param value value of the scalar field in Fourier space, real
       *  part
       *  @param i the i-th cell along the x-axis
       *  @param j the j-th cell along the y-axis
       *  @param k the k-th cell along the z-axis
       *  @param add 1 &rarr; add to the current value; 0 &rarr;
       *  overwrite the value
       *
       *  
       */
      void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the scalar field in Fourier space,
       * real part
       * 
       * @param value value of the scalar field in Fourier space, real
       * part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * 
       */
      void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0);

      ///@}


      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       * @brief get the value of the scalar field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field
       */
      double ScalarField (const int i, const int j, const int k) const;
            
      /**
       * @brief get the value of the scalar field
       * 
       * @param pos vector containing the point coordinates
       *
       * @return the value of the vector field
       */
      double ScalarField (const std::vector<double> pos) const;

      /**
       * @brief get the values of the scalar field
       *
       * @return the values of the scalar field
       */
      std::vector<double> ScalarField () const;

      /**
       * @brief get the value of the scalar field, Fourier space, real
       * part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, real part
       */
      double ScalarField_FourierSpace_real (const int i, const int j, const int k) const;

      /**
       * @brief get the value of the scalar field, Fourier space,
       * complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, complex part
       */
      double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const;
       
      ///@}
      

      /**
       *  @name Member functions to compute operations on the field
       */
      ///@{
      /**
       * @brief perform the Fourier transform on the field
       * 
       */
      void FourierTransformField ();

      /**
       * @brief perform the anti-Fourier transform on the field
       * 
       */
      void FourierAntiTransformField ();

      /**
       * @brief perform a smoothing of the field with a gaussian kernel
       * @param kernel_size size of the gaussian kernel
       * 
       */
      void GaussianConvolutionField (const double kernel_size);

      /**
       * @brief set to 0 the fields
       *
       * 
       */
      void reset();

      ///@}
    };


    // ==================================================================================================================
    // ==================================================================================================================
    // ==================================================================================================================

    
    /**
     *  @class VectorField3D Field3D.h "Headers/Field3D.h"
     *
     *  @brief The class VectorField3D
     *
     *  This class is used to handle objects of type <EM> VectorField3D
     *  </EM>
     */
    class VectorField3D : public Field3D{
    protected:

      /// vector field
      std::vector<double *> m_field;

      /// vector field in fourier space
      std::vector<fftw_complex *> m_field_FourierSpace;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      VectorField3D () = default;

      /**
       *  @brief constructor
       *
       *  @param deltaR size of the cubic cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       */
      VectorField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief constructor
       *
       *  @param nx number of x-axis cells
       *  @param ny number of y-axis cells
       *  @param nz number of z-axis cells
       *  @param minX lower x bound
       *  @param maxX upper x bound
       *  @param minY lower y bound
       *  @param maxY upper y bound
       *  @param minZ lower z bound
       *  @param maxZ upper z bound
       *
       */
      VectorField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);
      /**
       *  @brief default destructor
       *  
       */
      ~VectorField3D () = default;

      ///@}


      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{
      
      /**
       * @brief set the value of the vectorr field
       * 
       * @param value vector containing values of the vector field
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * 
       */
      void set_VectorField (const std::vector<double> value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the vector field, Fourier space, real
       * part
       * 
       * @param value vector containing values of the vector field,
       * Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       */
      void set_VectorField_FourierSpace_real (const std::vector<double> value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the vector field, Fourier space,
       * complex part
       * 
       * @param value vector containing values of the vector field, 
       * Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       */
      void set_VectorField_FourierSpace_complex (const std::vector<double> value, const int i, const int j, const int k, const bool add=0);

      ///@}


      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       * @brief get the value of the vector field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field
       */
      std::vector<double> VectorField (const int i, const int j, const int k) const;
      
      /**
       * @brief get the value of the vector field
       * 
       * @param pos vector containing the point coordinates
       *
       * @return vector containing the value of the vector field
       */
      std::vector<double> VectorField (const std::vector<double> pos) const;

      /**
       * @brief get the value of the vector field, Fourier space, real
       * part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, real part
       */
      std::vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const;

      /**
       * @brief get the value of the vector field, Fourier space,
       * complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, complex part
       */
      std::vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const;

      ///@}

      
      /**
       *  @name Member functions to compute operations on the field
       */
      ///@{
      
      /**
       * @brief perform the anti-Fourier transform on the field
       * 
       */
      void FourierTransformField ();

      /**
       * @brief perform the anti-Fourier transform on the field
       * 
       */
      void FourierAntiTransformField ();

      /**
       * @brief set to 0 the fields
       *
       * 
       */
      void reset();

      ///@}

    };

  }
}

#endif
