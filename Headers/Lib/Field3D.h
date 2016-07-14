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
 *  @file Headers/Lib/Field3D.h
 *
 *  @brief Implementation of the field3D data structure
 *
 *  This file defines the interface of the class Field3D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __Field3D__
#define __Field3D__

#include "Data.h"
#include <fftw3.h>


namespace cosmobl {

  namespace data {
  
    /**
     *  @class Field3D Field3D.h "Headers/Lib/Field3D.h"
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
       *  @return object of type Field3D
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
       *  @return object of type Field3D
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
       *  @return object of type Field3D
       */
      Field3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Field3D() = default;

      ///@}

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
       *  @return none
       */    
      void set_parameters (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

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
       *  @return none
       */
      void set_parameters (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);

      /**
       * @brief return the private member m_nX
       * @return the number of cells along the x-axis 
       */
      int nx () const { return m_nX; }

      /**
       * @brief return the private member m_nY
       * @return the number of cells along the y-axis 
       */   
      int ny() const {return m_nY;}

      /**
       * @brief return the private member m_nZ
       * @return the number of cells along the Z-axis 
       */
      int nz () const { return m_nZ; }

      /**
       * @brief return the private member m_nZF
       * @return the number of cells along the z-axis, Fourier space
       */
      int nzFourier () const { return m_nZF; }
    
      /**
       * @brief return the private member m_nCells
       * @return the number of cells
       */
      int nCells () const { return m_nCells; }

      /**
       * @brief return the private member m_nCells_Fourier
       * @return the number of cells, Fourier space
       */
      int nCellsFourier () const { return m_nCells_Fourier; }

      /**
       * @brief return the private member m_MinX
       * @return the lower x bound 
       */
      double MinX () const { return m_MinX; }

      /**
       * @brief return the private member m_MinY
       * @return the lower y bound 
       */
      double MinY () const { return m_MinY; }

      /**
       * @brief return the private member m_MinZ
       * @return the lower z bound 
       */   
      double MinZ () const { return m_MinZ; }

      /**
       * @brief return the private member m_MaxX
       * @return the upper x bound 
       */                                    
      double MaxX () const { return m_MaxX; }
   
      /**
       * @brief return the private member m_MaxY
       * @return the upper y bound 
       */ 
      double MaxY() const {return m_MaxY;}
   
      /**
       * @brief return the private member m_MaxZ
       * @return the upper z bound 
       */
      double MaxZ () const { return m_MaxZ; }
   
      /**
       * @brief return the private member m_deltaX
       * @return the X cell size 
       */
      double deltaX () const { return m_deltaX; }
    
      /**
       * @brief return the private member m_deltaY
       * @return the Y cell size 
       */   
      double deltaY () const { return m_deltaY; }
    
      /**
       * @brief return the private member m_deltaZ
       * @return the Z cell size
       */
      double deltaZ () const { return m_deltaZ; }

      /**
       * @brief return the private member m_Volume
       * @return the box volume 
       */
      double Volume () const { return m_Volume; }

      /**
       * @brief perform the Fourier transform on the field
       * @return none
       */
      virtual void FourierTransformField ()
      { ErrorMsg("Error in FourierTransformField of Field3D"); }

      /**
       * @brief perform the anti-Fourier transform on the field
       * @return none
       */
      virtual void FourierAntiTransformField ()
      { ErrorMsg("Error in FourierAntiTransformField of Field3D"); }

      /**
       * @brief perform a smoothing of the field with a gaussian kernel
       * @param kernel_size size of the gaussian kernel
       * @return none
       */
      virtual void GaussianConvolutionField (const double kernel_size)
      { (void)kernel_size; ErrorMsg("Error in GaussianConvolutionField of Field3D"); }

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
       * @return none
       */
      virtual void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_ScalarField of Field3D"); }
    
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
       * @return none
       */
      virtual void set_VectorField (const vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_Vectorield of Field3D"); }

      /**
       * @brief set the value of the scalar field in Fourier space, real part
       * 
       * @param value value of the scalar field in Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      virtual void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_ScalarField_FourierSpace_real of Field3D"); }
    
      /**
       * @brief set the value of the scalar field in Fourier space, complex part
       * 
       * @param value value of the scalar field in Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      virtual void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_ScalarField_FourierSpace_complex of Field3D"); }

      /**
       * @brief set the value of the vector field, Fourier space, real part
       * 
       * @param value vector containing values of the vector field,
       * Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      virtual void set_VectorField_FourierSpace_real (const vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_VectorField_FourierSpace_real of Field3D"); }

      /**
       * @brief set the value of the vector field, Fourier space, complex part
       * 
       * @param value vector containing values of the vector field, 
       * Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      virtual void set_VectorField_FourierSpace_complex (const vector<double> value, const int i, const int j, const int k, const bool add=0)
      { (void)value; (void)i; (void)j; (void)k; (void)add; ErrorMsg("Error in set_VectorField_FourierSpace_complex of Field3D"); }

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
      { (void)i; (void)j; (void)k; ErrorMsg("Error in Scalarield of Field3D"); double vv; return vv; }

      /**
       * @brief get the value of the vector field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field
       */
      virtual vector<double> VectorField (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorMsg("Error in VectorField of Field3D"); double vv; return {vv}; }

      /**
       * @brief get the value of the scalar field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, real part
       */
      virtual double ScalarField_FourierSpace_real (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorMsg("Error in ScalarField_FourierSpace of Field3D"); double vv; return vv; }

      /**
       * @brief get the value of the scalar field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, complex part
       */
      virtual double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorMsg("Error in ScalarField_FourierSpace of Field3D"); double vv; return vv; }

      /**
       * @brief get the value of the vector field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, real part
       */
      virtual vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorMsg("Error in VectorField_FourierSpace of Field3D"); double vv; return {vv}; }

      /**
       * @brief get the value of the vector field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, complex part
       */
      virtual vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const
      { (void)i; (void)j; (void)k; ErrorMsg("Error in VectorField_FourierSpace_complex of Field3D"); double vv; return {vv}; }
    };

    class ScalarField3D : public Field3D{
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
       *  @return object of type ScalarField3D
       */
      ScalarField3D() {}

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
       *  @return object of type ScalarField3D
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
       *  @return object of type ScalarField3D
       */
      ScalarField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~ScalarField3D() {}

      ///@}

      /**
       * @brief perform the Fourier transform on the field
       * @return none
       */
      void FourierTransformField ();

      /**
       * @brief perform the anti-Fourier transform on the field
       * @return none
       */
      void FourierAntiTransformField ();

      /**
       * @brief perform a smoothing of the field with a gaussian kernel
       * @param kernel_size size of the gaussian kernel
       * @return none
       */
      void GaussianConvolutionField (const double kernel_size);

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
       * @return none
       */
      void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the scalar field in Fourier space, real part
       * 
       * @param value value of the scalar field in Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the scalar field in Fourier space, real part
       * 
       * @param value value of the scalar field in Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0);

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
       * @brief get the value of the scalar field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, real part
       */
      double ScalarField_FourierSpace_real (const int i, const int j, const int k) const;

      /**
       * @brief get the value of the scalar field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return the value of the vector field, Fourier space, complex part
       */
      double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const;

    };

    class VectorField3D : public Field3D{
    protected:

      /// vector field
      vector<double *> m_field;

      /// vector field in fourier space
      vector<fftw_complex *> m_field_FourierSpace;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of type VectorField3D
       */
      VectorField3D () {}

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
       *  @return object of type VectorField3D
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
       *  @return object of type VectorField3D
       */
      VectorField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ);
      /**
       *  @brief default destructor
       *  @return none
       */
      ~VectorField3D () {}

      ///@}

      /**
       * @brief perform the anti-Fourier transform on the field
       * @return none
       */
      void FourierTransformField ();

      /**
       * @brief perform the anti-Fourier transform on the field
       * @return none
       */
      void FourierAntiTransformField ();

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
       * @return none
       */
      void set_VectorField (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the vector field, Fourier space, real part
       * 
       * @param value vector containing values of the vector field,
       * Fourier space, real part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      void set_VectorField_FourierSpace_real (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief set the value of the vector field, Fourier space, complex part
       * 
       * @param value vector containing values of the vector field, 
       * Fourier space, complex part
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       * @param add  1 &rarr; add to the current value; 0
       *  &rarr; overwrite the value
       *
       * @return none
       */
      void set_VectorField_FourierSpace_complex (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      /**
       * @brief get the value of the vector field
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field
       */
      vector<double> VectorField (const int i, const int j, const int k) const;

      /**
       * @brief get the value of the vector field, Fourier space, real part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, real part
       */
      vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const;

      /**
       * @brief get the value of the vector field, Fourier space, complex part
       * 
       * @param i the i-th cell along the x-axis
       * @param j the j-th cell along the y-axis
       * @param k the k-th cell along the z-axis
       *
       * @return vector containing the value of the vector field, Fourier space, complex part
       */
      vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const;

    };

  }

}

#endif
