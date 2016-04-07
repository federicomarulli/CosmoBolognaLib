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

#include "Cosmology.h"
#include <fftw3.h>


namespace cosmobl {

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

    int m_nX;
    int m_nY;
    int m_nZ;
    int m_nZF;
    
    int m_nCells;
    int m_nCells_Fourier;

    double m_deltaX;
    double m_deltaY;
    double m_deltaZ;

    double m_MinX;
    double m_MinY;
    double m_MinZ;

    double m_MaxX;
    double m_MaxY;
    double m_MaxZ;

    double m_Volume;

    long int inds_to_index(int i, int j, int k) const {return k+m_nZ*(j+m_nY*i);} 
    long int inds_to_index_Fourier(int i, int j, int k) const {return k+m_nZF*(j+m_nY*i);}

  public:
    
    Field3D () {}

    virtual ~Field3D () {}

    Field3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

    Field3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

    void set_parameters(const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

    void set_parameters(const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

    int nx() const {return m_nX;}
    
    int ny() const {return m_nY;}

    int nz() const {return m_nZ;}

    int nzFourier() const {return m_nZF;}

    int nCells() const {return m_nCells;}

    int nCellsFourier() const {return m_nCells_Fourier;}

    double MinX() const {return m_MinX;}

    double MinY() const {return m_MinY;}
    
    double MinZ() const {return m_MinZ;}
                                      
    double MaxX() const {return m_MaxX;}
    
    double MaxY() const {return m_MaxY;}

    double MaxZ() const {return m_MaxZ;}

    double deltaX() const {return m_deltaX;}
    
    double deltaY() const {return m_deltaY;}

    double deltaZ() const {return m_deltaZ;}

    double Volume() const {return m_Volume;}

    virtual void FourierTransformField ()
    { ErrorMsg("Error in FourierTransformField of Field3D"); }

    virtual void FourierAntiTransformField ()
    { ErrorMsg("Error in FourierAntiTransformField of Field3D"); }

    virtual void GaussianConvolutionField (const double kernel_size)
    { ErrorMsg("Error in GaussianConvolutionField of Field3D"); }

    virtual void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_ScalarField of Field3D"); }

    virtual void set_VectorField (const vector<double> value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_Vectorield of Field3D"); }

    virtual void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_ScalarField_FourierSpace_real of Field3D"); }

    virtual void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_ScalarField_FourierSpace_complex of Field3D"); }

    virtual void set_VectorField_FourierSpace_real (const vector<double> value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_VectorField_FourierSpace_real of Field3D"); }

    virtual void set_VectorField_FourierSpace_complex (const vector<double> value, const int i, const int j, const int k, const bool add=0)
    { ErrorMsg("Error in set_VectorField_FourierSpace_complex of Field3D"); }

    virtual double ScalarField (const int i, const int j, const int k) const
    { ErrorMsg("Error in Scalarield of Field3D"); double vv; return vv; }

    virtual vector<double> VectorField (const int i, const int j, const int k) const
    { ErrorMsg("Error in VectorField of Field3D"); double vv; return {vv}; }

    virtual double ScalarField_FourierSpace_real (const int i, const int j, const int k) const
    { ErrorMsg("Error in ScalarField_FourierSpace of Field3D"); double vv; return vv; }

    virtual double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const
    { ErrorMsg("Error in ScalarField_FourierSpace of Field3D"); double vv; return vv; }

    virtual vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const
    { ErrorMsg("Error in VectorField_FourierSpace of Field3D"); double vv; return {vv}; }

    virtual vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const
    { ErrorMsg("Error in VectorField_FourierSpace_complex of Field3D"); double vv; return {vv}; }
  };

  class ScalarField3D : public Field3D{
    protected:
      double *m_field;
      fftw_complex *m_field_FourierSpace;

    public:
      ScalarField3D() {}

      ~ScalarField3D() {}

      ScalarField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

      ScalarField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

      void FourierTransformField ();

      void FourierAntiTransformField ();

      void GaussianConvolutionField (const double kernel_size);

      void set_ScalarField (const double value, const int i, const int j, const int k, const bool add=0);

      void set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add=0);

      void set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add=0);

      double ScalarField (const int i, const int j, const int k) const;

      double ScalarField_FourierSpace_real (const int i, const int j, const int k) const;

      double ScalarField_FourierSpace_complex (const int i, const int j, const int k) const;

  };

  class VectorField3D : public Field3D{
    protected:
      vector<double *> m_field;
      vector<fftw_complex *> m_field_FourierSpace;

    public:
      VectorField3D () {}

      ~VectorField3D () {}

      VectorField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

      VectorField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double MaxZ);

      void FourierTransformField ();

      void FourierAntiTransformField ();

      void set_VectorField (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      void set_VectorField_FourierSpace_real (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      void set_VectorField_FourierSpace_complex (const vector<double> value, const int i, const int j, const int k, const bool add=0);

      vector<double> VectorField (const int i, const int j, const int k) const;

      vector<double> VectorField_FourierSpace_real (const int i, const int j, const int k) const;

      vector<double> VectorField_FourierSpace_complex (const int i, const int j, const int k) const;

  };
}

#endif
