/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Histogram.h
 *
 *  @brief Class used to handle binned variables
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __HIST__
#define __HIST__ 


#include "Kernel.h"


// =====================================================================================


namespace cbl {

  namespace glob {

    /**
     * @enum HistogramType
     * @brief the histogram type
     */
    enum class HistogramType {

      /// the binned counts, \f$N_V[i]\f$
      _N_V_,

      /// the normalised binned counts, i.e. \f$n_V[i]=N_V[i]/fact\f$, where the factor \f$fact\f$ is a number provided in input
      _n_V_,

      /// \f$n_V[i]/(edge[i+1]-edge[i])\f$, where \f$edge\f$ are the bin limits
      _dn_dV_,

      /// \f$n_V[i]/(\log_{10}(edge[i+1])-\log_{10}(edge[i]))\f$, where \f$edge\f$ are the bin limits
      _dn_dlogV_

    };

    /**
     * @brief return a vector containing the
     * HistogramType names
     * @return a vector containing the
     * HistogramType names
     */
    inline std::vector<std::string> HistogramTypeNames ()
    { return {"dn_dV", "dn_dlogV", "N_V", "n_V"}; }

    /**
     * @brief cast an enum of type HistogramType
     * from its index
     * @param histogramTypeIndex the histogramType index
     * @return object of class HistogramType
     */
    inline HistogramType HistogramTypeCast (const int histogramTypeIndex)
    { return castFromValue<HistogramType>(histogramTypeIndex); }

    /**
     * @brief cast an enum of type HistogramType
     * from its name
     * @param histogramTypeName the histogramType name
     * @return object of class HistogramType
     */
    inline HistogramType HistogramTypeCast (const std::string histogramTypeName)
    { return castFromName<HistogramType>(histogramTypeName, HistogramTypeNames()); }

    /**
     * @brief cast an enum of type HistogramType
     * from indeces
     * @param histogramTypeIndeces the histogramType indeces
     * @return object of class HistogramType
     */
    inline std::vector<HistogramType> HistogramTypeCast (const std::vector<int> histogramTypeIndeces)
    { return castFromValues<HistogramType>(histogramTypeIndeces); } 

    /**
     * @brief cast an enum of type HistogramType
     * from thier names
     * @param histogramTypeNames the histogramType names
     * @return vector of HistogramType enums
     */
    inline std::vector<HistogramType> HistogramTypeCast (const std::vector<std::string> histogramTypeNames)
    { return castFromNames<HistogramType>(histogramTypeNames, HistogramTypeNames()); }

    
    /**
     *  @class Histogram Histogram.h "Headers/Histogram.h"
     *
     */
    class Histogram
    {
      public:

	/**
	 *  @brief default constructor
	 *  @return object of class Histogram
	 */
	Histogram () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Histogram () = default;

	/**
	 *  @name Functions to set the private members of the class
	 */
	///@{
	/**
	 * @brief set the histogram variables
	 *
	 * @param nbins the number of bins
	 * @param minVar the variable minimum 
	 * @param maxVar the variable maximum 
	 * @param shift the shift of the bin
	 * @param bin_type the binning type
	 *
	 * @return none
	 */
	virtual void set (const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift=0.5, const BinType bin_type=BinType::_linear_)
	{ (void)nbins; (void)minVar; (void)maxVar; (void)shift; (void)bin_type; ErrorCBL("Error in set of Histogram.h!"); }

	/**
	 * @brief get the histogram index
	 *
	 * @param var value of the var
	 *
	 * @return the histogram index
	 */
	virtual int digitize (const double var)
	{ (void)var; ErrorCBL("Error in digitize of Histogram.h!"); return 0;}

	/**
	 * @brief get the histogram indeces
	 *
	 * @param var values of the var
	 *
	 * @return the histogram indeces
	 */
	virtual std::vector<int> digitize (const std::vector<double> var)
	{ (void)var; ErrorCBL("Error in digitize of Histogram.h!"); std::vector<int> vv; return vv;}

	/**
	 * @brief bin the data
	 *
	 * @param var value of the var
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	virtual void put (const double var, const double weight)
	{ (void)var; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param var values of the var
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	virtual void put (const std::vector<double> var, const std::vector<double> weight)
	{ (void)var; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param bin value of the bin
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	virtual void put (const int bin, const double weight)
	{ (void)bin; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param bins values of the bin
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	virtual void put (const std::vector<int> bins, const std::vector<double> weight)
	{ (void)bins; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief set the histogram variables
	 *
	 *  @param nbins1 the number of bins for the first variable
	 *  @param nbins2 the number of bins for the second variable
	 *  @param minVar1 minimum range  for the first variable
	 *  @param maxVar1 maximmum range for the first variable
	 *  @param minVar2 minimum range for the second variable 
	 *  @param maxVar2 maximmum range for the second variable
	 *  @param shift1 bin shift for the first variable
	 *  @param shift2 bin shift for the second variable
	 *  @param bin_type1 the binning type for the first variable
	 *  @param bin_type2 the binning type for the second variable
	 *
	 * @return none
	 */
	virtual void set (const size_t nbins1, const size_t nbins2, const double minVar1=par::defaultDouble, const double maxVar1=par::defaultDouble, const double minVar2=par::defaultDouble, const double maxVar2=par::defaultDouble, const double shift1=0.5, const double shift2=0.5, const BinType bin_type1=BinType::_linear_, const BinType bin_type2=BinType::_linear_)
	{ (void)nbins1; (void)minVar1; (void)maxVar1; (void)nbins2; (void)minVar2; (void)maxVar2; (void)shift1; (void)shift2; (void)bin_type1; (void)bin_type2; ErrorCBL("Error in set of Histogram.h!"); }

	/**
	 * @brief get the histogram index
	 *
	 * @param var1 value of the first var
	 * @param var2 value of the second var
	 *
	 * @return the histogram index
	 */
	virtual std::vector<int> digitize (const double var1, const double var2)
	{ (void)var1; (void)var2; ErrorCBL("Error in digitize of Histogram.h!"); std::vector<int> vv; return vv;}

	/**
	 * @brief get the histogram indeces
	 *
	 * @param var1 values of the first var
	 * @param var2 values of the second var
	 *
	 * @return the histogram indeces
	 */
	virtual std::vector<std::vector<int>> digitize (const std::vector<double> var1, const std::vector<double> var2)
	{ (void)var1; (void)var2; ErrorCBL("Error in digitize of Histogram.h!"); std::vector<std::vector<int>> vv; return vv;}

	/**
	 * @brief bin the data
	 *
	 * @param var1 value of the first var
	 *
	 * @param var2 value of the second var
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	virtual void put (const double var1, const double var2, const double weight)
	{ (void)var1; (void)var2; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param var1 values of the first var
	 *
	 * @param var2 values of the second var
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	virtual void put (const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> weight)
	{ (void)var1; (void)var2; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param bin1 value of the first bin
	 *
	 * @param bin2 value of the second bin
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	virtual void put (const int bin1, const int bin2, const double weight)
	{ (void)bin1; (void)bin2; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	/**
	 * @brief bin the data
	 *
	 * @param bins values of the bins
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	virtual void put (const std::vector<std::vector<int>> bins, const std::vector<double> weight)
	{ (void)bins; (void)weight; ErrorCBL("Error in put of Histogram.h!");}

	///@}

	/**
	 *  @name Functions to get the private members of the class
	 */
	///@{

	/**
	 * @brief return the number
	 * of bins
	 *
	 * @return the number of bins
	 */
	virtual size_t nbins () const 
	{ErrorCBL("Error in nbins of Histogram.h!"); return 1; }

	/**
	 * @brief return the bin size
	 *
	 * @return the bin size
	 */
	virtual double binSize () const
	{ErrorCBL("Error in binSize of Histogram.h!"); return 0.;}

	/**
	 * @brief return the lower limit of the
	 * histogram
	 *
	 * @return the lower limit of the histogram
	 */
	virtual double minVar() const
	{ErrorCBL("Error in minVar of Histogram.h!"); return 0.;}
	
	/**
	 * @brief return the upper limit of the
	 * histogram
	 *
	 * @return the lower upper of the histogram
	 */
	virtual double maxVar() const
	{ErrorCBL("Error in maxVar of Histogram.h!"); return 0.;}

	/**
	 * @brief return the bin shift
	 *
	 * @return the bin shift
	 */
	virtual double shift() const
	{ErrorCBL("Error in shift of Histogram.h!"); return 0.;}

	/**
	 * @brief return the bin type
	 *
	 * @return the bin type
	 */
	virtual BinType bin_type() const
	{ErrorCBL("Error in bin_type of Histogram.h!"); return BinType::_linear_;}

	/**
	 * @brief return the i-th bin
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin
	 */
	virtual double bin (const size_t i) const 
	{(void)i; ErrorCBL("Error in bin of Histogram.h!"); return 0.;}

	/**
	 * @brief return the bins
	 *
	 * @return the histogram bins
	 */
	virtual std::vector<double> bins () const
	{ ErrorCBL("Error in bins of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the i-th edge
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge
	 */
	virtual double edge (const size_t i) const
	{(void)i; ErrorCBL("Error in edge of Histogram.h!"); return 0.;}

	/**
	 * @brief return the histogram edges
	 *
	 * @return the histogram edges
	 */
	virtual std::vector<double> edges () const
	{ ErrorCBL("Error in edges of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the number
	 * of bins for the first variable
	 *
	 * @return the number of bins for the first variable
	 */
	virtual size_t nbins1 () const
	{ErrorCBL("Error in bins1 of Histogram.h!"); return 1; }

	/**
	 * @brief return the first variable bin size
	 *
	 * @return the first variable bin size
	 */
	virtual double binSize1 () const
	{ ErrorCBL("Error in binSize1 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the lower limit of the
	 * histogram for the first variable
	 *
	 * @return the lower limit of the histogram
	 * for the first variable
	 */
	virtual double minVar1() const
	{ ErrorCBL("Error in minVar1 of Histogram.h!"); return 1.; }
	
	/**
	 * @brief return the upper limit of the
	 * histogram for the first variable
	 *
	 * @return the upper limit of the histogram
	 * for the first variable
	 */
	virtual double maxVar1() const
	{ ErrorCBL("Error in maxVar1 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the bin shift for the
	 * first variable
	 *
	 * @return the bin shift for the first variable
	 */
	virtual double shift1() const
	{ ErrorCBL("Error in shift1 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the bin type for the first
	 * variable 
	 *
	 * @return the bin type for the first variable
	 */
	virtual BinType bin_type1() const
	{ ErrorCBL("Error in bin_type1 of Histogram.h!"); return BinType::_linear_; }
	
	/**
	 * @brief return the number
	 * of bins for the second variable
	 *
	 * @return the number of bins for the second variable
	 */
	virtual size_t nbins2 () const
	{ErrorCBL("Error in bins2 of Histogram.h!"); return 1; }

	/**
	 * @brief return the second variable bin size
	 *
	 * @return the second variable bin size
	 */
	virtual double binSize2 () const
	{ErrorCBL("Error in binSize2 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the lower limit of the
	 * histogram for the second variable
	 *
	 * @return the lower limit of the histogram
	 * for the second variable
	 */
	virtual double minVar2() const
	{ ErrorCBL("Error in minVar2 of Histogram.h!"); return 1.; }
	
	/**
	 * @brief return the upper limit of the
	 * histogram for the second variable
	 *
	 * @return the upper limit of the histogram
	 * for the second variable
	 */
	virtual double maxVar2() const
	{ ErrorCBL("Error in maxVar2 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the bin shift for the
	 * second variable
	 *
	 * @return the bin shift for the second variable
	 */
	virtual double shift2() const
	{ ErrorCBL("Error in shift2 of Histogram.h!"); return 1.; }

	/**
	 * @brief return the bin type for the second
	 * variable 
	 *
	 * @return the bin type for the second variable
	 */
	virtual BinType bin_type2() const
	{ ErrorCBL("Error in bin_type2 of Histogram.h!"); return BinType::_linear_; }

	/**
	 * @brief return the i-th bin of the first variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin of the first variable
	 */
	virtual double bin1 (const size_t i) const
	{ (void)i; ErrorCBL("Error in histogram_bin1 of Histogram.h!"); return 0.;}

	/**
	 * @brief return the first variable bins
	 *
	 * @return the first variable bins
	 */
	virtual std::vector<double> bins1 () const
	{ErrorCBL("Error in histogram_bins1 of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the i-th edge of the first variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge of the first variable
	 */
	virtual double edge1 (const size_t i) const
	{ (void)i; ErrorCBL("Error in histogram_edge1 of Histogram.h!"); return 0.;}

	/**
	 * @brief return the histogram edges of the first variable
	 *
	 * @return the histogram edges of the first variable
	 */
	virtual std::vector<double> edges1 () const
	{ErrorCBL("Error in histogram_edges1 of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the i-th bin of the second variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin of the second variable
	 */
	virtual double bin2 (const size_t i) const
	{ (void)i; ErrorCBL("Error in histogram_bin2 of Histogram.h!"); return 0.;}

	/**
	 * @brief return the second variable bins
	 *
	 * @return the second variable bins
	 */
	virtual std::vector<double> bins2 () const
	{ErrorCBL("Error in histogram_bins2 of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the i-th edge of the second variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge of the second variabl
	 */
	virtual double edge2 (const size_t i) const
	{ (void)i; ErrorCBL("Error in histogram_edge2 of Histogram.h!"); return 0.;}

	/**
	 * @brief return the histogram edges of the second variable
	 *
	 * @return the histogram edges of the second variabl
	 */
	virtual std::vector<double> edges2 () const
	{ErrorCBL("Error in histogram_edges2 of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the histogram
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	virtual double operator() ( const int i, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); return 0.;}

	/**
	 * @brief return the histogram at (i,j)
	 *
	 * @param i the i-th first variable bin
	 *
	 * @param j the j-th second variable bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram at (i,j)
	 */
	virtual double operator() ( const int i, const int j, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)j; (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); return 0.;}


	/**
	 * @brief return the histogram
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	virtual std::vector<double> operator() (const HistogramType hist_type, const double fact=1.) const
	{ (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the bin normalization
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the bin normalization
	 */
	virtual double normalization (const int i, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)hist_type; (void)fact; ErrorCBL("Error in normalization of Histogram.h!"); return 0.;}

	/**
	 * @brief return the bin normalization
	 *
	 * @param i i-th bin
	 *
	 * @param j j-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the bin normalization
	 */
	virtual double normalization (const int i, const int j, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)j; (void)hist_type; (void)fact; ErrorCBL("Error in normalization of Histogram.h!"); return 0.;}

	/**
	 * @brief return the bin weight
	 *
	 * @param i i-th bin
	 *
	 * @return the bin weight
	 */
	virtual double weight (const int i) const
	{ (void)i; ErrorCBL("Error in weight of Histogram.h!"); return 0.;}

	/**
	 * @brief return the weights
	 *
	 * @return the weights
	 */
	virtual std::vector<double> weights () const
	{ ErrorCBL("Error in weights of Histogram.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return the bin weight
	 *
	 * @param i i-th bin
	 *
	 * @param j j-th bin
	 *
	 * @return the bin normalization
	 */
	virtual double weight (const int i, const int j) const
	{ (void)i; (void)j; ErrorCBL("Error in weight of Histogram.h!"); return 0.;}

	/**
	 * @brief return the poisson error of
	 * the histogram
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	virtual double poisson_error ( const int i, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); return 0.;}

	/**
	 * @brief return the poisson error of
	 * the histogram at (i,j)
	 *
	 * @param i the i-th first variable bin
	 *
	 * @param j the j-th second variable bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram at (i,j)
	 */
	virtual double poisson_error ( const int i, const int j, const HistogramType hist_type, const double fact=1.) const
	{ (void)i; (void)j; (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); return 0.;}


	/**
	 * @brief return the poisosn error 
	 * of t he histogram
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	virtual std::vector<double> poisson_error ( const HistogramType hist_type, const double fact=1.) const
	{ (void)hist_type; (void)fact; ErrorCBL("Error in operator() of Histogram.h!"); std::vector<double> vv; return vv;}

	///@}

	/**
	 *  @name input/output member functions (customized in all the derived classes)
	 */
	///@{

	/**
	 * @brief write the histogram
	 *  
	 * @param dir output directory
	 *
	 * @param file output file
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return none
	 */
	virtual void write (const std::string dir, const std::string file, const HistogramType hist_type, const double fact=1.) const
	{ (void)dir; (void)file; (void)hist_type; (void)fact; ErrorCBL("Error in write of Histogram.h!");}


	///@}

    };

    /**
     *  @class Histogram Histogram.h "Headers/Histogram.h"
     *
     *  @brief The class Histogram
     *
     *  This class is used to bin 1D variable.
     */
    class Histogram1D : public Histogram
    {

      private:
	
	/// GSL histogram
	std::shared_ptr<gsl_histogram> m_histo;
	
	/// histogram weights
	std::vector<double> m_weight;

	/// variable bins
	std::vector<double> m_bins;

	/// variable edges
	std::vector<double> m_edges;

	/// the number of bins
	size_t m_nbins;

	/// the bin size
	double m_binSize;

	/// the shift of the bin
	double m_shift;

	/// minimum var value
	double m_minVar;

	/// maximum var value
	double m_maxVar;

	/// bin type
	BinType m_binType;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class Histogram1D
	 */
	Histogram1D () = default;

	/**
	 *  @brief constructor
	 *
	 *  @param var the variable values
	 *
	 *  @param weight the variable weights
	 *
	 *  @param nbins the number of bins
	 *
	 *  @param minVar minimum range 
	 *  
	 *  @param maxVar maximmum range
	 *  
	 *  @param shift the shift of the bin
	 *  
	 *  @param bin_type the binning type
	 *
	 *  @return object of class Histogram1D
	 */
	Histogram1D (const std::vector<double> var, const std::vector<double> weight, const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift = 0.5,  const BinType bin_type=BinType::_linear_);

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Histogram1D () = default;

	///@}

	/**
	 *  @name Functions to set the private members of the class
	 */
	///@{

	/**
	 * @brief set the histogram variables
	 *
	 * @param nbins the number of bins
	 * @param minVar the variable minimum 
	 * @param maxVar the variable maximum 
	 * @param shift the shift of the bin
	 * @param bin_type the binning type
	 *
	 * @return none
	 */
	void set (const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift=0.5, const BinType bin_type=BinType::_linear_) override;

	/**
	 * @brief get the histogram index
	 *
	 * @param var value of the var
	 *
	 * @return the histogram index
	 */
	int digitize (const double var) override;

	/**
	 * @brief get the histogram indeces
	 *
	 * @param var values of the var
	 *
	 * @return the histogram indeces
	 */
	std::vector<int> digitize (const std::vector<double> var) override;

	/**
	 * @brief bin the data
	 *
	 * @param var value of the var
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	void put (const double var, const double weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param var values of the var
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	void put (const std::vector<double> var, const std::vector<double> weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param bin value of the bin
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	void put (const int bin, const double weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param bin value of the bin
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	void put (const std::vector<int> bin, const std::vector<double> weight) override;

	///@}

	/**
	 *  @name Functions to get the private members of the class
	 */
	///@{

	/**
	 * @brief return the number
	 * of bins
	 *
	 * @return the number of bins
	 */
	size_t nbins () const override { return m_nbins; }

	/**
	 * @brief return the bin size
	 *
	 * @return the bin size
	 */
	double binSize () const override { return m_binSize; }

	/**
	 * @brief return the lower limit of the
	 * histogram
	 *
	 * @return the lower limit of the histogram
	 */
	double minVar() const override {return m_minVar;}
	
	/**
	 * @brief return the upper limit of the
	 * histogram
	 *
	 * @return the lower upper of the histogram
	 */
	double maxVar() const override {return m_maxVar;}

	/**
	 * @brief return the bin shift
	 *
	 * @return the bin shift
	 */
	double shift() const override {return m_shift;}

	/**
	 * @brief return the bin type
	 *
	 * @return the bin type
	 */
	BinType bin_type() const override {return m_binType;}

	/**
	 * @brief return the i-th bin
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin
	 */
	double bin (const size_t i) const override { return m_bins[i]; }

	/**
	 * @brief return the bins
	 *
	 * @return the histogram bins
	 */
	std::vector<double> bins () const override { return m_bins; }

	/**
	 * @brief return the i-th edge
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge
	 */
	double edge (const size_t i) const override { return m_edges[i]; }

	/**
	 * @brief return the histogram edges
	 *
	 * @return the histogram edges
	 */
	std::vector<double> edges () const override { return m_edges; }

	/**
	 * @brief return the bin normalization
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the bin normalization
	 */
	double normalization (const int i, const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the bin weight
	 *
	 * @param i i-th bin
	 *
	 * @return the bin weight
	 */
        double weight (const int i) const override {return m_weight[i]; }

	/**
	 * @brief return the bin weights
	 *
	 * @return the bin weights
	 */
        std::vector<double> weights () const override {return m_weight; }

	/**
	 * @brief return the histogram
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	double operator() ( const int i, const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the histogram
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	std::vector<double> operator() ( const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the histogram
	 *
	 * @param i i-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	double poisson_error ( const int i, const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the poisson error of the
	 * histogram
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram
	 */
	std::vector<double> poisson_error ( const HistogramType hist_type, const double fact=1.) const override;

	///@}
	
	/**
	 *  @name input/output member functions (customized in all the derived classes)
	 */
	///@{

	/**
	 * @brief write the histogram
	 *  
	 * @param dir output directory
	 *
	 * @param file output file
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return none
	 */
	void write (const std::string dir, const std::string file, const HistogramType hist_type, const double fact=1.) const override;

	///@}

    };

    /**
     *  @class Histogram2D Histogram.h "Headers/Histogram.h"
     *
     *  @brief The class Histogram2D
     *
     *  This class is used to bin 2D variables.
     */
    class Histogram2D : public Histogram
    {
      private:
	
	/// GSL histogram
	std::shared_ptr<gsl_histogram2d> m_histo;
	
	/// histogram weights
	std::vector<std::vector<double>> m_weight;

	/// first variable bins
	std::vector<double> m_bins1;

	/// first variable edges
	std::vector<double> m_edges1;

	/// the number of bins for the first variable
	size_t m_nbins1;

	/// the binSize for the first variable
	double m_binSize1;

	/// the bin shift for the first variable
	double m_shift1;

	/// minimum first variable value
	double m_minVar1;

	/// maximum first variable value
	double m_maxVar1;

	/// first variable bin type
	BinType m_binType1;

	/// variable bins
	std::vector<double> m_bins2;

	/// variable edges
	std::vector<double> m_edges2;

	/// the number of bins
	size_t m_nbins2;

	/// the binSize for the second variable
	double m_binSize2;

	/// the bin shift for the second variable
	double m_shift2;

	/// minimum second variable value
	double m_minVar2;

	/// maximum second variable value
	double m_maxVar2;

	/// second variable bin type
	BinType m_binType2;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class Histogram2D
	 */
	Histogram2D () = default;

	/**
	 * @brief constructor
	 *
	 * @param var1 values of the first var
	 *
	 * @param var2 values of the second var
	 * 
	 * @param weight weights of the var
	 * 
	 * @param nbins1 the number of bins for the first variable
	 *
	 * @param nbins2 the number of bins for the second variable
	 *
	 * @param minVar1 minimum range  for the first variable
	 * 
	 * @param maxVar1 maximmum range for the first variable
	 * 
	 * @param minVar2 minimum range for the second variable 
	 * 
	 * @param maxVar2 maximmum range for the second variable
	 * 
	 * @param shift1 bin shift for the first variable
	 * 
	 * @param shift2 bin shift for the second variable
	 * 
	 * @param bin_type1 the binning type for the first variable
	 * 
	 * @param bin_type2 the binning type for the second variable
	 *
	 * @return object of class Histogram1D
	 */
	Histogram2D (const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> weight, const size_t nbins1, const size_t nbins2, const double minVar1=par::defaultDouble, const double maxVar1=par::defaultDouble, const double minVar2=par::defaultDouble, const double maxVar2=par::defaultDouble, const double shift1=0.5, const double shift2=0.5, const BinType bin_type1=BinType::_linear_, const BinType bin_type2=BinType::_linear_);

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Histogram2D () = default;

	///@}

	/**
	 *  @name Functions to set the private members of the class
	 */
	///@{

	/**
	 * @brief set the histogram variables
	 *
	 * @param nbins1 the number of bins for the first variable
	 * 
	 * @param nbins2 the number of bins for the second variable
	 *
	 * @param minVar1 minimum range  for the first variable
	 *
	 * @param maxVar1 maximmum range for the first variable
	 * 
	 * @param minVar2 minimum range for the second variable 
	 *
	 * @param maxVar2 maximmum range for the second variable
	 * 
	 * @param shift1 bin shift for the first variable
	 * 
	 * @param shift2 bin shift for the second variable
	 * 
	 * @param bin_type1 the binning type for the first variable
	 * 
	 * @param bin_type2 the binning type for the second variable
	 *
	 * @return none
	 */
	void set (const size_t nbins1, const size_t nbins2, const double minVar1=par::defaultDouble, const double maxVar1=par::defaultDouble, const double minVar2=par::defaultDouble, const double maxVar2=par::defaultDouble, const double shift1=0.5, const double shift2=0.5, const BinType bin_type1=BinType::_linear_, const BinType bin_type2=BinType::_linear_) override;

	/**
	 * @brief get the histogram index
	 *
	 * @param var1 value of the first var
	 * @param var2 value of the second var
	 *
	 * @return the histogram index
	 */
	std::vector<int> digitize (const double var1, const double var2) override;

	/**
	 * @brief get the histogram indeces
	 *
	 * @param var1 values of the first var
	 * @param var2 values of the second var
	 *
	 * @return the histogram indeces
	 */
	std::vector<std::vector<int>> digitize (const std::vector<double> var1, const std::vector<double> var2) override;

	/**
	 * @brief bin the data
	 *
	 * @param var1 value of the first var
	 *
	 * @param var2 value of the second var
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	void put (const double var1, const double var2, const double weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param var1 values of the first var
	 *
	 * @param var2 values of the second var
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	void put (const std::vector<double> var1, const std::vector<double> var2, const std::vector<double> weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param bin1 value of the first bin
	 *
	 * @param bin2 value of the second bin
	 *
	 * @param weight weight of the var
	 *
	 * @return none
	 */
	void put (const int bin1, const int bin2, const double weight) override;

	/**
	 * @brief bin the data
	 *
	 * @param bins std::vector containing
	 * the values of the bins
	 *
	 * @param weight weights of the var
	 *
	 * @return none
	 */
	void put (const std::vector<std::vector<int>> bins, const std::vector<double> weight) override;

	///@}

	/**
	 *  @name Functions to get the private members of the class
	 */
	///@{

	/**
	 * @brief return the number
	 * of bins for the first variable
	 *
	 * @return the number of bins for the first variable
	 */
	size_t nbins1 () const override { return m_nbins1; }

	/**
	 * @brief return the bin size for the first 
	 * variable
	 *
	 * @return the bin size for the first variable
	 */
	double binSize1 () const override { return m_binSize1; }

	/**
	 * @brief return the lower limit of the
	 * histogram for the first variable
	 *
	 * @return the lower limit of the histogram
	 * for the first variable
	 */
	double minVar1 () const override { return m_minVar1; }
	
	/**
	 * @brief return the upper limit of the
	 * histogram for the first variable
	 *
	 * @return the upper limit of the histogram
	 * for the first variable
	 */
	double maxVar1 () const override { return m_maxVar1; }

	/**
	 * @brief return the bin shift for the
	 * first variable
	 *
	 * @return the bin shift for the first variable
	 */
	double shift1 () const override { return m_shift1; }

	/**
	 * @brief return the bin type for the first
	 * variable 
	 *
	 * @return the bin type for the first variable
	 */
	BinType bin_type1 () const override {return m_binType1;}

	/**
	 * @brief return the number
	 * of bins for the second variable
	 *
	 * @return the number of bins for the second variable
	 */
	size_t nbins2 () const override { return m_nbins2; }

	/**
	 * @brief return the bin size for the second 
	 * variable
	 *
	 * @return the bin size for the second variable
	 */
	double binSize2 () const override { return m_binSize2; }

	/**
	 * @brief return the lower limit of the
	 * histogram for the second variable
	 *
	 * @return the lower limit of the histogram
	 * for the second variable
	 */
	double minVar2 () const override {return m_minVar2;}
	
	/**
	 * @brief return the upper limit of the
	 * histogram for the second variable
	 *
	 * @return the upper limit of the histogram
	 * for the second variable
	 */
	double maxVar2 () const override {return m_maxVar2;}

	/**
	 * @brief return the bin shift for the
	 * second variable
	 *
	 * @return the bin shift for the second variable
	 */
	double shift2 () const override {return m_shift2;}

	/**
	 * @brief return the bin type for the second
	 * variable 
	 *
	 * @return the bin type for the second variable
	 */
	BinType bin_type2 () const override {return m_binType2;}

	/**
	 * @brief return the i-th bin of the first variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin of the first variable
	 */
	double bin1 (const size_t i) const override { return m_bins1[i]; }

	/**
	 * @brief return the first variable bins
	 *
	 * @return the first variable bins
	 */
	std::vector<double> bins1 () const override { return m_bins1; }

	/**
	 * @brief return the i-th edge of the first variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge of the first variable
	 */
	double edge1 (const size_t i) const override { return m_edges1[i]; }

	/**
	 * @brief return the histogram edges of the first variable
	 *
	 * @return the histogram edges of the first variable
	 */
	std::vector<double> edges1 () const override { return m_edges1; }

	/**
	 * @brief return the i-th bin of the second variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th bin of the second variable
	 */
	double bin2 (const size_t i) const override { return m_bins2[i]; }

	/**
	 * @brief return the second variable bins
	 *
	 * @return the second variable bins
	 */
	std::vector<double> bins2 () const override { return m_bins2; }

	/**
	 * @brief return the i-th edge of the second variable
	 *
	 * @param i the i-th index 
	 *
	 * @return the i-th edge of the second variabl
	 */
	double edge2 (const size_t i) const override { return m_edges2[i]; }

	/**
	 * @brief return the histogram edges of the second variable
	 *
	 * @return the histogram edges of the second variabl
	 */
	std::vector<double> edges2 () const override { return m_edges2; }

	/**
	 * @brief return the bin normalization
	 *
	 * @param i i-th bin
	 *
	 * @param j j-th bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the bin normalization
	 */
	double normalization (const int i, const int j, const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the bin weight
	 *
	 * @param i i-th bin
	 *
	 * @param j j-th bin
	 *
	 * @return the bin weight
	 */
        double weight (const int i, const int j) const override {return m_weight[i][j];}

	/**
	 * @brief return the histogram at (i,j)
	 *
	 * @param i the i-th first variable bin
	 *
	 * @param j the j-th second variable bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram at (i,j)
	 */
	double operator() ( const int i, const int j, const HistogramType hist_type, const double fact=1.) const override;

	/**
	 * @brief return the poisson error
	 * of the histogram at (i,j)
	 *
	 * @param i the i-th first variable bin
	 *
	 * @param j the j-th second variable bin
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return the histogram at (i,j)
	 */
	double poisson_error (const int i, const int j, const HistogramType hist_type, const double fact=1.) const override;

	///@}

	/**
	 *  @name input/output member functions (customized in all the derived classes)
	 */
	///@{

	/**
	 * @brief write the histogram
	 *  
	 * @param dir output directory
	 *
	 * @param file output file
	 *
	 * @param hist_type the type of histogram
	 *
	 * @param fact the factor used to normalized the
	 * histogram
	 *
	 * @return none
	 */
	void write (const std::string dir, const std::string file, const HistogramType hist_type, const double fact=1.) const override;

	///@}

    };

  }
}

#endif
