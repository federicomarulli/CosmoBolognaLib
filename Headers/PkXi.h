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
 *******************************************************************/

/**
 *  @file Headers/PkXi.h
 *
 *  @brief The class PkXi
 *
 *  This file defines the interface of the class PkXi, used to model
 *  the power spectrum and two-point correlation function of
 *  extra-galactic sources
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __PKXI__
#define __PKXI__

#include "Cosmology.h"


// ===================================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class PkXi PkXi.h "Headers/PkXi.h"
     *
     *  @brief The class PkXi
     *
     *  This class is used to handle objects of type <EM> PkXi
     *  </EM>. It is used to model the power spectrum and two-point
     *  correlation function of extra-galactic sources
     */
    class PkXi {

    private:

      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// the normalisation of the power spectrum for Eisenstein & Hu [http://background.uchicago.edu/~whu/transfer/transferpage.html]
      double m_Pk0_EH = 1.;

      /// the normalisation of the power spectrum for CAMB [http://camb.info/]
      double m_Pk0_CAMB = 1.;

      /// the normalisation of the power spectrum for MPTbreeze [http://arxiv.org/abs/1207.1465]
      double m_Pk0_MPTbreeze = 1.; 

      /// the normalisation of the power spectrum for CLASS [http://class-code.net/]
      double m_Pk0_CLASS = 1.; 
      
      /**
       *  @brief write and read the table where the dark matter power
       *  spectrum, computed with either CAMB or MPTbreeze, is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/] or MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]

       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param [in] redshift redshift
       *  
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       */
      void m_Table_Pk_CAMB_MPTbreeze (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100.) const;

      /**
       *  @brief write and read the table where the dark matter power
       *  spectrum, computed with either CAMB or MPTbreeze, is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/] or MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]

       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of vectors containing the log(k) at
       *  each redshift
       *
       *  @param [out] lgPk vector of vectors containing the log(P(k)) at
       *  each redshift
       *
       *  @param [in] redshift vector of redshifts
       *  
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       */
      void m_Table_Pk_CAMB_MPTbreeze (const std::string code, const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100.) const;

      /**
       *  @brief write and read the table where the dark matter power
       *  spectrum computed with CLASS is stored
       *
       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param [in] redshift redshift
       *  
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       */
      void m_Table_Pk_CLASS (const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100.) const;

      /**
       *  @brief write and read the table where the dark matter power
       *  spectrum computed with CLASS is stored
       *
       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of vectors containing the log(k) at
       *  each redshift
       *
       *  @param [out] lgPk vector of vectors containing the log(P(k)) at
       *  each redshift
       *
       *  @param [in] redshift vector of redshifts
       *  
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       */
      void m_Table_Pk_CLASS (const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100.) const;
      
      /**
       *  @brief write and read the table where the dark matter power
       *  spectrum is stored; it is used when a parameter file is
       *  provided in input
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       *
       *  @warning the input output_root parameter must be the same as
       *  the one in the parameter file
       */
      void m_Table_Pk_parameterFile (const std::string code, const std::string file_par, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const std::string output_root="test") const;

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      PkXi () = default;
      
      /**
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      PkXi (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}
      
      /**
       *  @brief default destructor
       */
      ~PkXi () = default;

      ///@}


      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }

      /**
       *  @brief Get the private member Cosmology::m_Pk0_EH
       *
       *  @return the normalisation of the power spectrum for
       *  Eisenstein & Hu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       */
      double Pk0_EH () const { return m_Pk0_EH; };

      /**
       *  @brief Get the private member Cosmology::m_Pk0_CAMB
       *
       *  @return the normalisation of the power spectrum for CAMB
       *  [http://camb.info/]
       */
      double Pk0_CAMB () const { return m_Pk0_CAMB; };

      /**
       *  @brief Get the private member Cosmology::m_Pk0_MPTbreeze 
       *
       *  @return the normalisation of the power spectrum for MPTbreeze
       *  [http://arxiv.org/abs/1207.1465]
       */
      double Pk0_MPTbreeze () const { return m_Pk0_MPTbreeze; };

      /**
       *  @brief Get the private member Cosmology::m_Pk0_CLASS
       *
       *  @return the normalisation of the power spectrum for CLASS [http://class-code.net/]
       */
      double Pk0_CLASS () const { return m_Pk0_CLASS; };

      ///@}

      
      /**
       *  @name Functions to estimate the power spectrum and two-point
       *  correlation function models
       */
      ///@{ 
      
      /**
       *  @brief The path to the power spectrum output
       *
       *  @param code method used to compute the power spectrum
       * 
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       * 
       *  @param redshift the redshift
       *
       *  @param run true \f$\rightarrow\f$ write or read the table
       *  where the dark matter power spectrum is stored
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param k_max the maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return the path to the power spectrum output
       */
      std::string Pk_output_file (const std::string code, const bool NL, const double redshift, const bool run=0, const bool store_output=true, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString);

      /**
       *  @brief run CAMB [http://camb.info/]
       *  
       *  this function runs CAMB [http://camb.info/], after editing the parameter file
       *  appropriately (if file_par=NULL)
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If NULL, the 
       *  output will be deleted after running CAMB
       *
       *  @param output_dir std::string containing the output directory
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       */
      void run_CAMB (const bool NL, const double redshift, const std::string output_root=par::defaultString, const std::string output_dir=par::defaultString, const double k_max=100., const std::string file_par=par::defaultString) const;

      /**
       *  @brief run CAMB [http://camb.info/] and read the matter
       *  power spectrum
       *  
       *  this function runs CAMB [http://camb.info/], after editing
       *  the parameter file appropriately (if file_par=NULL) and
       *  store the matter power spectrum in two vectors
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If NULL, the 
       *  output will be deleted after running CAMB
       *
       *  @param output_dir std::string containing the output directory
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       */
      void run_CAMB (std::vector<double> &lgkk, std::vector<double> &lgPk, const bool NL, const double redshift, const std::string output_root="test", const std::string output_dir=par::defaultString, const double k_max=100., const std::string file_par=par::defaultString) const;
       /**
       *  @brief write or read the table where the dark matter power
       *  spectrum is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]
       *
       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param [in] redshift redshift
       *  
       *  @param [in] store_output if true the output files created
       *  are stored; if false the output files created are removed
       *
       *  @param [in] output_root output_root of the parameter file
       *  used to compute the power spectrum; it can be any name. If
       *  this parameter is different from the default value it will
       *  be used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       */
      void Table_PkCodes (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const;

      
      /**
       *  @brief write or read the table where the dark matter power
       *  spectrum is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]
       *
       *  @param [in] NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of vectors containing the log(k) at
       *  each redshift
       *
       *  @param [out] lgPk vector of vectors containing the log(P(k))
       *  at each redshift
       *
       *  @param [in] redshift vector of redshifts
       *  
       *  @param [in] store_output if true the output files created
       *  are stored; if false the output files created are removed
       *
       *  @param [in] output_root output_root of the parameter file
       *  used to compute the power spectrum; it can be any name. If
       *  this parameter is different from the default value it will
       *  be used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       */
      void Table_PkCodes (const std::string code, const bool NL, std::vector<std::vector<double>> &lgkk, std::vector<std::vector<double>> &lgPk, const std::vector<double> redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const;
      
           
      /**
       *  @brief remove the output generated by the methods CAMB,
       *  MPTbreeze or CLASS
       *
       *  @param code method used to compute the power spectrum; valid
       *  codes are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465]

       *  @param NL false \f$\rightarrow\f$ linear power
       *  spectrum; true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param redshift redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       */
      void remove_output_Pk_tables (const std::string code, const bool NL, const double redshift, const std::string output_root="test") const;

      /**
       *  @brief write or read the table where the dark matter two-point
       *  correlation function is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] rr vector of comoving separations
       *
       *  @param [out] xi vector of the binned values of &xi;(r)
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file
       *  used to compute the power spectrum; it can be any name. If
       *  this parameter is different from the default value it will
       *  be used also in the output directory name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       */
      void Table_XiCodes (const std::string code, const bool NL, std::vector<double> &rr, std::vector<double> &xi, const double redshift, const bool store_output, const std::string output_root, const double k_max, std::string file_par) const;

      /**
       *  @brief normalisation of the power spectrum
       *
       *  this function sets the value of the private member m_Pk0_*,
       *  i.e. the normalisation of the power spectrum
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CAMB_wrapper (running CAMB wrapper), CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default
       *  value it will be used also in the output directory name
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation;
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation;
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       */
      void Pk_0 (const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      


      /**
       *  @brief the dark matter power spectrum
       *
       *  this function provides the dark matter power spectrum;
       *  it can use either CAMB, CLASS, MPTbreeze or the analytic
       *  approximation by Eisenstein & Hu
       *
       *  @param kk vector of wave vector modules
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CAMB_wrapper (running CAMB wrapper), CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default
       *  value it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the dark matter power spectrum 
       */
      std::vector<double> Pk_matter (const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      /**
       *  @brief The dark matter power spectrum
       *
       *  This function computes the dark matter power spectrum; it
       *  can use either CAMB, CLASS, MPTbreeze or the analytic
       *  approximation by Eisenstein & Hu
       *
       *  @param kk vector of wave vector modules
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param redshift vector of redshifts
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default
       *  value it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the dark matter power spectrum 
       */
      std::vector<std::vector<double>> Pk_matter (const std::vector<double> kk, const std::string method_Pk, const bool NL, const std::vector<double> redshift, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);

      /**
       *  @brief The dark matter two-point correlation function
       *
       *  This function computes the dark matter correlation function,
       *  obtained by Fourier transforming the matter power spectrum
       *
       *  @param rr the module of the comoving separation
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL false \f$\rightarrow\f$ linear DM two-point
       *  correlation function; true \f$\rightarrow\f$ non-linear DM
       *  two-point correlation function
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, const ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>DM</SUB>(r): the spherically
       *  averaged (monopole) of the two-point correlation function of
       *  dark matter
       */
      double xi_matter (const double rr, const std::string method_Pk, const bool NL, const double redshift, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the dark matter projected correlation function
       *
       *  this function provides the dark matter projected correlation
       *  functions, obtained by Fourier transforming the matter power
       *  spectrum
       *
       *  @param rp r<SUB>p</SUB>: projected separation
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL false \f$\rightarrow\f$ linear DM two-point
       *  correlation function; true \f$\rightarrow\f$ non-linear DM
       *  two-point correlation function
       *
       *  @param redshift the redshift
       *
       *  @param pimax the upper limit of the line-of-sight integration
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return w<SUB>p,DM</SUB>(&theta;): the projected correlation
       *  function of dark matter
       */
      double wp_DM (const double rp, const std::string method_Pk, const bool NL, const double redshift, const double pimax, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=cbl::par::defaultString);

      
      /**
       *  @brief the dark matter rms mass fluctuation
       *
       *  @param RR radius inside which the dark matter rms mass
       *  fluctuation is computed
       *
       *  @param corrType 0 \f$\rightarrow\f$ the projected correlation function,
       *  w(&theta;), is used; 1 \f$\rightarrow\f$ the spherically averaged
       *  correlation function, &xi;(r), is used
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param pimax the upper limit of the line-of-sight integration
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *  
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &sigma;<SUB>R</SUB>: the dark matter rms mass
       *  fluctuation
       */
      double sigmaR_DM (const double RR, const int corrType, const std::string method_Pk, const double redshift, const double pimax=40, const bool store_output=true, const std::string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      


      /**
       *  @brief the dark matter rms mass fluctuation within 8 Mpc/h
       *
       *  this function provides the rms mass fluctuation within 8
       *  Mpc/h, estimated directly from the power spectrum
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *  
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &sigma;<SUB>8</SUB>: the dark matter rms mass
       *  fluctuation within 8 Mpc/h
       */
      double sigma8_Pk (const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool NL=0, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString) const; 
      
      /**
       *  @brief the k<SUB>*</SUB> parameter 
       *
       *  this function provides the k<SUB>*</SUB> parameter used to
       *  model the BAO (see e.g. Chuang & Wang 2012, Crocce et
       *  al. 2006, Matsubara 2008)
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return k<SUB>*</SUB>
       */
      double k_star (const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const;

      
      /**
       *  @brief Get the dark matter two-point correlation function
       *
       *  this function provides the dark matter correlation function,
       *  obtained by Fourier transforming the matter power spectrum
       *
       *  @param [out] rr vector of r, the module of the comoving
       *  separation
       *
       *  @param [out] Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param [in] method_Pk method used to compute the power
       *  spectrum (i.e. the Boltzmann solver); valid choices for
       *  method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] store_output if true the output files created by
       *  the Boltzmann solver are stored; if false the output files
       *  are removed
       *
       *  @param [in] output_root output_root of the parameter file used
       *  to compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param [in] xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param [in] xiNL 0 \f$\rightarrow\f$ linear power spectrum;
       *  1 \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [in] norm 0 \f$\rightarrow\f$ don't normalise the
       *  power spectrum; 1 \f$\rightarrow\f$ normalise the power
       *  spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is
       *  set
       *
       *  @param [in] r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param [in] GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param [in] prec accuracy of the integration
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, const ignoring the cosmological parameters of the object
       */
      void get_xi (std::vector<double> &rr, std::vector<double> &Xi, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0.001, const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief Get the barred dark matter correlation functions
       *
       *  this function provides the dark matter \e barred correlation
       *  functions, used to model the two-point correlation function
       *  in redshift-space
       *
       *  @param [in] rr vector of r, the module of the comoving
       *  separation
       *
       *  @param [in] Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param [out] Xi_ vector of barred &xi;(r),
       *
       *  @param [out] Xi__ vector of double-barred &xi;(r)
       *
       *  @param [in] method_Pk method used to compute the power
       *  spectrum (i.e. the Boltzmann solver); valid choices for
       *  method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] xiType 0 \f$\rightarrow\f$ standard; 1 \f$\rightarrow\f$ Chuang & Wang
       *  model
       *
       *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param [in] xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param [in] norm 0 \f$\rightarrow\f$ don't normalise the
       *  power spectrum; 1 \f$\rightarrow\f$ normalise the power
       *  spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is
       *  set
       *
       *  @param [in] r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param [in] prec accuracy of the integration
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       */
      void get_barred_xi (std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string method_Pk, const double redshift, const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0.001, const double k_max=100., const double aa=0., const double prec=1.e-2, const std::string file_par=par::defaultString) const;

      /**
       *  @brief The dark matter power spectrum without BAO wiggles
       *
       *  This function computes the No Wiggles dark matter power
       *  spectrum. It follows the method proposed in Vlah et al. 2015
       *  (https://arxiv.org/abs/1509.02120, Appendix A).  The no
       *  wiggles power spectrum is obtained by smoothing the
       *  oscillatory part of the linear power spectrum with gaussian
       *  filter.
       * 
       *  It implements two methods:
       *
       *  - Gaussian 3D:
       *
       *  \f[ P_{\mathrm{nw}}(k) = \int d^{3}q \, P(q)
       *  \mathcal{F}_{G}(|\mathbf{k} - \mathbf{q}|) =
       *  \frac{\sqrt{2}}{\sqrt{\pi} \lambda} \int dq \, q^{2} P(q)
       *  \exp\left(-\frac{1}{2\lambda^{2}} \left(q^{2} +
       *  k^{2}\right)\right) \frac{\sinh\left(\frac{k
       *  q}{\lambda^{2}}\right)}{k q} \f]
       *
       *  - Gaussian 1D:
       *
       *  \f[
       *    P_{\mathrm{nw}}\left(10^{k_{\mathrm{log}}}\right)=\frac{1}{\sqrt{2
       *    \pi} \lambda} \int d q_{\mathrm{log}}
       *    P\left(10^{q_{\mathrm{log}}}\right) \exp \left(-\frac{1}{2
       *    \lambda^{2}}\left(k_{\mathrm{log}}-q_{\mathrm{log}}\right)^{2}\right)
       *    \f]
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param kk array containing the wave vector module
       *
       *  @param PkLin array that contains the linear power spectrum
       *  
       *  @param PkApprox array that contains the approximated no-wiggle power spectrum
       *
       *  @param lambda size of the kernel
       *
       *  @param method gaussian smoothing method; it can be 
       *  "gaussian_1d" or "gaussian_3d"
       *
       *  @return P;<SUB>NW</SUB>(k): the No-Wiggle par of the power
       *  spectrum of dark matter
       */
      std::vector<double> Pk_matter_NoWiggles_gaussian (const std::vector<double> kk, const std::vector<double> PkLin, const std::vector<double> PkApprox, const double lambda, const std::string method);

      /**
       *  @brief The dark matter power spectrum without BAO wiggles
       *
       *  This function computes the No Wiggles dark matter power
       *  spectrum. It follows the method proposed in Vlah et al. 2015
       *  (https://arxiv.org/abs/1509.02120, Appendix A).  The no
       *  wiggles power spectrum is obtained by interpolating the
       *  oscillatory part of the linear power spectrum with a basis
       *  spline of a given order and number of knots.
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param kk array containing the wave vector module
       *
       *  @param PkLin array that contains the linear power spectrum
       *  
       *  @param PkApprox array that contains the approximated no-wiggle power spectrum
       *
       *  @param order basis spline order
       *
       *  @param nknots number of knots
       *
       *  @return P;<SUB>NW</SUB>(k): the No-Wiggle par of the power
       *  spectrum of dark matter
       */
      std::vector<double> Pk_matter_NoWiggles_bspline (const std::vector<double> kk, const std::vector<double> PkLin, const std::vector<double> PkApprox, const int order, const int nknots);
      
      /**
       *  @brief The dark matter power spectrum without BAO wiggles
       *
       *  This function computes the no-wiggles dark matter power
       *  spectrum. It can be computed with: - EisensteinHu
       *  approximate formulas [http://arxiv.org/abs/1207.1465]: -
       *  bspline interpolation (see
       *  cbl::cosmology::PkXi::Pk_matter_NoWiggles_bspline).  -
       *  gaussian 3d smoothing (see
       *  cbl::cosmology::PkXi::Pk_matter_NoWiggles_gaussian).  -
       *  gaussian 1d smoothing (see
       *  cbl::cosmology::PkXi::Pk_matter_NoWiggles_gaussian).
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param method method to obtain power spectrum
       *  with no wiggles. It can be "EisensteinHu", "bspline"
       *  "gaussian_3d" or "gaussian_1d"
       *
       *  @param kk array containing the wave vector module
       *
       *  @param redshift the redshift
       *
       *  @param linear_method method to compute the linear
       *  power spectrum. It can be "CAMB" or "CLASS"
       *  This is used only when method=="bspline"
       *
       *  @param order basis spline order. 
       *  This is used only when method=="bspline"
       *
       *  @param nknots number of knots
       *  This is used only when method=="bspline"
       *
       *  @param lambda width of the gaussian filter
       *  This is used only when method=="gaussian_3d"
       *  or "gaussian_1d"
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @return P;<SUB>DW</SUB>(k): the De-Wiggled power
       *  spectrum of dark matter
       */
      std::vector<double> Pk_matter_NoWiggles (const std::string method, const std::vector<double> kk, const double redshift, const std::string linear_method="CAMB", const int order=4, const int nknots=10, const double lambda=0.25, const bool store_output=true, const std::string output_root="test", const bool norm=1, const double prec=1.e-4);

      /**
       *  @brief  The dark matter linear power spectrum
       *
       *  This function computes the linear dark matter power
       *  spectrum.  Valid choices are: CAMB [http://camb.info/],
       *  CLASS [http://class-code.net/]
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param method method to obtain power spectrum
       *  with no wiggles. It can be "CAMB" or "CLASS"
       *
       *  @param kk array containing the wave vector module
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @return P;<SUB>DW</SUB>(k): the De-Wiggled power
       *  spectrum of dark matter
       */
      std::vector<double> Pk_matter_Linear (const std::string method, const std::vector<double> kk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool norm=1, const double prec=1.e-4);

      
      
      /**
       *  @brief the interpolated matter power spectrum, used to compute
       *  \f$\sigma_{\rm R}\f$ and \f$\sigma_{\rm M}\f$ (see cbl::cosmology::Sigma)
       *
       *  this function computes the matter power spectrum and return its interpolation
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CAMB_wrapper (running CAMB wrapper), CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum; it can be any name. If this
       *  parameter is different from the default value it will be
       *  used also in the output directory name
       * 
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param kmax maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return the interpolated matter power spectrum
       */
      std::function<double(double)> Pk_interp_for_sigma (const std::string method_Pk, const double redshift, const bool store_output = true, const std::string output_root = "test", const std::string interpType = "Linear", const double kmax = 100., const std::string input_file = par::defaultString, const bool is_parameter_file = true) const;
      
      /**
       *  @brief The dark matter power spectrum, de-wiggled (see
       *  e.g. Anderson et al 2014)
       *
       *  this function provides the De-Wiggled dark matter power
       *  spectrum
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param linear_method method to obtain linear power spectrum
       *
       *  @param nowiggles_method method to obtain power spectrum
       *  with no wiggles
       *
       *  @param kk array containing the wave vector module
       *
       *  @param redshift the redshift
       *
       *  @param sigma_NL the non linear BAO damping
       *
       *  @param order basis spline order. 
       *  This is used only when method=="bspline"
       *
       *  @param nknots number of knots
       *  This is used only when method=="bspline"
       *
       *  @param lambda width of the gaussian filter
       *  This is used only when method=="gaussian_3d"
       *  or "gaussian_1d"
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @return P;<SUB>DW</SUB>(k): the De-Wiggled power
       *  spectrum of dark matter
       */
      std::vector<double> Pk_matter_DeWiggled (const std::string linear_method, const std::string nowiggles_method, const std::vector<double> kk, const double redshift, const double sigma_NL, const int order=4, const int nknots=10, const double lambda=0.25, const bool store_output=true, const std::string output_root="test", const bool norm=1, const double prec=1.e-4);

      /**
       *  @brief Monopole of the real-space two-point correlation
       *  function
       *
       *  \f[ \xi_0(r) =  b^2 \cdot \xi_{\rm DM}(r) \f]
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param bias the linear bias
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass) (i.e. the Boltzmann solver); valid choices
       *  for method_Pk are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>0</SUB>
       *
       */
      std::vector<double> xi0_bias_realSpace (const std::vector<double> rad, const double bias, const std::string method_Pk, const double redshift, const bool store_output=true, const std::string output_root="test", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      
      /**
       *  @brief the dark matter two-point correlation function,
       *  de-wiggled (see e.g. Anderson et al 2014)
       *
       *  this function provides the dark matter correlation function,
       *  obtained by Fourier transforming the De-Wiggled matter power
       *  spectrum
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param rr the module of the comoving separation
       *
       *  @param redshift the redshift
       *
       *  @param sigma_NL the non linear BAO damping
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param prec accuracy of the integration 
       *
       *  @return &xi;<SUB>DW</SUB>(r): the De-Wiggled spherically
       *  averaged (monopole) of the two-point correlation function of
       *  dark matter
       */
      double xi_matter_DeWiggle (const double rr, const double redshift, const double sigma_NL, const bool store_output=true, const std::string output_root="test", const bool norm=1, const double k_min=0.001, const double k_max=100., const double aa=1., const double prec=1.e-2);
 
       /**
       *  @brief The normalised amplitude of the growing mode at a
       *  given redshift, computed using the linear matter power
       *  spectrum.
       *
       *  This member function computes the square root of the ratio 
       *  between the linear matter power spectrum at given redshift
       *  and the one computed at given redshift_norm:
       *
       *  \f[D(z) = \sqrt{\frac{P_m(z,k)}{P_m(z_n,k)}} \ \ \f]
       *
       *  @param redshift the redshift
       *
       *  @param redshift_norm the redshift_norm
       *
       *  @param kk vector of wave vector modules
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CAMB_wrapper (running CAMB wrapper), CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If this parameter is different from the default
       *  value it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed to estimate the power spectrum
       *  normalisation; this parameter is used only if norm=1
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the amplitude of the growing mode
       */
      double DN_Pk (const double redshift, const double redshift_norm, const double kk, const std::string method_Pk, const bool NL, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);     

      ///@}
      

      /**
       *  @name Functions to estimate the two-point correlation function, bias and related quantities
       */
      ///@{

      /**
       *  @brief The dark matter angular two-point correlation
       *  function
       *
       *  this function provides the dark matter angular correlation function,
       *  obtained by integrating the 2PCF using Limber approximation 
       *  (see Salazar et al. 2014, eqs. 5-8):
       *
       *  \f[ w(\theta) = \int_{z_{min}}^{z_{max}} \mathrm{d}z_1
       *  \int_{z_{min}}^{z_{max}} \mathrm{d}z_2 \phi(z_1) \phi(z_2)
       *  \xi(r,z) \f]
       *
       *  where \f$ r = \sqrt{\chi(z_1)^2+\chi(z_2)^2
       *  -2\chi(z_1)\chi(z_2)\cos(\theta)} \f$, \f$\chi\f$ is the
       *  comoving distance and \f$z = (z_1+z_2)/2 \f$
       *
       *  @param theta the angular separation
       *
       *  @param zz the redshift range
       *
       *  @param phiz the number density
       *
       *  @param interpolationMethod the method in interpolation
       *
       *  @param coordUnits the angular separation units
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param method_Pk method used to compute the power spectrum
       *  (i.e. the Boltzmann solver); valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalize the power
       *  spectrum; 1 \f$\rightarrow\f$ normalize the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return \f$w_{DM}(\theta)\f$: the angular two point
       *  correlation function of dark matter
       */
      double wtheta_DM (const double theta, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const CoordinateUnits coordUnits = CoordinateUnits::_degrees_, const bool GSL=false, const std::string method_Pk="CAMB", const bool NL=false, const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief The dark matter angular two-point correlation
       *  function
       *
       *  This function computes the dark matter angular correlation
       *  function, obtained by integrating the 2PCF using Limber
       *  approximation (see Sawangwit et al. 2011, eqs. 14,15):
       *
       *  \f[ w(\theta) = \frac{\int_{z_{min}}^{z_{max}} \mathrm{d}z_1
       *  \int_{z_{min}}^{z_{max}} \mathrm{d}z_2 f(z_1) f(z_2)
       *  \xi(r,z)} { \left\{ \int_{z_{min}}^{z_{max}} f(z) dz
       *  \right\}^2} \f]
       *
       *  where \f$ r = \sqrt{\chi(z_1)^2+\chi(z_2)^2
       *  -2\chi(z_1)\chi(z_2)\cos(\theta)} \f$, \f$\chi\f$ is the
       *  comoving distance and \f$z = (z_1+z_2)/2 \f$.
       *
       *  The function \f$f(z)\f$ is the number of objects per unitar
       *  volume:
       *
       *  \f[ f(z) = \frac{\mathrm{d}V}{\mathrm{d}z \mathrm{d\Omega}}
       *     n(z) \phi(z) \f]
       *
       *  where \f$\frac{\mathrm{d}V}{\mathrm{d}z \mathrm{d\Omega}}\f$
       *  is the comoving volume element, \f$n(z)\f$ is the comoving
       *  number density and \f$\phi(z)\f$ is the selection function
       *
       *  @param theta the angular separation
       *
       *  @param kk the wave vector module
       *
       *  @param Pk linear power spectrum
       *
       *  @param zz the redshift range
       *
       *  @param nz the comoving number density
       *
       *  @param phiz the selection function
       *
       *  @param interpolationType the method in interpolation
       *
       *  @param coordUnits the angular separation units
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param redshift_Pk the redshift of the input power spectrum
       *
       *  @return \f$w_{DM}(\theta)\f$: the angular two point
       *  correlation function of dark matter
       */
      double wtheta_DM (const double theta, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> zz, const std::vector<double> nz, const std::vector<double> phiz, const std::string interpolationType="Spline", const CoordinateUnits coordUnits = CoordinateUnits::_degrees_, const bool GSL=false, const double redshift_Pk=0);
      
      /**
       * @brief The dark matter angular linear power spectrum
       * \f$C_l\f$
       *
       * This function computes the angular linear power spectrum
       * using the Limber approximation up to a given \f$l_{max}\f$:
       * \f[ C_l = \int_{z_{min}}^{z_{max}}\phi^2(z) D^2(z)
       * P(\frac{l}{D_C(z)}) \frac{H(z)}{c \cdot D_c(z)} dz \f]
       *
       * where l is the multipole order, \f$\phi(z)\f$ is the tracers
       * redshift distribution, \f$P\f$ is the \f$z=0\f$ power
       * spectrum and \f$D(z)\f$ is the growth factor. \f$D_C(z)\f$
       * and \f$H(z)\f$ are respectively the comoving function and the
       * Hubble parameter.
       *
       * @param lmax the maximum multipole order
       *
       * @param zz the redshift range
       *
       * @param phiz the number density
       *
       * @param interpolationMethod the method in interpolation
       *
       * @param method_Pk method used to compute the power spectrum
       * (i.e. the Boltzmann solver); valid choices for method_Pk are:
       * CAMB [http://camb.info/], CLASS
       * [http://class-code.net/], MPTbreeze-v1
       * [http://arxiv.org/abs/1207.1465], EisensteinHu
       * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       * @param store_output if true the output files created by the
       * Boltzmann solver are stored; if false the output files are
       * removed
       *
       * @param output_root output_root of the parameter file used to
       * compute the power spectrum and &sigma;(mass); it can be any
       * name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       * @param norm 0 \f$\rightarrow\f$ don't normalize the power
       * spectrum; 1 \f$\rightarrow\f$ normalize the power spectrum;
       * -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       * @param k_min minimum wave vector module up to which the power
       * spectrum is computed
       *
       * @param k_max maximum wave vector module up to which the power
       * spectrum is computed
       *
       * @param prec accuracy of the integration
       *
       * @param file_par name of the parameter file; if a parameter
       * file is provided (i.e. file_par!=NULL), it will be used,
       * ignoring the cosmological parameters of the object
       *
       * @return vector containing the angular linear power spectrum up to \f$l_{max}\f$
       *
       */
      std::vector<double> C_l_DM (const int lmax, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const std::string method_Pk="CAMB", const bool store_output=true, const std::string output_root="test", const int norm=-1, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const std::string file_par=par::defaultString);

      

      /**
       *  @brief the function &xi;<SUB>*</SUB> of the Chuang & Wang 2012
       *  model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rr comoving separation
       *
       *  @param redshift the redshift
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name. If this parameter is different from the default value
       *  it will be used also in the output directory name
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed in order to estimate the power
       *  spectrum normalisation; this parameter is used only if
       *  either norm=1, or norm=-1 and sigma8 is set
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>*</SUB>
       */
      double xi_star (const double rr, const double redshift, const bool store_output=true, const std::string output_root="test", const double k_star=-1., const double k_min=0.001, const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
      
      ///@}
      
      
      /**
       *  @name Functions to estimate the multipoles/wedges covariance matrix
       */
      ///@{
      
      /**
       *  @brief The covariance matrix of the first three non-null
       *  multipoles of the two-point correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param nn order of the moment
       *
       *  @param Volume the volume
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the covariance matrix of the first three non-null
       *  multipoles of the two-point correlation function
       */
      std::vector<std::vector<double>> XiMonopole_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const int IntegrationMethod=1);

      /**
       *  @brief The covariance matrix of the first three non-null
       *  multipole moments of the two-point correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param nn order of the moment
       *
       *  @param Volume the volume
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param Pk2 vector containing the quadrupole of the power
       *  spectrum
       * 
       *  @param Pk4 vector containing the hexadecapole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the covariance matrix of the first three non-null
       *  multipole moments of the two-point correlation function
       */
      std::vector<std::vector<double>> XiMultipoles_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod=1);

      /**
       *  @brief The first three non-null multipoles of the two-point
       *  correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param Pk2 vector containing the quadrupole of the power
       *  spectrum
       * 
       *  @param Pk4 vector containing the hexadecapole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the matrix containing the first three non-null
       *  multipoles of the two-point correlation function
       */
      std::vector<std::vector<double>> XiMultipoles (const int nbins, const double rMin, const double rMax, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod=1);

      ///@}
      
    };
  }
  

  // =====================================================================================
  
  
  namespace glob {
    
    /// @cond glob
    
    double func_xi_EH_GSL (double, void *);

    struct STR_xi_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double little_h;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      std::string method_Pk;
      double rr;
      double redshift;
      double aa;
    };

    struct STR_sigma2M_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double little_h;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      std::string method_Pk;
      double redshift;
      double mass;
      double rho;
    };

    /// @endcond

  }
  
}

#endif
