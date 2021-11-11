/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Alfonso Veropalumbo     *
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
 *  @file Headers/EisensteinHu.h
 *
 *  @brief The class EisensteinHu 
 *
 *  This file defines the class Eisenstein, used for fitting 
 *  formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. 
 *  See Daniel J. Eisenstein & Wayne Hu, 1997. 
 *
 *  @author Alfonso Veropalumbo, Federico Marulli 
 *
 *  @author alfonso.veropalumbo@unibo.it, federico.marulli3@unibo.it
 */

#ifndef __EH__
#define __EH__


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class EisensteinHu EisensteinHu.h
     *  "Headers/EisensteinHu.h"
     *
     *  @brief The class EisensteinHu
     *
     *  This class is used to handle objects of type <EM>
     *  EisensteinHu </EM>. It contains all the functions and
     *  fitting formulae for CDM + Baryon + Massive Neutrino (MDM)
     *  cosmologies
     */
    class EisensteinHu {

    protected:

      /// sqrt(alpha_nu) 
      double alpha_gamma;	

      /// the small-scale suppression 
      double alpha_nu;

      /// the correction to the log in the small-scale 
      double beta_c;		

      /// number of degenerate massive neutrino species
      double num_degen_hdm;	

      /// baryon fraction
      double f_baryon;

      /// baryon + Massive Neutrino fraction
      double f_bnu;		

      /// baryon + CDM fraction
      double f_cb;		

      /// CDM fraction 
      double f_cdm;		

      /// massive Neutrino fraction 
      double f_hdm;		

      /// D_1(z) -- the growth function as k->0 
      double growth_k0;	

      /// D_1(z)/D_1(0) -- the growth relative to z=0 
      double growth_to_z0;

      /// need to pass Hubble constant to TFmdm_onek_hmpc() 
      double hhubble;	

      /// the comoving wave number of the horizon at equality
      double k_equality;	

      /// Omega_baryon * hubble^2 
      double obhh;		

      /// = 1 - omega_matter - omega_lambda 
      double omega_curv;

      /// Omega_lambda at the given redshift 
      double omega_lambda_z;     

      /// Omega_matter at the given redshift 
      double omega_matter_z;	

      /// Omega_matter * hubble^2 
      double omhh;		

      /// Omega_hdm * hubble^2 
      double onhh;		

      /// the correction to the exponent before drag epoch 
      double p_c;		

      /// the correction to the exponent after drag epoch 
      double p_cb;		

      /// the sound horizon at the drag epoch 
      double sound_horizon_fit;  

      /// the temperature of the CMB; in units of 2.7 K 
      double theta_cmb;	

      /// ratio of z_equality to z_drag 
      double y_drag;	

      /// redshift of the drag epoch 
      double z_drag;		

      /// redshift of matter-radiation equality 
      double z_equality;	

      /// effective \f$gamma\f$ 
      double gamma_eff;	
	  
      /// growth factor for CDM+Baryon perturbations 
      double growth_cb;
	  
      /// growth factor for CDM+Baryon+Neutrino pert.
      double growth_cbnu;	
	  
      /// correction near maximal free streaming 
      double max_fs_correction; 
	  
      /// wavenumber rescaled by \f$\gamma\f$ 
      double qq;		
	  
      /// wavenumber rescaled by effective Gamma 
      double qq_eff;		
	  
      /// wavenumber compared to maximal free streaming 
      double qq_nu;		
	  
      /// master TF 
      double tf_master;	
	  
      /// suppressed TF 
      double tf_sup;		

      /// the epoch of free-streaming for a given scale 
      double y_freestream; 	

      /// the transfer function for density-weighted CDM + Baryon perturbations.
      double tf_cb;

      /// the transfer function for density-weighted CDM + Baryon + Massive Neutrino perturbations.
      double tf_cbnu;

      /// the power spectrum normalization
      double pk_normalization;

      /// the spectral index
      double ns;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      EisensteinHu () = default;

      /**
       *  @brief default destructor
       */
      ~EisensteinHu () = default;

      ///@}

	  
      /**
       * @brief set the cosmological parameters 
       *
       *  This routine takes cosmological parameters and a redshift
       *  and sets up all the internal scalar quantities needed to
       *  compute the transfer function.
       *
       *  @param omega_matter density of CDM, baryons, and massive
       *  neutrinos, in units of the critical density
       *
       *  @param omega_baryon density of baryons, in units of
       *  critical
       *
       *  @param omega_hdm density of massive neutrinos, in units of
       *  critical
       *
       *  @param degen_hdm number of degenerate massive neutrino
       *  species
       *
       *  @param omega_lambda cosmological constant 
       *
       *  @param hubble Hubble constant, in units of 100 km/s/Mpc 
       *
       *  @param redshift the redshift at which to evaluate
       *
       *  @param As the amplitude of the curvature perturbations
       *  
       *  @param k_pivot the scalar pivot k in \f$Mpc^{-1}\f$
       *
       *  @param n_spec the primordial spectral index
       *
       *  @return 0 if all is well, 1 if a warning was issued   
       */
      int TFmdm_set_cosm (double omega_matter, double omega_baryon, double omega_hdm, int degen_hdm, double omega_lambda, double hubble, double redshift, double As=2.56e-9, double k_pivot=0.05, double n_spec=0.96)
      {
	double z_drag_b1, z_drag_b2, omega_denom;
	int qwarn;
	qwarn = 0;

	theta_cmb = 2.728/2.7; // assuming T_cmb = 2.728 K 

	/* Look for strange input */

	if (omega_baryon<0.0) {
	  fprintf(stderr,
		  "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
	  qwarn = 1;
	}
	    
	if (omega_hdm<0.0) {
	  fprintf(stderr,
		  "TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
	  qwarn = 1;
	}

	if (hubble<=0.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
	  exit(1);  /* Can't recover */
	} else if (hubble>2.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
	  qwarn = 1;
	}

	if (redshift<=-1.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
	  exit(1);
	} else if (redshift>99.0) {
	  fprintf(stderr,
		  "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
	  qwarn = 1;
	}

	if (degen_hdm<1) degen_hdm = 1;
	num_degen_hdm = (float) degen_hdm;	
	/* Have to save this for TFmdm_onek_mpc() */
	/* This routine would crash if baryons or neutrinos were zero, 
	   so don't allow that */
	if (omega_baryon<=0) omega_baryon=1e-5;
	if (omega_hdm<=0) omega_hdm=1e-5;

	omega_curv = 1.0-omega_matter-omega_lambda;
	omhh = omega_matter*pow(hubble, 2);
	obhh = omega_baryon*pow(hubble, 2);
	onhh = omega_hdm*pow(hubble, 2);
	f_baryon = omega_baryon/omega_matter;
	f_hdm = omega_hdm/omega_matter;
	f_cdm = 1.0-f_baryon-f_hdm;
	f_cb = f_cdm+f_baryon;
	f_bnu = f_baryon+f_hdm;

	/* Compute the equality scale. */
	z_equality = 25000.0*omhh*pow(theta_cmb, -4);	/* Actually 1+z_eq */
	k_equality = 0.0746*omhh*pow(theta_cmb, -2);

	/* Compute the drag epoch and sound horizon */
	z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
	z_drag_b2 = 0.238*pow(omhh,0.223);
	z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2));
	y_drag = z_equality/(1.0+z_drag);

	sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

	/* Set up for the free-streaming & infall growth function */
	p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
	p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

	omega_denom = omega_lambda+pow(1.0+redshift, 2)*(omega_curv+omega_matter*(1.0+redshift));
	omega_lambda_z = omega_lambda/omega_denom;
	omega_matter_z = omega_matter*pow(1.0+redshift,2)*(1.0+redshift)/omega_denom;

	growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)-omega_lambda_z+(1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
	growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)-omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
	double growth_z = 1./(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)-omega_lambda_z + (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
	growth_to_z0 = growth_k0/growth_to_z0;	

	/*
	  growth_k0 = z_equality*5.*omega_matter_z/(2*(1+redshift))/(1./70.+209./140.*omega_matter_z-pow(omega_matter_z,2)/140.+pow(omega_matter_z,4./7.)); 
	  growth_to_z0 = z_equality*5.*omega_matter_z/(2*(1+redshift))/(1./70.+209./140.*omega_matter_z-pow(omega_matter_z,2)/140.+pow(omega_matter_z,4./7.)); 
	  growth_to_z0 = growth_k0/growth_to_z0;	
	*/

	/* Compute small-scale suppression */
	alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
	  pow(1+y_drag,p_cb-p_c)*
	  (1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
	  (1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*
	  (1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
	alpha_gamma = sqrt(alpha_nu);
	beta_c = 1/(1-0.949*f_bnu);

	/* Done setting scalar variables */
	hhubble = hubble;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */

	ns = n_spec;

	//double norm = 2.*par::pi*par::pi*m_scalar_amp*pow(2.*pow(par::cc/100.,2)/(5*m_Omega_matter)*DN(redshift),2)*pow(k/m_scalar_pivot, m_n_spec-1)*k_ov_h*pow(fact, -3);
	pk_normalization = 2.*par::pi*par::pi*As*pow(2.*pow(par::cc/100.,2)/(5*omega_matter)*growth_z,2)*pow(1./k_pivot, ns-1);

	return qwarn;
      }

      /** 
       * @brief compute the transfer function
       * \f$ T(k) \f$
       *
       * Given  a  wavenumber  in   Mpc^-1,  return  the  transfer
       * function   for  the   cosmology   held   in  the   global
       * variables. 
       * @param kk Wavenumber in Mpc^-1
       *
       * @return the transfer function for density-weighted
       * CDM+Baryon perturbations
       */
      double TFmdm_onek_mpc (double kk)
      {
	double tf_sup_L, tf_sup_C;
	double temp1, temp2;

	qq = kk/omhh*pow(theta_cmb, 2);

	/* Compute the scale-dependent growth functions */
	y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*pow(num_degen_hdm*qq/f_hdm, 2);
	temp1 = pow(growth_k0, 1.0-p_cb);
	temp2 = pow(growth_k0/(1+y_freestream),0.7);
	growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
	growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

	/* Compute the master function */
	gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/(1+pow(kk*sound_horizon_fit*0.43, 4)));
	qq_eff = qq*omhh/gamma_eff;

	tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
	tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
	tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*pow(qq_eff, 2));

	qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
	max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
	tf_master = tf_sup*max_fs_correction;

	/* Now compute the CDM+HDM+baryon transfer functions */
	tf_cb = tf_master*growth_cb/growth_k0;
	tf_cbnu = tf_master*growth_cbnu/growth_k0;
	  
	return tf_cb;
      }

      /** 
       * @brief compute the transfer function \f$ T(k) \f$
       *
       * Given  a  wavenumber  in h  Mpc^-1,  return  the  transfer
       * function   for  the   cosmology   held   in  the   global
       * variables. 
       * @param kk Wavenumber in h Mpc^-1
       *
       * @return the transfer function for density-weighted
       * CDM+Baryon perturbations
       */
      double TFmdm_onek_hmpc (double kk)
      {
	return TFmdm_onek_mpc(kk*hhubble);
      }

      /**
       * @brief get the power spectrum in unit of h^{-3}Mpc^{-3}
       *
       * @param kk wavenumber in h Mpc^-1 
       *
       * @return the power spectrum in unit of h^{-3}Mpc^{-3}
       */
      double Pk (double kk)
      {
	double k = kk*hhubble;

	double func = TFmdm_onek_mpc(k); // transfer function

	return pk_normalization*pow(func,2)*kk*pow(k, ns-1);
      }

      /**
       * @brief get the power spectrum in unit of h^{-3}Mpc^{-3}
       *
       * @param kk vector of wavenumbers in h Mpc^-1 
       *
       * @return vector with the power spectrum in unit of
       * h^{-3}Mpc^{-3}
       */
      std::vector<double> Pk(std::vector<double> kk) {
	    
	std::vector<double> _Pk(kk.size(), 0);
	    
	for (size_t i=0; i<kk.size(); i++)
	  _Pk[i] = Pk(kk[i]);
	    
	return _Pk;
      }
	  
    };
    
  }
}

#endif
