import os
import sys
include "exCosmology.pyd"
from classy import Class
import camb4py

cosmopar ={ "Omega_matter" : 0.27,
    "Omega_baryon" : 0.046,
    "Omega_neutrinos" : 0.,
    "massless_neutrinos" : 3.04,
    "massive_neutrinos" : 0,
    "Omega_DE" : 0.73,
    "Omega_radiation": 0.,
    "hh" : 0.7,
    "scalar_amp" : 2.46e-9,
    "n_spec" : 0.96,
    "w0" :-1.,
    "wa" : 0.,
    "f_NL": 0.,
    "type_NG": 1,
    "model" : "LCDM",
    "unit" : 1
    }


# ============================================================================================================


#cdef shared_ptr[Cosmology] get_cosmology():

cdef class pyCosmologyCBL:

#  cdef Cosmology *thisptr
  cdef shared_ptr[Cosmology] thisptr
  
  def __cinit__(self,dict):

    for kk in dict.keys():
      if kk in cosmopar.keys():
        cosmopar[kk] = dict[kk]

    parameters = cosmopar
#    self.thisptr = new Cosmology(cosmopar["Omega_matter"],cosmopar["Omega_baryon"],cosmopar["Omega_neutrinos"],cosmopar["massless_neutrinos"],cosmopar["massive_neutrinos"],cosmopar["Omega_DE"],cosmopar["Omega_radiation"],cosmopar["hh"],cosmopar["scalar_amp"], cosmopar["n_spec"],cosmopar["w0"],cosmopar["wa"],cosmopar["f_NL"],cosmopar["type_NG"],cosmopar["model"],cosmopar["unit"])
    self.thisptr = shared_ptr[Cosmology](new Cosmology(cosmopar["Omega_matter"],cosmopar["Omega_baryon"],cosmopar["Omega_neutrinos"],cosmopar["massless_neutrinos"],cosmopar["massive_neutrinos"],cosmopar["Omega_DE"],cosmopar["Omega_radiation"],cosmopar["hh"],cosmopar["scalar_amp"], cosmopar["n_spec"],cosmopar["w0"],cosmopar["wa"],cosmopar["f_NL"],cosmopar["type_NG"],cosmopar["model"],cosmopar["unit"]))

#  def __dealloc__(self):
#    del self.thisptr
#  cdef shared_ptr[Cosmology] get_cosmology(self):
#    return shared_ptr[Cosmology](self.thisptr)
  cdef shared_ptr[Cosmology] get_cosmology(self):
    return self.thisptr
  def Omega_matter(self):
    return deref(self.thisptr).Omega_matter()
  def Omega_baryon(self):
    return deref(self.thisptr).Omega_baryon()
  def Omega_neutrinos(self):
    return deref(self.thisptr).Omega_neutrinos()
  def massless_neutrinos(self):
    return deref(self.thisptr).massless_neutrinos()
  def massive_neutrinos(self):
    return deref(self.thisptr).massive_neutrinos()
  def Omega_DE(self):
    return deref(self.thisptr).Omega_DE()
  def Omega_radiation(self):
    return deref(self.thisptr).Omega_radiation()
  def Omega_k(self):
    return deref(self.thisptr).Omega_k()
  def Omega_CDM(self):
    return deref(self.thisptr).Omega_CDM()
  def H0(self):
    return deref(self.thisptr).H0()
  def hh(self):
    return deref(self.thisptr).hh()
  def t_H(self):
    return deref(self.thisptr).t_H()
  def D_H(self):
    return deref(self.thisptr).D_H()
  def sigma8(self):
    return deref(self.thisptr).sigma8()
  def scalar_amp(self):
    return deref(self.thisptr).scalar_amp()
  def n_spec(self):
    return deref(self.thisptr).n_spec()
  def w0(self):
    return deref(self.thisptr).w0()
  def wa(self):
    return deref(self.thisptr).wa()
  def RhoZero(self):
    return deref(self.thisptr).RhoZero()
  def fNL(self):
    return deref(self.thisptr).fNL()
  def type_NG(self):
    return deref(self.thisptr).type_NG()
  def Pk0_EH(self):
    return deref(self.thisptr).Pk0_EH()
  def Pk0_CAMB(self):
    return deref(self.thisptr).Pk0_CAMB()
  def PK0_MPTbreeze(self):
    return deref(self.thisptr).Pk0_MPTbreeze()
  def Pk0_CLASS(self):
    return deref(self.thisptr).Pk0_CLASS()
  def Pk_0 (self,author,redshift,Model="LCDM", k_min=0, k_max=100, GSL=1, prec=1.e-2, file_par="NULL"): 
    return deref(self.thisptr).Pk_0(author,redshift,Model,k_min,k_max,GSL,prec,file_par)

  def model(self):
    return deref(self.thisptr).model()
  def unit(self):
    return deref(self.thisptr).unit()
  def print_parameters(self):
    deref(self.thisptr).print_parameters()

  def set_Omega(self,Omega_matter):
    deref(self.thisptr).set_Omega(Omega_matter)
  def set_OmegaB(self,Omega_baryon):
    deref(self.thisptr).set_OmegaB(Omega_baryon)
  def set_OmegaM(self,Omega_matter):
    deref(self.thisptr).set_OmegaM(Omega_matter)
  def set_OmegaDE(self,Omega_DE):
    deref(self.thisptr).set_OmegaDE(Omega_DE) 
  def set_OmegaNu(self,Omega_neutrinos=0., massless_neutrinos=3.04, massive_neutrinos = 0.):
    deref(self.thisptr).set_OmegaNu(Omega_neutrinos,massless_neutrinos,massive_neutrinos)
  def set_H0(self,H0):
    deref(self.thisptr).set_H0(H0)
  def set_sigma8(self,sigma8):
    deref(self.thisptr).set_sigma8(sigma8)
  def set_w0(self,w0):
    deref(self.thisptr).set_w0(w0)
  def set_wa(self,wa):
    deref(self.thisptr).set_wa(wa)

  def OmegaM(self,redshift):
    return deref(self.thisptr).OmegaM(redshift)
  def OmegaDE(self,redshift):
    return deref(self.thisptr).OmegaDE(redshift)
  def OmegaR(self,redshift):
    return deref(self.thisptr).OmegaR(redshift)
  def OmegaK(self,redshift):
    return deref(self.thisptr).OmegaK(redshift)
  def Omega(self,redshift):
    return deref(self.thisptr).Omega(redshift)

  def Rho (self,Omega_matter,Omega_neutrinos,unit=1):
    deref(self.thisptr).Rho(Omega_matter,Omega_neutrinos,unit)
  def DeltaR(self,Delta_critical,redshift):
    deref(self.thisptr).DeltaR(Delta_critical,redshift)

  def w_CPL(self,redshift):
    return deref(self.thisptr).w_CPL(redshift)
  def f_DE(self,redshift):
    return deref(self.thisptr).f_DE(redshift)

  def EE(self,redshift):
    return deref(self.thisptr).EE(redshift)
  def HH(self,redshift):
    return deref(self.thisptr).HH(redshift)

  def gg(self,redshift):
    return deref(self.thisptr).gg(redshift)
  def DD(self,redshift):
    return deref(self.thisptr).DD(redshift) 
  def lookback_time(self,redshift):
    return deref(self.thisptr).lookback_time(redshift) 
  def cosmic_time(self,redshift):
    return deref(self.thisptr).cosmic_time(redshift) 
  def EE2(self,redshift):
    return deref(self.thisptr).EE2(redshift) 
  def qq(self,redshift):
    return deref(self.thisptr).qq(redshift) 
  def Hdot(self,redshift):
    return deref(self.thisptr).Hdot(redshift) 

  def z_acc(self):
    return deref(self.thisptr).z_acc()
  def z_eq(self):
    return deref(self.thisptr).z_eq()

  def Mag_Volume_limited(self,redshift_max,magnitude_limit):
    return deref(self.thisptr).Mag_Volume_limited(redshift_max,magnitude_limit)
  def Lum_bol(self,redshift,flux):
    return deref(self.thisptr).Lum_bol(redshift,flux)
 
  def Redshift(self,comoving_distance,redshift_guess=(0.,1.),method="generic",go_fast=True,prec=0.0001):
    method_list = ("LCDM","generic")
    if method not in method_list: 
      print "Error in Redshift, provided method does not exist, choose one of LCDM,generic"
      sys.exit(1)
    elif method=="LCDM":
      return deref(self.thisptr).Redshift_LCDM(comoving_distance,redshift_guess[0],redshift_guess[1],go_fast,prec)
    elif method == "generic":
      return deref(self.thisptr).Redshift(comoving_distance,redshift_guess[0],redshift_guess[1],prec)
   
  def Redshift_time(self,time,redshift_guess=(0.,1.)):
     return deref(self.thisptr).Redshift_time(time,redshift_guess[0],redshift_guess[1])

  def deltac(self,redshift):
     return deref(self.thisptr).deltac(redshift)
#  def Deltavir(self,redshift):
#     return deref(self.thisptr).Deltavir(redshift)

  def D_C(self,redshift,method="generic"):
    method_list = ("LCDM","generic")
    if method not in method_list: 
      print "Error in D_C, provided method does not exist, choose one of LCDM,generic"
      sys.exit(1)
    elif method=="generic":
      return deref(self.thisptr).D_C(redshift)
    elif method=="LCDM":
      return deref(self.thisptr).D_C_LCDM(redshift)

  def D_M(self,redshift):
    return deref(self.thisptr).D_M(redshift)
  def D_A(self,redshift):
    return deref(self.thisptr).D_A(redshift)
  def D_L(self,redshift):
    return deref(self.thisptr).D_L(redshift)
  def D_V(self,redshift):
    return deref(self.thisptr).D_V(redshift)  
  
  def Volume_z1z2Area(self,redshift_min,redshift_max,area):
    return deref(self.thisptr).Volume(redshift_min,redshift_max,area)
  def Volume(self,redshift):
    return deref(self.thisptr).Volume(redshift)
  def max_redshift(self,volume,area,redshift_min):
    return deref(self.thisptr).max_redshift(volume,area,redshift_min) 
  def dV_dZdOmega (self,redshift,AngleUnit=1):
    return self.thiptr.dV_dZdOmega(redshift,AngleUnit)  
  def Pk(self,kk,author,NL,redshift,Model="LCDM",norm=-1,k_min=0,k_max=100,GSL=1,prec=1.e-2,file_par="NULL"):
    return deref(self.thisptr).Pk(kk,author,NL,redshift,Model,norm,k_min,k_max,GSL,prec,file_par)

  def xi_DM (self, r,author,redshift,Model="LCDM" , NL=1, norm=-1, k_min=0., k_max=100., aa=0., GSL=1, prec=1.e-2, file_par="NULL"):
    return deref(self.thisptr).xi_DM(r,author,redshift,Model,NL,norm,k_min,k_max,aa,GSL,prec,file_par)
    #if (author == "EisensteinHu" or author == "MPTbreeze-v1"):
    #  return deref(self.thisptr).xi_DM(r,author,redshift,Model,NL,norm,k_min,k_max,aa,GSL,prec,file_par)
    #if (author == "CAMB"):

  def xi_DM_DeWiggle(self, r,redshift,sigma_NL,Model="LCDM" , norm=-1, k_min=0., k_max=100., aa=0., prec=1.e-2):
    xi = deref(self.thisptr).xi_DM_DeWiggle(r,redshift,sigma_NL,Model,norm,k_min,k_max,aa,prec)
    return xi

  def rs_CLASS(self):
    params = { 
	"A_s" : deref(self.thisptr).scalar_amp(),
	"n_s" : deref(self.thisptr).n_spec(),
	"h" : deref(self.thisptr).hh(),
	"omega_b": deref(self.thisptr).Omega_baryon()*deref(self.thisptr).hh()*deref(self.thisptr).hh(),
	"omega_cdm": deref(self.thisptr).Omega_CDM()*deref(self.thisptr).hh()*deref(self.thisptr).hh(),
	}
   # print params
    cosmoclass = Class()
    cosmoclass.set(params)
    rs= cosmoclass.rs_drag()
    cosmoclass.struct_cleanup()
    cosmoclass.empty()
    return rs

  def rs_EH(self,T_CMB=2.7255):
    return deref(self.thisptr).rs_EH(T_CMB)
  def rs(self,method,T_CMB=2.7255):
    if method=="EH":
      return deref(self.thisptr).rs_EH(T_CMB)
    elif method == "CLASS":
      return self.rs_CLASS()
    elif method == "CAMB":
      return deref(self.thisptr).rs_CAMB()
    else:
      print "Unknown method"
      sys.exit(1)
