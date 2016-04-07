"""
 Assume that the CLASS file given as input is:
 
 l corr[1]-[1] corr[1]-[2] corr[1]-[...] corr[2]-[2] corr[2]-[3] corr...
 
 where corr[i][j] are C_l's (not rescaled by [l(l+1)/2pi]) for l>=0

 
 From these Cl's, the following is calculated in sequence:
 - Cl(z'-z) at a fixed multipole around z_mean
 - xi(z'-z) at a fixed theta around z_mean
 - xi(theta) at a given redshift, and its variance
 - Reduced covariance Matrix for xi(theta) computed above (OPTIONAL)
 
 For further information, see Di Dio et al. (2013) [arXiv:XXXX.YYYY]
""" 
 

import numpy as np
import math
import scipy.special as sp
import matplotlib.pyplot as plt
from pylab import *
import multiprocessing as multi


"""
----------------------------------------------------------------
----> Input: MODIFY THIS SECTION ACCORDING TO THE DESIRED OUTPUT
----------------------------------------------------------------
"""

"Reference redshift at which compute the output"
selection_mean= 1.
"Number of bins calculated in CLASS, assumed to be odd"
Nbins= 51

"Cosmology"
h= 0.67
Omega_b= 0.022/h**2.
Omega_cdm= 0.1199/h**2.
Omega_m= Omega_b+Omega_cdm
Omega_Lambda= 1.-Omega_m

"CLASS input Cl's file"
InputFile="./output/test_cl.dat"

"** Cl(z-z') *"
"multipole at which evaluate Cl(Dz)"
multipole=0

"** xi(z-z') *"
"Dz sampling with which Cl's have been computed"
Dz=0.002
"angle, in degrees, at which compute xi(Dz)"
theta=0.

"** xi(z-z') and xi(theta) *"
"width of the gaussian smoothing high ell's when converting Cl's --> xi(theta), to reduce oscillations in xi(theta)"
ls=600.

"** Reduced covariance Matrix *"
"if 'yes', evaluate reduced covariance matrix, see arXiv:1004.4640"
evaluate_cov = 'no'
"accessible fraction of the sky"
fsky=1.
"shot_noise=1/nbar, where nbar is the number of galaxies per steradiant"
shot_noise=3.48*10.**(-7.)

"""
-----------------------------------------------------------------------
----> END of Input section. Modify what follows only to change outputs.
-----------------------------------------------------------------------
"""




"""
--------------------------------------------------------------------------------
********************************************************************************
--------------------------------------------------------------------------------
"""
"""
-----------------------
----> Compute Cl(z-z'):
-----------------------
"""
print("Computing Cl(Dz)...")

Fin_Cl = np.loadtxt(InputFile,unpack=True)

N_central_with_class_array = 1
for corr_i in range( (Nbins + 1)/2 ):
  N_central_with_class_array += (Nbins - corr_i)

"Take first column, all lines"
ell = Fin_Cl[0,:]
Nell=len(ell)

"Fill in half of a symmetric matrix from the input array"
Cl_matrix = np.zeros(shape=(Nbins,Nbins,Nell))
for ell_index in range(Nell):
  counter=0
  for z_row in range(Nbins):
    for z_column in range(Nbins):
      if (z_row>z_column):
        Cl_matrix[z_column,z_row,ell_index] = 0
      else:
        Cl_matrix[z_column,z_row,ell_index] = Fin_Cl[counter+1,ell_index]
        counter+=1

"Define Cl array containing only useful bins [central][central], [central-1][central+1], [central-2][central+2], ..."
N_corr_Dz = (Nbins+1)/2
Cl = np.zeros(shape=(N_corr_Dz,Nell))
for i_Dz in range(N_corr_Dz):
  Cl[i_Dz,:] = Cl_matrix[N_corr_Dz+i_Dz-1,N_corr_Dz-i_Dz-1,:]

"Export plot"
Dz_axis=[2*i*Dz for i in range(N_corr_Dz)]
ell_plot = multipole-ell[0]
Cl_plot =Cl[:,ell_plot]

"Export data file"
Fout_clDz = open('output/Cl_Dz.dat','w')
Fout_clDz.write('# Dz   Cl\n')
for i in range(N_corr_Dz):
  Fout_clDz.write("%e\t %e\n"%(Dz_axis[i], Cl_plot[i]))
Fout_clDz.close()
print("--> 'output/Cl_Dz.dat' created")

figure()
plot(Dz_axis, Cl_plot, linewidth=1.0, marker='o', linestyle='-', color='b', label='positive')
plot(Dz_axis, -Cl_plot, linewidth=1.0, marker='o', linestyle='--', color='b', label='negative')
xlabel(r'$\Delta z$')
ylabel(r'$C_{\ell}$')
yscale('log')
title(r'$\ell=$%g, $z_s=$%g.' %(multipole,selection_mean))
grid(True)
legend()
savefig("output/Cl_Dz.pdf")
print("--> 'output/Cl_Dz.pdf' created")
#show()


"""
-------------------------
----> Calculate xi(z-z'):
-------------------------
"""
print("Computing xi(Dz)...")

cosTh=math.cos(theta*math.pi/180.)
xi=[0]*N_corr_Dz
for row in range(len(ell)):
  xi += (2.*ell[row]+1.)*Cl[:,row]*sp.eval_legendre(ell[row],cosTh) * math.exp(-ell[row]*(ell[row]+1.)/ls**2.)
xi = xi/(4.*math.pi)


Dz_axis=[2.*i*Dz for i in range(N_corr_Dz)]

"#Export data file"
Fout_xiDz = open('output/xi_Dz.dat','w')
Fout_xiDz.write('# Dz   xi\n')
for i in range(N_corr_Dz):
  Fout_xiDz.write("%e\t %e\n"%(Dz_axis[i], xi[i]))
Fout_xiDz.close()
print("--> 'output/xi_Dz.dat' created")

"Export plot"
figure()
plot(Dz_axis, xi, linewidth=1.0, marker='o', linestyle='-', color='b', label='positive')
plot(Dz_axis, -xi, linewidth=1.0, marker='o', linestyle='--', color='b', label='negative')
xlabel(r'$\Delta z$')
ylabel(r'$\xi(\Delta z)$')
yscale('log')
title(r'Correlation at $\theta=$%g deg, around $z=$%g. Gaussian cutoff $\ell_s=$%g.' %(theta*180./math.pi,selection_mean,ls))
grid(True)
legend()
savefig("output/xi_Dz.pdf")
print("--> 'output/xi_Dz.pdf' created")
#show()


"""
-------------------------
----> Calculate xi(theta):
-------------------------
"""
print("Computing xi(theta)...")

Ntheta=100 # Must be larger than theta_max
theta_max=10. #in degrees, better en integer number (for plots)

Dz_bin=0. # bin number at which to compute xi(theta); 0 corresponds to corr[central]-[central]

xi=[0]*Ntheta
var=[0]*Ntheta

for th_line in range(Ntheta):
  cosTh=math.cos(th_line*theta_max/Ntheta*math.pi/180.)
  for row in range(Nell):
    xi[th_line] += (2.*ell[row]+1.)*Cl[Dz_bin,row] /(4.*math.pi) * math.exp(-ell[row]*(ell[row]+1.)/ls**2.) * sp.eval_legendre(ell[row],cosTh)
    var[th_line] += 2./fsky* (2.*ell[row]+1.) /(4.*math.pi)**2. *(Cl[Dz_bin,row]+shot_noise)**2.*sp.eval_legendre(ell[row],cosTh)**2.

var=sqrt(var)

theta_axis=[i*theta_max/Ntheta for i in range(Ntheta)]
#Export data file
Fout_xiTh = open('output/xi_theta.dat','w')
Fout_xiTh.write('# Correlation around z=%g.\n # Gaussian cutoff l_s=%g\n' %(selection_mean,ls))
Fout_xiTh.write('# Columns: theta [deg]   xi   variance\n')
for th_line in range(Ntheta):
  Fout_xiTh.write("%e\t %e\t %e\n" %(theta_axis[th_line],xi[th_line],var[th_line]))
Fout_xiTh.close()
print("--> 'output/xi_theta.dat' created")


figure()
plot(theta_axis, xi, linewidth=1.0, marker='o', linestyle='-', color='b', label='xi, f_{sky}=%g'%fsky)

#plot the negative value
for th_line in range(Ntheta):
  xi[th_line]*=-1.  
plot(theta_axis, xi, linewidth=1.0, marker='o', linestyle='--', color='b', label='-xi')
for th_line in range(Ntheta):
  xi[th_line]*=-1.  #Multiply back by -1, in case later use is needed
  
plot(theta_axis, var, linewidth=1.0, marker='x', linestyle='-', color='black', label='variance')
xlabel(r'$\theta$ [deg]')
ylabel(r'$\xi(\theta)$')
yscale('log')
title(r'Correlation at $z=$%g. Gaussian cutoff $\ell_s=$%g' %(selection_mean,ls))
grid(True)
legend()
savefig("output/xi_theta.pdf")
print("--> 'output/xi_theta.pdf' created")
#show()


"""
--------------------------------
----> Reduced covariance Matrix:
--------------------------------
"""
if (evaluate_cov != 'yes'):
  exit()
print("Computing reduced covariance matrix...")

Npixel=50

cov=np.zeros(shape=(Npixel,Npixel))
cov_line=multi.Manager().list([0]*Npixel)

def ComputeCov(th_line):
  for th_column in range(Npixel):
    if th_column>=th_line:
      for row in range(Nell):
        cosTh=math.cos(th_line*theta_max/Npixel*math.pi/180.)
        cosTh_prime=math.cos(th_column*theta_max/Npixel*math.pi/180.)
        """
         The following matrix is not handled by multi.Manager(), hence it will
         be forgotten outside this ComputeCov() function when running in parallel
        """
        cov[th_column,th_line] += 2./fsky* (2.*ell[row]+1.) /(4.*math.pi)**2. *(Cl[Dz_bin,row]+shot_noise)**2.*sp.eval_legendre(ell[row],cosTh)*sp.eval_legendre(ell[row],cosTh_prime)
    """
     Fill the th_line-element with all the columns in the line number=th_line.
     Hence cov_line is a vector where every element is a vector containing a
     line of cov
    """
    cov_line[th_line]=cov[:,th_line]


"--- Parallel computation"
multi.Pool(processes=multi.cpu_count()).map(ComputeCov,range(Npixel))

"Fill the covariance matrix line by line"
for th_line in range(Npixel):
  cov[:,th_line]=cov_line[th_line]
"Fill the symmetric part of the covariance matrix"
for th_line in range(Npixel):
  for th_column in range(Npixel):
    if th_column<th_line:
      cov[th_column,th_line]=cov[th_line,th_column]


"Normalization to obtain the reduced cov"
for th_line in range(Npixel):
  for th_column in range(Npixel):
    if th_column != th_line:
      cov[th_column,th_line] /= math.sqrt(cov[th_line,th_line]*cov[th_column,th_column])
for th_line in range(Npixel):
  cov[th_line,th_line] = 1.


"Export data file"
Fout_covXi = open('output/cov_ThThp.dat','w')

cov_line=[0]*Npixel
for th_line in range(Npixel):
  cov_line[th_line]=[cov[th_column,th_line] for th_column in range(Npixel)]
  Fout_covXi.write("\t".join(map(str, cov_line[th_line]))+"\n")
  
Fout_covXi.close()
print("--> 'output/cov_ThThp.dat' created")

from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

figure()
pcolor(cov)
cb=colorbar()
for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)
xlim([0,Npixel])
ylim([0,Npixel])
ticks=[i*Npixel/theta_max for i in range(int(theta_max)+1)]
labels=[i for i in range(int(theta_max)+1)]
xticks(ticks,labels,fontsize=20)
yticks(ticks,labels,fontsize=20)
xlabel(r"$\theta \left[{\rm deg}\right]$",fontsize=20)
ylabel(r"$\theta' \left[{\rm deg}\right]$",fontsize=20)

savefig("output/cov_ThThp.pdf")
print("--> 'output/xi_ThThp.pdf' created")
#show()
