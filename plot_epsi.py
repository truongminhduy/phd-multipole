import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as plt
import pandas as pd
import sys
import seaborn as sns
# print(sys.version) 
from function_e import *
import matplotlib.gridspec as gridspec

import os
import sys
import matplotlib.pyplot as plt
plt.rc('font', size=13,)
pi = np.pi
plt.style.use('seaborn-muted')
plt.rc('text', usetex=True)
plt.rc('font', **{'family' : "serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)
#plt.rcParams['text.latex.preamble'] = [r'\boldmath']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.titlesize'] = 18

pi    = np.pi
micro = 1e-6
femto = 1e-15
peta  = 1e15
cc    = 299792458
mu_0  = pi*4e-7
eps_0 = 1./(mu_0*cc**2)
## Incident plane wave parameters
scale_light = 1e8
micro       = 1e-6
scale_frequ = scale_light/micro # 1e14
cel = cc/scale_light
mu0 = mu_0  *scale_light
ep0 = eps_0 *scale_light
#------------------------------------------------------------------------------------------#
filetxt = "mau/Si_Green_data.txt"
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
Lam_dat, n_dat, k_dat = np.loadtxt(filetxt, usecols=(0,1,2), skiprows = 1, unpack = True )
w_BG    = 2*pi*cc/(Lam_dat*micro)/peta
Chi_dat = (n_dat + 1j*k_dat)**2-np.ones_like(w_BG)
#------------------------------------------------------------------------------------------#

N_num = 10
N_den =	10
w = np.linspace(-w_BG.max()*1.5, w_BG.max()*1.5,100)
p_norm, q_norm, Mat_pol_amp = Fitting_PQ(w, w_BG, Chi_dat, N_num, N_den)

##----change the unit from Mau to our problem-----##
Mat_pol_amp = eliminate_negative(Mat_pol_amp)
Mat_eps = Mat_pol_amp * (peta/scale_frequ)

# omegas = 2*pi*cel/Lam_dat
omegas = np.linspace(17,80,500)

N_tru = 1
num,den,num2 = getEpsrNumDen(Mat_eps,N_tru) 
eps1 = getrational(omegas,num,den)

N_tru = 3
num,den,num2 = getEpsrNumDen(Mat_eps,N_tru) 
eps3 = getrational(omegas,num,den)

N_tru = 4
num,den,num2 = getEpsrNumDen(Mat_eps,N_tru) 
eps4 = getrational(omegas,num,den)

# print (pole_root.imag)
# print (-pole_root.real)
# #------------------------------------------------------------------------------------------#
plt.figure()
gs = gridspec.GridSpec(2, 1)
ax1 = plt.subplot(gs[0,0])
ax1.plot(omegas, eps1.real, 'C1:',label='1 pole',lw=1.8)
ax1.plot(omegas, eps3.real, 'C3-.',label='3 poles',lw=1.8)
ax1.plot(omegas, eps4.real, 'C2-',label='4 poles',lw=1.8)
ax1.plot(w_BG*10,1+ Chi_dat.real, 'C0o',label='data', ms=4) 
plt.legend()
# plt.xlabel(r"$\omega$")
plt.ylabel(r"$\Re(\varepsilon)$")

ax2 = plt.subplot(gs[1,0])
ax2.plot(omegas, eps1.imag, 'C1:',label='1 pole', linewidth=1.8) 
ax2.plot(omegas, eps3.imag, 'C3-.',label='3 poles',lw=1.8) 
ax2.plot(omegas, eps4.imag, 'C2-',label='4 poles',lw=1.8) 
ax2.plot(w_BG*10, Chi_dat.imag, 'C0o',label='data', ms=4)
plt.legend() 
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\Im(\varepsilon)$")

## plt.show()
#------------------------------------------------##

plt.show()