################################################################################
#### PLOT COMPARE ################################################
####################################################################!!!#########
import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as plt
import pandas as pd
import sys
import seaborn as sns
from function_e import *
# print(sys.version) 
np.set_printoptions(precision=4)
import os
import sys
import matplotlib.pyplot as plt
plt.rc('font', size=13)
pi = np.pi
import matplotlib.gridspec as gridspec

plt.rcParams['font.family'] = 'serif'
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

# #### rational function (N/D)' ###############################################
# res = np.load('output.npz')['saveddic'].item()
res = np.load('pole4_PML_TE/output.npz')['saveddic'].item()
# res = np.load('pole4_closed_TM/output.npz')['saveddic'].item()
# res = np.load('pole4_closed_TE/output.npz')['saveddic'].item()
# res = np.load('pole1_closed_TE/output.npz')['saveddic'].item()


E2   = res['E2']
E0_0 = res['E0_0']
E0_1 = res['E0_1']
E0_2 = res['E0_2']
omegas = res['omegas']
num = res['num']
den = res['den']

pole_root = res['pole_root']
# epsi_root = res['epsi_root']

eigs = res['eigs']

# ##### PLOT 1 #################################################################################################################
plt.figure()
# gs = gridspec.GridSpec(5, 1)
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 5]) 
##### subplot 1 #####################################
ax1 = plt.subplot(gs[1])
# ax1 = plt.subplot(2, 1, 1)
ax1.plot(eigs.real, eigs.imag, 'C0o',markersize=3 ,mfc='none',label='Eigenmode')
# for k in range(len(eigs)):
#     plt.annotate(k, xy=(eigs.real[k],eigs.imag[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
omega_s = np.linspace(0,100,500)
aaa = -1*np.angle(3.5+6.5*(np.cos(pi/10)+1j*np.sin(pi/10)))
aaa = -1*np.angle(2.8+5.2*np.exp(1j*pi/10))
print (aaa)
omega_2 = omega_s*aaa
plt.plot(omega_s,omega_2 ,'C3-')
plt.ylabel(r"$\Im(\omega)$")
# plt.xlabel(r"$\omega.real$")
plt.plot(pole_root.imag,-pole_root.real,'C2+',label='Pole',  mew=3, ms=7)
# plt.plot(epsi_root.imag,-epsi_root.real,'C3+',label=r"$\omega_2$",  mew=3, ms=7)
plt.xlabel(r"$\Re(\omega)$")
plt.legend()
ax1.xaxis.grid()
plt.xlim(0, 100) 
plt.ylim(-23, 1) 
plt.legend(loc=3)
# plt.legend()
# plt.xlim(-20, 100)  


# ax2 = plt.subplot(gs[:2,0] ,sharex=ax1)
ax2 = plt.subplot(gs[0])
# ax2 = plt.subplot(2, 1, 2)
# omega_s = np.linspace(-max(eigs.real),max(eigs.real),500)
omega_s = np.linspace(0,100,500)
# omega_s = np.linspace(-6,15,500)
epsils = getrational(omega_s,num,den)
color = 'C0'
ax2.plot(omega_s, epsils.real, '-',markersize=3 ,color=color,label=r'$\Re(\varepsilon)$')
ax2.set_ylabel(r'$\Re(\varepsilon)$', color=color, size = 18) 
ax2.tick_params(axis='y', labelcolor=color)
ax2.xaxis.grid()
# plt.legend(loc=2)
ax3 = ax2.twinx()

color = 'C2'
ax3.plot(omega_s, epsils.imag, '-.',markersize=3 ,color=color, label=r'$\Im(\varepsilon)$')
ax3.set_ylabel(r'$\Im(\varepsilon)$', color=color , size = 18) 
ax3.tick_params(axis='y', labelcolor=color)
plt.xlim(0, 100) 
# plt.xlim(-20, 100)  
# plt.tight_layout()
# plt.legend(loc=1)

# ##### PLOT 2 #################################################################################################################
# plt.figure(figsize=(12,15))
plt.figure(figsize=(8,10))
# plt.figure()
###------------------------------------------------------------------------------##
ax1 = plt.subplot(3, 1, 1)
# ax1.plot(omegas,E2[:,0],'C1-',label='Direct' , fillstyle='none',markersize=5)
ax1.plot(omegas,E2[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
ax1.plot(omegas,E2[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax1.plot(omegas,E2[:,2],'C2:',label=r'$f_{\rho}=\lambda$',linewidth=1.8)
# plt.axvline(x=29.633,color='C4')
# plt.axvline(x=31.628,color='C4')
plt.ticklabel_format(axis='y', scilimits=(-4,-3))
plt.ylabel(r"$\vert E \vert$")
# plt.ylabel(r"$\vert H \vert$")
plt.xlabel(r"$\omega$")
# plt.yscale('symlog')
plt.yscale('log')
# ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#           ncol=3, fancybox=True, shadow=True)
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4)
# plt.legend()


####-----------------------------------######
ax2 = plt.subplot(3, 1, 2)
# ax2.plot(omegas,E0_1.real[:,0],'C1-',  fillstyle='none',label='Direct', linewidth=1.6)

ax2.plot(omegas,E0_1.real[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
ax2.plot(omegas,E0_1.real[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax2.plot(omegas,E0_1.real[:,2],'C2:' ,label=r'$f_{\rho}=\lambda$',linewidth=1.8)
plt.ylabel(r'$\Re(E_1)$')
# plt.ylabel(r'$\Re(H_1)$')
# pl.xlabel("omega")
# plt.yscale('symlog')
# plt.axvline(x=29.633,color='C4')
# plt.axvline(x=31.628,color='C4')
plt.ticklabel_format(axis='y', scilimits=(-4,-3))
# plt.legend()

ax3 = plt.subplot(3, 1, 3)
# ax3.plot(omegas,E0_1.imag[:,0],'C1-',  fillstyle='none',label='Direct', linewidth=1.6)

ax3.plot(omegas,E0_1.imag[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
ax3.plot(omegas,E0_1.imag[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
ax3.plot(omegas,E0_1.imag[:,2],'C2:' ,label=r'$f_{\rho}=\lambda$',linewidth=1.8)
plt.ylabel(r'$\Im(E_1)$')
# plt.ylabel(r'$\Im(H_1)$')
plt.xlabel(r"$\omega$")
plt.ticklabel_format(axis='y', scilimits=(-4,-3))
# # plt.legend(loc=2)


# # ##------------------------------------------------------------------------------##

plt.show()
