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

# ##### PLOT #################################################################################################################
res = np.load('pole1_closed_TE/output.npz')['saveddic'].item()
E2   = res['E2']
E0_0 = res['E0_0']
E0_1 = res['E0_1']
E0_2 = res['E0_2']
omegas = res['omegas']

# res = np.load('new_pole1_closed/output_D.npz')['saveddic'].item()
res = np.load('pole1_closed_TE/output_lambda.npz')['saveddic'].item()
# res = np.load('output.npz')['saveddic'].item()
E2a   = res['E2']
omegasa = res['omegas']

# tab_colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8']
# plt.figure(figsize=(12,15))
plt.figure(figsize=(9,10))
gs = gridspec.GridSpec(3, 1, height_ratios=[4,1,4]) 
# plt.figure()
###------------------------------------------------------------------------------##
ax1 = plt.subplot(gs[0])
# ax1.plot(omegas,E2[:,0],'C1-',label='Direct' , fillstyle='none',markersize=5)
ax1.plot(omegas,E2[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
# ax1.plot(omegas,E2[:,1],'C0-',label=r'$f_{\rho}=1$', linewidth=1.6)
# ax1.plot(omegas,E2[:,2],'C2:',label=r'$f_{\rho}=\lambda$',linewidth=1.8)
# ax1.plot(omegas,E2[:,2],'C2-.',label=r'$f_{\rho}=\lambda$',linewidth=1.8)
ax1.plot(omegas,E2[:,3],'C3-.',label=r'$f_{\rho}=\lambda-\lambda_0$',linewidth=1.8)
# ax1.plot(omegas,E2[:,4],'C4-.',label=r'$f_{\rho}=\lambda-\lambda_0$',linewidth=1.8)
ax1.plot(omegasa,E2a[:,1],'C0-',label=r'$f_{\rho}=\lambda^2$',linewidth=1.8)
# ax1.plot(omegasa,E2a[:,2],'C5-.',label=r'$f_{\rho}=\lambda^3$',linewidth=1.8)
# plt.axvline(x=29.633,color='C4')
plt.axvline(x=31.628,color='C4')
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

res = np.load('pole1_closed_TE/output_D.npz')['saveddic'].item()
# res = np.load('output.npz')['saveddic'].item()
E2   = res['E2']
omegas = res['omegas']
ax2 = plt.subplot(gs[2])
ax2.plot(omegas,E2[:,0],'C18',label='Direct' , fillstyle='none',markersize=5)
ax2.plot(omegas,E2[:,1],'C0-',label=r'$g_{\sigma}=1$', linewidth=1.6)
ax2.plot(omegas,E2[:,2],'C3:',label=r'$g_{\sigma}=\lambda$',linewidth=1.8)
ax2.plot(omegas,E2[:,3],'C2-.',label=r'$g_{\sigma}=\lambda^2$',linewidth=1.8)
plt.ticklabel_format(axis='y', scilimits=(-4,-3))
plt.ylabel(r"$\vert E \vert$")
plt.xlabel(r"$\omega$")
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4)

# # ##------------------------------------------------------------------------------##

plt.show()
