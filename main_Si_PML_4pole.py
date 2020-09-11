import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import scipy.signal as signal
import matplotlib.gridspec as gridspec
from function_e import *
import os
import sys
np.set_printoptions(precision=4)
# os.environ['OPENBLAS_NUM_THREADS'] = '4'
import time
start = time.time()
## Mau's constants
pi    = np.pi
micro = 1e-6
femto = 1e-15
peta  = 1e15
cc    = 299792458
mu_0  = pi*4e-7
eps_0 = 1./(mu_0*c**2)
#------------------------------------------------------------------------------------------#
# filetxt = "mau/Au_Johnson_data.txt"
filetxt = "mau/Si_Green_data.txt"
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
Lam_dat, n_dat, k_dat = np.loadtxt(filetxt, usecols=(0,1,2), skiprows = 1, unpack = True )
w_BG    = 2*pi*c/(Lam_dat*micro)/peta
Chi_dat = (n_dat + 1j*k_dat)**2-np.ones_like(w_BG)
#------------------------------------------------------------------------------------------#
N_tru = 4
N_num, N_den = 10, 10
# N_num, N_den = 6, 6
w = np.linspace(-w_BG.max()*1.5, w_BG.max()*1.5,100)
p_norm, q_norm, Mat_pol_amp = Fitting_PQ(w, w_BG, Chi_dat, N_num, N_den)

##############################################################################################################################
###### MAIN PROGRAM ##########################################################################################################
##############################################################################################################################
##--- PARAMETER ----------------#########################
## Incident plane wave parameters
scale_light = 1e8
micro       = 1e-6
scale_frequ = scale_light/micro # 1e14
cel = cc/scale_light
mu0 = mu_0  *scale_light
ep0 = eps_0 *scale_light

# e_scale = 1e-1
# e_scale = 1e-1
e_scale = 8e-2
r_ellipse_1, r_ellipse_2 = 2.5*e_scale, 1.5*e_scale
# r_ellipse_1, r_ellipse_2 = 2.5*e_scale, 2.5*e_scale
xS, yS   = -3  *e_scale, 1  *e_scale
xD1, yD1 =  2.5*e_scale, 1.5*e_scale
# xD1, yD1 =  1  *e_scale, 2  *e_scale
xD2, yD2 =  2  *e_scale,-2  *e_scale

lambda_m  = 10*e_scale
# lambda_m  = 8*e_scale
lambda_vp = 10*e_scale
lambda0   = 6.5*e_scale
# lambda0   = 6.5*e_scale
R_pml_in  = 3.5*e_scale
R_pml_out = R_pml_in+lambda0
r_source  = 1.5e-1*e_scale

E_or_H = 0 # 0 or 1

eps_mil = 1.
par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'w')
par_gmsh_getdp.write('lambda_m  = %3.3e;\n'%(lambda_m))
par_gmsh_getdp.write('lambda_vp = %3.3e;\n'%(lambda_vp))
# par_gmsh_getdp.write('A         = %3.3e;\n'%(A))
par_gmsh_getdp.write('R_pml_in  = %3.3e;\n'%(R_pml_in))
par_gmsh_getdp.write('R_pml_out = %3.3e;\n'%(R_pml_out))
par_gmsh_getdp.write('r_ellipse_1 = %3.3e;\n'%(r_ellipse_1))
par_gmsh_getdp.write('r_ellipse_2 = %3.3e;\n'%(r_ellipse_2))
par_gmsh_getdp.write('r_source  = %3.3e;\n'%(r_source))
par_gmsh_getdp.write('xS        = %3.3e;\n'%(xS))
par_gmsh_getdp.write('yS        = %3.3e;\n'%(yS))
par_gmsh_getdp.write('xD1       = %3.3e;\n'%(xD1))
par_gmsh_getdp.write('yD1       = %3.3e;\n'%(yD1))
par_gmsh_getdp.write('xD2       = %3.3e;\n'%(xD2))
par_gmsh_getdp.write('yD2       = %3.3e;\n'%(yD2))
par_gmsh_getdp.write('cel       = %3.3e;\n'%(cel))
par_gmsh_getdp.write('mu0       = %3.3e;\n'%(mu0))
par_gmsh_getdp.write('epsilon0   = %3.3e;\n'%(ep0))
par_gmsh_getdp.write('eps_mil_re = %3.3e;\n'%(eps_mil.real))
par_gmsh_getdp.write('eps_mil_im = %3.3e;\n'%(eps_mil.imag))
par_gmsh_getdp.write('angle      = Pi/10;\n')
# par_gmsh_getdp.write('angle      = Pi/6;\n')
par_gmsh_getdp.write('E_or_H     = %d;\n'%(E_or_H))
par_gmsh_getdp.close()

str_gmsh_path  = ''
str_getdp_path = ''
# mesh_filename  = 'mesh_ellipse.msh'
# mesh_geo       = 'mesh_ellipse.geo'
mesh_filename  = 'mesh_ellipse_PML2D.msh'
mesh_geo       = 'mesh_ellipse_PML2D.geo'
os.system(str_gmsh_path+'gmsh '+mesh_geo+' -2 -o '+mesh_filename)


N_tru = 4
##----change the unit from Mau to our problem-----##
Mat_pol_amp = eliminate_negative(Mat_pol_amp)
Mat_eps = Mat_pol_amp * (peta/scale_frequ)
os.system('rm resolution.pro')
num,den,num2 = getEpsrNumDen(Mat_eps,N_tru) # get function N(iw)/D(iw) and [N(iw)*w^2/D(iw)]'
dnum, dden = get_der_rational(num2,den)
solve_den = den[::-1]
pole_root = np.roots(solve_den)

omegas = 2*pi*cel/Lam_dat
Chis = getrational(omegas,num,den)
Chi_list = Chi_PP(omegas,Mat_eps,2*N_tru)

print (Mat_eps)
print (pole_root.imag)
print (-pole_root.real)
# # ------------------------------------------------------------------------------------------#
# plt.figure()
# gs = gridspec.GridSpec(2, 1)
# ax1 = plt.subplot(gs[0,0])
# ax1.plot(omegas, Chis.real, 'C0-.') 
# ax1.plot(w_BG*10, Chi_dat.real, 'C03') 
# ax1.plot(omegas, Chi_list.real+1, 'C1--') 
# ax2 = plt.subplot(gs[1,0])
# ax2.plot(omegas, Chis.imag, 'C0-.') 
# ax2.plot(w_BG*10, Chi_dat.imag, 'C03') 
# ax2.plot(omegas, Chi_list.imag, 'C1--') 
# # plt.show()
# # #------------------------------------------------##

######## SPECTRAL PROBLEM ########################
# n_eig = np.array([150])  # pi/10 8e-2 +6.5
n_eig = np.array([150,100,250,200,150,52])  # pi/10 8e-2 +6.5
# n_eig = np.array([150,70,200,150,150,40])  # pi/6 8e-2 +6.5
# n_eig = np.array([150,100,250,200,90,70])  # pi/20 8e-2 +6.5
# n_eig = np.array([130,70,150,100,50,32])  # pi/10 8e-2 +4.5
# n_eig = np.array([20,20,100,80,10,7])  # pi/10 3e-2
n_region = len(n_eig)
n_cum    = np.cumsum([0, *n_eig])
eigs     = np.zeros(sum(n_eig),dtype=complex)
dd       = np.zeros(sum(n_eig),dtype=complex)
####0
print_w2_chi_new(num2,den,0.3,43,n_eig[0])
# print_w2_chi_new(num2,den,2,32,n_eig[0]) # pi/10 3e-2
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints -50,50,0,50'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral0.res')
eigs[range(n_cum[0],n_cum[1])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
# ####1
print_w2_chi_new(num2,den,4.7,54.7,n_eig[1]) 
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 0,50,50,59.5'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral1.res')
eigs[range(n_cum[1],n_cum[2])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
####2
print_w2_chi_new(num2,den,4.5,62,n_eig[2]) 
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 0,50,59.5,64.6'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral2.res')
eigs[range(n_cum[2],n_cum[3])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
####3
print_w2_chi_new(num2,den,8,71,n_eig[3]) #pi/8
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints -50,50,64,71.7'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral3.res')
eigs[range(n_cum[3],n_cum[4])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
####4
##print_w2_chi_new(num2,den,15,73,n_eig[4]) 
print_w2_chi_new(num2,den,15,75,n_eig[4]) 
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 0,100,72.08,120'
# slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 15,100,72.08,100'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral4.res')
eigs[range(n_cum[4],n_cum[5])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])
####5
# print_w2_chi_new(num2,den,6.5,65,n_eig[5]) 
print_w2_chi_new(num2,den,12,55,n_eig[5]) 
# slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 4,50,50,72.08' # 3e-2
slepc_options_rational = ' -nep_max_it 100 -nep_type nleigs -nep_rational -nep_target_magnitude -rg_interval_endpoints 7,50,50,72.08'
eigen_getdp = str_getdp_path+'getdp spectral.pro -pre Projection -msh '+mesh_filename+' -cal -pos postop_e -slepc'# -v 0' 
os.system(eigen_getdp+ slepc_options_rational)
os.system('mv spectral.res spectral5.res')
eigs[range(n_cum[5],n_cum[6])] = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])

###----------------------------##
# eigs = eigs[range(n_cum[4],n_cum[5])]
# ### plot eigen mode -----------###
plt.figure()
plt.plot(eigs.real, eigs.imag, 'C0o',markersize=3 ,mfc='none')
# for k in range(len(eigs)):
#       plt.annotate(k, xy=(vpr[k],vpi[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
plt.plot(pole_root.imag,-pole_root.real,'C2o')
plt.title("Eigenvalues")
plt.xlabel("Re")
plt.ylabel("Im")

# # ########## NORM ##########################
for x in range(n_region):
	os.system('rm norm1.txt norm2.txt')
	norm_getdp = str_getdp_path+'getdp norm.pro -pre Projection -msh '+mesh_filename+' -res spectral'+str(x)+'.res -pos postop_norm -v 0'
	os.system(norm_getdp)
	norm1 = np.loadtxt('norm1.txt',usecols=[5]) + 1j*np.loadtxt('norm1.txt',usecols=[6] ) 
	norm2 = np.loadtxt('norm2.txt',usecols=[5]) + 1j*np.loadtxt('norm2.txt',usecols=[6] ) 
	deri  = 1j*getrational(eigs[range(n_cum[x],n_cum[x+1])],dnum,dden) 
	d1 = 2/cel**2 * eigs[range(n_cum[x],n_cum[x+1])] * norm1 
	d2 = 1/cel**2 * deri * norm2 
	dd[range(n_cum[x],n_cum[x+1])] = d1+d2

###################################################################################################################################
N_frequency = 150
E0_0 = np.zeros((N_frequency,4),dtype=complex)
E0_1 = np.zeros((N_frequency,4),dtype=complex)
E0_2 = np.zeros((N_frequency,4),dtype=complex)
E2   = np.zeros((N_frequency,4))
omegas = np.linspace(20,90,N_frequency)
for i,omega0 in enumerate(omegas) : 
	lambda0 = 2.*pi*cel/omega0	
	epsr    = getrational(omega0,num,den)
	# epsr    = 1
	print(epsr)
	par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'a')
	par_gmsh_getdp.write('lambda0 = %3.3e;\n'%(lambda0))
	par_gmsh_getdp.write('eps_diff_re = %3.3e;\n'%(epsr.real ))
	par_gmsh_getdp.write('eps_diff_im = %3.3e;\n'%(epsr.imag ))
	par_gmsh_getdp.close()
	############ direct ##################################################################
	direct_getdp = str_getdp_path+'getdp direct.pro -pre Scattering -msh '+mesh_filename+' -cal -pos postop_scat -v 0' 
	os.system(direct_getdp)
	E2[i,0] = np.loadtxt('E2.txt')[1]
	E0_0[i,0] = np.loadtxt('E0_0.txt',usecols=[5]) + 1j*np.loadtxt('E0_0.txt',usecols=[6] ) 
	E0_1[i,0] = np.loadtxt('E0_1.txt',usecols=[5]) + 1j*np.loadtxt('E0_1.txt',usecols=[6] )
	E0_2[i,0] = np.loadtxt('E0_2.txt',usecols=[5]) + 1j*np.loadtxt('E0_2.txt',usecols=[6] )

# 	############ Projection 0 ###############################################################################################################################
	for x in range(n_region):
		reso = eigs[range(n_cum[x],n_cum[x+1])]
		os.system('rm Jns.txt')
		Jn_getdp = str_getdp_path+'getdp Jn.pro -pre Projection -msh '+mesh_filename+' -res spectral'+str(x)+'.res -pos postop_Jn -v 0'
		os.system(Jn_getdp)
		Jn  = np.loadtxt('Jns.txt',usecols=[5]) + 1j*np.loadtxt('Jns.txt',usecols=[6] ) 
		Pns = Jn/(omega0-reso)/dd[range(n_cum[x],n_cum[x+1])]
		file_Pns = open('Pns.dat', 'w')
		file_Pns.write('neig = %d;\n'%(n_eig[x]))
		for k, Pn in enumerate(Pns):
			file_Pns.write('Pns_re_%d = %3.3e;\n'%(k,Pn.real))
			file_Pns.write('Pns_im_%d = %3.3e;\n'%(k,Pn.imag))
		file_Pns.close()
		project_getdp = str_getdp_path+'getdp projection.pro -pre Projection -msh  '+mesh_filename+' -res spectral'+str(x)+'.res -pos  postop_'+str(x)+' -bin -v 0'
		os.system(project_getdp)
	#### Combine ---------------------------
	projection_coeff_getdp = str_getdp_path+'getdp projection_coeff_6.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
	os.system(projection_coeff_getdp)
	E2[i,1] = np.loadtxt('E2.txt')[1]
	E0_0[i,1] = np.loadtxt('E0_0.txt',usecols=[5]) + 1j*np.loadtxt('E0_0.txt',usecols=[6] ) 
	E0_1[i,1] = np.loadtxt('E0_1.txt',usecols=[5]) + 1j*np.loadtxt('E0_1.txt',usecols=[6] )
	E0_2[i,1] = np.loadtxt('E0_2.txt',usecols=[5]) + 1j*np.loadtxt('E0_2.txt',usecols=[6] )

	############# Projection 1 ###############################################################################################################################
	for x in range(n_region):
		reso = eigs[range(n_cum[x],n_cum[x+1])]
		os.system('rm Jns.txt')
		Jn_getdp = str_getdp_path+'getdp Jn.pro -pre Projection -msh '+mesh_filename+' -res spectral'+str(x)+'.res -pos postop_Jn -v 0'
		os.system(Jn_getdp)
		Jn  = np.loadtxt('Jns.txt',usecols=[5]) + 1j*np.loadtxt('Jns.txt',usecols=[6] ) 
		Pns = Jn/(omega0-reso)/dd[range(n_cum[x],n_cum[x+1])]*(reso/omega0)
		file_Pns = open('Pns.dat', 'w')
		file_Pns.write('neig = %d;\n'%(n_eig[x]))
		for k, Pn in enumerate(Pns):
			file_Pns.write('Pns_re_%d = %3.3e;\n'%(k,Pn.real))
			file_Pns.write('Pns_im_%d = %3.3e;\n'%(k,Pn.imag))
		file_Pns.close()
		project_getdp = str_getdp_path+'getdp projection.pro -pre Projection -msh  '+mesh_filename+' -res spectral'+str(x)+'.res -pos  postop_'+str(x)+' -bin -v 0'
		os.system(project_getdp)
	#### Combine ---------------------------
	projection_coeff_getdp = str_getdp_path+'getdp projection_coeff_6.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
	os.system(projection_coeff_getdp)
	E2[i,2] = np.loadtxt('E2.txt')[1]
	E0_0[i,2] = np.loadtxt('E0_0.txt',usecols=[5]) + 1j*np.loadtxt('E0_0.txt',usecols=[6] ) 
	E0_1[i,2] = np.loadtxt('E0_1.txt',usecols=[5]) + 1j*np.loadtxt('E0_1.txt',usecols=[6] )
	E0_2[i,2] = np.loadtxt('E0_2.txt',usecols=[5]) + 1j*np.loadtxt('E0_2.txt',usecols=[6] )

# os.system('gmsh '+mesh_geo+' us.pos uss.pos ')
# # os.system('gmsh '+mesh_geo+' us.pos ')
dumpglob2npz('output.npz',globals()) 

# print (E2)
# print (E0_0)

end = time.time()
print(end - start)
plt.show()