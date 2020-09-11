import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as pl
import os
import sys
import scipy.signal as signal
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from function_e import *

np.set_printoptions(precision=3)
pi = np.pi
#### MAIN 1 #######################################################################
### define inpute files ###########################################################
# path = '../../../Downloads/gmsh-4.3.0-Linux64/bin/'
path = ''
# os.system(path+'gmsh mesh_ellipse_PML2D.geo us.pos uss.pos -merge compare.geo')
# os.system(path+'gmsh mesh_ellipse_PML2D.geo u.pos')
# os.system(path+'gmsh mesh_ellipse_PML2D_closed.geo us.pos uss.pos -merge compare.geo')
# os.system(path+'gmsh mesh_duy_simple.geo u.pos')
os.system(path+'gmsh mesh_ellipse_PML2D_closed.geo ')
# os.system(path+'gmsh mesh_triangle_closed.geo us.pos up1.pos -merge compare.geo')
# os.system('gmsh mesh_triangle.geo us.pos uss.pos up1.pos up2.pos up0.pos')
# os.system(path+'gmsh mesh_triangle_closed.geo us.pos up1.pos -merge compare.geo')
# os.system(path+'gmsh mesh_ellipse_PML2D.geo us.pos up0.pos -merge compare.geo')
# os.system(path+'gmsh mesh_triangle.geo us.pos')
# os.system(path+'gmsh mesh_triangle.msh ')


