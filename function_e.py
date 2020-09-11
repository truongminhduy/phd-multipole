import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import scipy.signal as signal
from numpy.polynomial import polynomial as poly
# np.set_printoptions(precision=2)
#Mathematical constants
pi    = np.pi
j     = complex(0,1)
micro = 1e-6
femto = 1e-15
nm    = 1e-9
peta  = 1e15
#Physical constants
c     = 299792458
mu_0  = pi*4e-7
eps_0 = 1./(mu_0*c**2)
###################################################
#-------------------------------------------------------------------------------------------------------#
def Fitting_PQ(w, w_BG, Chi_dat, N_num, N_den):
    #This function :
	#(i)   computes the fitting parameters P_k's and Q_l's
    #(ii)  computes the poles Omega_j and the associated amplitudes A_j
    #(iii) sorts them by decreasing order of amplitude modulus
    j = complex(0,1)
    def error_N2(A,B):
        return linalg.norm(A-B)/linalg.norm(B)*100.

    def Coefs_pq(w_BG, Chi_dat, N_num, N_den):
        #Normalizing the variable w for stability purposes
        w_scale = 0.5*(w_BG.max() + w_BG.min())
        w_BG    = w_BG/w_scale

        #Enforcing Hermitian Symmetry
        Chi_dat   = np.concatenate((np.conj(Chi_dat[::-1]),Chi_dat), axis = 0)
        w_BG      = np.concatenate((-w_BG[::-1]           , w_BG  ), axis = 0)

        #Definition of the matrix Xi
        Xi = np.asarray([(+w_BG*j)**n if n in range(N_num+1) 
                        else -Chi_dat*(+w_BG*j)**(n-N_num) 
                        for n in range(N_num+N_den +1)]).T

        #Computation of the vector r via scipy's least squares routine
        r  = linalg.lstsq(Xi, Chi_dat)[0]
        r = np.insert(r, N_num+1, 1.0)

        P_norm = r[:N_num+1]
        Q_norm = r[N_num+1:N_num+2+N_den]

        chi_aprox  = np.asarray(sum([(+w_BG*j)**n*P_norm[n] for n in range(len(P_norm))])) 
        chi_aprox /= np.asarray(sum([(+w_BG*j)**n*Q_norm[n] for n in range(len(Q_norm))]))
        error      = error_N2(chi_aprox, Chi_dat)

        # print('--'*10)
        # print('Fitting Error PQ(%)')
        # print(error)

        # Renormalizing the P_k's and Q_l's
        p_norm = np.asarray([P_norm[n]/(w_scale)**n for n in range(len(P_norm))])
        q_norm = np.asarray([Q_norm[n]/(w_scale)**n for n in range(len(Q_norm))])

        return  p_norm, q_norm
    #----------------------------------------------------------------------------------------------#
    #Setting the p's and q's as global variables
    p_norm, q_norm = Coefs_pq(w_BG, Chi_dat, N_num, N_den)

    def Chi_final(w):
        chi  = np.asarray(sum([(+j*w)**n*p_norm[n] for n in range(len(p_norm))]))
        chi /= np.asarray(sum([(+j*w)**n*q_norm[n] for n in range(len(q_norm))]))
        return chi
    #----------------------------------------------------------------------------------------------#

    Chi_f = Chi_final(w)

    #Getting the poles by a numpy routine
    poles = np.roots(q_norm[::-1])/j

    # The Amplitudes A_j's are computed by another least squares procedure
    Mat_part_frac = np.asarray([1/(w-poles[n]) for n in range(len(poles))]).T
    Amps = linalg.lstsq(Mat_part_frac, Chi_f)[0]

    #The results are sorted by the modulus of the A_j's
    Mat_pol_amp = np.array([abs(Amps),Amps, poles])
    Mat_pol_amp = Mat_pol_amp[:, np.argsort(Mat_pol_amp[0])[::-1]]
    return p_norm, q_norm, Mat_pol_amp
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def Chi_aprox(w, p_norm, q_norm):
    #Chi as a rational function
    chi  = np.asarray(sum([(1j*w)**n*p_norm[n] for n in range(len(p_norm))]))
    chi /= np.asarray(sum([(1j*w)**n*q_norm[n] for n in range(len(q_norm))]))
    return chi
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def Chi_PP(w,Mat_pol_amp,M):
    #Principal part of Chi, this can be truncated.
    if M == None:
        M = len(Mat_pol_amp[0])
    return np.asarray(sum([Mat_pol_amp[1][n]/(w-Mat_pol_amp[2][n]) for n in range(M) ]))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def error(A,B,order):
    return linalg.norm(A-B, ord = order)/linalg.norm(B, ord = order)*100.
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def eliminate_negative(Mat_pol):
    index = [i for i in range(len(Mat_pol[2,:])) if Mat_pol[2,i].imag > 0]
    Mat_new = np.delete(Mat_pol,index,1)
    return Mat_new
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def getEpsrNumDen(Mat_pol,N_true_poles):
    ## create num and den
    list_num,list_den = [],[]
    for k in range(2*N_true_poles):
        list_num.append((Mat_pol[1,k]))
        list_den.append((-Mat_pol[2,k],1.)) 
    den = (1)
    for k in range(2*N_true_poles):
        den = poly.polymul(den, list_den[k])
    inter  = []
    pinter = (1)
    for k1 in range(2*N_true_poles):
        pinter = list_num[k1]
        for k2 in range(2*N_true_poles):
            if (k1!=k2): 
                pinter = poly.polymul(pinter, list_den[k2])
        inter.append(pinter)  
    num = inter[0]
    for k in range(1,2*N_true_poles): num = poly.polyadd(num,inter[k])
    for k in range(len(num)): num[k] /= (1j)**k
    for k in range(len(den)): den[k] /= (1j)**k
    num,den = num.real,den.real   
    # np.savetxt('text_num.out', num, delimiter=',')
    # np.savetxt('text_den.out', den, delimiter=',')
    num = poly.polyadd(num,den)    
    w2 = (0., 0., -1.)
    num2 = poly.polymul(w2,num)
    # print_w2_chi(num2,den)
    # dnum2, dden = get_der_rational(num2,den)
    return (num,den,num2)  
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def eli_positive(Mat_eps,N_tru):
    ## eliminate positive-imag-value poles
    index_positive = np.where( Mat_eps[2,:].imag > 0 )
    Mat_eps = np.delete(Mat_eps, index_positive, 1)
    return (Mat_eps,int(len(Mat_eps[2,:])/2))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def getrational(w,num,den):
    fun_w  = np.asarray(sum([(1j*w)**n*num[n] for n in range(len(num))]))
    fun_w /= np.asarray(sum([(1j*w)**n*den[n] for n in range(len(den))]))
    return (fun_w) 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def getfunction(w,num):
    fun_w  = np.asarray(sum([(1j*w)**n*num[n] for n in range(len(num))]))
    return (fun_w)     
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#    
def get_der_rational(num,den):
    dnum     = poly.polyder(num)
    dden     = poly.polyder(den)
    dnum_den = poly.polymul(dnum,den)
    dden_num = poly.polymul(num,dden)
    num_new  = poly.polysub(dnum_den,dden_num)
    den_new  = poly.polymul(den,den)
    return (num_new,den_new)     
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def print_w2_chi(num,den):
    # num = poly.polyadd(num,den)
    str_num, str_den = '{', '{'  
    for k in range(len(num)): str_num += ' %.3e,'%(num[len(num)-1-k])  
    for k in range(len(den)): str_den += ' %.3e,'%(den[len(den)-1-k])   
    str_num, str_den = str_num[:-1], str_den[:-1]   
    str_num += '}'
    str_den += '}'
    str_main =  \
'''Resolution {
  { Name Projection;
    System{{ Name M1; NameOfFormulation Ez; Type ComplexValue; }}
    Operation{
      GenerateSeparate[M1];
      EigenSolve[M1,neig,0,0,EigFilter[],
          { {-1}, {-1,0,0}, '''+ str_num +'''  } ,
          { { 1}, { 1}    , '''+ str_den +'''  } 
          ];
      SaveSolutions[M1];
    }
  }
}'''
#     str_main =  \
# '''Resolution {
#   { Name Projection;
#     System{{ Name M1; NameOfFormulation Ez; Type ComplexValue; }}
#     Operation{
#       GenerateSeparate[M1];
#         EigenSolve[M1,neig,0,0,EigFilter[],
#           { {-1}, {-1,0,0} , {-5,0,0} } ,
#           { {1} , {1}      , {1}      } ];
#       SaveSolutions[M1];
#     }
#   }
# }'''
    print(str_main, file=open('resolution.pro', 'w'))
    # print(str_num, file=open('str_num.pro', 'w'))
    # print(str_den, file=open('str_den.pro', 'w'))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#    
def print_w2_chi_new(num,den,imag,real,n_eig):
    str_num, str_den = '{', '{'  
    for k in range(len(num)): str_num += ' %.3e,'%(num[len(num)-1-k])  
    for k in range(len(den)): str_den += ' %.3e,'%(den[len(den)-1-k])   
    str_num, str_den = str_num[:-1], str_den[:-1]   
    str_num += '}'
    str_den += '}'
    str_main =  \
'''Resolution {
  { Name Projection;
    System{{ Name M1; NameOfFormulation Ez; Type ComplexValue; }}
    Operation{
      GenerateSeparate[M1];
      EigenSolve[M1, '''
    str_main +=  '%d,%.2e,%.2e' %(n_eig,imag,real)
    str_main += ''',EigFilter[],
          { {-1}, {-1,0,0}, '''+ str_num +'''  } ,
          { { 1}, { 1}    , '''+ str_den +'''  } 
          ];
      SaveSolutions[M1];
    }
  }
}''' 
    print(str_main, file=open('resolution.pro', 'w'))
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
def dumpglob2npz(filename, current_globals):
    dic2save = {}
    for glob_key in sorted(current_globals.keys()):  
        glob_type = type(current_globals[glob_key])
        if glob_type in (int, bool, float, np.float64, complex, np.ndarray, str, list, dict) and glob_key[0]!='_':
            dic2save[glob_key] = current_globals[glob_key]
    np.savez_compressed(filename,saveddic = dic2save)