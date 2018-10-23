import os
import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import interp1d
from math import factorial as fac

def binomial(x, y):
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def cumulativa_exponencial(Lambda1,t):
    Lambda1 = Lambda1
    return (1.0 - np.exp(-Lambda1*t))


def samples_inversion(numero_de_muestras,cumulativa,dominio,*parametros):

    """ esta funcion genera N muestras de una v.a 
    mediante la inversion de la CFD """
    
    param = parametros

    cumul = cumulativa
   
    dominio = dominio 
    
    tc   =   np.max(float(dominio))
    t   =   np.linspace(0, tc, num = int(1e7))
    y0   =   cumul(param[0],t)

    x    =   np.random.random(numero_de_muestras)
    fint =   interp1d(y0,t,copy = False)
    X    =   fint(x)
    return X


def momento(numero_de_momento,muestras):

    """ Calculo de momentos """ 
    nn = numero_de_momento
    muestras = muestras

    Mn = np.sum((muestras**nn))/np.size(muestras)

    return Mn


def cumulante(numero_de_cumulante_max,muestras):

    """ Calculo de cumulantes """ 
    n = numero_de_cumulante_max
    muestras = muestras 

    C  = np.diag(np.ones(n-1),1)
    M = np.array([momento(i,muestras) for i in range(1,n+1)])
    
    C[0:,0] = M
    for i in range(1,n):
        C[i:,i] = M[:-i]
    
    for i in range(2,n):
        for j in range(2,n):
            C[i,j] = binomial(i,j-1)*C[i,j] 
    
    K = np.zeros(n)
    for k in range(0,len(K)):
        K[k] = (-1)**(k)*np.linalg.det(C[:k+1,:k+1])#

    return K


def distribucion(muestras,bines):

    Pi, bin_edges = np.histogram(muestras,bines)


    Db = np.diff(bin_edges)
    Pi = Pi/float(len(muestras))/Db
    x  = np.cumsum(Db) + bin_edges[0]/2.0
    return Pi, x 

#def samples_inversion(numero_de_muestras,a,b,p):
#    a = a
#    b = b
#    p = p
#    q = 1 - p
#    tc = 16/max(a,b)
#    Tf = 16/min(a,b)
#    
#    t1 = np.linspace(0,tc, num = int(1e7))
#    t2 = np.linspace(tc,Tf,num = int(1e7))
#    t0 = np.hstack([t1[:-1],t2])
#    y0 = p - p*np.exp(-a*t0) + q - q*np.exp(-b*t0)
#    
#    x  = np.random.random(numero_de_muestras)
#    
#    fint       = interp1d(y0,t0,copy = False)
#    Ti         = fint(x)
#    Ti[0]    = 0.0
#
#    return Ti, t0

#def ruido_intermitente(a,b,p,N):
##3 generacion de realizaciones de ruido dictomico 
#    N  = N
#    Ti, t0 = samples_inversion(N,a,b,p)
#    intermitente_x = np.zeros(N) # "posicion"
#    intermitente_t = np.zeros(N) # tiempo 
#    
#    intermitente_x[:]    =  np.ones(N)
#    intermitente_x[1::2] = -1*np.ones(N/2)
#    intermitente_t[:]    =  np.cumsum(Ti)
#
#    return Ti, intermitente_x, intermitente_t
#
#
#def extender_ruido(Dt,ruido):
#
#    Dt    = Dt
#    ruido_t = ruido
#    #gamma = (np.exp(gamma0*Dt) -1)/Dt
#    #Dt    = 1e-2
#    TI = (ruido_t//Dt).astype(int)
#    TI = np.unique(TI)
#    I = np.zeros(TI[-1]+1)
#    I[TI[0::2]] =  2
#    I[TI[1::2]] = -2
#    I[0] = 1
#    I[-1] = 0 
#    ruido = np.cumsum(I)
#    t = np.arange(0,round(ruido_t[-1], int(1/Dt)),Dt)
#    return ruido, t
#
#
#def integrar(gamma,theta,Dt,ruido):
#    Dt   =  Dt    
#    theta = theta 
#    gamma = gamma 
#    ruido = ruido   
#    V = np.zeros(len(ruido))
#    for i in range(0,len(ruido)-1):
#        V[i+1] =  V[i] + (-gamma*V[i] + theta*ruido[i])*Dt
#    return V#
##V = integrar(1,1e-2,ruido)
#
#def Pn(V):
#    DV  = np.linspace(min(V),max(V), num = 100)
#    PVerg, bins = np.histogram(V ,DV, normed = True) 
#    Vbins = DV[:-1] + (DV[-1] - DV[-2])/2
#    return PVerg, Vbins
#
#def P(X,a,gamma,theta):
#    X = X
#    a = a 
##def P(a,b,p,gamma,theta):
##    theta = theta
##    a2  = a
##    b2  = b 
##    p = p 
##    q = 1- p
##    gamma = gamma 
##    num = (a2*q + b2*p)*a2*b2
##    den = ((a2**2)*q*(1+p) + (b2**2)*p*(1+q) - 2*a2*b2*p*q)
##    A   = (num/den/gamma - 1)
##    B   = ((gamma/theta/(a2*q+b2*p))**2)*den
#    A  = a/gamma - 1
#    B  = (gamma/theta)**2
#    N   = np.sqrt(B)*fungamma(A + 1.5)/fungamma(A + 1)/np.sqrt(np.pi)
##    X = np.linspace(-1/np.sqrt(B),1/np.sqrt(B),100)
#    Pst = N*(1 - B*(X**2))**A
#    return Pst
#
#
################ cheque histograma de tiempos ################
#
##Ti, t0 = samples_inversion(n*R,a,b,p)
##bins = t0[0::50000]
##data = Ti
##hist, bin_edges = np.histogram(data,bins) # make the histogram
##
##Db = np.zeros(len(bins)-1)
##for i in range(0,len(bins)-1):
##    Db[i] = (bins[i+1] - bins[i])
##
##Pt2 = hist/float(n*R)/Db
##tbins = bins[0:-1] + Db*0.5
##Pt2a = distribucion_intermitente(q,p,a,b,tbins)
##
##fig3 =  plt.figure()
##plt.title('histograma vs. distribucion')
##plt.plot(tbins, Pt2a ,tbins, Pt2 )
###plt.axvline(tc)
###plt.axhline(distribucion_intermitente(q,p,a,b,tc))
##plt.xscale('log')
##plt.yscale('log')
###############################################################
#################################################################
############### Graficas #######################################
#
##fig0 = plt.figure()
##fint2  = interp1d(intermitente_t[:11] ,intermitente_x[:11], kind = 'zero', axis = 0, bounds_error = False, fill_value = 1)
##Dt0 = 0.01
#t1 = np.arange(0,intermitente_t[10],Dt0) 
#plt.plot(t1, fint2(t1), intermitente_t[0:10], intermitente_x[0:10], 'o')
#

#ruido[0,:] = np.sin(2*np.pi*t)
#plt.plot(intermitente_t[0:8], dicotomico_x[0:8], 'ro', t[0:500], fint2(t[0:500]))

#plt.xlabel('tiempo [T]')
#plt.ylabel(r'$\theta [L/T^{2}]$')
plt.show() # muestra graficas

################## guardar datos #################################
#path = os.getcwd()
#path = path+'/datos/'
#
#
#np.save(path+'a'     , np.array([a]))
#np.save(path+'p'     , np.array([p]))
#np.save(path+'q'     , np.array([q]))
#np.save(path+'b'     , np.array([b]))
##np.save(path+'R'     , np.array([R]))
#np.save(path+'n'     , np.array([n]))
#
################# realizaciones #####################################
#
#np.save(path+'Ti'      , Ti                 )
#np.save(path+'intermitente_t' , intermitente_t )
#np.save(path+'intermitente_x' , intermitente_x )
