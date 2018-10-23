import numpy as np

################## funciones q voy a usar #########################

def cumulativa_exponencial(Lambda1,t):
    Lambda1 = Lambda1
    return 1 - np.exp(-Lambda1*t) 

def cumulativa_marcos(q,p,a,b,t):
    b = b 
    a = a
    q = q
    p = p
    Lambda1 = a*b/(a + p*(b-a))
    return 1 - np.exp(-Lambda1*t) 

def cumulativa_intermitente(q,p,Lambda1,Lambda2,t):
    q = q
    p = p
    Lambda1 = Lambda1
    Lambda2 = Lambda2
    return (p - p*np.exp(-Lambda1*t) + q - q*np.exp(-Lambda2*t))

def distribucion_exponencial(Lambda1,t):
    Lambda1 = Lambda1
    return Lambda1*np.exp(-Lambda1*t)

def distribucion_intermitente(q,p,Lambda1,Lambda2,t):
    q = q
    p = p
    Lambda1 = Lambda1
    Lambda2 = Lambda2
    return p*Lambda1*np.exp(-Lambda1*t) + q*Lambda2*np.exp(-Lambda2*t)

################################################################


