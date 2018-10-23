from funciones2 import *

N       = 1000000
a = 1.4
dominio = 17/a

X = samples_inversion(N,cumulativa_exponencial,dominio,a)

M = momento(1,X)
P,x = distribucion(X,np.linspace(0,dominio, num = 100))


K  = cumulante(4,X)

plt.plot(P,x,label = "histograma Exp")
plt.axvline(K[0], label = "<<X>>")
plt.axvline(K[2], label = "<<X2>>")
plt.legend()
plt.show()


