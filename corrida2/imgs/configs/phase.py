import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import sys
import funcs as fun

#cargar datos
data=np.loadtxt(sys.argv[1])
#ordenar por radios
data=fun.qsort(data,0,len(data)-1,1)
#graficar para cada radio \rho \lambda
r=data[0,1]
#    #lambda r_cyl rho hexord r_layer packaging_factor #layers  <ncoord>
# separar por radios
rr=np.zeros(7)
c=0
for i in range(len(data)):
    if not(r==data[i,1]):
        r=data[i,1]
        rr[c]=i
        c+=1
rr=[int(i) for i in rr]
rc('text', usetex=True)

for j in range(len(rr)-1):
    sc=plt.scatter(data[rr[j]:rr[j+1],2],data[rr[j]:rr[j+1],0],c=data[rr[j]:rr[j+1],3],cmap='cool')
    plt.colorbar(sc)
    plt.ylabel(r'\lambda')
    plt.xlabel(r'\rho')
    plt.title("Radius of cylinder =={}".format(data[rr[j],1]))
    plt.show()
    
