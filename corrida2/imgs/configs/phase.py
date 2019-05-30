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
#       0      1    2     3       4       5              6        7
#    #lambda r_cyl rho hexord r_layer packaging_factor #layers  <ncoord>
# separar por radios
rr=np.zeros(7)
c=0
for i in range(len(data)):
    if not(r==data[i,1]):
        r=data[i,1]
        rr[c]=i
        c+=1
rr[-1]=len(data)
rr=[int(i) for i in rr]
rc('text', usetex=True)
data[:,5]=data[:,5]/np.pi
for j in range(len(rr)-1):
    # hexatic
    sc=plt.scatter(data[rr[j]:rr[j+1],2],data[rr[j]:rr[j+1],0],c=data[rr[j]:rr[j+1],3],cmap='hsv',vmin=0,vmax=1,s=125)
    plt.ylabel(r'$\lambda$')
    plt.xlabel(r'$\rho$')
    plt.title("Radius of cylinder = {}".format(data[rr[j],1]))
    cbar=plt.colorbar(sc)
    cbar.ax.set_ylabel('hexatic order parameter '+r'$\psi_6$',rotation=270,labelpad=15)
    plt.savefig("hexatic_{}_layer_r{}.png".format(str(sys.argv[2]),data[rr[j],1]))
    plt.clf()
    #ncoord
    sc=plt.scatter(data[rr[j]:rr[j+1],2],data[rr[j]:rr[j+1],0],c=data[rr[j]:rr[j+1],7],cmap='hsv',s=125)
    plt.ylabel(r'$\lambda$')
    plt.xlabel(r'$\rho$')
    plt.title("Radius of cylinder = {}".format(data[rr[j],1]))
    cbar=plt.colorbar(sc)
    cbar.ax.set_ylabel('Coordination number',rotation=270,labelpad=15)
    plt.savefig("ncoord_{}_layer_r{}.png".format(str(sys.argv[2]),data[rr[j],1]))
    plt.clf()
    #packaging
    sc=plt.scatter(data[rr[j]:rr[j+1],2],data[rr[j]:rr[j+1],0],c=data[rr[j]:rr[j+1],7],cmap='hsv',s=125)
    plt.ylabel(r'$\lambda$')
    plt.xlabel(r'$\rho$')
    plt.title("Radius of cylinder = {}".format(data[rr[j],1]))
    cbar=plt.colorbar(sc)
    cbar.ax.set_ylabel('Packaging factor',rotation=270,labelpad=15)
    plt.savefig("packaging_{}_layer_r{}.png".format(str(sys.argv[2]),data[rr[j],1]))
    plt.clf()
    
