# esto recibe un archivo xyz 
# todo en python 3
# python3 config.py *.xyz path
import numpy as np
import sys
from funcs import *
from qtree import *
from octree import *
import matplotlib.pyplot as plt
fname=str(sys.argv[1])
#path=str(sys.argv[2])

# obtener parametros del sistema
lamb,r,rho=getparam(str(sys.argv[1]))
r0=(r+1.5*lamb)

# leer archivo
xyz=np.loadtxt(fname,skiprows=2,usecols=(1,2,3))
npart=len(xyz)

# cambio coord y dist
rtz,dist=cart2cyl(xyz,r0)

#   rtz[i,:] -> [x,y,z,r,r*\theta,parametro]
# ordenar con respecto al radio
qsort(rtz,0,len(rtz)-1,3)

# separar las capas 
pk,npks=pf(dist,2)
lyrs=lyr(rtz[:,3],pk[:,0],.2)

# parametro de orden para las capas
#for i in range(len(pk)):
#    hexord(rtz[int(lyrs[i,0]):int(lyrs[i,1]),3:6],.5*lamb)
#xyzwrite(rtz,fname[:-4]+"hx"+".xyz",6) # escribe todo en un archivo

# calcular promedio de ncoord
#eta=ncoord(xyz,lamb)
#nbar=np.dot(eta[:,0],eta[:,1])/float(npart)

#for i in range(len(pk)):
#    f=open(path+"lyr"+str(i)+'.dat','a')
#    avg=np.average(rtz[int(lyrs[i,0]):int(lyrs[i,1]),5])
#    p_f=pk[i,1]/(pk[i,0]*2*40)
#    #lambda r_cyl rho hexord r_layer packaging_factor #layers  <ncoord>
#    f.write('{} {} {} {} {} {} {} {}\n'.format(lamb,r,rho,avg,pk[i,0],p_f,npks,nbar))
#f.close()

#plt.plot(ncoord(xyz,lamb))
#plt.xlabel('ncoord')
#plt.ylabel('npart')
#plt.show()
#plt.xticks(np.arange(0, 20,1))
#plt.savefig(fname[:-4]+".png")
