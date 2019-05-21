import matplotlib.pyplot as plt
import numpy as np
from qtree import *
from funcs import *

a=np.array([1,0,0])
b=np.array([np.cos(np.pi/3),np.sin(np.pi/3),0])
ii=np.array([1,1,0])
hexgrid=np.zeros((49,3))
sc=np.zeros((49,3))
bsc=np.array([0,1,0])
c=0
for i in range(7):
    for j in range(7):
        rr=np.random.rand(1,2)
        rr=np.append(rr,0)
        rr=(ii-2*rr)*.2
        hexgrid[c,:] = a*i+b*j+rr
        sc[c,:] = a*i+bsc*j
        c+=1

        
#plotp(hexgrid[:,0:2])
#tree = qtree(cell([0,0],[15,9]))
#tree.insertArray(hexgrid[:,0:2])
#tree.plotree()
#plotp(hexgrid[:,0:2])
hexord(hexgrid,1.1)
#hexord(sc,1.1)
print(hexgrid[:,2])
a=plt.scatter(hexgrid[:,0], hexgrid[:,1], c=hexgrid[:,2],cmap='cool')
#a=plt.scatter(sc[:,0], sc[:,1], c=sc[:,2],cmap='cool')
#print(sc[:,2])
plt.colorbar(a)
plt.clim(min(hexgrid[:,2]),max(hexgrid[:,2]))
plt.show()


