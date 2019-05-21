# python layer.py leer.dat 
from math import atan,sqrt,pi,floor
import sys
import numpy as np

# encontrar dos picos
def pf(dist,n):
    pk=np.zeros((n))
    j=0
    while j < n:
        midx = np.argmax(dist[:,1]) # indice del maximo
        rm=dist[midx,0]
        npart=dist[midx,1]
        #print rm, npart
        if dist[midx-1,1] < npart and \
                dist[midx+1,1] < npart:
                if not(any(np.abs(pk-rm)<.5)):
                    pk[j] = rm
                    dist[midx,1]=0
                    j=j+1
                else:
                    dist[midx,1]=0
                         
    return pk
                        


#separar las capas
def lyr(pos,pk,fname,r0): # r0 es el "origen"
    for i in range(len(pk)):
        f  = open("lyr"+str(i+1)+fname,"w") 
        f.write(str(len(pos))+"\n") # este no es el numero de part
        f.write("\n")
        for j in range(len(pos)):
            rxy=np.linalg.norm(pos[j,0:2]-r0)
            print pos[j,0:2],rxy
            if pk[i]-.5 < rxy < pk[i]+.5:
                f.write("%s %s %s\n"%(pos[j,0],pos[j,1],pos[j,2]))
        f.close()
            




#archivos a leer
dfile = "dist"+ str(sys.argv[1])
pfile = str(sys.argv[1])
dist=np.loadtxt(dfile)
pos=np.loadtxt(pfile,skiprows=2,usecols=(1,2,3))
# lx.xxrx.xxrhox.xx.xyz
lamb = float(str(sys.argv[1])[1:4])
rcyl = float(str(sys.argv[1])[6:9])
# calculamos r0
#r0 = 2*(rcyl+1.5*lamb)
r0 = (rcyl+1.5*lamb)

pks=pf(dist,2)
print pks
lyr(pos,pks,pfile,r0)

