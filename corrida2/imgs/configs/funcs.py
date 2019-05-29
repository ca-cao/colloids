from math import atan,sqrt,pi,floor,ceil
import numpy as np
import matplotlib.pyplot as plt
from qtree import *
from octree import *


def getparam(fname):
	# fname -> l*_*r*_*rho*_*.xyz	
	fname=[i.replace('_','.') for i in fname]
	fname=''.join(fname)
	lamb=float(fname[1:4])
	r=float(fname[5:8])
	rho=float(fname[11:14])
	return lamb,r,rho

def cart2cyl(xyz,r0):
    dr=.2
    dist=np.zeros((int(np.ceil(r0/dr)),2))
    dist[:,0]=np.arange(0,r0,dr)
    rtz=np.zeros((len(xyz),6))
    cnt=0
    rtz[0:len(xyz),0:3]=xyz
    for vec in xyz:
        x = vec[0]-r0
        y = vec[1]-r0
        r = sqrt(x**2+y**2)
        if x>0:
            if y>0:
                theta = atan(y/x)
            else:
                theta = 2*pi + atan(y/x)
        if x<0:
            theta = pi + atan(y/x)

        dist[int(floor(r/dr)),1] += 1.0
        rtz[cnt,3]=r
        rtz[cnt,4]=theta*r
        cnt +=1
    return rtz,dist
    
# encontrar n picos 
def pf(dist,n):
    # dist-> [r,npart]
    pk=np.zeros((1,2))
    for i in range(1,len(dist)-1):
        if (dist[i-1,1] < dist[i,1] and dist[i+1,1]< dist[i,1]):
            pk=np.append(pk,[[dist[i,0],dist[i,1]]],axis=0)
    pk=np.delete(pk,0,axis=0)# borra el primer elemento
    qsort(pk,0,len(pk)-1,1) # ordena respecto a npart
    npks=len(pk)
    if (len(pk)>=n):
        pk=pk[len(pk)-n:len(pk),:] # regresa los utlimos n elementos
    qsort(pk,0,len(pk)-1,0) # ordena respecto a r
    for i in range(len(pk)):
        if pk[i,1]==0:
            pk=np.delete(pk,i,axis=0)
    return pk,npks #el radio del pico y su npart 

# asume que las particulas estan ordenadas por su coord r
#regresa un arreglo con los indices de inicio y fin de las capas 
def lyr(r,pk,dr):
    lyrs=np.zeros((len(pk),2))
    for p in range(len(pk)):
        for i in range(len(r)-1):
            if r[i]>=pk[p]-dr and r[i-1]<pk[p]-dr:
                lyrs[p,0]=i
            if r[i]>=pk[p]+dr and r[i-1]<pk[p]+dr:
                lyrs[p,1]=i
                break
    for i in range(len(pk)):
        if lyrs[i,1]==0:
            lyrs[i,1]=len(r)-1
    return lyrs

# escribir un archivo xyz
def xyzwrite(xyz,fname,ncols):
    npart = len(xyz)
    f = open(fname,"w")
    f.write(str(npart)+"\n")
    f.write("\n")
    for i in range(npart):
        for j in range(ncols):
            f.write(str(xyz[i,j])+" ")
        f.write("\n")
    f.close()
            
def hexord(rtz,rcut): #elemento y radio de corte
    # rtz es un arreglo [r*\theta,z,0]
    # definir arbol con centro en z/2 2*pi*r
    z0=max(rtz[:,1])/2
    r2=abs(min(rtz[:,0])-max(rtz[:,0]))/2
    r0=min(rtz[:,0])+r2
    center=[r0,z0]
    hlen=[r2*1.1,z0*1.1]
    treebndry=cell(center,hlen) # inicializa la frontera del arbol
    tree=qtree(treebndry) # inicializa el arbol
    tree.insertArray(rtz[:,0:2])# crea el arbol para las posiciones
    for i in range(len(rtz)): 
        hx=0
        vec=[]
        acut=cell(rtz[i,0:2],[rcut,rcut])#define la region para encontrar vecinos
        tree.search(acut,vec)# busca vecinos en el arbol
        removePoint(rtz[i,0:2],vec) # quita el punto central
        if len(vec)>1:
            vec=np.array(vec)-np.array(rtz[i,0:2]) # traslada al punto central
            for x in range(len(vec)):
                v1=vec[x%len(vec)]/np.linalg.norm(vec[x%len(vec)])
                v2=vec[(x+1)%len(vec)]/np.linalg.norm(vec[(x+1)%len(vec)])
                c = np.dot(v1,v2)
                s=np.sqrt(abs(1-c**2)) # esto deberia ser siempre positivo pero hay fallos en la precision
                hx=hx+(c+s*1j)**6  # \sum e^{6i\theta}
            rtz[i,2] = np.absolute(hx)/len(vec) # esto es el parametero de orden para uno
        else:
            rtz[i,2] = 0

# ordena de menor a mayor
def qsort(a,lo,hi,col):
    if lo<hi:
        p=partition(a,lo,hi,col)
        qsort(a,lo,p-1,col)
        qsort(a,p+1,hi,col)
    return a

def partition(a,lo,hi,col):
    pivot=a[hi,col]
    i=lo
    for j in range(lo,hi):
        if a[j,col]<pivot:
            a[j,:],a[i,:]=swap(a[j,:],a[i,:])
            i +=1 
    a[i,:],a[hi,:]=swap(a[i,:],a[hi,:])
    return i

def swap(a,b):
     a=a+b
     b=a-b
     a=a-b
     return a,b

def writeLayers(rtz,lyrs,fname):
    for i in range(len(lyrs)):
        l=rtz[int(lyrs[i,0]):int(lyrs[i,1])+1,:]
        xyzwrite(l,'lyr'+str(i)+fname,6)

def removePoint(point,array):
    c=0
    for element in array:
        if all(point==element):
            del(array[c])
        c+=1

# recibe la configuracion y regresa la distribucion de \eta

def ncoord(xyz,rc):
    eta=np.zeros((200,2))
    eta[:,0]=np.arange(0,200,1)
    hlen=[(max(xyz[:,0])-min(xyz[:,0]))/2,(max(xyz[:,1])-min(xyz[:,1]))/2,(max(xyz[:,2])-min(xyz[:,2]))/2]
    cntr=[(max(xyz[:,0])+min(xyz[:,0]))/2,(max(xyz[:,1])+min(xyz[:,1]))/2,(max(xyz[:,2])+min(xyz[:,2]))/2]
    tree3d=octree(cell3d(cntr,hlen))
    tree3d.insertArray3d(xyz)
    nc=[]
    for part in xyz:
        nc=[]
        rcut=cell3d(part,[rc,rc,rc])
        tree3d.search3d(rcut,nc)
        eta[len(nc)-1,1]+=1
    return eta 
        


