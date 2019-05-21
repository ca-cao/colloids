import numpy as np
#import matplotlib.pyplot as plt
import sys

def hexord(ps,p,r,dr): #posiciones, elemento y distancias para vecinos
    #definir vecinos cercanos
    # hacer una lista de vecinos cercanos vec
    vec=[ps[0]]
    for k in ps:
        if any(k!=p):
            if r < np.linalg.norm(p-k) < r*dr:
                vec=np.append(vec,[k],axis=0) # agrega elementos
    vec=np.delete(vec,0,axis=0)
    vec=np.add(vec,np.negative(p)) # traslacion al punto elegido
    hx=0+0j
    #print vec
    if len(vec)==0:
        return 0
    else:
        v0=vec[0]/(np.linalg.norm(vec[0]) # v0 es unitario
        vec=np.delete(vec,0,axis=0) # para no hacer dot consigo mismo
        if len(vec)>0:
            for x in vec:
                c=np.dot(v0,x)/np.linalg.norm(x))
                s=np.sqrt(1-c**2)
                hx=hx+(c+s*1j)**6
            return np.absolute(hx)/len(vec) # esto es el parametero de orden para uno
        else:
            return 0
#dmin=10
#for x in range(len(pos)):
#    for i in range(len(pos)-x-1):
#        if 0<np.linalg.norm(pos[x,:]-pos[i+1])<dmin:
#            dmin=np.linalg.norm(pos[x,:]-pos[i+1])
#print dmin


def cart2cyl(pos,fname):
    # lyrxlx.xxrx.xxrhox.xx.xyz
    lamb = float(fname[5:9])
    rcyl = float(fname[11:14])
    cartpos=np.zeros((len(pos),3))
    # calculamos r0
    r0 = (rcyl+1.5*lamb)
    rbar=0
    for i in pos: # radio promedio
        rbar=rbar+np.linalg.norm(pos[0:2]-r0)
    rbar=rbar/len(pos)
    for i in range(len(pos)):
        cartpos[i,0]=rbar*np.arctan2(pos[i,1]-r0,pos[i,0]-r0)
        cartpos[i,1]=pos[i,2]
    return cartpos
    
        
pos=np.loadtxt(str(sys.argv[1]),skiprows=2)


#cartpos=cart2cyl(pos,str(sys.argv[1]))
hexo=np.zeros((len(pos)))
r=1
dr=float(sys.argv[5:9])
idx=0
f=open("hex"+str(sys.argv[1]),"w")
for p in pos:
    hexo[idx]=hexord(cartpos,p,r,dr)
    idx=idx+1
    # al final de esto hay [l,z,hex]
f=open("hex"+str(sys.argv[1]),"w")
f.write("%s\n"%(len(cartpos)))
f.write("\n")

for j in range(len(pos)):
    f.write("%s %s %s\n"%(pos[j,0],pos[j,1],hexo[j]))
f.close()


#with open("hex"+str(sys.argv[1]),"w") as f:
#    for i in range(len(pos)):
#        print pos[i,0],pos[i,1],pos[i,2],t[i]

    
#a=plt.scatter(pos[:,0], pos[:,1], c=t,cmap='cool')
#plt.colorbar(a)
#plt.clim(0,1)
#plt.show()


# hacer capas mas delgadas
