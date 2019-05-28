import numpy as np
import matplotlib.pyplot as plt
import sys

# todo esto es en 2D tal vez algun dia lo haga para dimensiones superiores pero no hoy
# hoy es ese dia 

class cell3d:
    # constructor o algo asi
    # c = cell3d(arr,arr)
    def __init__(self,center,halflength):
        self.cntr=center
        self.hlen=halflength

    def containsPoint3d(self,point):
        res=1
        for i in range(3):
            res *= (abs(point[i]-self.cntr[i])<=self.hlen[i])
        if res == 0:
            return False
        else:
            return True
    
    def intersects3d(self,other):
        # se puede hacer tmbn para una region circular pero que flojera
        res=1
        for i in range(3):
            res *= (abs(self.cntr[i]-other.cntr[i])<(self.hlen[i]+other.hlen[i]))
        if res == 0:
            return False
        else:
            return True

class octree:
    # cap es un entero igual para todos
    cap=4
    # boundary es tipo cell3d
    # children es un array de octrees
    # points es un array
    def __init__(self,boundary,children=None,points=None,full=None,divided=None):
        if children is None:
            self.node=[]
        else:
            self.node=children
        if points is None:
            self.pts=[]
        else:
            self.pts=points
        self.bndry=boundary

    # conform or be cast out
    def subdivide3d(self):
       # crear nuevos centros
       l2 = [ i/2. for i in self.bndry.hlen]
       c = self.bndry.cntr
       for i in range(2):
           for j in range(2):
               for k in range(2):
                    self.node.append(octree(cell3d([c[0]+l2[0]*(-1)**i,c[1]+l2[1]*(-1)**j,c[2]+l2[2]*(-1)**k],l2)))
    
    def passp3d(self):
        if (len(self.node)>0 and len(self.pts)>0):
            for p in self.pts:
                for i in range(8):
                    if self.node[i].bndry.containsPoint3d(p):
                        self.node[i].pts.append(p)
                        break
            self.pts=[]

   # esta es la func recursiva 
    def insertp3d(self,pt):
        if not(self.bndry.containsPoint3d(pt)):
            return False
        if(len(self.pts) <self.cap and len(self.node)==0):
            self.pts.append(pt)
            return True
        if(len(self.node)==0):
            self.subdivide3d()
            self.passp3d()
        for i in range(8):
            if (self.node[i].insertp3d(pt)):
                return True
        return False

    def search3d(self,rcut,points):
        # rcut es un tipo cell3d tiene que ser un cubo
        rc=rcut.hlen[0]
        if not(self.bndry.intersects3d(rcut)):
            return points 
        if len(self.node)==0:
            for i in range(len(self.pts)):
                if (rcut.containsPoint3d(self.pts[i]) and np.linalg.norm(self.pts[i]-rcut.cntr)<rc):
                    points.append(self.pts[i])
        else:
            for i in range(8):
                self.node[i].search3d(rcut,points) 
        return points


    def insertArray3d(self,array):
        for i in array:
            self.insertp3d(i)

#    def plotree(self):
#        rect(self.bndry.cntr,self.bndry.hlen)
#        if len(self.node)>0:
#            for i in range(8):
#                self.node[i].plotree()
#def rect(c,l,fmt='k-'):
#    cx,cy=c[0],c[1]
#    lx,ly=l[0],l[1]
#    e1=[cx-lx,cy+ly]
#    e2=[cx-lx,cy-ly]
#    e3=[cx+lx,cy-ly]
#    e4=[cx+lx,cy+ly]
#    plt.plot((e1[0],e2[0]),(e1[1],e2[1]),fmt,lw=.5)
#    plt.plot((e2[0],e3[0]),(e2[1],e3[1]),fmt,lw=.5)
#    plt.plot((e3[0],e4[0]),(e3[1],e4[1]),fmt,lw=.5)
#    plt.plot((e4[0],e1[0]),(e4[1],e1[1]),fmt,lw=.5)
#
#def plotp(arr,fmt='ko',size=1):
#    for e in arr:
#        plt.plot(e[0],e[1],fmt,markersize=size)
#
#np.random.seed(1)
#pts=np.zeros((64,3))
#c=0
#for i in range(4):
#    for j in range(4):
#        for k in range(4):
#            pts[c]=[.25*i,.25*k,.25*j]
#            c+=1
#
#for i in range(len(pts)):
#    print("{} {} {}".format(pts[i,0],pts[i,1],pts[i,2]))
#c0=cell3d([.5,.5,.5],[.5,.5,.5])
#tree3d=octree(c0)
#tree3d.insertArray(pts)
#ccut=cell3d(pts[22],[1,1,1])
#vec=[]
#tree3d.search3d(ccut,vec,.26)
#print(len(vec))
#print("")
#for i in range(len(vec)):
#    print("{} {} {}".format(vec[i][0],vec[i][1],vec[i][2]))
#
#
#