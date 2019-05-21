import numpy as np
import matplotlib.pyplot as plt
import sys

# todo esto es en 2D tal vez algun dia lo haga para dimensiones superiores pero no hoy
class cell:
    # constructor o algo asi
    # c = cell(arr,arr)
    def __init__(self,center,halflength):
        self.cntr=center
        self.hlen=halflength

    def containsPoint(self,point):
        res=1
        for i in range(2):
            res *= (abs(point[i]-self.cntr[i])<=self.hlen[i])
        if res == 0:
            return False
        else:
            return True
    
    def intersects(self,other):
        # se puede hacer tmbn para una region circular pero que flojera
        res=1
        for i in range(2):
            res *= (abs(self.cntr[i]-other.cntr[i])<(self.hlen[i]+other.hlen[i]))
        if res == 0:
            return False
        else:
            return True

class qtree:
    # cap es un entero igual para todos
    cap=4
    # boundary es tipo cell
    # children es un array de qtrees
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
    def subdivide(self):
       # crear nuevos centros
       x=self.bndry.hlen[0]/2.    # 0 | 2
       y=self.bndry.hlen[1]/2.    # -----
       cx=self.bndry.cntr[0]      # 1 | 3
       cy=self.bndry.cntr[1]
       for i in range(1,3):
           for j in range(2):
               self.node.append(qtree(cell([cx+x*(-1)**i,cy+y*(-1)**j],[x,y])))
    
    def passp(self):
        if (len(self.node)>0 and len(self.pts)>0):
            for p in self.pts:
                for i in range(4):
                    if self.node[i].bndry.containsPoint(p):
                        self.node[i].pts.append(p)
                        break
            self.pts=[]

   # esta es la func recursiva 
    def insertp(self,pt):
        if not(self.bndry.containsPoint(pt)):
            return False
        if(len(self.pts) <self.cap and len(self.node)==0):
            self.pts.append(pt)
            return True
        if(len(self.node)==0):
            self.subdivide()
            self.passp()
        for i in range(4):
            if (self.node[i].insertp(pt)):
                return True
        return False

    def search(self,rcut,points):
        # rcut es un tipo cell
        #points=[]
        if not(self.bndry.intersects(rcut)):
            return points
        #for i in range(len(self.pts)):
        #    if rcut.containsPoint(self.pts[i]):
        #        points.append(self.pts[i])
        #if len(self.node)==0:
        #    return points
        if len(self.node)==0:
            for i in range(len(self.pts)):
                if rcut.containsPoint(self.pts[i]):
                    points.append(self.pts[i])
        else:
            for i in range(4):
                self.node[i].search(rcut,points)
        return points

    def plotree(self):
        rect(self.bndry.cntr,self.bndry.hlen)
        if len(self.node)>0:
            for i in range(4):
                self.node[i].plotree()

    def insertArray(self,array):
        for i in array:
            self.insertp(i)


def rect(c,l,fmt='k-'):
    cx,cy=c[0],c[1]
    lx,ly=l[0],l[1]
    e1=[cx-lx,cy+ly]
    e2=[cx-lx,cy-ly]
    e3=[cx+lx,cy-ly]
    e4=[cx+lx,cy+ly]
    plt.plot((e1[0],e2[0]),(e1[1],e2[1]),fmt,lw=.5)
    plt.plot((e2[0],e3[0]),(e2[1],e3[1]),fmt,lw=.5)
    plt.plot((e3[0],e4[0]),(e3[1],e4[1]),fmt,lw=.5)
    plt.plot((e4[0],e1[0]),(e4[1],e1[1]),fmt,lw=.5)

def plotp(arr,fmt='ko',size=1):
    for e in arr:
        plt.plot(e[0],e[1],fmt,markersize=size)
#
## esto genera puntos aleatorios los pone en un qtree y lo muestra
##    python3 qtree.py num_puntos
#points= np.random.normal([.3,.3],.3,[int(sys.argv[1]),2])
#tree=qtree(cell([0,0],[1,1]))
#tree.insertArray(points)
#rcut=cell([-.1,.1],[.3,.3])
#
#nb=[]
#tree.search(rcut,nb)
#
#tree.plotree()
#rect(rcut.cntr,rcut.hlen,'r-')
#plotp(points)
#plotp(nb)
#plotp(nb,'rx',2)
#plt.show()
