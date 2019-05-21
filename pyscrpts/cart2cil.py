# python cart2cil.py leer.dat escribir.dat r0
#  para hacerlo mas general calculamos el radio r0 a aprtir de parametros de la simulacion
#  parametros que si tenemos
# python cart2cil.py leer.dat escribir.dat lambda rcyl
from math import atan,sqrt,pi,floor
import sys

# archivos a escribir
f  = open(str(sys.argv[1]),"r") 
#f2 = open(str(sys.argv[2]),"w")
dfile ="dist"+ str(sys.argv[1])
f3 = open(dfile,"w")
# variables para calcular r0
#lamb = float(sys.argv[3])
#rcyl = float(sys.argv[4])
# lx.xxrx.xxrhox.xx.xyz
lamb = float(str(sys.argv[1])[1:4])
rcyl = float(str(sys.argv[1])[6:9])
# calculamos r0
#r0 = 2*(rcyl+1.5*lamb)
r0 = (rcyl+1.5*lamb)
#para la distribucion radial dividimos r0 en 500
n=100.0
dr= r0/n
#definimos un arreglo para la distribucion
dist=[0.0]*int(n)
#print(rcyl)

#saltamos dos lineas por el formato del archivo .xyz
# en la primera linea esta el numero de part
line = f.readline()
w = line.split()
npart= w[0]
#f2.write("%s\n"%npart)
line = f.readline()

# comenzamos a leer las posiciones
while line:
    line = f.readline()
    if line=="":
        break
    w = line.split() 
    x = float(w[1])-r0
    y = float(w[2])-r0
    z = float(w[3])
    r = sqrt(x**2+y**2)
    if x>0:
        if y>0:
            theta = atan(y/x)
        else:
            theta = 2*pi + atan(y/x)
    if x<0:
        theta = pi + atan(y/x)
#    l = theta*r
#    f2.write("%s %s\n"%(x,y))
# r = dr*i la particula esta en algun intervalo i
    dist[int(floor(r/dr))]=dist[int(floor(r/dr))]+1.0
    #if .8<r<1.6:    # esto selecciona la capa
     #   f2.write("%s %s %s\n"%(x,y,z))

# escribimos al archivo de distribuciones
for j in range(int(n)):
    xx=dr*(j+1)
    f3.write("%s %s\n"%(xx,dist[j]))

#cerramos archivos
f.close()
#f2.close()
f3.close()
