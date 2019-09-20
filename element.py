import numpy
import math
A=None
B=None
def homogenousneumann(mesh,f,h):
    n_nodes=len(mesh.nodes)
    global A,B
    A=numpy.zeros((n_nodes,n_nodes))
    B=numpy.zeros(n_nodes);
    for triangle in mesh.triangles:
        gradients=[]
        solute=[]
        for i in range(0,3):
            solute.append(mesh.nodes[triangle[i]]+[1])
        for i in range(0,3):
            basis=[0]*3
            basis[i]=1
            gradients.append(numpy.linalg.solve(numpy.matrix(solute),numpy.array(basis))[:2])
        for i in range(0,3):
            for j in range(i,3):
                A[triangle[i],triangle[j]]+=numpy.dot(gradients[i],gradients[j]);
                A[triangle[j],triangle[i]]+=A[triangle[i],triangle[j]]
        area=abs(numpy.cross(numpy.subtract(solute[1][:2],solute[0][:2]),numpy.subtract(solute[2][:2],solute[0][:2])))/2
        centroid=[(solute[0][0]+solute[1][0]+solute[2][0])/3,(solute[0][1]+solute[1][1]+solute[2][1])/3]
        for vertex in triangle:
            B[vertex]+=area*f(centroid[0],centroid[1])/3
        for i in range(0,3):
            for j in range(i+1,3):
                if (h(triangle[i])!=0) and h(triangle[j]!=0):
                    B[triangle[i]]+=h(triangle[i])*numpy.linalg.norm(numpy.subtract(mesh.nodes[triangle[i]],mesh.nodes[triangle[j]]))
                    B[triangle[j]]+=h(triangle[j])*numpy.linalg.norm(numpy.subtract(mesh.nodes[triangle[i]],mesh.nodes[triangle[j]]))
    return numpy.linalg.solve(A,B)
    
class mesh:
    nodes=[]
    triangles=[]

test_mesh=mesh()
L1=20
L2=20
for i in range(0,L1):
    for j in range(0,L2):
        test_mesh.nodes.append([i,j])
        if (i!=0) and (j!=0):
            v00=L2*(i-1)+j-1
            v01=L2*(i-1)+j
            v11=L2*i+j
            v10=L2*i+j-1
            test_mesh.triangles+=[[v00,v01,v11],[v00,v11,v10]]
def test_h(i):
    if i<L2:
        return (i-L2/2)*0.00001
    else:
        return 0
def test_f(x,y):
    return 0
    if x<y:
        return 1
    else:
        return -1
    #return (x-L1/2)*1+(y-L2/2)*4
    
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

x=arange(5,L1,1)
y=arange(1,L2-1,1)
X,Y=meshgrid(x,y)
#Z=test_f(X,Y)
#im=imshow(Z)

solve=homogenousneumann(test_mesh,test_f,test_h)
def solve_f(x,y):
    return (solve[x*L2+y])
Z2=solve_f(X,Y)
minz=numpy.min(Z2)
if minz<0:
    minz=abs(minz)+0.001
else:
    minz=0.001
print(Z2)
for i in range(0,len(Z2)):
    for j in range(0,len(Z2[0])):
        #if Z2[i][j]>0:
        Z2[i][j]=Z2[i][j]#math.log(minz+(Z2[i][j]))
im2=imshow(Z2)
show()