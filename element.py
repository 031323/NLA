import numpy
import math
A=None
B=None
def homogenousdirichlet(mesh,f,fixed):
    node_map=[]
    n_free_nodes=0
    for i in range(0,len(mesh.nodes)):
        if fixed(i):
            node_map.append(-1)
        else:
            node_map.append(n_free_nodes)
            n_free_nodes+=1
    global A,B
    A=numpy.zeros((n_free_nodes,n_free_nodes))
    B=numpy.zeros(n_free_nodes);
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
                if fixed(triangle[i]) or fixed(triangle[j]):
                    continue
                addition=numpy.dot(gradients[i],gradients[j]);
                A[node_map[triangle[i]],node_map[triangle[j]]]+=addition
                if i!=j:
                    A[node_map[triangle[j]],node_map[triangle[i]]]+=addition
        area=abs(numpy.cross(numpy.subtract(solute[1][:2],solute[0][:2]),numpy.subtract(solute[2][:2],solute[0][:2])))/2
        centroid=[(solute[0][0]+solute[1][0]+solute[2][0])/3,(solute[0][1]+solute[1][1]+solute[2][1])/3]
        for vertex in triangle:
            if fixed(vertex):
                continue
            B[node_map[vertex]]+=area*f(centroid[0],centroid[1])/3
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
    if i<L2 or i>L2*(L1-1) or i%L2==0 or (i+1)%L2==0:
        return 1
    else:
        return 0
def test_f(x,y):
    #return 1
    #if x<y:
    #    return 1
    #else:
    #    return 0
    #return (x-L1/2)*1+(y-L2/2)*4
    if x>L1/2:
        return 1
    else:
        return 0

solve=homogenousdirichlet(test_mesh,test_f,test_h)

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

bafar=0
b = numpy.arange(bafar, L2-bafar, 1)
d = numpy.arange(bafar, L1-bafar, 1)
nu = numpy.zeros( (b.size, d.size) )
counter= 0

for i in range(bafar,L1-bafar):
    for j in range(bafar,L2-bafar):
        if not test_h(i*L2+j):
            nu[i-bafar][j-bafar]=solve[counter]
            counter+=1

X, Y = numpy.meshgrid(d, b)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('persp')
ax.plot_surface(X, Y, nu,cmap=cm.RdBu)
plt.show()

