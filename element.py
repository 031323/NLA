import numpy
import math
A=None
B=None
def homogenousdirichlet(mesh,boundary,f,fixed,g=lambda x:0,h=lambda x:0):
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
        area=abs(numpy.cross(numpy.subtract(solute[1][:2],solute[0][:2]),numpy.subtract(solute[2][:2],solute[0][:2])))/2
        for i in range(0,3):
            for j in range(i,3):
                if fixed(triangle[i]) or fixed(triangle[j]):
                    continue
                addition=numpy.dot(gradients[i],gradients[j])*area;
                A[node_map[triangle[i]],node_map[triangle[j]]]+=addition
                if i!=j:
                    A[node_map[triangle[j]],node_map[triangle[i]]]+=addition
                
                if boundary(triangle[i]) and boundary(triangle[j]):
                    dL=numpy.linalg.norm(numpy.subtract(solute[i][:2],solute[j][:2]))
                    addition=(h(triangle[i])+h(triangle[j]))*dL/2
                    B[node_map[triangle[i]]]+=addition
                    B[node_map[triangle[j]]]+=addition
                
        area=abs(numpy.cross(numpy.subtract(solute[1][:2],solute[0][:2]),numpy.subtract(solute[2][:2],solute[0][:2])))/2
        centroid=[(solute[0][0]+solute[1][0]+solute[2][0])/3,(solute[0][1]+solute[1][1]+solute[2][1])/3]
        #print("centroid:",centroid)
        #print("solute:",solute)
        for vertex in triangle:
            if fixed(vertex):
                continue
            B[node_map[vertex]]+=area*f(centroid[0],centroid[1])/3
        G_gradient=numpy.zeros(2)
        for i in range(0,3):
            if fixed(triangle[i]):
                G_gradient=numpy.add(G_gradient,g(triangle[i])*gradients[i])
        for i in range(0,3):
            if not fixed(triangle[i]):
                B[node_map[triangle[i]]]-=area*numpy.dot(gradients[i],G_gradient)
    #print("A:\n",A,"\nB:\n",B)
    return numpy.linalg.solve(A,B)
    
class mesh:
    nodes=[]
    triangles=[]

test_mesh=mesh()
L1=20
L2=20
for i in range(0,L1):
    for j in range(0,L2):
        test_mesh.nodes.append([i/(L1-1),j/(L2-1)])
        if (i!=0) and (j!=0):
            v00=L2*(i-1)+j-1
            v01=L2*(i-1)+j
            v11=L2*i+j
            v10=L2*i+j-1
            test_mesh.triangles+=[[v00,v01,v11],[v00,v11,v10]]
def test_boundary(i):
    return i<L2 or i>L2*(L1-1) or i%L2==0 or (i+1)%L2==0
def test_fixed(i):
    #return 0
    if i>L2*(L1-1) or i%L2==0 or (i+1)%L2==0:
        return 1
    else:
        return 0
def test_g(i):
    if i<L2:
        return i/L2
    elif (i+1)%L2==0:
        #return (1-i/L2/L1)
        return (L2-1)/L2
    elif i>L2*(L1-1):
        #return (1-(i-L2*(L1-1))/L2)
        return (i-L2*(L1-1))/L2
    elif i%L2==0:
        #return i/L2/L1
        return 0
def test_f(x,y):
    #return 10
    if x>1/2:
        return 20
    else:
        return -1
    #return (x-L1/2)*1+(y-L2/2)*4
    if x>L1/2:
        return 1
    else:
        return 0
def test_h(i):
    return 1
    return (L2*L1/2-i)*2/L2/L1

solve=homogenousdirichlet(test_mesh,test_boundary,test_f,test_fixed,test_g,test_h)

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

bafar=0
b = numpy.arange(bafar, L2-bafar, 1)
d = numpy.arange(bafar, L1-bafar, 1)
nu = numpy.zeros( (d.size, b.size) )
counter= 0

for i in range(bafar,L1-bafar):
    for j in range(bafar,L2-bafar):
        if not test_fixed(i*L2+j):
            nu[i-bafar][j-bafar]=solve[counter]
            counter+=1
        else:
            nu[i-bafar][j-bafar]=test_g(i*L2+j)

X, Y = numpy.meshgrid(d, b)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('persp')
ax.plot_surface(Y, X, nu,cmap=cm.RdBu)
plt.show()

