import numpy
import math
A=None
B=None
def pxy(degree,c,x,y):
    sum=0
    i=0
    for j in range(0,degree+1):
        for k in range(0,j+1):
            sum+=c[i]*pow(x,j-k)*pow(y,k)
            i+=1
    return sum
def dx_pxy(degree,c,x,y):
    sum=0
    for j in range(0,degree+1):
        for k in range(0,j+1):
            sum+=c[i]*(j-k)*pow(x,j-k-1)*pow(y,k)
            i+=1
    return sum

def dy_pxy(degree,c,x,y):
    sum=0
    for j in range(0,degree+1):
        for k in range(0,j+1):
            sum+=c[i]*(k)*pow(x,j-k)*pow(y,k-1)
            i+=1
    return sum

def homogenousdirichlet(mesh,degree,boundary,f,fixed,g=lambda x:0,h=lambda x:0):
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
        M=[]
        for node in triangle:
            [x,y]=mesh.nodes[node]
            Mi=[]
            for j in range(0,degree+1):
                for k in range(0,j+1):
                    Mi.append(pow(x,j-k)*pow(y,k))
            M.append(Mi)
    #print("A:\n",A,"\nB:\n",B)
    if n_free_nodes==len(mesh.nodes):
        print("lstsq")
        return numpy.linalg.lstsq(A,B)[0]
    else:
        print("solve")
        return numpy.linalg.solve(A,B)
    #A=numpy.delete(A,0,0)
    #A=numpy.delete(A,0,1)
    #B=numpy.delete(B,0,0)
    #x=numpy.linalg.solve(A,B)
    #return numpy.concatenate((numpy.zeros(1),x))

class mesh:
    nodes=[]
    triangles=[]

test_mesh=mesh()
degree=2
l1=1
l2=1
L1=1*degree*2+1
L2=1*degree*2+1
n_elements=0
for i in range(0,L1):
    for j in range(0,L2):
        test_mesh.nodes.append([i/(L1-1)*l1,j/(L2-1)*l2])
        if (i!=0) and (j!=0) and (i%degree==0) and (j%degree==0):
            n_elements+=1
            t1=[]
            t2=[]
            for i_ in range(i-degree,i+1):
                for j_ in range(j-degree,j+1):
                    node=L2*i_+j_
                    if (n_elements+int((n_elements-1)/(L2-1)))%2==0:
                        if (i-i_+j-j_>=degree):
                            t1+=[node]
                        if (i-i_+j-j_<=degree):
                            t2+=[node]
                    else:
                        if (i-i_<=j-j_):
                            t1+=[node]
                        if (i-i_>=j-j_):
                            t2+=[node]
            mesh.triangles+=[t1,t2]
def test_boundary(i):
    return i<L2 or i>=L2*(L1-1) or i%L2==0 or (i+1)%L2==0
def test_fixed(i):
    return 0
    if i<L2:return 1
    else:return 0
    #if i==L2/2 or i==L2/2-1 :return 1
    #else:return 0
    return 0
    if test_boundary(i):
        if i<L2 or i%L2==0:return 1
    else:return 0
    if i==L2-1:
        return 1
    else:
        return 0
    if i>L2*(L1-1) or i%L2==0:
        return 1
    else:
        return 0
def test_g(i):
    return 0#i/(L2-1)
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
    return 4
    #return -2*x*x-2*y*y
    if x>y:
        return 1
    else:
        return -1
    #return (x-L1/2)*1+(y-L2/2)*4
    if x>L1/2:
        return 1
    else:
        return 0
def test_xy(i):
    return test_mesh.nodes[i][0],test_mesh.nodes[i][1]
def test_h(i,towards):
    return -1
    x,y=test_xy(i)
    #if (x==0 or x==1) and (y==0 or y==1):
    #    return 0
    #else:return -1
    #return x*(1-x)+y*(1-y)
    if i<L2:
        return 0
    elif (i+1)%L2==0:
        return 2*x*x+2
    elif i%L2==0:
        return -2
    elif i>=L2*(L1-1):
        return 2*y*y
    else:raise ValueError
    multiplier=10
    if i>0 and i<L2-1:
        return -1
    elif i==0 and towards!=L2:
        return -1
    elif i==L2-1 and towards==L2-2:
        return -1
    else:return 0
        #return -i/L2*(L2-1-i)/L2*multiplier
    if (i+1)%L2==0:
        return -(L1-(i+1)/L2)/L1*(i/L2)/L1*multiplier
    else:return 0
    if i>L2*(L1-1):
        return -(i-L2*(L1-1))/L1*(L2*L1-i)/L1*multiplier
    if i%L2==0:
        return -i/L2/L1*(L1-i/L2)/L1*multiplier
    
    if i<L2*L1/2:
        return 1
    else:
        return -1
    return i*(L2-1-i)/L2/L2*10

solve=homogenousdirichlet(test_mesh,degree,test_boundary,test_f,test_fixed,test_g,test_h)

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

bafar=0
b = numpy.arange(0, l2+l2/(L2), l2/(L2-1))
d = numpy.arange(0, l1+l1/(L1), l1/(L1-1))

nu = numpy.zeros( (d.size, b.size) )
counter= 0


for i in range(bafar,L1-bafar):
    for j in range(bafar,L2-bafar):
        if not test_fixed(i*L2+j):
            nu[i-bafar][j-bafar]=solve[counter]
            counter+=1
        else:
            nu[i-bafar][j-bafar]=test_g(i*L2+j)
        
off=nu[0][0]
for i in range(0,L1-0):
    for j in range(0,L2-0):
        x,y=test_xy(i*L2+j)
        #nu[i][j]-=off+x*x*y*y+2*y

X, Y = numpy.meshgrid(b, d)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('persp')
ax.plot_surface(Y, X, nu,cmap=cm.RdBu)
plt.show()

def XC(x1,x2):
    for i in range(len(x1)):
        for j in range(len(x1[0])):
            if abs(x1[i][j]-x2[i*2][j*2])<0.0001:
                print(i,j)

