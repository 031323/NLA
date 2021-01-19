import numpy
import math
from matplotlib import cm
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
	i=0
	for j in range(0,degree+1):
		for k in range(0,j+1):
			sum+=c[i]*(j-k)*pow(x,j-k-1)*pow(y,k)
			i+=1
	return sum

def dy_pxy(degree,c,x,y):
	sum=0
	i=0
	for j in range(0,degree+1):
		for k in range(0,j+1):
			sum+=c[i]*(k)*pow(x,j-k)*pow(y,k-1)
			i+=1
	return sum

gq1d=[

[[0,2]],

[[0.57735,1],
[-0.57735,1]],

[[0,0.888888],
[0.774597,0.555555],
[-0.774597,0.555555]]

]

gq2d=[

[[[0.333333,0.333333,0.333333],1]],

[[[0.666666666666666,0.166666666666666,0.166666666666666],0.333333333333333],
[[0.166666666666666,0.666666666666666,0.166666666666666],0.333333333333333],
[[0.166666666666666,0.166666666666666,0.666666666666666],0.333333333333333]],

[],

[[[0.108103018168070,0.445948490915965,0.445948490915965],0.223381589678011],
[[0.445948490915965,0.108103018168070,0.445948490915965],0.223381589678011],
[[0.445948490915965,0.445948490915965,0.108103018168070],0.223381589678011],
[[0.816847572980459,0.091576213509771,0.091576213509771],0.109951743655322],
[[0.091576213509771,0.816847572980459,0.091576213509771],0.109951743655322],
[[0.091576213509771,0.091576213509771,0.816847572980459],0.109951743655322]
],

[],

[[[0.501426509658179,0.249286745170910,0.249286745170910],0.116786275726379],
[[0.249286745170910,0.501426509658179,0.249286745170910],0.116786275726379],
[[0.249286745170910,0.249286745170910,0.501426509658179],0.116786275726379],
[[0.873821971016996,0.063089014491502,0.063089014491502],0.050844906370207],
[[0.063089014491502,0.873821971016996,0.063089014491502],0.050844906370207],
[[0.063089014491502,0.063089014491502,0.873821971016996],0.050844906370207],
[[0.053145049844817,0.310352451033784,0.636502499121399],0.082851075618374],
[[0.053145049844817,0.636502499121399,0.310352451033784],0.082851075618374],
[[0.310352451033784,0.053145049844817,0.636502499121399],0.082851075618374],
[[0.636502499121399,0.053145049844817,0.310352451033784],0.082851075618374],
[[0.310352451033784,0.636502499121399,0.053145049844817],0.082851075618374],
[[0.636502499121399,0.310352451033784,0.053145049844817],0.082851075618374]
]

]

error=float('inf')
def solve(mesh,degree,boundary,f,fixed,element_type,g=lambda x:0,h=lambda x:0,P1_skip=[]):
	node_map=[]
	#global n_free_nodes
	n_free_nodes=0
	for i in range(0,len(mesh.nodes)):
		if fixed(i) or (element_type=='P1' and P1_skip[i]):
			node_map.append(-1)
		else:
			node_map.append(n_free_nodes)
			n_free_nodes+=1
	global A,B
	A=numpy.zeros((n_free_nodes,n_free_nodes))
	B=numpy.zeros(n_free_nodes);
	Ms=[]
	for tn,triangle in enumerate(mesh.triangles):
		v=[mesh.nodes[mesh.triangle_vertices[tn][0]],mesh.nodes[mesh.triangle_vertices[tn][1]],mesh.nodes[mesh.triangle_vertices[tn][2]]]
		area=abs(numpy.cross(numpy.subtract(v[1],v[0]),numpy.subtract(v[2],v[0])))*0.5
		M=[]
		#print(triangle)
		for node in triangle:
			[x,y]=mesh.nodes[node]
			Mi=[]
			for j in range(0,degree+1):
				for k in range(0,j+1):
					if mesh.node_type[node]==' ':
						Mi.append(pow(x,j-k)*pow(y,k))
					elif mesh.node_type[node]=='x':
						Mi.append((j-k)*pow(x,max(j-k-1,0))*pow(y,k))
					elif mesh.node_type[node]=='y':
						Mi.append((k)*pow(x,j-k)*pow(y,max(k-1,0)))
					'''
					elif mesh.node_type[node]=='l' or mesh.node_type[node]=='r':
						d=0
						if mesh.node_type[node]=='l':d=1
						else:d=-1
						for i in range(0,3):
							if [x,y]==v[i]:v_node=i;break;
						node2=v[(i+d)%3]
						l=math.sqrt(pow(node2[0]-x,2)+pow(node2[1]-y,2))
						Dx=(node2[0]-x)/l
						Dy=(node2[1]-y)/l
						Mi.append((j-k)*pow(x,max(j-k-1,0))*pow(y,k)*Dx+Dy*(k)*pow(x,j-k)*pow(y,max(k-1,0)))
					'''
			#print(M)
			M.append(Mi)
		ce=numpy.linalg.inv(M)
		Ms.append(ce)
		
		#  a(ei,ej) in omega
		
		for i in range(0,len(triangle)):
			for j in range(i,len(triangle)):
				if fixed(triangle[i]) or fixed(triangle[j]):
					continue
				integration=0
				p=2*degree-2
				quadratures=gq2d[p-1]
				for k in range(0,len(quadratures)):
					x=0
					y=0
					for l in range(0,3):
						x+=v[l][0]*quadratures[k][0][l]
						y+=v[l][1]*quadratures[k][0][l]
					integration+=area*quadratures[k][1]*\
					(dx_pxy(degree,ce[:,i],x,y)*dx_pxy(degree,ce[:,j],x,y)+dy_pxy(degree,ce[:,i],x,y)*dy_pxy(degree,ce[:,j],x,y))
					
				A[node_map[triangle[i]],node_map[triangle[j]]]+=integration
				if i!=j:
					A[node_map[triangle[j]],node_map[triangle[i]]]+=integration
		
		for i in range(0,len(triangle)):
		
		# f(ei) in omega
		
			if fixed(triangle[i]):continue
			integration=0
			p=2*degree
			quadratures=gq2d[p-1]
			for k in range(0,len(quadratures)):
				x=0
				y=0
				for l in range(0,3):
					x+=v[l][0]*quadratures[k][0][l]
					y+=v[l][1]*quadratures[k][0][l]
				integration+=area*quadratures[k][1]*pxy(degree,ce[:,i],x,y)*f(x,y)
			
			#if triangle[i]==2 or triangle[i]==22:print([triangle,triangle[i],integration])
			
			# a(g,ei) in omega
			
			p=2*degree-2
			quadratures=gq2d[p-1]
			for j in range(0,len(triangle)):
				if not fixed(triangle[j]):continue
				for k in range(0,len(quadratures)):
					x=0
					y=0
					for l in range(0,3):
						x+=v[l][0]*quadratures[k][0][l]
						y+=v[l][1]*quadratures[k][0][l]
					integration-={' ':g(triangle[j]), 'x':0, 'y':0}[mesh.node_type[triangle[j]]]*area*quadratures[k][1]*\
					(dx_pxy(degree,ce[:,i],x,y)*dx_pxy(degree,ce[:,j],x,y)+dy_pxy(degree,ce[:,i],x,y)*dy_pxy(degree,ce[:,j],x,y))
			
			#if triangle[i]==2 or triangle[i]==22:print([triangle,triangle[i],integration])
			
			# h*ei on (boundary of omega - Gamma). Note that g0 is zero on this part.
			
			if(boundary(triangle[i])):
				for a in range(0,3):
					for b in range(a,3):
						if a==b or (not boundary(mesh.triangle_vertices[tn][a])) or (not boundary(mesh.triangle_vertices[tn][b])):continue
						if  fixed(mesh.triangle_vertices[tn][a]) and fixed(mesh.triangle_vertices[tn][b]):continue
						length=math.sqrt(pow(v[a][0]-v[b][0],2)+pow(v[a][1]-v[b][1],2))
						n=math.ceil((degree+1)/2)
						quadratures=gq1d[n-1]
						for k in range(0,len(quadratures)):
							x=(v[a][0]*(1-quadratures[k][0])+v[b][0]*(1+quadratures[k][0]))/2
							y=(v[a][1]*(1-quadratures[k][0])+v[b][1]*(1+quadratures[k][0]))/2
							integration+=length*quadratures[k][1]/2*pxy(degree,ce[:,i],x,y)*h(x,y)
			
			#if triangle[i]==2 or triangle[i]==22:print([triangle,triangle[i],integration])
			B[node_map[triangle[i]]]+=integration
					
	#print("A:\n",A,"\nB:\n",B)
	if n_free_nodes==len(mesh.nodes):
		print("lstsq")
		node_values=numpy.linalg.lstsq(A,B)[0]
	else:
		print("solve")
		node_values=numpy.linalg.solve(A,B)
	
	#Eriksson-Johnson error indicator
	def coeff(i):
		node_z=[0,0,0]
		for j in range(0,3):
			if fixed(mesh.triangles[i][j]):
				node_z[j]=g(mesh.triangles[i][j])
			else:
				node_z[j]=node_values[node_map[mesh.triangles[i][j]]]
		return numpy.matmul(Ms[i],numpy.transpose([node_z]))
	if degree==1 and element_type=='lagrange':
		maxerr=0
		for n in range(0,len(mesh.nodes)):
			if fixed(n):continue
			maxerr=max(maxerr,abs(2*mesh.nodes[n][0]*(1-mesh.nodes[n][0])-node_values[node_map[n]]))
		print('Linf '+str(maxerr))
		triangle_error=numpy.full(len(mesh.triangles),numpy.inf)
		du_x=numpy.full(len(mesh.triangles),0.0)
		du_y=numpy.full(len(mesh.triangles),0.0)
		c=numpy.full((len(mesh.triangles),2),0.0)
		for tn,t in enumerate(mesh.triangles):
			for j in range(0,3):
				c[tn]+=numpy.array(mesh.nodes[mesh.triangle_vertices[tn][j]])
			c[tn]/=3
			coeffs=coeff(tn)
			du_x[tn]=dx_pxy(degree,coeffs,c[tn][0],c[tn][1])
			du_y[tn]=dy_pxy(degree,coeffs,c[tn][0],c[tn][1])
		maxerr=0
		maxtn=[0,0]
		for tn1 in range(0,len(mesh.triangles)):
			for tn2 in range(tn1+1,len(mesh.triangles)):
				n_common=0
				for i in range(0,3):
					for j in range(0,3):
						if mesh.triangles[tn1][i]==mesh.triangles[tn2][j]:n_common+=1
				if n_common!=2:continue
				distance=numpy.linalg.norm(c[tn1]-c[tn2])
				if max(abs(du_x[tn1]-du_x[tn2])/distance,abs(du_y[tn1]-du_y[tn2])/distance)>maxerr:maxtn=[tn1,tn2]
				maxerr=max(maxerr,abs(du_x[tn1]-du_x[tn2])/distance,abs(du_y[tn1]-du_y[tn2])/distance)
		print('error: '+str(maxerr))
		print('maxtnc')
		print(c[maxtn[0]],c[maxtn[1]])
					
				
	
	if element_type=='P1':
		tri_points=[]
		for i in range(0,len(Ms)):
			node_z=[0,0,0]
			for j in range(0,3):
				if fixed(mesh.triangles[i][j]):
					node_z[j]=g(mesh.triangles[i][j])
				else:
					node_z[j]=node_values[node_map[mesh.triangles[i][j]]]
			coeffs=numpy.matmul(Ms[i],numpy.transpose([node_z]))
			print(coeffs)
			print ("")
			tn=i
			v=[mesh.nodes[mesh.triangle_vertices[tn][0]],mesh.nodes[mesh.triangle_vertices[tn][1]],mesh.nodes[mesh.triangle_vertices[tn][2]]]
			for j in range(0,3):
				v_z=coeffs[0]+coeffs[1]*v[j][0]+coeffs[2]*v[j][1]
				tri_points+=[(v[j][0],v[j][1],v_z[0])]
		
		from itertools import chain
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D
# Get the X, Y and Z coordinates of each point
		print(tri_points[0])
		x, y, z = zip(*tri_points)
# Make list of triangle indices ([(0, 1, 2), (3, 4, 5), ...])
		tri_idx = [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(int(len(tri_points)/3))]
# Make 3D axes
		ax = plt.figure().gca(projection='3d')
# Plot triangles
		ax.plot_trisurf(x, y, z, triangles=tri_idx,cmap=cm.RdBu)
		plt.show()

	return node_values	
	#A=numpy.delete(A,0,0)
	#A=numpy.delete(A,0,1)
	#B=numpy.delete(B,0,0)
	#x=numpy.linalg.solve(A,B)
	#return numpy.concatenate((numpy.zeros(1),x))


class mesh:
	nodes=[]
	node_type=[]
	triangles=[]
	triangle_vertices=[]

test_mesh=mesh()
element_type='lagrange'
# hermite triangle pair
# 0 1 1 0
# 1 0 1 1
# 1 1 0 1
# 0 1 1 0
degree=1
if element_type=='hermite':degree=3
if element_type=='P1':degree=1
l1=1
l2=1
n_block_per_side=10
if element_type=='lagrange' or element_type=='P1':
	if element_type=='lagrange':spe=degree # sections per edge
	else: spe=2
	L1=n_block_per_side*spe*2+1
	L2=n_block_per_side*spe*2+1
else:
	L1=n_block_per_side*1*2+1
	L2=n_block_per_side*1*2+1
	spe=1
n_elements2=0
hermite_map=[]
hermite_extra_nodes=[]
hermite_extra_node_type=[]
P1_skip=[]
for i in range(0,L1):
	for j in range(0,L2):
		test_mesh.nodes.append([i/(L1-1)*l1,j/(L2-1)*l2])
		test_mesh.node_type.append(' ')
		if element_type=='P1':
			if i%spe==0 and j%spe==0:P1_skip+=[True]
			else: P1_skip+=[False]
		if element_type=='hermite':
			hermite_map+=[len(hermite_extra_nodes)]
			hermite_extra_nodes.append([i/(L1-1)*l1,j/(L2-1)*l2])
			hermite_extra_node_type.append('x')
			hermite_extra_nodes.append([i/(L1-1)*l1,j/(L2-1)*l2])
			hermite_extra_node_type.append('y')
		if (i!=0) and (j!=0) and (i%spe==0) and (j%spe==0):
			n_elements2+=1
			F2=((n_elements2+int((n_elements2-1)*spe/(L2-1)))%2==0)
			#print([n_elements2,(n_elements2-1)*degree,L2-1,int((n_elements2-1)*degree/(L2-1)),(n_elements2+int((n_elements2-1)*degree/(L2-1))),F2])
			t1=[]
			t2=[]
			for i_ in range(i-spe,i+1):
				for j_ in range(j-spe,j+1):
					if element_type=='P1' and (i_==i-spe or i_==i ) and (j_==j-spe or j_==j):continue
					node=L2*i_+j_
					if F2:
						if (i-i_+j-j_>=spe):
							t1+=[node]
							if element_type=='hermite':
								t1+=[L1*L2+hermite_map[node],L1*L2+hermite_map[node]+1]
						if (i-i_+j-j_<=spe):
							t2+=[node]
							if element_type=='hermite':
								t2+=[L1*L2+hermite_map[node],L1*L2+hermite_map[node]+1]
					else:
						if (i-i_<=j-j_):
							t1+=[node]
							if element_type=='hermite':
								t1+=[L1*L2+hermite_map[node],L1*L2+hermite_map[node]+1]
						if (i-i_>=j-j_):
							t2+=[node]
							if element_type=='hermite':
								t2+=[L1*L2+hermite_map[node],L1*L2+hermite_map[node]+1]
							
			v00=L2*(i-spe)+j-spe
			v01=L2*(i-spe)+j
			v11=L2*i+j
			v10=L2*i+j-spe
			
			if F2:
				test_mesh.triangle_vertices+=[[v00,v01,v10],[v10,v11,v01]]
				c1=[(test_mesh.nodes[v00][0]+test_mesh.nodes[v01][0]+test_mesh.nodes[v10][0])/3,(test_mesh.nodes[v00][1]+test_mesh.nodes[v01][1]+test_mesh.nodes[v10][1])/3]
				c2=[(test_mesh.nodes[v10][0]+test_mesh.nodes[v11][0]+test_mesh.nodes[v01][0])/3,(test_mesh.nodes[v10][1]+test_mesh.nodes[v11][1]+test_mesh.nodes[v01][1])/3]
			else:
				test_mesh.triangle_vertices+=[[v00,v11,v10],[v00,v01,v11]]
				c1=[(test_mesh.nodes[v00][0]+test_mesh.nodes[v11][0]+test_mesh.nodes[v10][0])/3,(test_mesh.nodes[v00][1]+test_mesh.nodes[v11][1]+test_mesh.nodes[v10][1])/3]
				c2=[(test_mesh.nodes[v00][0]+test_mesh.nodes[v01][0]+test_mesh.nodes[v11][0])/3,(test_mesh.nodes[v00][1]+test_mesh.nodes[v01][1]+test_mesh.nodes[v11][1])/3]
				
			if element_type=='hermite':
				hermite_extra_nodes.append(c1)
				hermite_extra_node_type.append(' ')
				t1+=[L1*L2+len(hermite_extra_nodes)-1]
				
				hermite_extra_nodes.append(c2)
				hermite_extra_node_type.append(' ')
				t2+=[L1*L2+len(hermite_extra_nodes)-1]
			
			test_mesh.triangles+=[t1,t2]

test_mesh.nodes+=hermite_extra_nodes
test_mesh.node_type+=hermite_extra_node_type
def test_boundary(i):
	return test_mesh.nodes[i][0]==0 or test_mesh.nodes[i][0]==1 or test_mesh.nodes[i][1]==0 or test_mesh.nodes[i][1]==1
	#return i<L2 or i>=L2*(L1-1) or i%L2==0 or (i+1)%L2==0
def test_fixed(i):
	#return test_mesh.nodes[i][0]==0 or test_mesh.nodes[i][0]==1
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
	#return 0
	return i/(L2-1)
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
	return -2*x*x-2*y*y
	if x>y:
		return 40
	else:
		return -40
	#return (x-L1/2)*1+(y-L2/2)*4
	if x>L1/2:
		return 1
	else:
		return 0
def test_xy(i):
	return test_mesh.nodes[i][0],test_mesh.nodes[i][1]
def test_h(i,towards):
	#return 0
	return -1
	x,y=test_xy(i)
	#if (x==0 or x==1) and (y==0 or y==1):
	#	return 0
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

solve=solve(test_mesh,degree,test_boundary,test_f,test_fixed,element_type,test_g,test_h,P1_skip)


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
		if not (test_fixed(i*L2+j) or (element_type=='P1' and P1_skip[i*L2+j])):
			nu[i-bafar][j-bafar]=solve[counter]
			counter+=1
		elif element_type != 'P1':
			nu[i-bafar][j-bafar]=test_g(i*L2+j)
		else:
			nu[i-bafar][j-bafar]=0
		
off=nu[0][0]
for i in range(0,L1-0):
	for j in range(0,L2-0):
		x,y=test_xy(i*L2+j)
		#nu[i][j]-=off+x*x*y*y+2*y

X, Y = numpy.meshgrid(b, d)

if element_type!='P1':fig = plt.figure()
if element_type!='P1':ax = fig.add_subplot(111, projection='3d')
if element_type!='P1':ax.set_proj_type('persp')
if element_type!='P1':ax.plot_surface(Y, X, nu,cmap=cm.RdBu)
if element_type!='P1':plt.show()
	

def XC(x1,x2):
	for i in range(len(x1)):
		for j in range(len(x1[0])):
			if abs(x1[i][j]-x2[i*2][j*2])<0.0001:
				print(i,j)


