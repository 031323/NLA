import numpy
def homogenousdirichlet(mesh,f):
    A=numpy.zeros((mesh.n_nodes,mesh.n_nodes))
    B=numpy.zeros(mesh.n_nodes)
    for triangle in mesh.triangles:
        gradients=[]
        solute=[]
        for i in range(0,3):
            solute.push(mesh.nodes[triangle.vertices[i]]+[1])
        for i in range(0,3):
            basis=[0]*3
            basis[i]=1
            gradients.push(numpy.linalg.solve(numpy.matrix(solute),numpy.array(basis))[:2])
        for i in range(0,3):
            for j in range(i,3):
                A[triangles.vertices[i],triangle.vertices[j]]=numpy.dot(gradients[i].gradients[j]);
                A[triangles.vertices[j],triangle.vertices[i]]=A[triangles.vertices[i],triangle.vertices[j]]
        area=abs(numpy.cross(numpy.subtract(solute[1],solute[0]),numpy.subtract(solute[2],solute[0])))/2
        centroid=[(solute[0][0]+solute[1][0]+solute[2][0])/3,(solute[0][1]+solute[1][1]+solute[2][1])/3]
        for vertex in triangle.vertices:
            B[vertex]+=triangle_area*f(centroid)/3
        