import copy
import math
import numpy
def innerproduct(a,b):
	if len(a)!=len(b):
		raise ValueError
	sum=0
	for i in range(0,len(a)):
		sum+=(a[i]*b[i])
	return sum

def transpose(M):
    i=len(M)
    j=len(M[0])
    M_=[]
    for i_ in range(0,j):
        row=[]
        for j_ in range(0,i):
            row.append(M[j_][i_])
        M_.append(row)
    return M_
    
def product(A,B):
    if len(A[0])!=len(B):
        raise ValueError
    l1=len(A);l2=len(B);l3=len(B[0])
    C=[]
    for l1_ in range(0,l1):
        row=[]
        for l3_ in range(0,l3):
            sum=0
            for l2_ in range(0,l2):
                sum+=A[l1_][l2_]*B[l2_][l3_]
            row.append(sum)
        C.append(row)
    return C

def thomasalgorithm(a,b,c,r):
    for i in range(0,len(r)):
        if i>0:
            b[i]-=c[i-1]*a[i]
            r[i]-=r[i-1]*a[i]
            a[i-1]=0
        if i<len(c):c[i]/=b[i]
        r[i]/=b[i]
        b[i]=1
    x=[0]*len(r)
    x[len(r)-1]=r[len(r-1)]
    for i in range(len(r)-2,-1,-1):
        x[i]=r[i]-c[i]*x[i+1]
    return x


def cholesky(A):
    r=[[0 for i in range(len(A))] for j in range(len(A))]
    for i in range(0,len(A)):
        for j in range(i,len(A)):
            if i==j:
                s=0
                for k in range(0,i):
                    s+=r[i][k]*r[i][k]
                r[i][i]=math.sqrt(A[i][i]-s)
            else:
                s=0
                for k in range(0,i):
                    s+=r[i][k]*r[j][k]
                r[j][i]=(A[i][j]-s)/r[i][i]
    return r

def gaussianelimination(A_,b_,scaledpartialpivot):
    A=copy.deepcopy(A_)
    b=copy.deepcopy(b_)
    if len(A)!=len(A[0]):
        raise ValueError
    if len(A[0])!=len(b):
        raise ValueError
    n=len(A)
    for j in range(0,n-1):
        if scaledpartialpivot:
            Max=0
            Maxi=-1
            for i in range(j,n):
                max=0
                for j_ in range(j,n):
                    if abs(A[i][j_])>max:
                        max=abs(A[i][j_])
                if max==0:
                    raise ValueError
                for j_ in range(j,n):
                    A[i][j_]/=max
                b[i]/=max
                if abs(A[i][j])>Max:
                    Max=abs(A[i][j])
                    Maxi=i
            if Maxi==-1:
                raise ValueError
            for j_ in range(j,n):
                A[j][j_],A[Maxi][j_]=A[Maxi][j_],A[j][j_]
            b[j],b[Maxi]=b[Maxi],b[j]
        if A[j][j]==0:
            for i in range(j+1,n):
                if A[i][j]!=0:
                    for j_ in range(j,n):
                        A[j][j_],A[i][j_]=A[i][j_],A[j][j_]
                    b[j],b[i]=b[i],b[j]
            if A[j][j]==0:
                raise ValueError
        for i in range(j+1,n):
            m=A[i][j]/A[j][j]
            for j_ in range(0,n):
                A[i][j_]-=m*A[j][j_]
            b[i]-=m*b[j]
    if A[n-1][n-1]==0:
        raise ValueError
    x=[0]*n
    for i in range(n-1,-1,-1):
        s=0
        for j in range(i+1,n):
            s+=x[j]*A[i][j]
        x[i]=(b[i]-s)/A[i][i]
    return x

def givens_parameters(a,b):
    c,s=0,0
    if a==0:c,s=1,0
    elif abs(a)>abs(b):
        t=b/a
        s=1/math.sqrt(1+t*t)
        c=s*t
    else:
        t=a/b
        c=1/math.sqrt(1+t*t)
        s=c*t
    return c,s
    
def givens_multiplier(A,i,j,c,s):
    A_=numpy.array(A,dtype=numpy.float)
    a=A_[i,:].copy()
    b=A_[j,:].copy()
    A_[i,:]=c*a+s*b
    A_[j,:]=-s*a+c*b
    return A_.tolist()
def givens(A_):
    A=copy.deepcopy(A_)
    m=len(A)
    n=len(A[0])
    Q=[[0 for i in range(m)] for j in range(m)]
    for i in range(m):Q[i][i]=1
    for i in range(0,min(m-1,n)):
        for j in range(i+1,m):
            c,s=givens_parameters(A[j][i],A[i][i])
            A=givens_multiplier(A,i,j,c,s)
            Q=givens_multiplier(Q,i,j,c,s)
    return transpose(Q),A

def jacobi(A_,b_,iterations,threshold):
    A=copy.deepcopy(A_)
    b=copy.deepcopy(b_)
    n=len(b)
    x=[0]*n
    x2=[0]*n
    K=0
    t=threshold
    for k in range(iterations):
        K+=1
        for i in range(n):
            sum=b[i]
            for j in range(n):
                if j!=i:
                    sum-=A[i][j]*x[j]
            x2[i]=sum/A[i][i]
        x=x2.copy()
        t=numpy.linalg.norm(numpy.array(b)-numpy.array(A)@numpy.array(x))/numpy.linalg.norm(numpy.array(b))
        if t<threshold:break
       
    return x,K,t

def gauss_siedel(A_,b_,iterations,threshold):
    A=copy.deepcopy(A_)
    b=copy.deepcopy(b_)
    n=len(b)
    x=[0]*n
    K=0
    t=threshold
    for k in range(iterations):
        K+=1
        for i in range(n):
            sum=b[i]
            for j in range(n):
                if j!=i:
                    sum-=A[i][j]*x[j]
            x[i]=sum/A[i][i]
        t=numpy.linalg.norm(numpy.array(b)-numpy.array(A)@numpy.array(x))/numpy.linalg.norm(numpy.array(b))
        if t<threshold:break
    return x,K,t
    
def SOR(A_,b_,w,iterations,threshold):
    A=copy.deepcopy(A_)
    b=copy.deepcopy(b_)
    n=len(b)
    x=[0]*n
    x2=[0]*n
    K=0
    t=threshold
    for k in range(iterations):
        K+=1
        for i in range(n):
            sum=b[i]*w
            for j in range(n):
                if j!=i:
                    sum-=w*A[i][j]*x2[j]
            x2[i]=sum/A[i][i]+(1-w)*x[i]
        x=x2.copy()
        t=numpy.linalg.norm(numpy.array(b)-numpy.array(A)@numpy.array(x))/numpy.linalg.norm(numpy.array(b))
        if t<threshold:break
       
    return x,K,t

def pseudo_inverse(A):
    return product(numpy.linalg.inv(product(transpose(A),A)).tolist(),transpose(A))

def ro(A):
    eigenvalues=numpy.linalg.eigvals(A)
    for i in range(len(eigenvalues)):
        eigenvalues[i]=abs(eigenvalues[i])
    return abs(max(eigenvalues))

def condition_number(A):
    A1=pseudo_inverse(A)
    return math.sqrt(ro(product(A,transpose(A)))*ro(product(A1,transpose(A1))))

def main():
    A=[[2,0],[3,4]]
    print("\nA=\n",numpy.array(A))
    print("\ncholesky of A*A'  \n",numpy.array(cholesky(product(A,transpose(A)))))
    
    B=[5,8]
    print("\nB=",B)
    x=gaussianelimination(A,B,True)
    print("\ngaussian elimination(spp) of Ax=B  x=",x)
    print("Ax=\n",numpy.array(product(A,transpose([x]))))

    print("\nJacobi of Ax=B  x,iterations,error=\n",(jacobi(A,B,10,0.0001)))
    print("\nGauss Siedel of Ax=B  x,iterations,error=\n",(gauss_siedel(A,B,10,0.0001)))
    print("\nSOR of Ax=B  x,iterations,error=\n",(SOR(A,B,1.2,100,0.0001)))
    
    print("\nA=\n",numpy.array(A))
    print("\nGivens of A:\n")
    Q,R=givens(A)
    print("Q:\n",numpy.array(Q))
    print("R:\n",numpy.array(R))
    print("QR:\n",numpy.array(product(Q,R)))

    A=[[1,6,1],[4,2,1],[0,3,5],[2,6,9],[-1,5,-8]]
    print("\nA=:\n",numpy.array(A))
    b=transpose([[2,12,5,4,7]])
    print("\nb=:\n",numpy.array(b))
    x=product(pseudo_inverse(A),b)
    print("\nsolution of Ax=b using pseudo_inverse: x=:\n",numpy.array(x))
    print("\nAx:\n",numpy.array(product(A,x)))
    
    A=[[2,4,8,16],[3,8,27,81],[1,5,25,125],[4,16,64,256],[1,6,36,215],[1,7,49,343]]
    print("\nLet A be a 6x4 Vandermonde matrix. A:\n",numpy.array(A))
    print("\nIts condition_number is: ",condition_number(A))
    print("As per matlab, 1.1733e+03\n")

if __name__ == "__main__":
    main()