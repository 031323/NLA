def innerproduct(a,b):
	if len(a)!=len(b):
		return
	sum=0
	for i in range(0,len(a)):
		sum+=(a[i]*b[i])
	return sum

def gaussianelimination(A,b,scaledpartialpivot):
    if len(A)!=len(A[0]):
        return
    if len(A[0])!=len(b):
        return
    n=len(A)
    for j in range(0,n-1):
        if scaledpartialpivot:
            Max=0
            Maxi=-1
            for i in range(j,n):
                max=0
                for j_ in range(j,n):
                    if A[i][j_]>max:
                        max=A[i][j_]
                if max==0:
                    return
                for j_ in range(j,n):
                    A[i][j_]/=max
                b[i]/=max
                if A[i][j]>Max:
                    Max=A[i][j]
                    Maxi=i
            if Maxi==-1:
                return
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
                return
        for i in range(j+1,n):
            m=A[i][j]/A[j][j]
            for j_ in range(0,n):
                A[i][j_]-=m*A[j][j_]
            b[i]-=m*b[j]
    if A[i][i]==0:
        return
    x=[0]*n
    for i in range(n-1,-1,-1):
        s=0
        for j in range(i+1,n):
            s+=x[j]*A[i][j]
        x[i]=(b[i]-s)/A[i][i]
    return x