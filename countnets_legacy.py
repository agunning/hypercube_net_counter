import math as math
import numpy as np
import scipy.special
import sympy
def gen_partitions(n,m=None,xt=None,add=None):
    if m is None:
        m=n
    if xt is None:
        xt=n
    if add is None:
        add = np.zeros(xt).astype(int)
    l=[]
    if n==0:
        return [add]
    for i in range(1,min(m,n)+1):
        nadd=add.copy()
        nadd[i-1]+=1
        l+=gen_partitions(n-i,i,xt,nadd)
    return l
def gen_partitions_st(n):
    l=[]
    for i in range(0,n+1):
        l+=[np.vstack((v,u)) for u in gen_partitions(i,n,n) for v in gen_partitions(n-i,n,n)]
    return  l
def count_centralizer_st(partst):
    n=partst.shape[1]
    t=1
    for i in range(0,n):
        t*=(2*(i+1))**(int(partst[0,i]+partst[1,i]))*math.factorial(int(partst[0,i]))*math.factorial(int(partst[1,i]))
    return t    
def build_count_spanning_trees_st(partst):
    if partst[0,0]+partst[0,1]>0:
        if partst[0,0]==0:
            for i in range(0,(partst.shape[1]+1)//2):
                if partst[0,2*i]!=0:
                    return 0
        if partst[0,0]==1:
            return 0
        l=[]
        for i in range(0,partst.shape[1]):
            l+=[i+1,-i-1]*partst[0,i]
            l+=[2*(i+1)]*partst[1,i]
        A=np.zeros((len(l),len(l))).astype(int)
        for i in range(0,len(l)):
            for j in range(0,len(l)):
                #i->j iff i|j
                if l[j]%l[i]!=0:
                    continue
                if (l[j]+l[i]==0) and (abs(i-j)==1) and ((j-i)*l[j]<0):
                    A[i][j]-=abs(l[i])-1
                else:
                    A[i][j]-=abs(l[i])
        for i in range(0,len(l)):
            A[i,i]-=A[:,i].sum()
        A[0,0]+=1
        if sympy.Matrix(A).det()==0:
            print(partst)
        if partst[0,0]>0:
            return sympy.Matrix(A).det()
        
        else:
            return 2*partst[0,1]*sympy.Matrix(A).det()
    else:
        return 0
def count_nets(n):
    parts=gen_partitions_st(n)
    bignum=math.factorial(n)*2**n
    if n<2:
        return 0
    for part in parts:
        a=build_count_spanning_trees_st(part)
        if a!=0:
            n+=a*(bignum//count_centralizer_st(part))
    return n//bignum
    

    
    
