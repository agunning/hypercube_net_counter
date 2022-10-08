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
        l+=[(v,u) for u in gen_partitions(i,n,n) for v in gen_partitions(n-i,n,n)]
    return  l
def count_centralizer_st(partst):
    
    n=partst[0].shape[0]
    t=1
    for i in range(0,n):
        t*=(2*(i+1))**(int(partst[0][i])+int(partst[1][i]))*math.factorial(int(partst[0][i]))*math.factorial(int(partst[1][i]))
    return t
def count_trees_cycle_size(sp,unsp,u,v,l=0):
    n=2*sp+unsp
    if l==0:
        return v*(n*u+v)**(n-1-sp)*(n*u+v-2)**sp
    else:
        return u*(n*u)**(n-2-sp)*(n*u-2)**sp
def count_trees2(partst):
    n=partst[0].shape[0]
    if partst[0][0]+partst[0][1]==0:
            return 0
    if partst[0][0]==1:
        return 0
    if partst[0][0]==0:
        for i in range((n+1)//2):
            if partst[0][2*i]!=0:
                return 0
    L=2*n
    divsum=[0]*L
    for (i,t) in enumerate(partst[0][:]):
        if t!=0:
            j=2
            while j*(i+1)-1<L:
                divsum[j*(i+1)-1]+=2*t*(i+1)
                j+=1
    for (i,t) in enumerate(partst[1][:]):
        if t!=0:
            j=2
            while 2*j*(i+1)-1<L:
                divsum[2*j*(i+1)-1]+=2*t*(i+1)
                j+=1
    prod=1
    for i in range(0,L):
        if i<n:
            sp = int(partst[0][i])
        else:
            sp=0
        if i%2!=0:
            unsp = int(partst[1][(i-1)//2])
        else:
            unsp = 0
        if sp+unsp==0:
            continue
        if i==0:
            prod*=count_trees_cycle_size(sp,unsp,i+1,0,1)
        elif i==1 and divsum[i]==0:
            if sp ==1 and unsp ==0:
                prod*=2
            else:
                prod*=count_trees_cycle_size(sp,unsp,i+1,0,1)*2*sp
        else:
            prod*=count_trees_cycle_size(sp,unsp,i+1,int(divsum[i]))
    return prod
def count_nets2(n):
    tot=0
    bignum = 2**n*math.factorial(n)
    for x in gen_partitions_st(n):
        t=count_trees2(x)
        if t!=0:
            tot+=t*(bignum//count_centralizer_st(x))
    return tot//bignum
    

    
    
