def gen_paired_partition(n):
    return flatten([[(u,v) for u in Partitions(n-i) for v in Partitions(i)] for i in range(0,n+1)],max_level=1)

def count_cent(pairpart):
    return 2**(len(pairpart[0])+len(pairpart[1]))*pairpart[0].centralizer_size()*pairpart[1].centralizer_size()

def count_trees(pairpart):
    if len(pairpart[0])==0:
        return 0
    if len(pairpart[0].to_exp())>1:
        if pairpart[0].to_exp()[0]+pairpart[0].to_exp()[1]==0:
            return 0
    if pairpart[0].to_exp()[0]==1:
        return 0
    if pairpart[0].to_exp()[0]==0:
        for i in range((len(pairpart[0].to_exp())+1)//2):
            if pairpart[0].to_exp()[2*i]!=0:
                return 0
    L=[2*k for k in pairpart[1]]+flatten([[k,-k] for k in pairpart[0]])
    A = matrix(ZZ,len(L)-1,len(L)-1,0)
    for (i,P) in enumerate(L):
        for (j,Q) in enumerate(L):
            if j==len(L)-1:
                continue
            if Q%P!=0:
                continue
            if (Q+P==0) and (abs(i-j)==1) and ((j-i)*Q<0):
                if i!=len(L)-1:
                    A[i,j]-=abs(P)-1
                A[j,j]+=abs(P)-1
            else:
                if i!=len(L)-1:
                    A[i,j]-=abs(P)
                A[j,j]+=abs(P)
    if pairpart[0].to_exp()[0]>0:
        return A.det()
    else:
        return 2*pairpart[0].to_exp()[1]*A.det()

def count_nets(n):
    tot=0
    for x in gen_paired_partition(n):
        t=count_trees(x)
        if t!=0:
            tot+=t/count_cent(x)
    return tot
    
