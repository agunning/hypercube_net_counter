import scipy.special
from Crypto.Util import number
import math as math

def getPrimeMod(N,n):
    '''
Generates a prime of the from 2**n*k+1 and a primitive 2**nth root
'''
    num = (number.getRandomNBitInteger(N)&((1<<N)-(1<<n)))|(1+(1<<N))
    while (not number.isPrime(num)):
        num = num+(1<<n)
    isnotgen=True
    while isnotgen:
        q=number.getRandomNBitInteger(N)%num
        g = pow(q,num>>n,num)
        r = pow(g,1<<(n-1),num)
        if r==num-1:
            isnotgen=False
    return num,g

def fftMod(values,p,n,g):
    '''
    fast fourier transform mod p, over 2**n elements and using g as a generator. Completely unnecessary optimization.
    '''
    output=[0]*(2**n)
    if len(values)==0:
        return output
    if n==0:
        return [values[0]]
    else:
        split1=fftMod(values[::2],p,n-1,pow(g,2,p))
        split2=fftMod(values[1::2],p,n-1,pow(g,2,p))
        for i in range(0,2**n):
            output[i]+=split1[i%(2**(n-1))]+pow(g,i,p)*split2[i%(2**(n-1))]
    return output


def ifftMod(values,p,n,g):
    '''
invesrse 

    '''
    ginv=pow(g,-1,p)
    ninv=pow(2**n,-1,p)
    noutput = fftMod(values,p,n,ginv)
    return [x*ninv%p for x in noutput]


def count_orthoplex_nets(d):
    ''' Counts the number of nets of a d dimensional orthoplex (exactly). Main issue is that the answer grows as exp(exp(d)) so that for d>20 or so this starts taking up a substantial chunk of memory'''
    
    
    rsum = math.prod((2*i)**scipy.special.comb(d,i,exact=True) for i in range(1,d+1))>>d
    N=d+2
    n=int(math.log2(d+1))+1
    p,g=getPrimeMod(N,n)
    jfft=fftMod([1,-1],p,n,g)

    kfftb=fftMod([1,0,1],p,n,g)
    lfftb=fftMod([1,1],p,n,g)
    x=[scipy.special.comb(d,i,exact=True) for i in range(0,d+1)]
    for k in range(0,(d-1)//2+1):
            l=d-1-2*k
            kfft=[pow(kfftb[i],k,p) for i in range(0,2**n)]
            lfft=[pow(lfftb[i],l,p) for i in range(0,2**n)]
            toifft=[(jfft[i]*kfft[i]*lfft[i]%p) for i in range(0,2**n)]
            adj=ifftMod(toifft,p,n,g)
            print(adj)
            #print(k,l,p,adj,x)
            trees= math.prod([(2*i)**(((x[i]+adj[i])%p)//2) for i in range(1,d+1)])
            mult=d*scipy.special.comb(d-1,2*k,exact=True)
            mult*=math.prod([2*i+1 for i in range(0,k)])
            #print(trees)
            #print(mult)
            rsum+=trees*mult
    return rsum//(2**d*math.factorial(d))
        
