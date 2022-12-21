import math as math
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import connected_components
A =np.arange(0,4**7).reshape((-1,1))
B =np.array([
    [0, 2, 2, 0, 0, 2, 0],
    [3, 1, 1, 3, 3, 1, 1],
    [0, 2, 2, 0, 0, 2, 2],
    [3, 1, 1, 3, 3, 1, 3],
    [0, 2, 2, 0, 4, 2, 4],
    [3, 1, 1, 3, 3, 5, 5],
    [0, 1, 2, 3, 4, 5, 6]]) ##Multplication table: elements are bab,aba,ba,ab,b,a,1
X=np.repeat(np.arange(0,6),np.arange(6,0,-1))
Y=np.repeat(np.cumsum(np.arange(7,1,-1)),np.arange(6,0,-1))-np.arange(21)-1
M=4**X+4**Y
N=4**B[X,Y]+4**B[Y,X]
U=A+N
V=~A&~A>>1&M<1
U-=A&N&10922|(A|N)&(A&N&5461)<<1 #Replacing 4x with 2x etc.
edges=np.array(np.nonzero(V))
einit=edges[0]
eend=U[edges[0],edges[1]]
G=sp.sparse.csr_matrix(sp.sparse.coo_matrix((np.ones_like(einit),(einit,eend)),shape=(4**7,4**7)))
elts,labels=connected_components(G,directed=False,return_labels=True)
sizes=np.bincount(labels)
sorter=np.argsort(labels,kind='stable')
reprs=sorter[np.searchsorted(labels,np.arange(elts),sorter=sorter)]
def get_poly(num):
    symbols=['','a','b','ab','ab','aba','bab']
    output=''
    nnum=np.zeros(7,dtype=int)
    X=np.array([int(t) for t in np.base_repr(num,4)])
    nnum[-X.shape[0]:]+=X
    for i in range(7):
        if (nnum[i]!=0):
            nstr=symbols[i]
            if (nnum[i]>1) or (i==0):
                nstr=str(nnum[i])+nstr
            else:
                nstr=' '+nstr
            if (not output.isspace()) and (output!=''):
                nstr='+'+nstr
            elif i!=0:
                nstr=' '+nstr
            output+=nstr
        else:
            output+=' '*(2+len(str(symbols[i]))-(i==0))
    if num==0:
        output='0'+output[:-1]
    return output
def get_res(path):
    ## Run this to list elements in our free monoid
    lines=[]
    lines.append('Element'+' '*17+'\t'+'Size of Isomorphism Class'+'\n')
    _=[lines.append(get_poly(reprs[i])+'\t'+str(sizes[i])+'\n') for i in range(elts)]
    lines.append('\n')
    lines.append(str(elts)+'elements')
    print(lines)
    with open(path,'a') as f:
        f.writelines(lines)
m=np.arange(4**7)
ms=[m[labels==i] for i in range(elts)]
mprimes=[(ms[i].reshape((-1,1))//(4**np.arange(7)))%4 for i in range(elts)]
Bprime=np.zeros((7,7,7),dtype=int)
Bprime[np.arange(7),np.arange(7).reshape((-1,1)),B]+=1
def test():
    ##because i haven't done relations of the form e.g x(y+z)^2=x(y+z) this checks if our relations are actually well defined (and all our elements are idempotent)
    for i in range(elts):
        mprimei =mprimes[i]
        for j in range(elts):
            print(i,j)
            mprimej =mprimes[j]
            T  = np.einsum('ai,bj,ijk->abk',mprimei,mprimej,Bprime)
    
            multtable = labels[((2*(T>1)+T%2)*4**np.arange(7)).sum(axis=-1)]
            S= mprimei.reshape((-1,1,7))+mprimej
            addtable = labels[((2*(S>1)+S%2)*4**np.arange(7)).sum(axis=-1)]
            
            if (multtable.min()<multtable.max()) or (addtable.min()<addtable.max()):
                return 1
            if (i==j) and (i!=multtable[0][0]):
                print(multtable)
                print(addtable)
                print(T)
                
                return 1
    return 0
                            
'''
mt=np.zeros((elts,elts,elts))
mt[np.arange(elts),np.arange(elts),multtable]+=1

termmult = np.nonzero(mt.sum(axis=1).prod(axis=0)*mt.sum(axis=0).prod(axis=0))

at=np.zeros((elts,elts,elts))
at[np.arange(elts),np.arange(elts),addtable]+=1

termadd = np.nonzero(at.sum(axis=0).prod(axis=0))

r=np.zeros((elts,elts,elts),dtype=int)

i = (r+np.arange(elts).reshape(-1,1,1)).reshape(-1)
j=(r+np.arange(elts).reshape(1,-1,1)).reshape(-1)
k=(r+np.arange(elts).reshape(-1)).reshape(-1)

ends = np.concatenate([multtable[i,k]*elts+multtable[j,k],multtable[k,i]*elts+multtable[j,i],addtable[i,k]*elts+addtable[j,k]])
starts = np.concatenate([elts*i+j]*3)

H=sp.sparse.csr_matrix(sp.sparse.coo_matrix((np.ones_like(einit),(einit,eend)),shape=(elts**2,elts**2)))
e,labels2=connected_components(G,connection='strong',return_labels=True)

endsl=ends[labels2]
startsl=starts[labels2]
'''
