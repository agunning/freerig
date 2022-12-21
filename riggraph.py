import math as math
import numpy as np
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
csize=np.ones(4**7,dtype=int)
cindex=np.arange(4**7)
for i in range(einit.shape[0]):
    cinit=einit[i]
    while cindex[cinit]!=cinit:
        cinit=cindex[cinit]
    cend=eend[i]
    while cindex[cend]!=cend:
        cend=cindex[cend]
    if cinit==cend:
        continue
    elif csize[cinit]>=csize[cend]:
        csize[cinit]+=csize[cend]
        cindex[cend]=cinit
    else:
        csize[cend]+=csize[cinit]
        cindex[cinit]=cend
for i in range(4**7):
    while cindex[i]!=cindex[cindex[i]]:
        cindex[i]=cindex[cindex[i]]
res=(cindex[(cindex-np.arange(4**7)==0)],csize[(cindex-np.arange(4**7)==0)].astype(int))

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
        output=output[:-1]+'0'
    return output
def get_res():
    ## Run this to list elements in our free monoid
    print('Element'+' '*17,'\t','Size of Isomorphism Class')
    _=[print(get_poly(res[0][i]),'\t', res[1][i]) for i in range(res[0].shape[0])]
    print('\n')
    print(res[0].shape[0],'elements')

    
