# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:18:21 2021

@author: kst
"""

from decimal import Decimal
from fractions import Fraction

import numpy as np
import matplotlib.pylab as plt

def floating(x, l):
    if x == 0:
        return 0,1,0,0
    zero = [1 if x == 0 else 0]
    sign = [1 if x < 0 else 0]
    x = abs(x)
    exponent = int( l - np.log2(x) )
    significand ,r = divmod(2**exponent, 1/x)
    
    if r > 1/x * 1/2:
        significand += 1
        
    return  int(significand), exponent, zero[0], sign[0]

def add(s1,s2,l):
#    print(s2,s1)
#    print(rec(s1))
#    print(rec(s2))
 
    ar = rec(s1)
    br = rec(s2)
    
    c = ar+br
    return floating(c,l)
    


def mult(s1,s2,l):
    v1,p1,z1,s1 = s1
    v2,p2,z2,s2 = s2
    print(v1,v2)
    v = v1*v2 #is now 2*l bits long
    b1 = len(bin(v)[2:]) 
    b = b1-l #Tjeck how many bits we have truncated
    print(v)
    v = int(bin(v)[2:l+2],2)
    z = z1 or z2
    s = s1 ^ s2
    p = (p1+p2-b)*(1-z)
    return v,p,z,s

def div(s1,s2,l):
    
    v1,p1,z1,s1 = s1
    v2,p2,z2,s2 = s2
    
    ar= (1- 2*s1)*(1-z1) *v1*2**-p1
    br =  (1- 2*s2)*(1-z2) *v2*2**-p2
    c = ar/br
    v3,p,z3,s3 = floating(c,l)
    
    y = v1
    x = v2+z2
    
    phi=5
    for i in range(1,phi):
       y = y*(2**(l+1)-x)
       y = int(bin(y)[2:l+2],2)
       x = x*(2**(l+1)-x)
       x = int(bin(x)[2:l+2],2)
    y = y*(2**(l+1)-x)

    v = int(bin(y)[2:l+2],2)

    z=z1
    s=s1^s2

    
    return v,p,z,s

def MatVecProd(A,v):
    '''
    Computes securely the matrix vector product.
    Input: A, v is secret shared matrix and vector respectively.
    Outout: secret shared vector.
    '''
    I,J = np.shape(A)[:2]
    li = []
    for i in range(I):
        s = floating(0,l)
        for j in range(J):
            s = add(s, mult(v[j], A[i,j],l),l)
        li.append(s)
    return np.array(li)

def VecMatProd(v,A):
    '''
    Computes securely the vector matrix product.
    Input: A, v is secret shared matrix and vector respectively.
    Outout: secret shared vector.
    '''
    I,J = np.shape(A)[:2]
    li = []
    for j in range(J):
        s = floating(0,l)
        for i in range(I):
            s = add(s, mult(v[i], A[i,j],l),l ) 
        li.append(s)
    return np.array(li)

def VecVectProd(v1,v2):
    '''
    Computes securely the vector vector product.
    Input: v1 ( n x 1 ), and v2 ( 1 X n ) are secret shared vectors.
    Outout: secret shared matrix.
    '''
    I = len(v1)
    li = []
    for i in range(I):
        s = []
        for j in range(I):
            s.append(mult(v2[i], v1[j],l))
        li.append(s)
    return np.array(li)

def VectVecProd(v1,v2):
    '''
    Computes securely the vector vector product.
    Input: v1 ( 1 x n ), and v2 ( n X 1 ) are secret shared vectors.
    Outout: secret shared scalar.
    '''
    I = len(v1)
    s=floating(0,l)
    for i in range(I):
        s= add(s, mult(v1[i], v2[i],l),l)
    return s

def scalMatProd(s1,A):
    '''
    Computes securely the elementwise multiplication of a scalar and a matrix.
    Input: s1 is a secret shared scalar and A is a secret shared matrix.
    Outout: secret shared matrix.
    '''
    I,J = np.shape(A)[:2]
    li = []
    for i in range(I):
        s = []
        for j in range(J):
            s.append(mult(s1,A[i,j],l))
        li.append(s)
    return np.array(li)

def scalVecProd(s1, v):
    '''
    Computes securely the elementwise multiplication of a scalar and a vector.
    Input: s1 is a secret shared scalar and v is a secret shared vector.
    Outout: secret shared vector.
    '''
    R = [mult(s1, i,l) for i in v]
    return np.array(R)


#def recMat(A):
#    '''
#    Computes the plaint text of the matrix A.
#    '''
#    I,J = np.shape(A)[:2]
#    R = np.zeros((I,J))
#    for i in range(I):
#        for j in range(J):
#            R[i,j] = rns.rec(A[i,j])
#    return R
def rec(s1):
    v1,p1,z1,s1 = s1
    p1 = int(p1)
    return (1- 2*s1)*(1-z1) *v1*2**-p1
    
#
def recVec(v):
    '''
    Computes the plaint text of the vector v.
    '''
    return np.array([rec(i) for i in v]).reshape(len(v),1)




def RLS(dx, X, Y):
    '''
    Computes securely the recursive least squares equations.
    Input: dx is dimension of the system. X and Y are observations.
    Outout: secret shared estimate of parameters of the system.
    '''
    obs = len(Y)                    #number of observations
    Po = np.identity(dx)            #P0 matrix
    li = []
    #Secret share the P0 matrix
    for i in range(dx):
        li.append([floating(j,l) for j in Po[i,:] ])
    
    P = np.array(li)
    
    w = np.array([floating(0,l) for i in range(dx)]) #W0 parameter estimate
    w_open = []
    #secret share w0
    for i in range(obs):
        x = [floating(j,l) for j in X[i,:] ]          #Secret shared version of observations
        y = floating(Y[i],l)                   
        
        # Start with all parts of P
        d1 = add(floating(1,l), VectVecProd(x, MatVecProd(P, x)),l)
        d1_inv = div( floating(1,l), d1,l)
        d2 = VecVectProd( MatVecProd(P,x) , VecMatProd(x,P) )
        D = scalMatProd(d1_inv, d2)
        P = P - D
        
        # g
        g = MatVecProd(P,x)
        
        # e
        e = y - VectVecProd(x, w)
        
        #w
        w = w + scalVecProd(e,g)
        w_open.append(recVec(w))
    return recVec(w), w_open
        
l = 10
#Main stuff:
np.random.seed(1)    
    
n=3                      #Number of parties
t = 1                    #privacy treshold
x = np.linspace(0.5,2,n) #Index numbers of parties

dx = 6                   #Dimension of system

#rns = RNS(x,n,t)         #Class for the real number secret sharing scheme (all operations)

obs = 10                 #number of observations

X = np.random.normal(0,50,(obs, dx))       #Observations of input
#X = np.random.randint(0,5,(obs,dx))
Beta = np.random.normal(0,30, dx)          #True Parameters
#Beta = np.array([3.5,1.2,2.8,4.1,2.9,3.3])

y = np.array([np.dot(Beta, X[i,:]) for i in range(obs)])        #Observations of output
y = y + np.random.normal(0,5, obs)      #Add some noise to the output


w, w_open = RLS(dx, X, y)

print('True             ', '   estimat')
for i in range(dx):
    print('{}   {}'.format(Beta[i], w[i][0]))

mse = np.zeros(obs)
for i in range(obs):
    mse[i] = np.mean((Beta.reshape(dx,1) - np.array(w_open[i]))**2)

plt.plot(mse)











    
#l = 10
#a = 7.565
#b = 4.7
#c = a+b
#af = floating(a,l)
#bf = floating(b,l)
#cf = floating(c,l)
#v1,p1,z1,s1 = af
#v2,p2,z2,s2 = bf
#v3,p3,z3,s3 = cf
##print(v,p)
#
#cf1 = add(af,bf,l)
#v31,p31,z31,s31 = cf1
#
#ar= (1- 2*s1)*(1-z1) *v1*2**-p1
#br =  (1- 2*s2)*(1-z2) *v2*2**-p2
#cr = (1- 2*s3)*(1-z3) *v3*2**-p3
#
#c1r = (1- 2*s31)*(1-z31) *v31*2**-p31
#
##print(a , ar)
##print(b , br)
##print(c , cr)
#
#print(c, ar+br)
#
#print(c, c1r )
#
#
#print(v3,v31,p3,p31)

#print(len(bin(v1*v2)[2:l+2]))