# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:25:42 2020

@author: kst
"""

import numpy as np
from RealNumberSecretSharing import RNS
import matplotlib.pyplot as plt

def MatVecProd(A,v):
    I,J = np.shape(A)[:2]
    l = []
    for i in range(I):
        s = 0
        for j in range(J):
            s+= rns.mult(v[j], A[i,j])
        l.append(s)
    return np.array(l)

def VecMatProd(v,A):
    I,J = np.shape(A)[:2]
    l = []
    for j in range(J):
        s = 0
        for i in range(I):
            s+= rns.mult(v[i], A[i,j])
        l.append(s)
    return np.array(l)

def VecVectProd(v1,v2):
    I = len(v1)
    l = []
    for i in range(I):
        s = []
        for j in range(I):
            s.append(rns.mult(v2[i], v1[j]))
        l.append(s)
    return np.array(l)

def VectVecProd(v1,v2):
    I = len(v1)
    s=0
    for i in range(I):
        s+= rns.mult(v1[i], v2[i])
    return s

def recMat(A):
    I,J = np.shape(A)[:2]
    R = np.zeros((I,J))
    for i in range(I):
        for j in range(J):
            R[i,j] = rns.rec(A[i,j])
    return R

def recVec(v):
    return np.array([rns.rec(i) for i in v]).reshape(len(v),1)

def scalMatProd(s1,A):
    I,J = np.shape(A)[:2]
    l = []
    for i in range(I):
        s = []
        for j in range(J):
            s.append(rns.mult(s1,A[i,j]))
        l.append(s)
    return np.array(l)

def scalVecProd(s1, v):
    R = [rns.mult(s1, i) for i in v]
    return np.array(R)


def RLS(dx, X, Y):
    xx = np.linspace(0.5,2,n)
    obs = len(Y)
    Po = np.identity(dx)
    l = []
    for i in range(dx):
        l.append([rns.sharing(j) for j in Po[i,:] ])
    
    P = np.array(l)
    
    w = np.array([rns.sharing(0) for i in range(dx)])
    w_open = []
    col = ['ro','bo','go','yo','mo','ko']
    for i in range(obs):
        x = [rns.sharing(j) for j in X[i,:] ]          #Consider the secret shared version of observations
        for j,jj in enumerate(x):
            plt.plot(xx,jj, col[j])
        y = rns.sharing(Y[i])
        
        # Start with all parts of P
        d1 = 1 + VectVecProd(x, MatVecProd(P, x))
        d1_inv = rns.div(d1)
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
        
 
n=3
t = 1
x = np.linspace(0.5,2,n)

dx = 6

rns = RNS(x,n,t)


np.random.seed(1)
obs = 1

#X = np.random.normal(0,5,(obs, dx))       #Observations of input
X = np.random.randint(0,5,(obs,dx))
#Beta = np.random.normal(5,10, dx)         #Parameter
Beta = np.array([3.5,1.2,2.8,4.1,2.9,3.3])

y = np.array([np.dot(Beta, X[i,:]) for i in range(obs)])        #Observations of output
y = y + np.random.normal(0,1, obs)


w, w_open = RLS(dx, X, y)

#print(recVec(w))
#print(Beta.reshape(dx,1))

#print('True             ', '   estimat')
#for i in range(dx):
#    print('{}   {}'.format(Beta[i], w[i][0]))
print(w)

mse = np.zeros(obs)
for i in range(obs):
    mse[i] = np.mean((Beta.reshape(dx,1) - np.array(w_open[i]))**2)

plt.plot(mse)

    