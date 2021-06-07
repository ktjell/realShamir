# -*- coding: utf-8 -*-
"""
Created on Thu May 20 10:26:49 2021

@author: kst
"""

import numpy as np
from RealNumberSecretSharing import RNS
import matplotlib.pyplot as plt

def MatVecProd(A,v):
    '''
    Computes securely the matrix vector product.
    Input: A, v is secret shared matrix and vector respectively.
    Outout: secret shared vector.
    '''
    I,J = np.shape(A)[:2]
    l = []
    for i in range(I):
        s = 0
        for j in range(J):
            s+= rns.mult(v[j], A[i,j])
        l.append(s)
    return np.array(l)

def VecMatProd(v,A):
    '''
    Computes securely the vector matrix product.
    Input: A, v is secret shared matrix and vector respectively.
    Outout: secret shared vector.
    '''
    I,J = np.shape(A)[:2]
    l = []
    for j in range(J):
        s = 0
        for i in range(I):
            s+= rns.mult(v[i], A[i,j])
        l.append(s)
    return np.array(l)

def VecVectProd(v1,v2):
    '''
    Computes securely the vector vector product.
    Input: v1 ( n x 1 ), and v2 ( 1 X n ) are secret shared vectors.
    Outout: secret shared matrix.
    '''
    I = len(v1)
    l = []
    for i in range(I):
        s = []
        for j in range(I):
            s.append(rns.mult(v2[i], v1[j]))
        l.append(s)
    return np.array(l)

def VectVecProd(v1,v2):
    '''
    Computes securely the vector vector product.
    Input: v1 ( 1 x n ), and v2 ( n X 1 ) are secret shared vectors.
    Outout: secret shared scalar.
    '''
    I = len(v1)
    s=0
    for i in range(I):
        s+= rns.mult(v1[i], v2[i])
    return s

def scalMatProd(s1,A):
    '''
    Computes securely the elementwise multiplication of a scalar and a matrix.
    Input: s1 is a secret shared scalar and A is a secret shared matrix.
    Outout: secret shared matrix.
    '''
    I,J = np.shape(A)[:2]
    l = []
    for i in range(I):
        s = []
        for j in range(J):
            s.append(rns.mult(s1,A[i,j]))
        l.append(s)
    return np.array(l)

def scalVecProd(s1, v):
    '''
    Computes securely the elementwise multiplication of a scalar and a vector.
    Input: s1 is a secret shared scalar and v is a secret shared vector.
    Outout: secret shared vector.
    '''
    R = [rns.mult(s1, i) for i in v]
    return np.array(R)


def recMat(A):
    '''
    Computes the plaint text of the matrix A.
    '''
    I,J = np.shape(A)[:2]
    R = np.zeros((I,J))
    for i in range(I):
        for j in range(J):
            R[i,j] = rns.rec(A[i,j])
    return R

def recVec(v):
    '''
    Computes the plaint text of the vector v.
    '''
    return np.array([rns.rec(i) for i in v]).reshape(len(v),1)

np.random.seed(1)    
    
n=3                      #Number of parties
t = 1                    #privacy treshold
x = np.linspace(0.5,2,n) #Index numbers of parties


rns = RNS(x,n,t)         #Class for the real number secret sharing scheme (all operations)