# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:02:02 2020

Implementation of the Real Number Secret Sharing scheme. 
The implementation consists of the sharing and recon algorithm of the scheme as well as multiplication and division operations.
@author: kst
"""

import random
import numpy as np
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt


class RNS():
    
    def __init__(self, x, n, t):
        self.x = x
        self.n = n
        self.t = t

    def sharing(self, secret):
        xt = np.insert(np.sort(random.sample(list(self.x),self.t)), 0, 0)
        y = [np.random.normal(0, 1000) for i in range(self.t)]
        y.insert(0,secret)
        poly = lagrange(xt,y)
        shares = np.polyval(poly.coef, self.x)
        return shares.reshape(self.n,1)
    
    def matrixsharing(self, A):
        I,J = np.shape(A)
        sh = [np.zeros((I,J)) for i in range(self.n)]
        for i in range(I):
            for j in range(J):
                shares = self.sharing(A[i,j])
                for N in range(self.n):
                    sh[N][i,j] = shares[N]
                    
        return sh
    
    def recmatrix(self, A):
        I,J = np.shape(A[0])
        B = np.zeros((I,J))
        
        for i in range(I):
            for j in range(J):
                S = []
                for N in range(self.n):
                    S.append(A[N][i,j])
                B[i,j] = self.rec(S)
        
        return B       
             
    
    def rec(self,y):
        s = 0
        for j in range(self.n):
            p = 1
            for i in [k for k in range(self.n) if k != j]:
                p *= -self.x[i]/(self.x[j]-self.x[i])
            s+=y[j]*p
        return s
    
    def triplet(self,i=1):
        
        a = np.random.normal(0, 1000,i)
        b = np.random.normal(0, 1000,i)
        c = [j*k for j,k in zip(a,b)]
        
        if i == 1:
            return [self.sharing(j) for j in a][0], [self.sharing(j) for j in b][0], [self.sharing(j) for j in c][0]
        
        return [self.sharing(j) for j in a], [self.sharing(j) for j in b], [self.sharing(j) for j in c]
    
    def mult(self,s1,s2):
        a,b,c = self.triplet()
        d = self.rec(s1-a)
        e = self.rec(s2-b)
        return d*e + d*b + a*e + c
        
    def div(self, s1):
        a,b,c = self.triplet()
        d = self.rec(self.mult(s1,a))
        return 1/d * a
    
    def multimult(self,s1):
        n = len(s1)
        a,b,c = self.triplet(n)
        d = [self.rec(i-j) for i,j in zip(s1,a)]
        e = [self.rec(i-j) for i,j in zip(s1,b)]
        return [k*l + k*j + i*l + h for i,j,k,l,h in zip(a,b,d,e,c)]
        
    
#
#n=3
#t = 1
#
#np.random.seed(1)
##
##secret = 5
####
#x = np.linspace(1,n,n)
##
#rns = RNS(x,n,t)
#
#s1 = rns.sharing(5)
#s2 = rns.sharing(3)
#
#s3 = s1*s2
#
#s4 = rns.mult(s3,s3)
#
#print(rns.rec(s4))

#
#A = np.random.normal(0,100,(3,1))
#
##print(rns.triplet(3))
##
###
#s1 = rns.matrixsharing(A)
###for i in s1:
###    print(i)
##
#print(rns.recmatrix(s1))
#print(A)
###s2 = rns.sharing(1.4)
#print(rns.rec(s2))
#s3 = rns.sharing(5.56)
#
#print(rns.multimult([s1,s2,s3]))

#plt.plot(x[:3],s1[:3], 'go')
#plt.plot(x[:3],s2[:3], 'bo')
#K, d,e= rns.mult(s1,s2)
#plt.plot(4,d, 'ro')
#plt.plot(5,e, 'yo')
#plt.plot(0,34.5, 'ro')
#print(rns.rec(K))
#print(34.5*3.42)
##print(rns.rec(rns.div(s1)))