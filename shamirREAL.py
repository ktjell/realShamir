# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 10:41:01 2020

@author: kst
"""

import random
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial
import itertools

kk = np.random.randint(0,10000)
np.random.seed(kk)
random.seed(kk)
# Creates shares of secrets using Shamir's secret sharing scheme.
# This is the "BAD" way to do it!
def share(secret, t, n, X):
    shares = []
    c = [np.random.normal(0, 1000) for i in range(t)]#random coefficients
#    print(c)
    c.insert(0,secret)  #
#    print(np.roots(c)) #Roots of the polynomial
    pol = np.poly1d(c[::-1]) #Make the polynomial from coefficients
    shares = pol(X) #Evaluate the polynomial in the X values (indices from participants)
    
    xlong = np.linspace(0.5,2, 1000) #For ploting, evaluate the pol in x values from -1,2
    plt.figure('The bad way')
    plt.plot(xlong, pol(xlong), label ='polynomial') 
    plt.plot(X, shares, 'rD', label = 'shares')
    plt.legend()
    
    
    return shares

#This is the better way to do it (we think!)
def create_shares(secret,t,n, x):

    xt = np.insert(np.sort(random.sample(list(x),t)), 0, 0) #Choose randomly x-value
    y = [np.random.normal(0, 1000) for i in range(t)] #choose the y-value

    y.insert(0,secret)  #Insert the secret
    poly = lagrange(xt,y)#Use lagrange interpolation to find the polynomial
#    print(poly) #get the polynomial
#    print(np.roots(y)) #get the roots of the polynomial
    shares = np.polyval(poly.coef, x) #Evaluate the polynomial in all x-values to get the shares

    xlong = np.linspace(0,2.1, 1000)
    plt.figure('The better way')
    plt.plot(xlong, np.polyval(poly.coef, xlong), label = 'polynomial')
    plt.plot(x,shares, 'bD', label = 'shares')
    plt.plot(xt,y, 'rD', label = 'random generated shares')
    plt.legend()
    
        
    return shares


# Creates the "recombination"-vector used to reconstruct a secret from its shares.
def basispoly(n):
    r = []
    C = range(1, n+1)

    for i in range(1, n+1):
        c = [k for k in C if k != i]
        p = 1
        for j in range(n - 1):
            p *= -c[j] / ((i) - c[j])
        r.append(p)
    return r


# reconstruct secret.
#def rec(x):
#    res = 0
#    n = len(x)
#    y = basispoly(n)
#    for i in range(len(x)):
#        res += x[i] * y[i]
#    return res

def rec2(x,y):
    s = 0
    for j in range(len(x)):
        p = 1
        for i in [k for k in range(len(x)) if k != j]:
            p *= -x[i]/(x[j]-x[i])
        s+=y[j]*p
    return s




#np.random.seed(1)
    
n=11 #Participants
t = 5 #Degree of polynomial
secret = 50  
X = np.linspace(0.5,2,n) # The x values for the polynomial (the index of the participants)

#Create shares the "bad" way:
SHARES1 = share(secret,t,n, X)


#Create shares the better way:
SHARES = create_shares(secret,t,n, X)
#plt.plot(X,SHARES, 'g')

#Evaluating the reconstruction of secret based on all different combinations of shares:
XX = set(itertools.combinations(set(np.arange(n)), t+1))
error = []
for i in XX:
    honest = list(i) #random.sample(range(n), t+1)
    x_h = list(map(X.__getitem__, honest))
    shares_h = list(map(SHARES.__getitem__, honest))
    e = abs(rec2(x_h,shares_h) - secret)
#    print(rec2(x_h,shares_h))
    error.append( e )
print(min(error), max(error))
    


#print(rec2(x,shares))



#
#print(rec2( list(map(X.__getitem__, honest)), list(map(shares.__getitem__, honest)), t+1))
#print(shares/(max(shares)))
#plt.plot(X,shares)
#for i in range(20):
#    shares = share(i,t,n)
#    print(rec2(X,shares,n))
##    print(rec(shares))
#    plt.plot(X,shares)
#    plt.plot(X[:t],shares[:t], 'bo')