# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 10:47:08 2021

@author: kst
"""

import numpy as np
np.random.seed(2)

A = np.random.randint(0,10,(1,3))
B = np.random.randint(0,10,(3,3))

X = np.random.randint(0,10, A.shape) 
Y = np.random.randint(0,10,B.shape)
Z = np.dot(X,Y)

E = A-X
D = B-Y

res = np.dot(E,D) + np.dot(E,Y) + np.dot(X,D) + Z

print(res)
print(np.dot(A,B))

a = np.random.randint(0,10)
b = np.random.randint(0,10)

x = np.random.randint(0,10)
y = np.random.randint(0,10)
z = x*y

d = a-x
e = b-y

res = e*d + d*y + x*e + z

print(res)
print(a*b)