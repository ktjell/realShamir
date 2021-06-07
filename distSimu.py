# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 12:53:57 2021

@author: kst
"""

import numpy as np
import matplotlib.pyplot as plt

n = 7
s = 1
P = np.linspace(0.5,2)
t = 4

Alfa_in = [np.arange(i,n-t+i+1) for i in range(t)]

A = [ [P[i] for i in Alfa_in[j]] for j in range(t) ]


Lb = 0

T = [i for i in range(t) if i != Lb]

x = P[2]

x_a = [ np.sort([x - i for i in A[j]]) for j in T ]

alf_alf = [ np.sort( [i - j for i in A[Lb] for j in A[k] if i-j != 0  ] ) for k in T  ]

div = [ [i /j for i in g1 for j in g2 ] for g1,g2 in zip(x_a,alf_alf) ]

p1 = np.outer(np.array(div[0]), np.array(div[1])).flatten()

for i in range(2,len(div)):
    p1 = np.outer(p1, np.array(div[i])).flatten()

y = np.random.normal(0,1,size=len(p1))


plt.hist(p1*y)



