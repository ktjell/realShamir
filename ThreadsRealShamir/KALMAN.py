# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:17:03 2021

@author: kst
"""

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import time

def kf_predict(X, P, A, Q, B, U):
    X = np.dot(A, X) + np.dot(B, U)
    P = np.dot(A, np.dot(P, A.T)) + Q
#    print(P)
#    print(np.dot(P, A.T))
    
    return(X,P) 
 
def kf_update(X, P, Y, H, R, K):
    IM = np.dot(H, X)

    IS = R + np.dot(H, np.dot(P, H.T))
    
    K = np.dot(P, np.dot(H.T, linalg.inv(IS)))

    X = X + np.dot(K, (Y-IM))
#    print(np.dot(H, P))
    P = P - np.dot(K, np.dot(H, P))
#    print(X)
#    print(np.dot(P, H.T))
    
    
    return (X,P,K)

def sys(A,X, B,u, w, C, v):
    x = np.dot(A,X) + B*u + B*w
    y = np.dot(C,X) + v
    yp = np.dot(C,X)
    
    return x,y,yp


n =3
np.random.seed(1)

A = np.array([[1.1269 ,  -0.4940,    0.1129], 
              [1.0000 ,   0 ,  0 ], 
              [0 ,   1.00,    0] ])
    

B = np.array([ -0.3832,  0.5919, 0.5191]).reshape(3,1)

C = np.array([1, 0 ,0]).reshape(1,3)

H=C


K = np.eye(n)
I=20
D = 0;

Q = 0.1
R= 1

t = np.arange(I)
u = np.sin(t/5)

w = np.sqrt(Q)*np.random.randn(len(t),1)
v = np.sqrt(R)*np.random.randn(len(t),1)
x = np.array([1,2,3])
X = x.reshape(n,1)
true = []
meas = []
STATE = []
for i,j in enumerate(u):
    state, measure, TRUE = sys(A,X,B,j,w[i],C,v[i].reshape(1,1))
    X = state
    STATE.append(state.reshape(3,))
    true.append(TRUE)
    meas.append(measure)
    
M = 0
ST = np.array(STATE)
#ST = ST[:,M]
#plt.plot(ST, label = 'STATE')
#    
ME = np.array(meas)[:,0]
plt.plot(ME, label = 'MEASURE')


TRUE = np.array(true)[:,0]
plt.plot(TRUE, label = 'TRUE')


X = x.reshape(n,1)
P = 10*np.eye(n)

x_est = []
st = time.time()
for i in range(I):
    X, P = kf_predict(X,P,A,Q,B,u[i])
#    print(P)

    X,P,K = kf_update(X,P,meas[i],C,R,K)
    x_est.append(X.reshape(3,))
#    print(X)
    

sl = time.time()

x = np.array(x_est)[:,M]
x = x.reshape(I,1)
plt.plot(x, label= 'KALMAN')  
##plt.legend()
#
##plt.plot(TRUE-ME, label='MEASURE')
##plt.plot(TRUE- x, label = 'KALMAN' )
plt.legend()
print('Exe time: ', sl-st)
