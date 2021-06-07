# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:39:36 2021

@author: kst
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:22:50 2021

@author: kst
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 14:47:14 2018

@author: kst
"""
import numpy as np
from threading import Thread
import FFArithmetic as field
import shamir_scheme as ss
import proc
import time
import queue as que
from ipcon import ipconfigs as ips
from ipcon import network as nw
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import random
from RealNumberSecretSharing import RNS
from numpy import linalg

class server:
    securecom = {}
    broadcasts = {}
    def __init__(self,F, n, t, numTrip, l = 7):
        self.b = ss.share(F,np.random.choice([-1,1]), t, n)
        self.triplets = [proc.triplet(F,n,t) for i in range(numTrip)]
        self.r, self.rb = proc.randomBitsDealer(F,n,t,l)

class communicationSimulation:
    def __init__(self,q):
        self.q = q
        
    def com(self,add, val):
        self.q[add].put(val)
    
class party(Thread):
    party_addr = ips.party_addr
    def __init__(self, i, n, shares, com, q):
        Thread.__init__(self)
        self.i = i
        self.n = n
        self.c = 0
        self.comr = 0
        self.shares = shares
        self.com = com
        self.q = q
        self.comr = 0        
        self.comtime = 0
        self.recv = {}
        
    def readQueue(self):
        while not self.q.empty():
            b = self.q.get()
            self.recv[b[0]] = b[1]
            
    def broadcast(self, name, s):
        for i in range(n):
            self.com.com(i, [name + str(self.i), s])
            
                    
    
    def get_shares(self, name):
        res = []
        for i in range(self.n):
            while name + str(i) not in self.recv:
                self.readQueue()    
            res.append(self.recv[name+str(i)])
            del self.recv[name + str(i)]
        return res


    
    def kf_predict(self, X, P, A, Q, B, U):
#        print(np.dot(A, X))
        X = np.dot(A, X) + np.dot(B, U)
        P = np.dot(A, np.dot(P, A.T)) + Q
        return(X,P) 
     
    def kf_update(self, X, P, Y, H, R, K):
        IM = np.dot(H, X)
        IS = R + np.dot(H, np.dot(P, H.T))
        K = np.dot(P, np.dot(H.T, linalg.inv(IS)))
        X = X + np.dot(K, (Y-IM))
    
        P = P - np.dot(K, np.dot(H, P))
        
        
        return (X,P,K)
    
        
        
    
    
    
        
    def run(self):
        
        A,B,C,X,P,Q,R,K,u,meas = self.shares
        
        
        
        x_est = []
        
        st = time.time()        
        for i in range(I):
            
            X, P = self.kf_predict(X,P,A,Q,B,u[i])
#            print(P)
#        
            X,P,K = self.kf_update(X,P,meas[i],C,R,K)
            x_est.append(X.reshape(3,))
#            print(X)
        sl = time.time()
        self.exe = sl-st
        self.x_est = x_est
    
   

def sys(A,X, B,u, w, C, v):
    x = np.dot(A,X) + B*u + B*w
    y = np.dot(C,x) + v
    yp = np.dot(C,X)
    
    return x,y,yp

np.random.seed(10)
I = 50

            
n = 3                    

#A = np.random.normal(0,1,(3,3))        
A = np.array([[1.1269 ,  -0.4940,    0.1129], 
              [1.0000 ,   0 ,  0 ], 
              [0 ,   1.00,    0] ])
    

B = np.array([ -0.3832,  0.5919, 0.5191]).reshape(n,1)

C = np.array([1, 0 ,0]).reshape(1,3)

H=C

K = np.eye(n)

Q = 0.1
R= 1

t = np.arange(I)
u = np.sin(t/5)

w = np.sqrt(Q)*np.random.randn(len(t),1)
v = np.sqrt(R)*np.random.randn(len(t),1)
x = np.array([1,2,3])
X = x.reshape(n,1)

X = x.reshape(n,1)
P = 10*np.eye(n)


true = []
meas = []
STATE = []

for i,j in enumerate(u):
    state, measure, TRUE = sys(A,X,B,j,w[i],C,v[i])
    X = state
    STATE.append(state.reshape(3,))
    true.append(TRUE)
    meas.append(measure)

#meas = np.array(meas).reshape(I,1)
    
n = 3
t = 1
x = np.linspace(1,n,n)
m=3

ques = [que.Queue() for k in range(m)]

com = communicationSimulation(ques)
 
#A,B,C,X,P,Q,R,K,u,meas
x = np.array([1,2,3])
X = x.reshape(n,1)



threads = []
for k in range(m):
    threads.append(party(k,m, [A, B, C, X, P, Q, R, K, u, meas], com, ques[k], ))

start = time.time()
    
for t in threads:
    t.start()


for t in threads:
    t.join()
    
end = time.time()
ex = end-start


times = []
for t in threads:
    times.append(t.exe)
times = np.array(times)
#    


ME = np.array(meas)[:,0].reshape(I,1)
plt.plot(ME[:,0], label = 'MEASURE')
#

TRUE = np.array(true)[:,0]
plt.plot(TRUE[:,0], label = 'TRUE')



x = np.array(threads[0].x_est)[:,0]
x = x.reshape(I,1)
plt.plot(x[:,0], label= 'KALMAN')  
#plt.legend()

#plt.plot(TRUE-ME, label='MEASURE')
#plt.plot(TRUE- x, label = 'KALMAN' )
plt.legend()
#            

plt.figure()

plt.plot(np.sqrt((x_sec - x[:,0])**2))
            
            
print('Avr. execution time: ', np.average(times)/I)
            


print('Execution time: ', ex, 's')   


