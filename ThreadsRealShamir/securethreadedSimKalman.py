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

#class server:
#    securecom = {}
#    broadcasts = {}
#    def __init__(self,F, n, t, numTrip, l = 7):
#        self.b = ss.share(F,np.random.choice([-1,1]), t, n)
#        self.triplets = [proc.triplet(F,n,t) for i in range(numTrip)]
#        self.r, self.rb = proc.randomBitsDealer(F,n,t,l)

class communicationSimulation:
    def __init__(self,q):
        self.q = q
        
    def com(self,add, val):
        self.q[add].put(val)
    
class party(Thread):
    party_addr = ips.party_addr
    def __init__(self, i, n, shares, com, q, triplet):
        Thread.__init__(self)
        self.i = i
        self.n = n
        self.c = 0
        self.comr = 0
        self.shares = shares
        self.triplet = triplet
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
        st = time.time()
        res = []
        for i in range(self.n):
            while name + str(i) not in self.recv:
                self.readQueue()    
            res.append(self.recv[name+str(i)])
            del self.recv[name + str(i)]
        sl = time.time()
        self.comtime += sl-st
        return res

    def mult(self,s1,s2):
        a,b,c = self.triplet
        a = a[0]
        b= b[0]
        c = c[0]
        ds = s1 - a
        es = s2 - b
        self.broadcast('s'+ str(self.comr), [ds,es])
        shares = self.get_shares('s'+ str(self.comr))
        self.comr+=1
        dss = [i[0] for i in shares]
        ess = [i[1] for i in shares]
        d = rns.rec(dss)
        e = rns.rec(ess)
        return d*e + d*b + a*e + c
    
    def multimult(self,s1,s2):
        n = len(s1)
        a,b,c = self.triplet
        ds = [i-j for i,j in zip(s1,a)]
        es = [i-j for i,j in zip(s2,b)]
        self.broadcast('s'+ str(self.comr), [ds,es])
        shares = self.get_shares('s'+ str(self.comr))
        self.comr+=1
        DS = np.array([i[0] for i in shares])
        ES = np.array([i[1] for i in shares])

        S = []
        for i in range(n):
            d = rns.rec(DS[:,i])
            e = rns.rec(ES[:,i])
            S.append(d*e + d*b[i] + a[i]*e + c[i])
                        
        return np.array(S) 

    def inv(self, s1):
        a,b,c = self.triplet
        a = a[0]
        ds = self.mult(s1,a)
        self.broadcast('is'+ str(self.comr), ds)
        dss = self.get_shares('is'+ str(self.comr))
        self.comr +=1
        d = rns.rec(dss)
        return 1/d * a        
    
    
    def MatMatProd(self,A,B):
        '''
        Computes securely the matrix matrix product.
        Input: A, B is secret shared matrices, respectively.
        Outout: secret shared matrix.
        '''
        I,J = np.shape(A)
        L,K = np.shape(B)
        
        s2 = np.tile(A.reshape(I*J,), K).reshape(1,I*J*K)
        
        s1 = np.repeat(B,I,axis=1).transpose().reshape(1,L*K*I)

        mul = self.multimult(s1,s2)

        S = np.sum(mul.reshape(J**2,I),axis=1).reshape(K,I).transpose()

        return S 
    
    
    def MatVecProd(self,A,v):
        '''
        Computes securely the matrix vector product.
        Input: A, v is secret shared matrix and vector respectively.
        Outout: secret shared vector.
        '''
        I,J = np.shape(A)
        
        s1 = np.tile(v.reshape(J,),I).reshape(I*J,1)
        s2 = A.reshape(I*J,1)
                
        mul = self.multimult(s1,s2)
        
        S = np.sum(mul.reshape(I,J),axis=1).reshape(I,1)

        return S

    def VecMatProd(self,v,A):
        '''
        Computes securely the vector matrix product.
        Input: A, v is secret shared matrix and vector respectively.
        Outout: secret shared vector.
        '''
        I,J = np.shape(A)
   
        s1 = np.tile(v.reshape(J,),J).reshape(I*J,1)
        s2 = A.transpose().reshape(I*J,1)
        mul = self.multimult(s1,s2)
        
        S = np.sum(mul.reshape(I,J),axis=1).reshape(1,J)
        

       
        return S
    
    def VecVectProd(self,v1,v2):
        '''
        Computes securely the vector vector product.
        Input: v1 ( n x 1 ), and v2 ( 1 X n ) are secret shared vectors.
        Outout: secret shared matrix.
        '''
        I = len(v1)
        s1 = np.repeat(v1, I)
        s2 = np.repeat(v2,I).reshape(I,I).transpose().reshape(I*I,1)
        mul = self.multimult(s1,s2)
        S = mul.reshape(I,I)
        
        return S
    
    def VectVecProd(self,v1,v2):
        '''
        Computes securely the vector vector product.
        Input: v1 ( 1 x n ), and v2 ( n X 1 ) are secret shared vectors.
        Outout: secret shared scalar.
        '''

        mul = self.multimult(v1,v2)
        S = np.sum(mul)
        return S
    
    def scalMatProd(self,s1,A):
        '''
        Computes securely the elementwise multiplication of a scalar and a matrix.
        Input: s1 is a secret shared scalar and A is a secret shared matrix.
        Outout: secret shared matrix.
        '''
        I,J = np.shape(A)
        s1 = np.repeat(s1,I*J).reshape(1,I*J)
        s2 = A.reshape(1,I*J)
        mul = self.multimult(s1,s2)
        
        S = mul.reshape(I,J)
        
        
        return S
    
    def scalVecProd(self,s1, v):
        '''
        Computes securely the elementwise multiplication of a scalar and a vector.
        Input: s1 is a secret shared scalar and v is a secret shared vector.
        Outout: secret shared vector.
        '''
        I,J = np.shape(v)
 
        i = max(I,J)

        s1 = np.repeat(s1,i).reshape(1,i)
        s2 = v.reshape(1,i)
        mul = self.multimult(s1,s2)

        S = mul.reshape(I,J)
        
        return S
    
    
    def kf_predict(self, x, P, A, Q, b, u):

        X1 = self.MatVecProd(A, x) 
        
        X2 = self.scalVecProd(u, b)
        X = X1+X2
        P1 = self.MatMatProd(P, A.T)
        P2 = self.MatMatProd(A, P1) 
        P = P2 + Q
      
        return X,P 
 
    def kf_update(self, x, P, y, H, R, K):

        IM = self.VectVecProd(H, x)
        IS1 = self.MatVecProd(P, H.T)
        IS = R + self.VectVecProd(H,  IS1)
        IS_inv = self.inv(IS)
        K = self.MatVecProd(P, self.scalVecProd(IS_inv, H.T))
        X = x + self.scalVecProd((y-IM), K)
        
        H1 = self.VecMatProd(H, P)
    
        P = P - self.VecVectProd(K, H1 )

        
        return (X,P,K)
    
        
        
    
    
    
        
    def run(self):
        
        A,B,C,X,P,Q,R,K,u,meas = self.shares
        
        
#        
        x_est = []
        st = time.time()
        for i in range(I):
            X, P = self.kf_predict(X,P,A,Q,B,u[i])
                
            X,P,K = self.kf_update(X,P,meas[i],C,R,K)
            x_est.append(X.reshape(3,))
            
            
        sl = time.time()
        self.x_est = x_est
##    
#        print(sl-st - self.comtime)
        self.exe = sl-st - self.comtime
        
#        #shares of secrets
#        As,bs, s1 = self.shares
#        
#        ## ADD
##        self.shares = s1+s2
#        
##        ## Mult
##        self.shares = self.multimult( [s1,s1,s2], [s2,s1,s2])
#        self.shares = self.VectVecProd(H,X) #+ 
#        inv = self.inv(u[0])
#        print(inv)
#        st = time.time()
#        self.shares = self.VecMatProd(C,P)
#        sl = time.time()
#        print(sl-st)
#        self.shares = self.scalVecProd(u[0],B)
##        self.shares = self.shares1 + self.shares2  
#        
###        ## Inv
###        self.shares = self.inv(s1)
##        
##        
##        #reconstruct
###        

##            print(x_est)
#        


def sys(A,X, B,u, w, C, v):
    x = np.dot(A,X) + B*u + B*w
    y = np.dot(C,X) + v
    yp = np.dot(C,X)
    
    return x,y,yp

np.random.seed(1)

I = 50

            
n = 3                    

A = np.random.normal(0,1,(3,3))      
#A = np.array([[1.1269 ,  -0.4940,    0.1129], 
#              [1.0000 ,   0 ,  0 ], 
#              [0 ,   1.00,    0] ])
    

B = np.random.normal(0,1,(n,1)) #np.array([ -0.3832,  0.5919, 0.5191]).reshape(n,1)

C = np.array([1, 0 ,0]).reshape(1,n)

H=C

K = np.eye(n)

Q = 0.1
R= 1

t = np.arange(I)
u = np.sin(t/5).reshape(I,1)

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

meas = np.array(meas).reshape(I,1)
    
n = 3
t = 1
x = np.linspace(1,n,n)
rns = RNS(x,n,t)

x = np.array([1,2,3])
X = x.reshape(n,1)


As = rns.matrixsharing(A)
Bs = rns.matrixsharing(B)
Cs = rns.matrixsharing(C)
Xs = rns.matrixsharing(X)
Ps = rns.matrixsharing(P)
Ks = rns.matrixsharing(K)
Qs = rns.sharing(Q)
Rs = rns.sharing(R)
us = rns.matrixsharing(u)
meass = rns.matrixsharing(meas)

atri, btri, ctri = rns.triplet(20)

ques = [que.Queue() for k in range(n)]

com = communicationSimulation(ques)
 
#A,B,C,X,P,Q,R,K,u,meas



threads = []
for k in range(n):
    threads.append(party(k,n, [As[k], Bs[k], Cs[k], Xs[k], Ps[k], Qs[k], Rs[k], Ks[k], us[k], meass[k]], com, ques[k],
                         [ [i[k] for i in atri], [i[k] for i in btri], [i[k] for i in ctri] ] ))

start = time.time()
    
for t in threads:
    t.start()


for t in threads:
    t.join()
    
end = time.time()
ex = end-start


x_shares = []
for t in threads:
    x_shares.append(t.x_est)
    
times = []
for t in threads:
    times.append(t.exe)
times = np.array(times)
    
x_est = []
for i in range(I):
    sh = []
    for t in x_shares:
        sh.append(t[i].reshape(n,1))
    x_est.append(rns.recmatrix(sh))


ME = np.array(meas)[:,0].reshape(I,1)
plt.plot(ME[:,0], label = 'MEASURE')
#

TRUE = np.array(true)[:,0]
plt.plot(TRUE[:,0], label = 'TRUE')



x = np.array(x_est)[:,0]
x = x.reshape(I,1)
plt.plot(x[:,0], label= 'KALMAN')  
#plt.legend()

#plt.plot(TRUE-ME, label='MEASURE')
#plt.plot(TRUE- x, label = 'KALMAN' )
plt.legend()
#            
            
            
            
            
            
print('Avr. execution time: ', np.average(times)/I)

print('Execution time: ', ex, 's')   


