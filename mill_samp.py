import numpy as np
import random
import copy
from matplotlib import pyplot as plt
import pickle

class conf:
    def __init__(self, rho_in, n):
        self.state = np.zeros(l)
        #self.cur = np.zeros(l)
        for i in range(l):
            if(np.random.uniform(0,1) <= rho_in[i]):
                self.state[i] = 1
            else:
                self.state[i] = 0
        self.ft = np.fft.ifft(self.state)
        #self.cur_ft = np.fft.ifft(self.cur)
        
    def update(self, n, par):
        self.cur = np.zeros(l)
        for p in range(n):
            no = 4*self.state[(3*p+par)%(3*n)] + 2*self.state[(3*p+1+par)%(3*n)] + self.state[(3*p+2+par)%(3*n)]
            if (no == 2):
                b = np.random.randint(2)
                #if (b==0):
                    #self.cur[(3*p+1+par)%(3*n)]=1
                #else:
                    #self.cur[(3*p+par)%(3*n)]=-1
                self.state[(3*p+par)%(3*n)] = b
                self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 1-b
            if (no == 1):
                #self.cur[(3*p+par)%(3*n)]=-1
                #self.cur[(3*p+1+par)%(3*n)]=-1
                self.state[(3*p+par)%(3*n)] = 1
                self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 0
            if (no == 4):
                #self.cur[(3*p+par)%(3*n)]=1
                #self.cur[(3*p+1+par)%(3*n)]=1
                self.state[(3*p+par)%(3*n)] = 0
                self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 1
        self.ft = np.fft.ifft(self.state)
        #self.cur_ft = np.fft.ifft(self.cur)

n = 8
l = 3*n
t = 100
A = 0.1
samples =50000000
freq = 2*np.pi*np.fft.fftfreq(l)
Q = freq[2]
#m_rho = np.zeros((l,t+1))
#m_j = np.zeros((l,t))
m_rho_ft = np.zeros((l,t+1),dtype='complex')
#m_j_ft = np.zeros((l,t),dtype='complex')
#qb2 = np.zeros(t+1,dtype='complex')
ini_conf = A*np.cos(Q*np.arange(l))+0.5
#ini_conf = (0.02*np.cos(freq[20]*np.arange(l)) + 0.02*np.cos(freq[40]*np.arange(l)))+0.6

for j in range(samples):
    x = conf(ini_conf,n)
    #m_rho[:,0] = copy.deepcopy(x.state) + m_rho[:,0]
    m_rho_ft[:,0] = copy.deepcopy(x.ft) + m_rho_ft[:,0]
    #qb2[0] = copy.deepcopy(x.ft[15])**2 + qb2[0]
    for i in range(t):
        x.update(n,np.random.randint(3))
        #qb2[i+1] = copy.deepcopy(x.ft[15])**2 + qb2[i+1]
        #m_rho[:,i+1] = copy.deepcopy(x.state) + m_rho[:,i+1]
        m_rho_ft[:,i+1] = copy.deepcopy(x.ft) + m_rho_ft[:,i+1]
        #m_j[:,i] = copy.deepcopy(x.cur) + m_j[:,i]
        #m_j_ft[:,i] = copy.deepcopy(x.cur_ft) + m_j_ft[:,i]

with open(f'../data/fms_0.5_A{A}_t{t}_k2_L{l}.pickle', 'wb') as data:
        #pickle.dump(qb2/samples,data)
	#pickle.dump([m_rho/samples, m_rho_ft/samples, m_j/samples, m_j_ft/samples], data)
	pickle.dump(m_rho_ft[2]/samples, data)
