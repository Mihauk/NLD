import numpy as np
import random
import copy
from joblib import Parallel , delayed
import h5py

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
        #self.cur = np.zeros(l)
        for p in range(n):
            no = 4*self.state[(3*p+par)%(3*n)] + 2*self.state[(3*p+1+par)%(3*n)] + self.state[(3*p+2+par)%(3*n)]
            if (no == 2):
                b = np.random.randint(2)
                '''if (b==0):
                    self.cur[(3*p+1+par)%(3*n)]=1
                else:
                    self.cur[(3*p+par)%(3*n)]=-1'''
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



def realisation(ic):
    rho_ft = np.zeros((l,t+1),dtype='complex')
    x = conf(ic,n)
    rho_ft[:,0] = copy.deepcopy(x.ft)
    for i in range(t):
        x.update(n,np.random.randint(3))
        rho_ft[:,i+1] = copy.deepcopy(x.ft)
    return rho_ft[q,:]



n = 100
n_jobs = 8
l = 3*n
t = 300
A = 0.1
samples =100000
freq = 2*np.pi*np.fft.fftfreq(l)
q = 30
Q = freq[q]
ini_conf = A*np.cos(Q*np.arange(l))+0.5
#ini_conf = (0.02*np.cos(freq[20]*np.arange(l)) + 0.02*np.cos(freq[40]*np.arange(l)))+0.6

d_ft = np . asarray ( Parallel ( n_jobs = n_jobs ) ( delayed (realisation ) ( ini_conf ) for j in  range ( samples ) ) )

with h5py.File("../data/hts_0.5_A.1_t300_k13_L300_rhoQ.hdf5", "w") as f:
    dset = f.create_dataset("init", data=d_ft )
