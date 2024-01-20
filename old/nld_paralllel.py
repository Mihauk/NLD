# importing all the required Libraries
import numpy as np
#import random
import copy
from joblib import Parallel , delayed
#from numba import njit, prange
import h5py


# Creating a class that stores the configration of the system and had a method to update it

class conf:
    # initiates the class when created
    def __init__(self, rho_in, n):
        self.state = np.zeros(l,dtype=int) # the configation(density or state) of the sytem 
        #self.cur = np.zeros(l) # the current configration.
        # creating(or picking) a state randomly with probility distribution given by rho_in
        for i in range(l):
            if(np.random.uniform(0,1) <= rho_in[i]):
                self.state[i] = 1
            else:
                self.state[i] = 0
        self.ft = np.fft.ifft(self.state) # Fourier Transform of the state configration
        #self.cur_ft = np.fft.ifft(self.cur) # Fourier Transform of the current configration
    
    # This method updates the state and current configrations according to specified rules 
    def update(self, n, par): # par vaiable here defines which partition you want
        #self.cur = np.zeros(l)
        for p in range(n):
            trip = f"{self.state[(3*p+par)%(3*n)]}{self.state[(3*p+1+par)%(3*n)]}{self.state[(3*p+2+par)%(3*n)]}"
            if (trip == "010"):
                b = np.random.randint(2)
                '''if (b==0):
                    self.cur[(3*p+1+par)%(3*n)]=1
                else:
                    self.cur[(3*p+par)%(3*n)]=-1'''
                self.state[(3*p+par)%(3*n)] = b
                self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 1-b
            elif (trip == "001"):
                #self.cur[(3*p+par)%(3*n)]=-1
                #self.cur[(3*p+1+par)%(3*n)]=-1
                self.state[(3*p+par)%(3*n)] = 1
                #self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 0
            elif (trip == "100"):
                #self.cur[(3*p+par)%(3*n)]=1
                #self.cur[(3*p+1+par)%(3*n)]=1
                self.state[(3*p+par)%(3*n)] = 0
                #self.state[(3*p+1+par)%(3*n)] = 0
                self.state[(3*p+2+par)%(3*n)] = 1
        self.ft = np.fft.ifft(self.state)
        #self.cur_ft = np.fft.ifft(self.cur)




# This function takes a realization, updates it and returns it
def realisation(ic):
    rho_ft = np.zeros((l,t+1),dtype='complex')
    x = conf(ic,n)
    rho_ft[:,0] = copy.deepcopy(x.ft)
    for i in range(t):
        x.update(n,np.random.randint(3))
        rho_ft[:,i+1] = copy.deepcopy(x.ft)
    return rho_ft[q,:]



n = 100
n_jobs = 8 # number of parallel jobs
l = 3*n # total number of sites
t = 300 # final time
A = 0.1 # Amplitude of Disturbances
samples =2 #total number of samples
freq = 2*np.pi*np.fft.fftfreq(l)
q = 30
Q = freq[q] # The frequency of the initial disturbance
ini_conf = A*np.cos(freq[Q]*np.arange(l))+0.5 # The initial density configration
#ini_conf = (0.02*np.cos(freq[20]*np.arange(l)) + 0.02*np.cos(freq[40]*np.arange(l)))+0.6

d_ft =  np.asarray ( Parallel ( n_jobs = n_jobs ) ( delayed (realisation)  ( ini_conf )  for j in  range ( samples ) ) )


with h5py.File("./hts_0.5_A.1_t300_k13_L300_rhoQ.hdf5", "w") as f:
    dset = f.create_dataset("init", data=d_ft )
