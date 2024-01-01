import numpy as np
import copy
from tqdm import tqdm
from joblib import Parallel , delayed
from matplotlib import pyplot as plt
#from scipy import integrate
import toolz.itertoolz

class conf:
    # initiates the class when created
    def __init__(self, n, rho_in=None):
        self.n  = n
        self.l = 3*n
        #self.state = np.random.randint(2,size=self.l,dtype=int)
        #self.cur = np.zeros(l) # the current configration.
        # creating(or picking) a state randomly with probility distribution given by rho_in
        self.state = np.zeros(self.l, dtype=int) # the configation(density or state) of the sytem 
        for i in range(self.l):
            self.state[i] = 1 if np.random.uniform(0,1) <= rho_in[i] else 0
        #self.ft = np.fft.ifft(self.state) # Fourier Transform of the state configration
        #self.cur_ft = np.fft.ifft(self.cur) # Fourier Transform of the current configration

    def trimer_updates(self, trip, dp):
        #b = np.random.randint(2)
        #toss = np.random.uniform(0,1)
        a = (1,0,0), (0,0,1)
        b = (1,0,0), (0,1,0)
        c = (0,0,1), (0,1,0)
        if (trip == (0, 1, 0)):
            return a[np.random.choice(len(a), size = 1, p = [0.5+dp, 0.5-dp])[0]]
        elif (trip == (0, 0, 1)):
            return b[np.random.choice(len(b), size = 1, p = [0.5+dp, 0.5-dp])[0]]
        elif (trip == (1, 0, 0)):
            return c[np.random.choice(len(c), size = 1, p = [0.5-dp, 0.5+dp])[0]]
        else:
            return trip

    def part(self, st, off):
        return toolz.itertoolz.partition(3,np.roll(st,-off))

    # This method updates the state and current configrations according to specified rules 
    def update(self, st, off):
        self.state = np.roll(list(map(self.trimer_updates, self.part(self.state, off), 
                    np.array(list(self.part(st, off)))[:,0] - np.array(list(self.part(st, off)))[:,2]) ), off).flatten()
        #self.ft = np.fft.ifft(self.state)

n = 160
l = 3*n
t = 1000
A = 0.1
n0 = 0.5
n_jobs = 8
samples = 5
freq = 2*np.pi*np.fft.fftfreq(l)
q = 40
Q = freq[q]
rho_in = n0 + A*np.cos(Q*np.arange(l))

def realisation(t,ic=None):
    rho_ft = np.zeros((l,t+1),dtype=complex)
    x = conf(n,rho_in=ic)
    rho_ft[:,0] = copy.deepcopy((np.fft.ifft(x.state)))
    for i in range(t):
        x.update(ic ,np.random.randint(3))
        rho_ft[:,i+1] = copy.deepcopy((np.fft.ifft(x.state)))
        #rho_ft[:,i+1] = copy.deepcopy(np.abs(np.fft.ifft(x.state))**2)
    return rho_ft[q//2]*rho_ft[q//2]

d_ft =  np.asarray ( Parallel ( n_jobs = n_jobs ) ( delayed (realisation) (t,ic = rho_in) for j in  range ( samples )  ) )
