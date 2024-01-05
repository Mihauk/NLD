 #import sys
import copy
import time
import argparse
import numpy as np
#from tqdm import tqdm
#import itertools as its
import toolz.itertoolz as itz
#from matplotlib import pyplot as plt
from joblib import Parallel , delayed

class conf:
    # initiates the class when created
    def __init__(self, l: int, rho_in: np.ndarray):
        self.l = l #system size
        # creating(or picking) a state randomly with probility distribution given by rho_in
        self.state = (np.random.random(size=self.l) <= rho_in).astype(int)

    def trimer_updates(self, trip: tuple) -> tuple:
        if (trip == (0, 1, 0)):
            return [(1,0,0), (0,0,1)][np.random.randint(2)]
        elif (trip == (0, 0, 1)):
            return [(1,0,0), (0,1,0)][np.random.randint(2)]
        elif (trip == (1, 0, 0)):
            return [(0,0,1), (0,1,0)][np.random.randint(2)]
        return trip

    # This method updates the state and current configrations according to specified rules 
    def update(self, off: int):
        #self.state = np.roll(list(map(self.trimer_updates, itz.partition(3,np.roll(self.state,-off)))), off).flatten()
        self.state = np.roll(list(itz.mapcat(self.trimer_updates, itz.partition(3,np.roll(self.state,-off)))),off)

def realisation(t: int, l: int,ic: np.ndarray) -> np.ndarray:
    rho = np.zeros((l,t+1),dtype=int)
    x = conf(l,rho_in=ic)
    rho[:,0] = copy.deepcopy(x.state)
    for i in range(t):
        x.update(np.random.randint(3))
        rho[:,i+1] = copy.deepcopy((x.state))
        #rho_ft[:,i+1] = copy.deepcopy(np.abs(np.fft.ifft(x.state))**2)
    return rho

parser = argparse.ArgumentParser()
parser.add_argument('--l', default='999999', type=int, help='System size')
parser.add_argument('--t', default='100', type=int, help='tmax')
parser.add_argument('--s', default='10', type=int, help='sample size')

args = parser.parse_args()

l = args.l
t = args.t
A = 0.1
n0 = 0.5
n_jobs = -1
samples = args.s
freq = 2*np.pi*np.fft.fftfreq(l)
q = 6
Q = freq[q]
rho_in = n0 + A*np.cos(Q*np.arange(l))

s = time.time()
d_ft =  np.mean(np.asarray ( Parallel ( n_jobs = n_jobs ) ( delayed (realisation) (t,l,ic = rho_in) for j in  range ( samples )  ) ), axis=0 )
e = time.time()

print ((e-s)/(t*l*samples))
