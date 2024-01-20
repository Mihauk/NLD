 #import sys
import copy
import time
import argparse
import numpy as np
#from tqdm import tqdm
#import itertools as its
import toolz.itertoolz as itz
from matplotlib import pyplot as plt
from joblib import Parallel , delayed

class conf:
    # initiates the class when created
    def __init__(self, lx: int, ly: int, rho_in: np.ndarray):
        self.lx, self.ly = lx, ly #system size
        # creating(or picking) a state randomly with probility distribution given by rho_in
        self.state = (np.random.random(size=(lx,ly)) <= rho_in).astype(int)

    def partition(self, a):    
        part = []
        for i in range(0,self.lx,2):
            for j in range(0,self.ly,2):
                part.append([a[i,j],a[i,(j+1)],a[(i+1),j],a[(i+1),(j+1)]])
        return part
    
    def combine(self, a):
        c = np.zeros((self.lx,self.ly),dtype=int)
        p = 0
        for k in range(int(self.lx/2)):
            for j in range(0,4,2):
                com = []
                [com.append([a[(i+k*int(self.ly/2)),j],a[(i+k*int(self.ly/2)),j+1]])  for i in range(0,int(self.ly/2))]
                c[p,:] = np.asarray(com).flatten()
                p+=1
        return c
    
    def tetra_update(self, tetramer):
        if tetramer == [1,0,0,0]:
            return [[0,1,0,0],[0,0,1,0], [0,0,0,1]][np.random.randint(3)]
        elif tetramer == [0,1,0,0]:
            return [[1,0,0,0],[0,0,1,0], [0,0,0,1]][np.random.randint(3)]
        elif tetramer == [0,0,1,0]:
            return [[1,0,0,0],[0,1,0,0], [0,0,0,1]][np.random.randint(3)]
        elif tetramer == [0,0,0,1]:
            return [[1,0,0,0],[0,1,0,0], [0,0,1,0]][np.random.randint(3)]
        return(tetramer)
    
    def update(self,bias_x,bias_y):
        self.state = np.roll(np.roll(self.combine(np.asarray(list(map(self.tetra_update,self.partition(np.roll(np.roll(self.state,bias_x,axis=0),bias_y,axis=1)))))),-bias_y,axis=1),-bias_x,axis=0)

def realisation(t: int, lx: int, ly: int, ic: np.ndarray) -> np.ndarray:
    rho = np.zeros((lx,ly,t+1),dtype=int)
    x = conf(lx,ly,rho_in=ic)
    rho[:,:,0] = copy.deepcopy(x.state)
    for i in range(t):
        x.update(np.random.randint(lx),np.random.randint(ly))
        rho[:,:,i+1] = copy.deepcopy((x.state))
        #rho_ft[:,i+1] = copy.deepcopy(np.abs(np.fft.ifft(x.state))**2)
    return rho

lx, ly = 100, 100
tmax = 200
n_jobs = -1
samples = 1000
freq_x = 2*np.pi*np.fft.fftfreq(lx)
freq_y = 2*np.pi*np.fft.fftfreq(ly)
q_x, q_y = 30, 30
Q_x, Q_y = freq_x[q_x], freq_y[q_y]

A_x, A_y = 0.1, 0.1
n0 = 0.5
xx, yy = np.meshgrid(np.arange(lx), np.arange(ly))
rho_in = n0 + A_x*np.cos(Q_x*xx) +  A_y*np.cos(Q_y*yy)

b_x, b_y = np.random.randint(lx),np.random.randint(ly)

s = time.time()
d_ft =  np.mean( np.asarray ( Parallel ( n_jobs = n_jobs, prefer= 'threads' ) ( delayed (realisation) (tmax,lx,ly,ic = rho_in) for j in  range ( samples )  ) ), axis=0)
e = time.time()
print ((e-s)/(tmax*lx*ly*samples))