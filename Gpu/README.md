# Julia code using CUDA GPU's to time evolve a lattice gas diffusion model(with non-linearity).

- The code takes in command line arguments
  -l chain_lenght
  -t time
  -s samples
  -a amplitude of the intital perturbation
  -q wave number mode of initial perturbation
  -n equilibrium density
  
And saves the output as .h5 file with the first four Cumulants stored with dataset names m1,m2,m3,m4. The file is saved in the current folder.


## An Example to run it 

```diff
@@ `nohup julia nld_1d_cuda.jl -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &` @@
@@ `setsid julia nld_1d_cuda.jl -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &` @@
```

- Here setsid and nohup keeps the process running in the background even if the terminal session has ended or disconnected and & in the end gives the cursor back to the terminal.
- The argumens are positional so the position matters and is in the specific order l,tmax,samples, total_samples.

## Threading
- set the number of CPU threads
- `nohup julia nld_1d_cuda.jl --threads 32 -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &`
- `setsid julia nld_1d_cuda.jl --threads 32 -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 -n 0.5 &`
- here note that the samples > number of threads.
