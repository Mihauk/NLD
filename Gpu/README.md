# $\textcolor{orange} {Julia \ code \ using  \ CUDA  \space GPU's  \space to  \space time  \space evolve  \space a  \space lattice \space  gas  \space diffusion  \space model(with  \space non-linearity).}$

- The code takes in command line arguments
  -l chain_lenght
  -t time
  -s samples
  -a amplitude of the intital perturbation
  -q wave number mode of initial perturbation
  -n equilibrium density
  
And saves the output as .h5 file with the first four Cumulants stored with dataset names m1,m2,m3,m4. The file is saved in the current folder.

$${\color{red}Welcome \space \color{lightblue}To \space \color{orange}Stackoverflow}$$

## An Example to run it 

```diff
nohup julia nld_1d_cuda.jl -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &
```
```diff
setsid julia nld_1d_cuda.jl -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &
```

- Here setsid and nohup keeps the process running in the background even if the terminal session has ended or disconnected and & in the end gives the cursor back to the terminal.
- The argumens are positional so the position matters and is in the specific order l,tmax,samples, total_samples.

## Threading
- set the number of CPU threads
```diff
nohup julia nld_1d_cuda.jl --threads 32 -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &
```
```diff
setsid julia nld_1d_cuda.jl --threads 32 -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 -n 0.5 &
```
- here note that the samples > number of threads.
