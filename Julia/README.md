# $\textcolor{lightblue} {Julia\ code\ to\ time\ evolve\ a\ lattice\ gas\ diffusion\ model(with\ non-linearity).}$

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
nohup julia nld.jl -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &
```
- The argumens are optional so the position does not matter.

## Threading
- set the number of CPU threads
```diff
nohup julia nld.jl --threads 32 -l 99999 -t 625 -s 100000 -a 0.1 -q 9999 &
```
- here note that the samples > number of threads.
