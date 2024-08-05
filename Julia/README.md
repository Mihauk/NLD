# $\textcolor{lightblue} {Julia code to time evolve a lattice gas diffusion model(with non-linearity).}

- The code takes in command line arguments
  -l chain_lenght
  -t time
  -s samples
  -a amplitude of the intital perturbation
  -q wave number mode of initial perturbation
  -n equilibrium density
  
And saves the output as .h5 file with the first four Cumulants stored with dataset names m1,m2,m3,m4. The file is saved in the current folder.

## An Example to run it 
### nohup julia nnld.jl 99999 625 1000 10000 &
- The argumens are optional so the position does not matter.

## Threading
- set the number of CPU threads
### setsid julia nld.jl --threads 32 99999 625 1000 10000 &
- here note that the samples > number of threads.
