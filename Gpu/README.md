# Julia code using CUDA GPU's to time evolve a lattice gas diffusion model(with non-linearity).

- The code takes in command line arguments l,tmax,samples and total_samples. And saves the output in two different files as avereage over all the samples the density configration and the current configration.

- Here is an example to run it 
### setsid julia nld.jl 99999 625 1000 10000 &
### nohup julia nld.jl 99999 625 1000 10000 &
- Here setsid and nohup keeps the process running in the background even if the terminal session has ended or disconnected and & in the end gives the cursor back to the terminal.
- The argumens are positional so the position matters and is in the specific order l,tmax,samples, total_samples.
