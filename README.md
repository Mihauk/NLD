# $\textcolor{orange}{NLD\ -\ Non\ Linear\ Diffusion}$ 
- Project on studying Celluar automata(Kinetically Constrained lattice gas) like simulations to study the behaviour of correlation functions at late time in the presence of non-linearity(Non Linear Diffusion).
- Here we have implemented the simulattion in 3 different languages [python](./old) [Julia](./Julia) and [Rust](./nld_rust) in respective order in search for perrforrmance gains.
- Rust multi-threaded code seems to be the most stable and as fast as julia multi-threaded one.
- the fastest of all is [Julia GPU](./Gpu) code with multi-threading.

| Language (Parallelism) | O() Time per one standard operation |
| :---: | :---: |
| Python (without CPU threading) | $1\ \mu s$ |
| Julia (without CPU threading) | $10^{-1}\ \mu s $ |
| Julia (with CPU threading) | $1\ ns$ |
| Rust (without CPU threading) | $10\ ns$ |
| Rust (with CPU threading) | $1\ ns$ |
| Julia (with CPU threading + GPU) | $10^{-1}\ ns$ |
