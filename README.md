# $\textcolor{orange}{NLD\ -\ Non\ Linear\ Diffusion}$ 
- Project on studying Celluar automata(Kinetically Constrained lattice gas) like simulations to study the behaviour of correlation functions at late time in the presence of non-linearity(Non Linear Diffusion).
- Here we have implemented the simulattion in 3 different languages [python](./old) [Julia](./Julia) and [Rust](./nld_rust) in respective order in search for performance gains.
- Rust multi-threaded code seems to be the most stable and as fast as julia multi-threaded one.
- the fastest of all is [Julia GPU](./Gpu) code with multi-threading.

| Language | CPU threading | GPU | O() Time per one standard operation |
| :---: | :---: | :---: | :---: |
| [Python](./old) | ❌ | ❌ | 1 $\mu$ s |
| Python | ✅ | ❌ | 0.1 $\mu$ s |
| [Julia](./Julia/nld.jl) | ❌ | ❌ | 0.1 $\mu$ s |
| [Rust](./nld_rust) | ❌ | ❌ | 10 $n$ s |
| [Julia](./Julia/nld.jl) | ✅ | ❌ | 1 $n$ s |
| [Rust](./nld_rust) | ✅ | ❌ | 1 $n$ s|
| [Julia](./Gpu) | ✅ | ✅ | 0.1 $n$ s|
