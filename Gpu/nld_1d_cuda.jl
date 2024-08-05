using CUDA
CUDA.allowscalar(false)
using FFTW
using HDF5
using Random
using ArgParse
using Statistics
using LinearAlgebra

mutable struct Total_data
    m1::Vector{Float64}
    m2::Vector{Float64}
    m3::Vector{Float64}
    m4::Vector{Float64}
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--chain_length", "-l"
            help = "Length of the 1D chain"
            arg_type = Int
            default = 999
        "--t_max", "-t"
            help = "Number of time steps"
            arg_type = Int
            default = 100
        "--samples", "-s"
            help = "Number of samples for each saved file"
            arg_type = Int
            default = 1000
        "--amplitude", "-a"
            help = "Amplitude of the intital perturbation from equilibrium density at n. n is constrained to be 0<n<1"
            arg_type = Float16
            default = 0.1
        "--eq_den", "-n"
            help = "Equilibrium density"
            arg_type = Float16
            default = 0.5
        "--wave_number", "-q"
            help = "wave number mode of initial perturbation. It is constrained to be 0<q<l/2"
            arg_type = Int
            default = 10
        "--t_samples", "-S"
            help = "Total number of samples for the entire run"
            arg_type = Int
            default = 10
        #"--flag1"
        #    help = "an option without argument, i.e. a flag"
        #    action = :store_true
        #"arg1"
        #    help = "a positional argument"
        #    required = true
    end

    return parse_args(s)
end

function update_d!(state,p)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    str = gridDim().x * blockDim().x
    for j in idx:str:length(p)
        idx1 = 3 * (j - 1) + 1
        idx2 = 3 * (j - 1) + 2
        idx3 = 3 * (j - 1) + 3
        if (state[idx1]==1 && state[idx2]==0 && state[idx3]==0)
            @inbounds state[idx1] = 0
            @inbounds state[idx2] = p[j]
            @inbounds state[idx3] = !p[j]
        elseif (state[idx1]==0 && state[idx2]==1 && state[idx3]==0)
            @inbounds state[idx1] = p[j]
            @inbounds state[idx2] = 0
            @inbounds state[idx3] = !p[j]
        elseif (state[idx1]==0 && state[idx2]==0 && state[idx3]==1)
            @inbounds state[idx1] = p[j]
            @inbounds state[idx2] = !p[j]
            @inbounds state[idx3] = 0
        end
    end
    return nothing
end

function with_threads(samples, tmax, l, A, q, n0)
    my_data = Total_data(zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1))
    
    freq = 2π * rfftfreq(l)
    Q = freq[q+1]
    Ax = cos.(Q * (1:l))
    rho_in = n0 .+ A * Ax
    
    for _ in samples
        CUDA.device!(rand(0:1))
        conf = rand(Float64, l) .<= rho_in
        rho_0 = sum(conf .* Ax)/l
        @inbounds my_data.m1[1] += rho_0
        @inbounds my_data.m2[1] += rho_0^2
        @inbounds my_data.m3[1] += rho_0^3
        @inbounds my_data.m4[1] += rho_0^4
        for t in 1:tmax
            lb3 = div(l,3)
            off = cu(rand(0:2))
            state_d = cu(circshift(conf, -off))
            p_d = CUDA.rand(Bool,lb3)
            kernel = @cuda launch=false update_d!(state_d,p_d)
            config = launch_configuration(kernel.fun)
            threads = min(convert(Int,l/3), config.threads) #MAX_THREADS_PER_BLOCK
            blocks = cld(convert(Int,l/3), threads)
            CUDA.@sync begin
                kernel(state_d,p_d; threads, blocks)
            end
            conf = Array(circshift(state_d, off))
            rho_t = sum(conf .* Ax)/l
            @inbounds my_data.m1[t+1] += rho_t
            @inbounds my_data.m2[t+1] += rho_t^2
            @inbounds my_data.m3[t+1] += rho_t^3
            @inbounds my_data.m4[t+1] += rho_t^4
        end
        CUDA.reclaim()
    end
    return my_data
end

function sum_multi_thread(samples, tmax, l, A, q, n0)
    chunks = Iterators.partition(1:samples, samples ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn with_threads(chunk,tmax,l, A, q, n0)
    end
    chunk_sums = fetch.(tasks)
    return [sum(map(p->p.m1, chunk_sums))/samples, sum(map(p->p.m2, chunk_sums))/samples, sum(map(p->p.m3, chunk_sums))/samples, sum(map(p->p.m4, chunk_sums))/samples]
end

parsed_args = parse_commandline()
println("Start!")
println("Running 1D non-linear simlation with the following parameters:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end


n0 = parsed_args["eq_den"]
A = parsed_args["amplitude"]
q = parsed_args["wave_number"]
l = parsed_args["chain_length"]
samples = parsed_args["samples"]
tmax = parsed_args["t_max"]
t_samples = parsed_args["t_samples"]

filename = "./data/rho_dotp-n_$n0-A_$A-q_$q-t_samples_$t_samples-samples_each_run_$samples-tmax_$tmax-l-$l.h5"
dataset_names = ["m1", "m2", "m3", "m4"]

h5open(filename, "w") do file
        for i in 1:4
            write(file, dataset_names[i], data[i])
        end
end
println("Done!")
