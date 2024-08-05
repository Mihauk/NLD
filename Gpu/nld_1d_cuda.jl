using CUDA
CUDA.allowscalar(false)
using FFTW
using HDF5
using Random
using ArgParse
using Statistics
using TimerOutputs
using LinearAlgebra

# Define the mutable struct that stores the data
mutable struct Total_data
    m1::Vector{Float64}
    m2::Vector{Float64}
    m3::Vector{Float64}
    m4::Vector{Float64}
end

# Function to parse command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--chain_length", "-l"
            help = "Length of the 1D chain"
            arg_type = Int
            default = 99999
        "--t_max", "-t"
            help = "Number of time steps"
            arg_type = Int
            default = 100
        "--samples", "-s"
            help = "Total number of samples"
            arg_type = Int
            default = 10000
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
            default = 9999
        #"--flag1"
        #    help = "an option without argument, i.e. a flag"
        #    action = :store_true
        #"arg1"
        #    help = "a positional argument"
        #    required = true
    end

    return parse_args(s)
end

# Function to calculate the cumulants from the moments
function moment_cumulant(data)
    μ = data[1]
    σ = sqrt.(data[2] .- μ.^2)
    γ = (data[3] .- 3*(μ.*(σ.^2)) .- μ.^3)./(σ.^3)
    κ = (data[4] .- 4*(μ.*data[3]) .+ 6*(μ.*σ).^2 .+ 3*μ.^4)./(σ.^4)
    return [μ, σ, γ, κ]
end

# Function to update the state based on certain conditions also the CUDA kernel
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

# Function to perform computations with multi-threading
function with_threads(samples, tmax, l, A, q, n0)
    
    # Calculate frequency and wave number and intital distribution to sample from
    freq = 2π * rfftfreq(l)
    Q = freq[q+1]
    Ax = cos.(Q * (1:l))
    rho_in = n0 .+ A * Ax
    
    # Initialize data structure to store results
    my_data = Total_data(zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1))

    # Loop over samples
    for _ in samples

        CUDA.device!(rand(0:1)) # Set the CUDA device randomly

        conf = rand(Float64, l) .<= rho_in  # Generate initial configuration
        rho_0 = sum(conf .* Ax)/l

        # Update data structure with initial configuration
        @inbounds my_data.m1[1] += rho_0
        @inbounds my_data.m2[1] += rho_0^2
        @inbounds my_data.m3[1] += rho_0^3
        @inbounds my_data.m4[1] += rho_0^4

        # Loop over time steps
        for t in 1:tmax
            lb3 = div(l,3)
            off = cu(rand(0:2))
            state_d = cu(circshift(conf, -off))
            p_d = CUDA.rand(Bool,lb3)

            # Define and launch CUDA kernel
            kernel = @cuda launch=false update_d!(state_d,p_d)
            config = launch_configuration(kernel.fun)
            threads = min(convert(Int,l/3), config.threads) #MAX_THREADS_PER_BLOCK
            blocks = cld(convert(Int,l/3), threads)
            CUDA.@sync begin
                kernel(state_d,p_d; threads, blocks)
            end

             # Update configuration and calculate new density
            conf = Array(circshift(state_d, off))
            rho_t = sum(conf .* Ax)/l

            # Update data structure with new configuration
            @inbounds my_data.m1[t+1] += rho_t
            @inbounds my_data.m2[t+1] += rho_t^2
            @inbounds my_data.m3[t+1] += rho_t^3
            @inbounds my_data.m4[t+1] += rho_t^4
        end
        CUDA.reclaim()
    end
    return my_data
end

# Function to sum the results from all threads
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

filename = "./data/rho_dotp-n_$n0-A_$A-q_$q-samples_$samples-tmax_$tmax-l-$l.h5"
dataset_names = ["m1", "m2", "m3", "m4"]

# Start timer
to = TimerOutput();

# Calculate and write data to file
@timeit to "Run time" data = moment_cumulant(sum_multi_thread(samples, tmax, l, A, q, n0))

println("Writing to file... - ", filename)

# Write data to file
h5open(filename, "w") do file
    for i in 1:4
        write(file, dataset_names[i], data[i])
    end
end
println("Done!")

print_timer(to)
