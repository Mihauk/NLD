using FFTW
using HDF5
using Plots
using Random
using ArgParse
using Statistics
using TimerOutputs
using LinearAlgebra

mutable struct Total_data
    m_rho::Vector{Float64}
    v_rho::Vector{Float64}
    s_rho::Vector{Float64}
    k_rho::Vector{Float64}
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
            help = "Total number of samples"
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
        "--update_type", "-u"
            help = "Type of update to use (1: AB-AB, 2: ABC-ABC)"
            arg_type = Int
            default = 2
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

function with_threads(samples, tmax, l, A, q, n0, u)
    my_data = Total_data(zeros(Float64, (u+1)*tmax+1), zeros(Float64, (u+1)*tmax+1), zeros(Float64, (u+1)*tmax+1), zeros(Float64, (u+1)*tmax+1))
    
    freq = 2π * rfftfreq(l)
    Q = freq[q+1]
    Ax = cos.(Q * (1:l))
    rho_in = n0 .+ A * Ax

    lb3 = div(l,3)

    for _ in samples
        conf = rand(Float64, l) .<= rho_in
        rho_0 = sum(conf .* Ax)/l
        @inbounds my_data.m_rho[1] += rho_0
        @inbounds my_data.v_rho[1] += rho_0^2
        @inbounds my_data.s_rho[1] += rho_0^3
        @inbounds my_data.k_rho[1] += rho_0^4
        for t in 1:tmax
            for k in 0:u
                circshift!(conf, -k)
                p = rand(Bool,lb3) 
                for j in eachindex(p)
                    if (conf[3(j-1)+1]==1 && conf[3(j-1)+2]==0 && conf[3(j-1)+3]==0)
                        @inbounds conf[3(j-1)+1] = 0
                        @inbounds conf[3(j-1)+2] = p[j]
                        @inbounds conf[3(j-1)+3] = !p[j]
                    elseif (conf[3(j-1)+1]==0 && conf[3(j-1)+2]==1 && conf[3(j-1)+3]==0)
                        @inbounds conf[3(j-1)+1] = p[j]
                        @inbounds conf[3(j-1)+2] = 0
                        @inbounds conf[3(j-1)+3] = !p[j]
                    elseif (conf[3(j-1)+1]==0 && conf[3(j-1)+2]==0 && conf[3(j-1)+3]==1)
                        @inbounds conf[3(j-1)+1] = p[j]
                        @inbounds conf[3(j-1)+2] = !p[j]
                        @inbounds conf[3(j-1)+3] = 0
                    end
                end
                circshift!(conf, k)
                rho_t = sum(conf .* Ax)/l
                @inbounds my_data.m_rho[(u+1)*t+k+1-u] += rho_t
                @inbounds my_data.v_rho[(u+1)*t+k+1-u] += rho_t^2
                @inbounds my_data.s_rho[(u+1)*t+k+1-u] += rho_t^3
                @inbounds my_data.k_rho[(u+1)*t+k+1-u] += rho_t^4
            end
        end
    end
    return my_data
end

function sum_multi_thread(samples, tmax, l, A, q, n0, u)
    chunks = Iterators.partition(1:samples, samples ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn with_threads(chunk,tmax,l, A, q, n0, u)
    end
    chunk_sums = fetch.(tasks)
    return [sum(map(p->p.m_rho, chunk_sums))/samples, sum(map(p->p.v_rho, chunk_sums))/samples, sum(map(p->p.s_rho, chunk_sums))/samples, sum(map(p->p.k_rho, chunk_sums))/samples]
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
u = parsed_args["update_type"]

filename = "./data/non_random_rho_dotp-n_$n0-A_$A-q_$q-samples_$samples-tmax_$tmax-l-$l-u-$u.h5"
dataset_names = ["m1", "m2", "m3", "m4"]

# Start timer
to = TimerOutput();

# Calculate and write data to file
@timeit to "Run time" data = moment_cumulant(sum_multi_thread(samples, tmax, l, A, q, n0, u))

println("Writing to file... - ", filename)

# Write data to file
h5open(filename, "w") do file
    for i in 1:4
        write(file, dataset_names[i], data[i])
    end
end
println("Done!")

print_timer(to)