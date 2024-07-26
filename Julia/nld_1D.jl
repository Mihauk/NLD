using FFTW
using HDF5
using Bumper
using Random
using ArgParse
using Statistics
using LinearAlgebra
#using BenchmarkTools

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
        "--t_samples", "-T"
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


function with_threads(samples, tmax, l, A, q, n0, buf)
    my_data = Total_data(zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1))
    @no_escape buf begin
        freq = @alloc(Float64, div(l,2)+1)
        Q = @alloc(Float64, 1)
        Ax = @alloc(Float64, l)
        rho_in = @alloc(Float64, l)
        freq = 2π * rfftfreq(l)
        Q = freq[q]
        Ax = cos.(Q * (1:l))
        rho_in = n0 .+ A * Ax
        for _ in samples
            @no_escape begin
                conf = @alloc(Bool, l)
                rho_0 = @alloc(Float64, 1)
                conf = rand(Float64, l) .<= rho_in
                rho_0 = sum(conf .* Ax)/l
                @inbounds my_data.m_rho[1] += rho_0
                @inbounds my_data.v_rho[1] += rho_0^2
                @inbounds my_data.s_rho[1] += rho_0^3
                @inbounds my_data.k_rho[1] += rho_0^4
                for t in 1:tmax
                    @no_escape begin
                        lb3 = @alloc(Int, 1)
                        lb3 = convert(Int,l/3)
                        off = @alloc(Int, 1)
                        p = @alloc(Bool, lb3)
                        off = rand(0:2)
                        circshift!(conf, -off)
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
                        circshift!(conf, off)
                        nothing
                    end
                    rho_t = sum(conf .* Ax)/l
                    @inbounds my_data.m_rho[t+1] += rho_t
                    @inbounds my_data.v_rho[t+1] += rho_t^2
                    @inbounds my_data.s_rho[t+1] += rho_t^3
                    @inbounds my_data.k_rho[t+1] += rho_t^4
                end
                nothing
            end
        end
        nothing
    end
    return my_data
end

function sum_multi_thread(samples, tmax, l, A, q, n0)
    chunks = Iterators.partition(1:samples, samples ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn with_threads(chunk,tmax,l, A, q, n0, default_buffer())
    end
    chunk_sums = fetch.(tasks)
    return [sum(map(p->p.m_rho, chunk_sums))/samples, sum(map(p->p.v_rho, chunk_sums))/samples, sum(map(p->p.s_rho, chunk_sums))/samples, sum(map(p->p.k_rho, chunk_sums))/samples]
end

parsed_args = parse_commandline()
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
dataset_names = ["m_rho", "v_rho", "s_rho", "k_rho"]

println("Start!")

for outer in 1:t_samples
    # Initialize or read previous data
    if outer == 1
        previous_data = [zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1), zeros(Float64, tmax+1)] # or some initial values
    else
        h5open(filename, "r") do file
            previous_data = [read(file[dataset_names[i]]) for i in 1:4]
        end
    end

    data = sum_multi_thread(samples, tmax, l, A, q, n0) # your computations here   
    # Save the data at the end of the inner loop
    h5open(filename, "w") do file
        for i in 1:4
            write(file, dataset_names[i], previous_data[i] + data[i])
        end
    end
end

println("Done!")
