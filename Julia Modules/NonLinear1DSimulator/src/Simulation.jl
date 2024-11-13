module Simulation

export with_threads, sum_multi_thread, run_simulation

using .Config
using .Update
using .Utils
using Random
using FFTW
using TimerOutputs
using Base.Threads

"""
    with_threads(samples::UnitRange{Int}, tmax::Int, l::Int, A::Float64, q::Int, n0::Float64, u::Int, mpl::Float64, mph::Float64) -> TotalData

Run simulations with threading over the given `samples`, for `tmax` time steps, chain length `l`, amplitude `A`, wave number `q`, equilibrium density `n0`, update type `u`, and waiting probabilities `mpl` and `mph`.
"""
function with_threads(samples::UnitRange{Int}, tmax::Int, l::Int, A::Float64, q::Int, n0::Float64, u::Int, mpl::Float64, mph::Float64)
    qb2 = Int(ceil(q / 2))
    num_timesteps = (u + 1) * tmax + 1
    my_data = TotalData(
        zeros(Float64, qb2 - 1, num_timesteps),
        zeros(Float64, qb2 - 1, num_timesteps),
        zeros(Float64, qb2 - 1, num_timesteps),
        zeros(Float64, qb2 - 1, num_timesteps)
    )

    for _ in samples
        my_conf = Conf(initQ(l, q - 1, n0, A))
        tmp_0 = fft(my_conf.st) / sqrt(l)
        # Initial data collection
        for k in 1:qb2 - 1
            rho_0 = real(tmp_0[k + 1] * tmp_0[q - k + 1])
            @inbounds my_data.m1[k, 1] += rho_0
            @inbounds my_data.m2[k, 1] += rho_0^2
            @inbounds my_data.m3[k, 1] += rho_0^3
            @inbounds my_data.m4[k, 1] += rho_0^4
        end
        # Time evolution
        for t in 1:tmax
            for d in 0:u
                single_upd!(my_conf, d, mpl, mph)
                tmp_t = fft(my_conf.st) / sqrt(l)
                idx = (u + 1) * t + d + 1 - u
                for k in 1:qb2 - 1
                    rho_t = real(tmp_t[k + 1] * tmp_t[q - k + 1])
                    @inbounds my_data.m1[k, idx] += rho_t
                    @inbounds my_data.m2[k, idx] += rho_t^2
                    @inbounds my_data.m3[k, idx] += rho_t^3
                    @inbounds my_data.m4[k, idx] += rho_t^4
                end
            end
        end
    end
    return my_data
end

"""
    sum_multi_thread(samples::Int, tmax::Int, l::Int, A::Float64, q::Int, n0::Float64, u::Int, mpl::Float64, mph::Float64) -> Vector{Matrix{Float64}}

Run the simulations using multithreading, splitting the `samples` among available threads.
"""
function sum_multi_thread(samples::Int, tmax::Int, l::Int, A::Float64, q::Int, n0::Float64, u::Int, mpl::Float64, mph::Float64)
    num_threads = nthreads()
    chunk_size = samples ÷ num_threads
    chunks = [((i - 1) * chunk_size + 1):min(i * chunk_size, samples) for i in 1:num_threads]
    tasks = map(chunks) do chunk
        Threads.@spawn with_threads(chunk, tmax, l, A, q, n0, u, mpl, mph)
    end
    chunk_sums = fetch.(tasks)
    total_samples = samples
    m1_sum = sum(map(p -> p.m1, chunk_sums))
    m2_sum = sum(map(p -> p.m2, chunk_sums))
    m3_sum = sum(map(p -> p.m3, chunk_sums))
    m4_sum = sum(map(p -> p.m4, chunk_sums))
    return [m1_sum / total_samples, m2_sum / total_samples, m3_sum / total_samples, m4_sum / total_samples]
end

"""
    run_simulation(params::Dict{String, Any}) -> Nothing

Run the simulation with the given parameters and save the output to an HDF5 file.
"""
function run_simulation(params::Dict{String, Any})
    n0 = params["eq_den"]
    A = params["amplitude"]
    q = params["wave_number"]
    l = params["chain_length"]
    samples = params["samples"]
    tmax = params["t_max"]
    u = params["update_type"]
    mpl = params["wait_probability_low_density"]
    mph = params["wait_probability_high_density"]

    filename = params["output_file"]
    dataset_names = ["m1", "m2", "m3", "m4"]

    # Start timer
    to = TimerOutput()

    # Run simulation and compute moments
    @timeit to "Run time" data = moment_cumulant(sum_multi_thread(samples, tmax, l, A, q, n0, u, mpl, mph))

    println("Writing to file... - ", filename)

    # Write data to HDF5 file
    h5open(filename, "w") do file
        for i in 1:4
            write(file, dataset_names[i], data[i])
        end
    end
    println("Done!")

    print_timer(to)
end

end # module
