using CUDA
CUDA.allowscalar(false)
using FFTW
using HDF5
using Random
using Statistics

mutable struct data
    t_st::Matrix{Float64}
    t_curr::Matrix{Float64}
end

mutable struct Conf
    st::BitVector
    curr::Vector{Int}
end

function update_d!(state,p,curr)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    str = gridDim().x * blockDim().x
    for j in idx:str:length(p)
        if (state[3(j-1)+1]==1 && state[3(j-1)+2]==0 && state[3(j-1)+3]==0)
            @inbounds state[3(j-1)+1] = 0
            @inbounds state[3(j-1)+2] = p[j]
            @inbounds state[3(j-1)+3] = !p[j]
            @inbounds curr[3(j-1)+1] = 1
            @inbounds curr[3(j-1)+2] = 1 - p[j]
        elseif (state[3(j-1)+1]==0 && state[3(j-1)+2]==1 && state[3(j-1)+3]==0)
            @inbounds state[3(j-1)+1] = p[j]
            @inbounds state[3(j-1)+2] = 0
            @inbounds state[3(j-1)+3] = !p[j]
            @inbounds curr[3(j-1)+1] = -p[j]
            @inbounds curr[3(j-1)+2] = 1 - p[j]
        elseif (state[3(j-1)+1]==0 && state[3(j-1)+2]==0 && state[3(j-1)+3]==1)
            @inbounds state[3(j-1)+1] = p[j]
            @inbounds state[3(j-1)+2] = !p[j]
            @inbounds state[3(j-1)+3] = 0
            @inbounds curr[3(j-1)+1] = 0 - p[j]
            @inbounds curr[3(j-1)+2] = -1
        else
            @inbounds state[3(j-1)+1] = state[3(j-1)+1]
            @inbounds state[3(j-1)+2] = state[3(j-1)+2]
            @inbounds state[3(j-1)+3] = state[3(j-1)+3]
            @inbounds curr[3(j-1)+1] = 0
            @inbounds curr[3(j-1)+2] = 0
        end
    end
    return nothing
end
function single_upd!(conf)
    l = length(conf.st)
    off = rand(0:2)
    state_d = cu(circshift(conf.st, -off))
    p_d = CUDA.rand(Bool,convert(Int,l/3))
    curr_d = cu(circshift(conf.curr, -off))
    kernel = @cuda launch=false update_d!(state_d,p_d,curr_d)
    config = launch_configuration(kernel.fun)
    threads = min(convert(Int,l/3), config.threads)
    blocks = cld(convert(Int,l/3), threads)
    CUDA.@sync begin
        kernel(state_d,p_d,curr_d; threads, blocks)
    end
    conf.st = circshift(Array(state_d), off)
    conf.curr = circshift(Array(curr_d), off)
    return
end

function with_threads(samples, tmax, l)
    n0 = 0.5
    A = 0.1
    q = convert(Int,floor(l/10))
    freq = 2π * rfftfreq(l)
    Q = freq[q]
    rho_in = n0 .+ A * cos.(Q * collect(1:l))
    norm = sum(rho_in.*rho_in)
    my_data = data(zeros(Float64, l, tmax+1), zeros(Float64, l, tmax))
    for i in samples
        my_conf = Conf(rand(Float64, l) .<= rho_in, zeros(Int, l))
        @inbounds my_data.t_st[:,1] .+= my_conf.st/length(samples)
        for t in 1:tmax
            single_upd!(my_conf)
            @inbounds my_data.t_st[:,t+1] .+= my_conf.st/length(samples)
            @inbounds my_data.t_curr[:,t] .+= my_conf.curr/length(samples)
        end
    end
    return my_data
end

function sum_multi_thread(samples, tmax, l)
    chunks = Iterators.partition(1:samples, samples ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn with_threads(chunk,tmax,l)
    end
    chunk_sums = fetch.(tasks)
    return [mean(map(p->p.t_st, chunk_sums), dims = 1)[1], mean(map(p->p.t_curr, chunk_sums), dims = 1)[1]]
end

l = parse(Int,ARGS[1])
tmax = parse(Int,ARGS[2])
samples = parse(Int,ARGS[3])
total_samples = parse(Int,ARGS[4])

for j in 1:total_samples
    k = j>=2 ? 2 : 1
    x,y = sum_multi_thread(samples, tmax, l)

    if isfile("../data/rho_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5")
        x .+= h5open("../data/rho_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5", "r") do file
                read(file, "A")
        end
    end
    h5open("../data/rho_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5", "w") do file
        write(file, "A", x/k)
    end

    if isfile("../data/curr_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5")
        y .+= h5open("../data/curr_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5", "r") do file
                read(file, "A")
        end
    end
    h5open("../data/curr_samples-($total_samples*$samples)_tmax-($tmax)_l-($l).h5", "w") do file
        write(file, "A", y/k)
    end
end
