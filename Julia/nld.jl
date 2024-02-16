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

function update!(config,p)
    for j in eachindex(p)
        if (config.st[3(j-1)+1]==1 && config.st[3(j-1)+2]==0 && config.st[3(j-1)+3]==0)
            @inbounds config.st[3(j-1)+1] = 0
            @inbounds config.st[3(j-1)+2] = p[j]
            @inbounds config.st[3(j-1)+3] = !p[j]
            @inbounds config.curr[3(j-1)+1] = 1
            @inbounds config.curr[3(j-1)+2] = 1 - p[j]
        elseif (config.st[3(j-1)+1]==0 && config.st[3(j-1)+2]==1 && config.st[3(j-1)+3]==0)
            @inbounds config.st[3(j-1)+1] = p[j]
            @inbounds config.st[3(j-1)+2] = 0
            @inbounds config.st[3(j-1)+3] = !p[j]
            @inbounds config.curr[3(j-1)+1] = -p[j]
            @inbounds config.curr[3(j-1)+2] = 1 - p[j]
        elseif (config.st[3(j-1)+1]==0 && config.st[3(j-1)+2]==0 && config.st[3(j-1)+3]==1)
            @inbounds config.st[3(j-1)+1] = p[j]
            @inbounds config.st[3(j-1)+2] = !p[j]
            @inbounds config.st[3(j-1)+3] = 0
            @inbounds config.curr[3(j-1)+1] = 0 - p[j]
            @inbounds config.curr[3(j-1)+2] = -1
        end
    end
    return
end

function single_upd!(conf)
    l = length(conf.st)
    off = rand(0:2)
    roll_conf = Conf(circshift(conf.st, -off), circshift(conf.curr, -off))
    p = rand(Bool,convert(Int,l/3)) 
    update!(roll_conf, p)
    conf.st = circshift(roll_conf.st, off)
    conf.curr = circshift(roll_conf.curr, off)
    return
end

function with_threads(samples, tmax, l)
    n0 = 0.5
    A = 0.1
    q = convert(Int,floor(l/10))
    freq = 2π * rfftfreq(l)
    Q = freq[q]
    rho_in = n0 .+ A * cos.(Q * collect(1:l))
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
