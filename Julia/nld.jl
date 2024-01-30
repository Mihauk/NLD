using FFTW
using Plots
#using Folds
using Random
using Statistics
#using Distributed
#addprocs(8)
using BenchmarkTools


function trimer_updates(trip)
    if trip == [0, 1, 0]
        return convert(Vector{Bool}, rand([[1, 0, 0], [0, 0, 1]]))
    elseif trip == [0, 0, 1]
        return convert(Vector{Bool}, rand([[1, 0, 0], [0, 1, 0]]))
    elseif trip == [1, 0, 0]
        return convert(Vector{Bool}, rand([[0, 0, 1], [0, 1, 0]]))
    end
    return convert(Vector{Bool},trip)
end

function update(state, off)
    #return circshift(collect(Iterators.flatten(map(trimer_updates, Iterators.partition(circshift(state, -off), 3)))), off)
    return circshift(collect(reduce(vcat,map(trimer_updates, Iterators.partition(circshift(state, -off), 3)))), off)
end


function my_func(samples, tmax, l)
    n0 = 0.5
    A = 0.1
    q = 6
    freq = 2π * rfftfreq(l)
    Q = freq[q]
    rho_in = n0 .+ A * cos.(Q * collect(1:l))
    a = zeros(Bool, samples, l, tmax+1)
    for i in 1:samples
        state = rand(Float64, l) .<= rho_in
        a[i,:,1] = state
        for t in 1:tmax
            state = update(state, rand(0:2))
            a[i,:,t+1] = state
        end
    end
    return a
end

trange = 10 * range(1, length=10)
srange = range(10, length=10, step=10)
lrange = 99 * range(1,stop=3)

x = [[[(@elapsed my_func(sr,tr,lr) )/(lr*tr*sr) for sr in srange] for tr in trange] for lr in lrange]

plot(srange, [x[3][i] for i in 1:length(trange)], label=[10 20 30 40 50 60 70 80 90 100])
xlabel!("samples")
ylabel!("time")
savefig("benchmark_jl.png")