module Config

export Conf, TotalData, initQ

using Random

# Structure representing the configuration state
mutable struct Conf
    st::BitVector
    curr::Union{Vector{Int}, Nothing}  # Optional current field
end

# Structure for storing simulation data
mutable struct TotalData
    m1::Matrix{Float64}
    m2::Matrix{Float64}
    m3::Matrix{Float64}
    m4::Matrix{Float64}
end

"""
    initQ(L::Int, nQ::Int, n0::Float64, A::Float64, compute_current::Bool) -> Conf

Initialize the configuration `st` as a BitVector of length `L`, where each bit represents the presence (1) or absence (0) of a particle. The initial configuration is perturbed from the equilibrium density `n0` with amplitude `A` and wave number `nQ`. If `compute_current` is `true`, initializes `curr` as a zero vector.
"""
function initQ(L::Int, nQ::Int, n0::Float64, A::Float64, compute_current::Bool)
    Q = 2 * π * nQ / L
    st = [x >= 0 ? true : false for x in n0 .+ A .* cos.(Q .* (0:L-1)) .+ (rand(L) .- 1)]
    curr = compute_current ? zeros(Int, L) : nothing
    return Conf(st, curr)
end

end # module
