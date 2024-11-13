module Utils

export moment_cumulant, parse_commandline

using ArgParse
using Statistics

"""
    moment_cumulant(data::Vector{Matrix{Float64}}) -> Vector{Matrix{Float64}}

Compute the mean (μ), standard deviation (σ), skewness (γ), and kurtosis (κ) from the data.
"""
function moment_cumulant(data::Vector{Matrix{Float64}})
    μ = data[1]
    σ = sqrt.(data[2] .- μ.^2)
    γ = (data[3] .- 3 .* (μ .* σ.^2) .- μ.^3) ./ (σ .^ 3)
    κ = (data[4] .- 4 .* μ .* data[3] .+ 6 .* (μ .* σ) .^ 2 .+ 3 .* μ .^ 4) ./ (σ .^ 4)
    return [μ, σ, γ, κ]
end

"""
    parse_commandline() -> Dict{String, Any}

Parse command-line arguments and return them as a dictionary.
"""
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
            help = "Amplitude of the initial perturbation from equilibrium density at n. n is constrained to be 0 < n < 1"
            arg_type = Float64
            default = 0.1
        "--eq_den", "-n"
            help = "Equilibrium density"
            arg_type = Float64
            default = 0.5
        "--wave_number", "-q"
            help = "Wave number mode of initial perturbation. It is constrained to be 0 < q < l/2"
            arg_type = Int
            default = 9999
        "--update_type", "-u"
            help = "Type of update to use (1: AB-AB, 2: ABC-ABC)"
            arg_type = Int
            default = 2
        "--wait_probability_low_density", "-x"
            help = "Probability of not moving low density triplets"
            arg_type = Float64
            default = 0.2
        "--wait_probability_high_density", "-y"
            help = "Probability of not moving high density triplets"
            arg_type = Float64
            default = 0.8
        "--output_file", "-o"
            help = "Output filename"
            arg_type = String
            default = "output.h5"
    end

    return parse_args(s)
end

end # module
