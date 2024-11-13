module NonLinear1DSimulator

using FFTW
using HDF5
using ArgParse
using TimerOutputs

include("Config.jl")
include("Update.jl")
include("Simulation.jl")
include("Utils.jl")

using .Config
using .Update
using .Simulation
using .Utils

export run_simulation, parse_commandline

end # module
