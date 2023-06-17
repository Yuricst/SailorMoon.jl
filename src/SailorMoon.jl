module SailorMoon

using LinearAlgebra
using Plots
using DifferentialEquations
using Roots
using JSON
using BSON
using CSV
using Distributed
using DataFrames
using Printf
using AstrodynamicsBase

import ForwardDiff
import DiffResults
import FiniteDiff

include("plot_func.jl")
include("lpo/build_lpo.jl")
include("integrator.jl")
include("dynamics.jl")

global param3b = dynamics_parameters()

include("deltaV_transcription.jl")
include("transformations.jl")
include("initialize_ode.jl")

include("propagator.jl")
include("integrator.jl")
include("multiple-shooting.jl")
include("fitness.jl")
include("ig_and_bounds.jl")
include("newton_method.jl")

end # module
