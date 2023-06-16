using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
using Printf
using JSON
using CSV
using DataFrames
using BenchmarkTools

include("../../julia-r3bp/R3BP/src/R3BP.jl")
include("../src/SailorMoon.jl")   # relative path to main file of module

### PARAMETERS ###################################
# csv file to load the initial solution
filename = "../run/data/grid_search_Tsit5_0525_EMrotThrust.csv"
# dv_dir function corresponding to the csv file 
dir_func = SailorMoon.dv_EMrotdir_sb1frame 

# parameters
param3b = SailorMoon.dynamics_parameters()
paramMulti = SailorMoon.multi_shoot_parameters(param3b)

n_arc = 5
##################################################

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess
df = DataFrame(CSV.File(filename))

# maybe want to use "for row in eachrow(df)" to automate the process...? 
row = df[21,:]

x0, lx, ux = SailorMoon.make_ig_bounds2(row, τ_ig, paramMulti.n_arc)
fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2(dir_func, paramMulti, x0)

xs, fs, converged = SailorMoon.differential_correction2(x0, eval_sft)

