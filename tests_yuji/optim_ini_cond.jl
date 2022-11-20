"""
Given initial guess, earth θf and ϕ0 that matches the initial condition at the LEO.
x = [ϕ0, θf]
"""

#push!(LOAD_PATH,"../src/")
using LinearAlgebra
using DifferentialEquations
using Plots
import ForwardDiff
using ForwardDiff: Dual
import DiffResults
using AstrodynamicsBase
using joptimise
using Printf
using DataFrames
using JSON
using CSV

include("../../julia-r3bp/R3BP/src/R3BP.jl")
include("../src/SailorMoon.jl")   # relative path to main file of module

param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

## set up of initial condition (Lyapunov orbit)
lp = 2
Az_km = 1200.0
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
# for Halo propagation, keep the tol as tight as possible 
sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12) #, saveat=LinRange(0, period, n+1))
monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
ys0 = R3BP.get_eigenvector(monodromy, true, 1) # monodromy eigenvector

ϵr = 1e-6
ϵv = 1e-6

LPOArrival = SailorMoon.CR3BPLPO2(
    res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12, 0.005
);

# initialize ODE 
svf_ = zeros(Float64, 1, 7)
params = [param3b.mu2, param3b.mus, 0.0, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0]
tspan = [0, -10.0]
prob = ODEProblem(R3BP.rhs_bcr4bp_thrust!, svf_, tspan, params)

function objective(x::AbstractVector{T}) where T
    # using x, generate the trajectory 
    # ϕ0  = ForwardDiff.value(x[1])
    # θf  = ForwardDiff.value(x[2])
    # tof = ForwardDiff.value(x[3])
    ϕ0, θf, tof = x[1], x[2], x[3]

    # println("test test test: ", ϕ0, " ", θf)
    xf = vcat(SailorMoon.set_terminal_state2(ϕ0, θf, param3b, LPOArrival), 1.0)

    println("xf: ", xf)

    ## make ensemble problems
    remake(prob, u0=xf, p=[param3b.mu2, param3b.mus, θf, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0],
                tspan=[0, -tof])

    cbs = CallbackSet()  # empty as of now, maybe we should add something? (periapsis?)
    sol = solve(prob, Tsit5(), callback=cbs, reltol=1e-12, abstol=1e-12);

    # this will allow the trajectories s.t. h = 200km but with really bad flight angle, need to enforce some condition?
    rp_target = (200 + 6375) / param3b.lstar  # altitude of LEO 
    rp = sqrt(sol.u[end][1]^2 + sol.u[end][1]^2 + sol.u[end][1]^2)
    f = (rp - rp_target)^2 * 10000
    return f
end


function obj_con(g::AbstractVector{T}, x::AbstractVector{T}) where T
    # compute objective
    f = objective(x)
    g[1] = 0.0
    return f
end

# initial guess: ϕ0[rad], θf[rad], tof[T]
x0 = [5.27787565803085; 0.125663706; 18.8977945731596]

# bounds on variables
lx = [5.0; -0.1; 17.5]
ux = [5.6; 0.4; 20.5]

# bounds on constriants
lg = [0.0]
ug = [0.0]
# number of constraints
ng = 1

## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "tol" => 1e-12
)

# currently using automatic differentiation
xopt, fopt, info = minimize(obj_con, x0, ng; lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt", options=ip_options);

println("Done with IPOPT!")
println(info)
println(xopt)