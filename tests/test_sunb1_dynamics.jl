"""
Minimal working test for Sun-B1 dynamics
"""

using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
import joptimise
using Printf

gr()

using R3BP
include("../src/SailorMoon.jl")   # relative path to main file of module


param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

# construct params
t_sidereal = 27.3217   * 86400      # sec
t_synodic  = 29.530588 * 86400      # sec
ωb = 2*pi*(t_synodic - t_sidereal)/t_synodic/t_sidereal*param3b.tstar
ωM = 2*pi / (t_synodic / param3b.tstar)
θM = π  # placeholder
params = [param3b.mu2, param3b.mus, param3b.as, θM, ωM, ωb, 0.0, 0.0, 0.0, 0.0, 0.0]

# initial state
u0_test = [
     387.6985504060269
  -2.7938137447835705e-7
  -5.620290959516347e-23
   1.0294967930424787e-6
  -1.203809978563231
   1.3395041586755617e-31
   1.0
]
_prob_base = ODEProblem(SailorMoon.rhs_bcr4bp_sb1frame2!, u0_test, [0,-23.5], params)

# propagate
sol = SailorMoon.integrate_rk4(_prob_base, 0.005);

# plot
function plot_circle(radius, x, y, n=100)
    circle = zeros(2,n)
    thetas = LinRange(0.0, 2π, n)
    for i = 1:n
        circle[1,i] = radius*cos(thetas[i]) + x
        circle[2,i] = radius*sin(thetas[i]) + y
    end
    return circle
end
moon = plot_circle(1-param3b.mu2, param3b.as , 0.0)
earth = plot_circle(param3b.mu2, param3b.as, 0.0)
moon_soi_outer = plot_circle(1-param3b.mu2+66000/param3b.lstar, param3b.as, 0.0)

pcart = plot(hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], set_aspect=:equal, size=(700,700))
plot!(pcart, earth[1,:], earth[2,:], c=:green, lw=1.0, label="earth")
plot!(pcart, moon[1,:], moon[2,:], c=:orange, lw=1.0, label="moon")
plot!(pcart, moon_soi_outer[1,:], moon_soi_outer[2,:], c=:black, lw=1.0, label="moon_soi_outer")
display(pcart)
