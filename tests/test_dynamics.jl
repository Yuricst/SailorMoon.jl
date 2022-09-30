"""
Test for dynamics
"""

using DifferentialEquations
using Plots

include("../src/SailorMoon.jl")   # relative path to main file of module

mu = 0.012150585609624
state0 = [8.3203900980615830e-1, 0.0, 1.2603694134707036e-1, 0.0, 2.4012675449219933e-1, 0.0]
tspan = [0, 2.7824079054366879]
params = [mu,]

# construct ODE problem and solve
ode_prob = ODEProblem(SailorMoon.rhs_cr3bp!, state0, tspan, params)
sol = solve(ode_prob, Tsit5(), reltol=1e-12, abstol=1e-12)

# plot 
ptraj = plot(frame_style=:box, set_aspect=:equal)
plot!(ptraj, sol; vars=(1,2,3))

display(ptraj)