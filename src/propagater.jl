"""
dynamics propagation (Sims-Flanagan Transcription)
"""

using DifferentialEquations
using Plots

include("dynamics.jl")
include("initialize_ode.jl")

μ1 = 0.012150585609624
μ2 = 0.00000303951
as = 388.709677419
ωM = 1  #???
n = 20  # segmentation of the SF transcription per arc

# initial condition
state0 = [8.3203900980615830e-1, 0.0, 1.2603694134707036e-1, 0.0, 2.4012675449219933e-1, 0.0]
t0 = 0

# extract optimization var
c_launch = x[1 : 6]   # Δt1, m0, v∞1, α1, δ1, θ1 
c_arr    = x[7 : 10]  # Δt2, m2, ϕ, θ2
tau1     = x[11 : (length(x)-10)/2+10]  # discretization numbers are the same for the first & second arc
tau2     = x[(length(x)-10)/2+11 : end]
tf       = t0 + c_launch[1] + c_arr[1]

# forward propagation
# set up the initial state
state0 = initialize_ode.set_initial_state(c_launch[3], c_launch[4], c_launch[5], c_launch[6])
for i in 1:n 
    τ, β, γ = tau1[8+3*n:10+3*n]
    params = [μ1, μ2, as, ωM, τ, β, γ]
    tspan = [t0, t0 + c_launch[1]/n]

    # construct ODE problem and solve
    ode_prob = ODEProblem(dynamics.rhs_bcr4bp!, state0, tspan, params)
    sol = solve(ode_prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    push!(t0, tspan[2])
    push!(state0, sol)
end

# backward propagation 
# set up the terminal state
state0 = initialize_ode.set_terminal_state(c_arr[3], c_arr[4])
for i in 1:n
    τ, β, γ = tau1[end-3*n+1:end-3*n+3]
    params = [μ1, μ2, as, ωM, τ, β, γ]
    tspan = [tf, tf - c_arr[1]/n]

    # construct ODE problem and solve
    ode_prob = ODEProblem(dynamics.rhs_bcr4bp!, state0, tspan, params)
    sol = solve(ode_prob, Tsit5(), reltol=1e-12, abstol=1e-12)
    push!(tf, tspan[2])
    push!(state0, sol)
end 


