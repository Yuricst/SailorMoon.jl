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


function sf_propagate(
    x::AbstractVector{Float64},
    p::Vector{Float64},
    n::Int,
    Propagator::ODEPropagator,
    SFset::SimsFlanaganSettings,
    param3b::Vector{Float64}
)

    # FIXME unpack parameters
    method = Propagator.method
    reltol = Propagator.reltol
    abstol = Propagator.abstol
    dt = Propagator.dt

    μ2, μS, as, θ0, ωM, mdot = p

    # initial condition
    t0 = 0

    # extract optimization var
    tof      = x[1]
    c_launch = x[2 : 6]  # m0, v∞1, α1, δ1, θ1 
    c_arr    = x[7 : 9]  # m2, ϕ, θ2
    tau1     = x[10 : (length(x)-9)/2 + 9]  # discretization numbers are the same for the first & second arc
    tau2     = x[(length(x)-9)/2 + 10 : end]
    tf       = t0 + c_launch[1] + c_arr[1]

    # forward propagation
    # set up the initial state
    state_fwd = zeros(7, n+1)
    param3b = [μ2,r0]
    state0 = initialize_ode.set_initial_state(param3b, c_launch[3:5], c_launch[6])
    state_fwd[:,1] = vcat(state0, c_launch[0])
    sol_fwd = []
    for i in 1:n 
        τ, γ, β = tau1[8+3*n:10+3*n]
        params = [μ2, μS, as, θ0, ωM, τ, γ, β, mdot]
        tspan = [t0, t0 + tof/2/n]

        # construct ODE problem and solve
        ode_prob = ODEProblem(dynamics.rhs_bcr4bp_with_mass!, state0, tspan, params)
        sol = solve(ode_prob, method, reltol=reltol, abstol=abstol)
        t0 = tspan[2]
        state0 = sol.u[end]
        state_fwd[:,i+1] = state0
        push!(sol_fwd, sol)
    end

    # backward propagation 
    # set up the terminal state
    state_bkwd = zeros(7, n+1)
    statef = initialize_ode.set_terminal_state(c_arr[3], c_arr[4], ωM)
    state_bkwd[:,1] = vcat(statef, c_arr[0])
    sol_bkwd = []

    for i in 1:n
        τ, γ, β = tau2[end-3*n+1:end-3*n+3]
        params = [μ2, μS, as, θ0, ωM, τ, γ, β, mdot]
        tspan = [tf, tf - tof/2/n]

        # construct ODE problem and solve
        ode_prob = ODEProblem(dynamics.rhs_bcr4bp_with_mass!, statef, tspan, params)
        sol = solve(ode_prob, method, reltol=reltol, abstol=abstol)
        tf = tspan[2]
        statef = sol.u[end]
        state_bkwd[:,i+1] = statef
        push!(sol_bkwd, sol)
    end 

    # take the residual 
    res = get_residual(state0, statef)
    
    return res, state_fwd, state_bkwd, sol_fwd, sol_bkwd
end
