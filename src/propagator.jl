"""
dynamics propagation (Sims-Flanagan Transcription)
"""

using DifferentialEquations
using Plots

abstract type AbstractPropagatorType end

Base.@kwdef struct ODEPropagator <: AbstractPropagatorType
    method=Tsit5()
    reltol::Real=1e-12
    abstol::Real=1e-12
    prob_base::ODEProblem  # base problem should have dynamics.rhs_bcr4bp_with_mass!
end


"""
Propagate trajectory and compute continuity residual
"""
function sf_propagate(
    x::AbstractVector{T},
    n::Int,
    SCparam::Vector,
    Propagator::ODEPropagator,
    param3b::AbstractParameterType,
    LPOArrival::AbstractTerminalType
) where T
    # ... but half of these are defined in r3bp-param FIXME
    mdot, tmax = SCparam  # what is θ0 ?? FIXME also probably need tmax!!!

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
    θ0 = c_launch[5]

    state0 = set_initial_state(param3b, c_launch[2:4], θ0)  
    state_fwd[:,1] = vcat(state0, c_launch[1])  
    sol_fwd = []
    for i in 1:n
        τ, γ, β = tau1[8+3*n:10+3*n]
        params = [param3b.mu2, param3b.mus, param3b.as, θ0, param3b.oml, τ, γ, β, mdot, tmax]
        tspan = [0, tof/2/n]

        # construct ODE problem and solve
        #ode_prob = ODEProblem(dynamics.rhs_bcr4bp_with_mass!, state0, tspan, params)
        _prob = remake(
            Propagator.prob_traj_design;
            tspan = tspan,
            u0 = state0,
            p = params,
        )
        sol = DifferentialEquations.solve(
            _prob, Propagator.method,
            reltol = Propagator.reltol, abstol = Propagator.abstol
        )
        t0 = tspan[2]
        state0 = sol.u[end]
        state_fwd[:,i+1] = state0
        push!(sol_fwd, sol)
        θ0 = θ0 + tof/2/n * oml  # FIXME: oml or oms? I believe it's sidereal ω 
    end

    # backward propagation
    # set up the terminal state
    state_bkwd = zeros(7, n+1)
    θf = c_arr[3]

    statef = set_terminal_state(c_arr[2], θf, param3b, LPOArrival)
    state_bkwd[:,1] = vcat(statef, c_arr[1])
    sol_bkwd = []

    for i in 1:n
        τ, γ, β = tau2[end-3*n+1:end-3*n+3]
        params = [param3b.mu2, param3b.mus, param3b.as, θf, param3b.oml, τ, γ, β, mdot, tmax]  
        tspan = [0, - tof/2/n]

        # construct ODE problem and solve
        #ode_prob = ODEProblem(dynamics.rhs_bcr4bp_with_mass!, statef, tspan, params)
        #sol = solve(ode_prob, method, reltol=reltol, abstol=abstol)
        _prob = remake(
            Propagator.prob_traj_design;
            tspan = tspan,
            u0 = statef,
            p = params,
        )
        sol = DifferentialEquations.solve(
            _prob, Propagator.method,
            reltol = Propagator.reltol, abstol = Propagator.abstol
        )
        tf = tspan[2]
        statef = sol.u[end]
        state_bkwd[:,i+1] = statef
        push!(sol_bkwd, sol)
        θf = θf - tof/2/n * oml  # FIXME: oml or oms? I believe it's sidereal ω
    end

    # take the residual
    res = statef - state0
    #res = matchpoint_residual(state0, statef)   # FIXED: in matchpoint.jl, just in case we need oher res

    return res, state_fwd, state_bkwd, sol_fwd, sol_bkwd
end
