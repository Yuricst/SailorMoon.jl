"""
Initial & Terminal condition
"""

"""

Providing the initial state of the SC

# Arguments
    - `param3b::dynamics_params`: 3bp parameters
    - `param_vinf' : parameters for defining the Vinf vector
        - `dv1` : delta-v at the initial state (w.r.t. the Earth),
        - `long` : initial longtitude,
        - `lat`  : initial latitude ,
    - `θ1`  : initial E-M line's angle w.r.t. Sun-B1 line
"""

function set_initial_state(param3b::dynamics_params, param_vinf, θ1)
    dv1, long, lat = param_vinf[1], param_vinf[2], param_vinf[3]

    # get the initial state of the earth
    ωE = param3b.oml
    earth0_in = [
        param3b.mu2*cos(π + θ1), param3b.mu2*sin(π+θ1), 0,
        ωE*param3b.mu2*sin(θ1), -ωE*param3b.mu2*cos(θ1), 0
    ]      # planar

    # SC position
    r_E_sc = param3b.r_park * [cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)]
    # SC circular velocity
    v_E_sc = param3b.v_park * [sin(long), cos(long), 0]
    Δv_in  = dv1 * [sin(long), cos(long), 0]  #[cos(long)*cos(lat), cos(long)*sin(lat), sin(lat)]

    # take sum of Earth position + parking orbit + delta-V
    return earth0_in + vcat(r_E_sc, v_E_sc)[:] + vcat([0,0,0], Δv_in)[:]
end


abstract type AbstractTerminalType end

Base.@kwdef struct CR3BPLPO <: AbstractTerminalType
    x0::Vector        # initial state
    period::Real    # period
    ys0::Vector       # stable eigenvector at x0
    prob_cr3bp_stm::ODEProblem
    ϵ::Real
    method=Tsit5()
    reltol::Real=1e-12
    abstol::Real=1e-12
end

Base.@kwdef struct CR3BPLPO2 <: AbstractTerminalType
    x0::Vector        # initial state
    period::Real    # period
    ys0::Vector       # stable eigenvector at x0
    prob_cr3bp_stm::ODEProblem
    ϵr::Real   # separate epsilon
    ϵv::Real   # separate epsilon
    method=Tsit5()
    reltol::Real=1e-12
    abstol::Real=1e-12
end


"""
    set_terminal_state(ϕ, θm, param3b::AbstractParameterType, LPOArrival::CR3BPLPO)

Providing the terminal state of the SC based on arrival to manifold.

# Arguments
    - `ϕ`: "angle" at the LPO, based on its periodic orbit
    - `θm`: terminal E-M line's angle w.r.t. Sun-B1 line
    - `param3b::AbstractParameterType`: angular velocity of E-M line w.r.t. Sun-B1 line
    - `LPOArrival::CR3BPLPO`: arrival periodic orbit object

# assumtion
    The initial velocity direction: the directions s.t. the SC is on the CR3BP invariant manifold...?
"""
function set_terminal_state(ϕ, param3b::AbstractParameterType, LPOArrival::CR3BPLPO)
    # propagate the periodic orbit until ϕT.
    x0_stm = vcat(LPOArrival.x0, reshape(I(6), (36,)))[:]
    _prob = remake(
        LPOArrival.prob_cr3bp_stm;
        tspan = (0.0, ϕ * LPOArrival.period),
        u0 = x0_stm,
        p=[param3b.mu2]
    )
    sol = DifferentialEquations.solve(
        _prob, LPOArrival.method,
        reltol = LPOArrival.reltol, abstol = LPOArrival.abstol
    )
    x_tf = sol.u[end][1:6]
    stm = transpose(reshape(sol.u[end][7:end], (6, 6)))

    # translate stable eigenvector and perturb final state
    ys = stm * LPOArrival.ys0
    if ys[1] > 0  # always perturb outward
        state_f = x_tf + LPOArrival.ϵ * ys/norm(ys)
    else
        state_f = x_tf - LPOArrival.ϵ * ys/norm(ys)
    end
    return state_f
end


"""
    set_terminal_state(ϕ, θm, param3b::AbstractParameterType, LPOArrival::CR3BPLPO)

Providing the terminal state of the SC based on arrival to manifold.

# Arguments
    - `ϕ`: "angle" at the LPO, based on its periodic orbit
    - `θm`: terminal E-M line's angle w.r.t. Sun-B1 line
    - `param3b::AbstractParameterType`: angular velocity of E-M line w.r.t. Sun-B1 line
    - `LPOArrival::CR3BPLPO`: arrival periodic orbit object

# assumtion
    The initial velocity direction: the directions s.t. the SC is on the CR3BP invariant manifold...?
"""
function set_terminal_state(ϕ, θm, param3b::AbstractParameterType, LPOArrival::CR3BPLPO)
    # propagate the periodic orbit until ϕT.
    x0_stm = vcat(LPOArrival.x0, reshape(I(6), (36,)))[:]
    _prob = remake(
        LPOArrival.prob_cr3bp_stm;
        tspan = (0.0, sign(ϕ)*mod(abs(ϕ),1) * LPOArrival.period),
        u0 = x0_stm,
        p=[param3b.mu2]
    )
    sol = DifferentialEquations.solve(
        _prob, LPOArrival.method,
        reltol = LPOArrival.reltol, abstol = LPOArrival.abstol
    )
    x_tf = sol.u[end][1:6]
    stm = transpose(reshape(sol.u[end][7:end], (6, 6)))

    # translate stable eigenvector and perturb final state
    ys = stm * LPOArrival.ys0
    if ys[1] > 0  # always perturb outward
        state_f = x_tf + LPOArrival.ϵ * ys/norm(ys)
    else
        state_f = x_tf - LPOArrival.ϵ * ys/norm(ys)
    end

    # coodinate transformation
    state_f_SunB1 = transform_EMrot_to_SunB1(state_f, π-θm, param3b.oms)  # FIXME is θs appropriate?

    return state_f_SunB1
end


"""
    set_terminal_state(ϕ, θm, param3b::AbstractParameterType, LPOArrival::CR3BPLPO)

Providing the terminal state of the SC based on arrival to manifold.

# Arguments
    - `ϕ`: "angle" at the LPO, based on its periodic orbit
    - `θm`: terminal E-M line's angle w.r.t. Sun-B1 line
    - `param3b::AbstractParameterType`: angular velocity of E-M line w.r.t. Sun-B1 line
    - `LPOArrival::CR3BPLPO`: arrival periodic orbit object

# assumtion
    The initial velocity direction: the directions s.t. the SC is on the CR3BP invariant manifold...?
"""
function set_terminal_state2(ϕ, θm, param3b::AbstractParameterType, LPOArrival::CR3BPLPO2)
    # propagate the periodic orbit until ϕT.
    x0_stm = vcat(LPOArrival.x0, reshape(I(6), (36,)))[:]
    _prob = remake(
        LPOArrival.prob_cr3bp_stm;
        tspan = (0.0, sign(ϕ)*mod(abs(ϕ),1) * LPOArrival.period),
        u0 = x0_stm,
        p=[param3b.mu2]
    )
    sol = DifferentialEquations.solve(
        _prob, LPOArrival.method,
        reltol = LPOArrival.reltol, abstol = LPOArrival.abstol
    )
    x_tf = sol.u[end][1:6]
    stm = transpose(reshape(sol.u[end][7:end], (6, 6)))

    # translate stable eigenvector and perturb final state
    ys = stm * LPOArrival.ys0
    y = ys/norm(ys)
    # println(vcat(LPOArrival.ϵr * y[1:3], LPOArrival.ϵv * y[4:6]))

    if ys[1] > 0  # always perturb outward
        state_f = x_tf + vcat(LPOArrival.ϵr * y[1:3], LPOArrival.ϵv * y[4:6])
    else
        state_f = x_tf - vcat(LPOArrival.ϵr * y[1:3], LPOArrival.ϵv * y[4:6])
    end

    return state_f
end
