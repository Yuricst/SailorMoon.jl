"""
Initial & Terminal condition
"""

"""

Providing the initial state of the SC

# Arguments
    - `param3b`: 3bp parameters
        - `μ2` : m_L / (m_L + m_E)
        - `r0` : distance of Earth center and SC (normalized)
    - `param_vinf' : parameters for defining the Vinf vector
        - `v∞1` : v-infinity at the initial state (w.r.t. the Earth),
        - `long`  : initial longtitude,
        - `lat`  : initial latitude ,
    - `θ1`  : initial E-M line's angle w.r.t. Sun-B1 line
"""

function set_initial_state(param3b, param_vinf, θ1)
    μ2, r0 = param3b[1], param3b[2]
    dv1, long, lat = param_vinf[1], param_vinf[2], param_vinf[3]
    
    # get the initial state of the earth
    ωE = 1
    earth0_in = [μ2*cos(π + θ1), μ2*sin(π+θ1), 0, ωE*μ2*sin(θ1), -ωE*μ2*cos(θ1), 0]

    # SC 
    esc0_in = r0 * [cos(long)*cos(lat), cos(long)*sin(lat), sin(lat)]
    Δv_in  = dv1 * [cos(long)*cos(lat), cos(long)*sin(lat), sin(lat)]
    push!(esc0_in, Δv_in)

    # take a sum
    return earth0_in + esc0_in
end



"""

Providing the terminal state of the SC

# Arguments
    - `ϕ`  : "angle" at the LPO, based on its periodic orbit
    - `θm` : terminal E-M line's angle w.r.t. Sun-B1 line
    - `ω`  : angular velocity of E-M line w.r.t. Sun-B1 line

# assumtion
    The initial velocity direction: the directions s.t. the SC is on the CR3BP invariant manifold...? 
"""

function set_terminal_state(ϕ, θm, ω)
    x0_lpo, period_lpo, monodromy = duhduhduh_function()

    # propagate the periodic orbit until ϕT.
    prob_ = ODEProblem(rhs_pcr3bp_sv!, x0_lpo, (0, ϕ*period_lpo), (μ));
    sol = solve(prob_, Tsit5(), reltol=1e-11, abstol=1e-11);
    x_t = sol.u[end][1:6]

    # add varation based on the eigenvector 
    ϵ = 1e-6  # koshiki_nanndakke()

    # get eigenvector for the stable manifold
    v_stb = get_eigenvector(monodromy, stable::True)

    state_f = x_t + ϵ * v_stb * norm(v_stb)

    # coodinate transformation
    state_f_SunB1 = transform_EMrot_to_SunB1(state_f, θm, -ω)

    return state_f_SunB1
end

    
