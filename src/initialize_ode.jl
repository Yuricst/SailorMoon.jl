"""
Initial & Terminal condition
"""

"""

Providing the initial state of the SC

# Arguments
    - `param3b`: 3bp parameters
        - `μ` : m_L / (m_L + m_E)
        - `r0` : distance of Earth center and SC (normalized)
    - `param_vinf' : parameters for defining the Vinf vector
        - `v∞1` : v-infinity at the initial state (w.r.t. the Earth),
        - `long`  : initial longtitude,
        - `lat`  : initial latitude ,
    - `θ1`  : initial E-M line's angle w.r.t. Sun-B1 line
"""

function set_initial_state(param3b, param_vinf, θ1)
    μ, r0 = param3b[1], param3b[2]
    dv1, long, lat = param_vinf[1], param_vinf[2], param_vinf[3]
    
    # get the initial state of the earth
    ωE = 1
    earth0_in = [μ*cos(π + θ1), μ*sin(π+θ1), 0, ωE*μ*sin(θ1), -ωE*μ*cos(θ1), 0]

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
    - `θ2` : terminal E-M line's angle w.r.t. Sun-B1 line

# assumtion
    The initial velocity direction: the directions s.t. the SC is on the CR3BP invariant manifold...? 
"""

function set_terminal_state(ϕ, θ2)
    x0_lpo, period_lpo, monodromy = duh_function()

    # propagate the periodic orbit until ϕT.
    x_t, stm_t = propagate(x0_lpo, ϕ*period_lpo)

    # add varation based on the eigenvector 
    ϵ = 1e-6   # koshiki_nanndakke()

    # get eigenvector for the stable manifold
    v_stb = get_eigenvector(monodromy, stable::True)

    state_f1 = x_t + ϵ * v_stb * norm(v_stb)
    state_f2 = x_t - ϵ * v_stb * norm(v_stb)

    # coodinate transformation

    return state_f
end
