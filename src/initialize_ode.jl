"""
Initial & Terminal condition
"""

include("../../julia-R3BP/R3BP/src/librationpointorbit/lpo_family.jl")

"""

Providing the initial state of the SC

# Arguments
    - `v∞1` : v-infinity at the initial state (w.r.t. the Earth),
    - `α1`  : initial right ascension,
    - `δ1`  : initial declination,
    - `θ1`  : initial E-M line's angle w.r.t. Sun-B1 line
"""

function set_initial_state(v∞1, α1, δ1, θ1)
    
    return state_0 
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
    # what is moon and Earth ID? 
    construct_arrival_condition(ACDict="lpo2", arrival_ID::Int, center_body_ID::Int)

    return state_f
end
