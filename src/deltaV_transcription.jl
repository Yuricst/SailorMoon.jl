"""
Delta-V transcriptions
"""

"""
    dv_lvlh2inertial(mu::Float64, state0::Vector{Float64}, vinf_params)

Convert velocity vector in LVLH frame to inertial frame

Args:
    `mu::Float64`: gravitational parameter
    `state0::Vector{Float}`: state before delta v is applied, should have at least length 6 or 7 (rx,ry,rz,vx,vy,vz,m)
    `vinf_params::Vector{Float64}`: [v, α, δ]
        v : V-inf magnitude
        α : right ascension
        δ : declination

Returns:
   `::Vector{Float}`: delta-v vector direction scaled by tau, in heliocentric frame
"""

function dv_lvlh2inertial(mu::Float64, state0, vinf_params)
    # get delta-V in LVLH frame
    v, α, δ = vinf_params[1], vinf_params[2], vinf_params[3]
    v_lvlh = v * [cos(α) * cos(δ), sin(α) * cos(δ), sin(δ)]
    # conversion matrix following Markley and Crassidis 2014 pg.36
    r, v = state0[1:3], state0[4:6]
    o3I = -r / norm(r)
    o2I = -cross(r, v) / norm(cross(r, v))
    o1I = cross(o2I, o3I)
    A_IO = reshape(vcat(o1I, o2I, o3I), 3, 3)
    dv_inertial = A_IO * dv_lvlh
    return dv_inertial
end

"""
    dv_inertial_angles(mu::Float64, state0::Vector{Float64}, vinf_params)

Construct delta-V vector based on angles w.r.t. inertial frame
"""
function dv_inertial_angles(state0, vinf_params)
    # unpack dv-parameters
    τ, γ, β = vinf_params[1], vinf_params[2], vinf_params[3]
    dv_vec = τ * [cos(γ) * cos(β), sin(γ) * cos(β), sin(β)]
    return dv_vec
end


"""
    dv_sun_dir_angles(μS::Float64, as::Float64, state0::Vector{Float64}, τ::Float64)

Construct delta-V vector, directing towards B2, Sun-(E-M barycenter) barycenter
"""
function dv_sun_dir_angles(μS::Float64, as::Float64, state0::Vector{Float64}, τ::Float64)
    b2 = [-μS/(μS+1)*as,0,0]
    sc = state0

    dir = (b2-sc)/ norm(b2-sc)
    return dir * τ
end