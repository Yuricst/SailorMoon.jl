"""
Delta-V transcriptions

INPUT 
    μS          Sun GM (normalized)
    as          Sun-B1 distance (normalized)
    θ           moon angle (in SB1 rotating frame)
    state0      state vector [x,y,z,vx,vy,vz] (do not include mass)
    p           control parameter [τ, γ, β]
"""

"""
dummy function for the no thrust mode
"""
function dv_no_thrust(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0::Vector{Float64}, p::Vector)
    return [0.0, 0.0, 0.0]
end

"""
    dv_sun_dir_angles(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)

Construct delta-V vector, directing towards B2, Sun-(E-M barycenter) barycenter
Frame: Sun-B1 rotating frame. origin = E-M barycenter (B1) 
"""
function dv_sun_dir_sb1frame(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0::Vector{Float64}, p::Vector)
    τ, γ, β = p[1], p[2], p[3]

    b2 = [-μS/(μS+1)*as, 0, 0]
    sc = state0[1:3]

    dir = (b2-sc) / norm(b2-sc)

    # add furhter rotation in gamma and beta
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)

    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]      

    return τ *  rot1 * rot2 * dir  

end


"""
    dv_tidal_dir_angles(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)

Construct delta-V vector, directing along with tidal force vector in the Sun-B1 frame
Frame: Sun-B1 rotating frame. origin = B1 (E-M barycenter)
"""
function dv_tidal_dir_sb1frame(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0::Vector{Float64}, p::Vector)

    τ, γ, β = p[1], p[2], p[3]
    # sun-B1 direction unit vector
    r = [-as, 0, 0]
    r = r / norm(r)

    # first, obtain the direction in sun-B1 rotating frame
    phi = μS / as^3 * (3 * r*transpose(r) - Matrix{Float64}(I, 3, 3)) * (state0[1:3] - [as, 0,0])
    dir = phi / norm(phi)
    # println("phi:", dir, "  pos:", state0[1:3] - [as, 0,0])

    # add furhter rotation in gamma and beta   
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]

    return τ * rot1 * rot2 * dir 
end


"""
    Earth-Moon rotating direction in SB1 rot frame
    Here, S-B1 frame refers to the reduced S-B1 frame where B2 = Sun 
"""
function dv_EMrotdir_sb1frame(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0::Vector{Float64}, p::Vector)
    τ, γ, β = p[1], p[2], p[3]

    # vector B1 -> SC 
    x = state0[1] - as
    y = state0[2]

    dir = [
        -y/sqrt(x^2+y^2) 
        x/sqrt(x^2+y^2)
        0
        ]
    
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    # rot about y-axis 
    rot1 = [
        cos_β  0 sin_β
        0      1 0
        -sin_β 0 cos_β
    ]

    # rot about z-axis 
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ  0
        0     0      1
    ]

    return τ * rot1 * rot2 * dir 

end



"""
    anti-Earth-Moon rotating direction in EM rot. frame
    According to Scheuerle et al. (2023 AAS), thrusting in the anti-Earth-Moon rotating direction 
    in CR3BP scenario (i.e., EMrot frame) is the optimal in order to enhance the Earth-Moon instantaneous 
    Jacobi constant. This is the implementation of his idea. 
    However,  we are propagating the dynamics in SB1 frame so we need to rotate the coordination at every timestep.
    Here, S-B1 frame refers to the reduced S-B1 frame where B2 = Sun. 
"""
function dv_maxJC_dir_sb1frame(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0, p::Vector)
    τ, γ, β = p[1], p[2], p[3]

    # change into EMrot fram 
    state_EMrot = transform_SunB1_to_EMrot(state0[1:6], θ, ωm, as)

    dir_EMrot = -state_EMrot[4:6]
    vec = vcat(dir_EMrot, [0.0,0.0,0.0])
    vec_sb1 = transform_EMrot_to_SunB1(vec, θ, ωm, as)
    dir = vec_sb1[1:3] - [as, 0, 0]

    # dir  = - state0[4:6] / norm(state0[4:6])  # opposite of the SC velocity
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    # rot about y-axis 
    rot1 = [
        cos_β  0 sin_β
        0      1 0
        -sin_β 0 cos_β
    ]

    # rot about z-axis 
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ  0
        0     0      1
    ]

    return τ * rot1 * rot2 * dir 
end


function dv_vel_dir_sb1frame(μS::Float64, as::Float64, θ::Float64, ωm::Float64, state0, p::Vector)
    if typeof(p) == Vector{Int64}
        p = float(p)
    end

    τ, γ, β = p[1], p[2], p[3]

    dir  = state0[4:6] / norm(state0[4:6])  # opposite of the SC velocity
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    # rot about y-axis 
    rot1 = [
        cos_β  0 sin_β
        0      1 0
        -sin_β 0 cos_β
    ]

    # rot about z-axis 
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ  0
        0     0      1
    ]

    return τ * rot1 * rot2 * dir 
end

function dv_max_drpdt_dir_sb1frame(μS::Float64, as::Float64, θ, ωm::Float64, state0, p::Vector)
    τ, γ, β = p[1], p[2], p[3]

    # sb1 -> Earth inertial 
    # state_Eine = transform_sb1_to_EearthIne(state0, θ, ωm, param3b.mu2, as)
    state_Eine = transform_sb1_to_EearthIne(state0, 0.0, ωm, param3b.mu2, as)


    koe = cart2kep(state0, param3b.mu1)
    sma, ecc, inc, OMEGA, omega, theta = koe
    r = norm(state0[1:3])
    p = sma*(1-ecc^2)
    h = norm(cross(state0[1:3], state0[4:6]))
    
    k1 = 2*sma^2*(1-ecc)*ecc *sin(theta) / h - sma*p/h * sin(theta)
    k2 = (1-ecc)*p - sma * ((p + r) *cos(theta) + r*ecc)

    gamma = atan(k1/k2)

    if k1*cos(gamma) + k2*sin(gamma) < 0
        gamma = gamma + pi
    end

    dir_RTN = [sin(gamma), cos(gamma), 0]

    # RTN -> xyz Earth-inertial
    dir_Eine = transform_RTN_to_EarthIne(dir_RTN, state_Eine)
    # dir_Eine = vcat(dir_Eine, [0.0,0.0,0.0])
    
    # Eine frame is parallel to SB1 rotating frame, so we do not need to change the direction from here 
    dir = dir_Eine    
    
    # xyz Earth-inertial -> rot coord & frame
    # dir_sb1 = transform_EearthIne_to_sb1(dir_Eine, θ, ωm, param3b.mu2, as) 
    # dir = dir_sb1[1:3]

    # println("dir_sb1:", norm(dir), " ", dir)

    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    # rot about y-axis 
    rot1 = [
        cos_β  0 sin_β
        0      1 0
        -sin_β 0 cos_β
    ]

    # rot about z-axis 
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ  0
        0     0      1
    ]

    return τ * rot1 * rot2 * dir
end


"""
    dv_inertial_angles(mu::Float64, state0::Vector{Float64}, vinf_params)

Construct delta-V vector based on angles w.r.t. inertial frame (Sun-B1 frame)
"""
function dv_inertial_angles(vinf_params)
    # unpack dv-parameters
    τ, γ, β = vinf_params[1], vinf_params[2], vinf_params[3]
    dv_vec = τ * [cos(γ) * cos(β), sin(γ) * cos(β), sin(β)]
    return dv_vec
end



### ====== NOT IN USE (06/09/2023) ==============================================================


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
    dv_tidal_dir_angles(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)

Construct delta-V vector, directing along with tidal force vector in the Sun-B1 frame
Frame: Sun-B1 rotating frame. origin = B2 
"""
function dv_tidal_dir_angles2(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)
    τ = control[1]
    # sun-B1 direction unit vector
    r = [-as, 0, 0]
    r = r / norm(r)

    # first, obtain the direction in sun-B1 rotating frame
    phi = μS / as^3 * (3 * r*transpose(r) - Matrix{Float64}(I, 3, 3)) * (state0[1:3] - [as,0,0])
    dir = phi / norm(phi)

    # add furhter rotation in gamma and beta
    γ = control[2]
    β = control[3]
    
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]

    return τ * rot1 * rot2 * dir 
end

"""
    dv_sun_dir_angles_emframe(μS::Float64, as::Float64, θ::Float64, state0::Vector{Float64}, p::Vector{Float64})

Construct delta-V vector, directing towards B2, Sun-(E-M barycenter) barycenter
"""
function dv_sun_dir_angles_emframe(μS::Float64, as::Float64, θ::Float64, state0, p::Vector{Float64})
    # sun-b1 distance vector
    τ = p[1]
    rs = [as * cos(θ), as * sin(θ), 0]
    sc = state0[1:3]
    dir = (rs - sc) / norm(rs - sc)

    # add furhter rotation in gamma and beta
    γ = p[2]
    β = p[3]
    
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]

    return τ * rot1 * rot2 * dir 
end

"""
    dv_tidal_dir_angles_emframe(μS::Float64, as::Float64, state0::Vector{Float64}, τ::Float64)

Construct delta-V vector, directing along with tidal force vector
"""

function dv_tidal_dir_angles_emframe(μS::Float64, as::Float64, θ::Float64, state0, p::Vector{Float64})
    τ = p[1]
    # sun-B1 direction unit vector
    r = [-as, 0, 0]
    r = r / norm(r)

    # first, obtain the direction in sun-B1 rotating frame
    phi = μS / as^3 * (3 * r*transpose(r) - Matrix{Float64}(I, 3, 3)) * state0[1:3]
    dir = phi / norm(phi)

    # convert it to the earth-moon rotating frame
    θ2 = pi - θ
    cos_θ2 = cos(-θ2)
    sin_θ2 = sin(-θ2)

    C = [
        cos_θ2 -sin_θ2 0
        sin_θ2  cos_θ2 0
        0       0      1
    ]

    dir = C * dir

    # add furhter rotation in gamma and beta
    γ = p[2]
    β = p[3]
    
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]

    return τ * rot1 * rot2 * dir 
end


"""
    dv_sun_dir_angles(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)

Construct delta-V vector, directing towards B2, Sun-(E-M barycenter) barycenter
Frame: Sun-B1 rotating frame. origin = B2 
"""
function dv_sun_dir2(μS::Real, as::Real, θ::Real, state0::Vector, control::Vector)
    b2 = [0.0, 0, 0]
    sc = state0[1:3]

    dir = (b2-sc) / norm(b2-sc)

    # add furhter rotation in gamma and beta
    γ = control[2]
    β = control[3]
    
    sin_β = sin(β)
    cos_β = cos(β)
    sin_γ = sin(γ)
    cos_γ = cos(γ)
    
    rot1 = [
        cos_β 0 sin_β
        0     1 0
        -sin_β 0 cos_β
    ]
    rot2 = [
        cos_γ -sin_γ 0 
        sin_γ cos_γ 0
        0 0 1
    ]

    return control[1] *  rot1 * rot2 * dir  

end

