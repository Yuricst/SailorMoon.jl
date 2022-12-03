"""
Transformations
"""


"""
    transform_EMrot_to_SunB1(state::Vector, θs::Real, ωs::Real)

Transform state from Earth-Moon rotating frame (origin: E-M barycenter B1) 
to Sun-B1 rotating frame (origin: S-B1 barycenter B2).
Careful with sign of ωs!! (should be negative)
"""
function transform_EMrot_to_SunB1(state::Vector, θs::Real, ωs::Real, as::Real)
    ωm = -ωs
    θm = π - θs
    cos_θm = cos(θm)
    sin_θm = sin(θm)
    C = [
        cos_θm      -sin_θm   0 0       0      0
        sin_θm      cos_θm    0 0       0      0
        0           0         1 0       0      0
        -ωm*sin_θm -ωm*cos_θm 0 cos_θm -sin_θm 0 
         ωm*cos_θm -ωm*sin_θm 0 sin_θm  cos_θm 0
         0          0         0 0       0      1
    ]
    state_conv = C * state
    state_conv = state_conv + [as, 0,0,0,0,0]
    
    return state_conv
end


"""
Transform state from Earth-Moon rotating frame to Sun-B1 rotating frame.
Careful with sign of ωm!! (should be positive)
"""
function transform_SunB1_to_EMrot(state, θm::Real, ωm::Real)
end


"""
Transform state from Earth inertial frame to Earth-Moon rotating frame.
Careful with sign of ωm!! (should be positive)
    θ : angle of E-M line from earth inertial frame's x-axis (CW: -, CCW: +)
"""
function transform_earthIne_to_EMrot(state::Vector, θ::Real, ωm::Real, μ2::Real)
    # move the orign from earth to B1
    state = state - [μ2*cos(θ), μ2*sin(θ), 0, 0, 0, 0]
    
    # construct transformation matrix
    cos_θ = cos(θ)
    sin_θ = sin(θ)
    rotmat = inv(
        [
            cos_θ -sin_θ 0.0 0.0 0.0 0.0
            sin_θ cos_θ 0.0 0.0 0.0 0.0
            0.0 0.0 1.0 0.0 0.0 0.0
            -ωm*sin_θ -ωm*cos_θ 0.0 cos_θ -sin_θ 0.0
            ωm*cos_θ -ωm*sin_θ 0.0 sin_θ cos_θ 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ],
    )
    # transform
    state_r = rotmat * state
    return state_r
end


function transform_sb1_to_EearthIne(state_sb1::Vector, θm::Real, ωm::Real, μ2::Real, as::Real)
    cos_θm = cos(θm)
    sin_θm = sin(θm)
    
    # step 1: SunB1 -> EM rot
    state_ = state_sb1 + [-as, 0, 0, 0, 0, 0]

    rot_mat1 = inv([
        cos_θm      -sin_θm   0 0       0      0
        sin_θm      cos_θm    0 0       0      0
        0           0         1 0       0      0
        -ωm*sin_θm -ωm*cos_θm 0 cos_θm -sin_θm 0 
        ωm*cos_θm -ωm*sin_θm 0 sin_θm  cos_θm 0
        0          0         0 0       0      1
    ])
    state_emrot = rot_mat1 * state_

    # step 2: EM rot to Earth Inertial

    rot_mat2 = [
        cos_θm -sin_θm 0.0 0.0 0.0 0.0
        sin_θm cos_θm 0.0 0.0 0.0 0.0
        0.0 0.0 1.0 0.0 0.0 0.0
        -ωm*sin_θm -ωm*cos_θm 0.0 cos_θm -sin_θm 0.0
        ωm*cos_θm -ωm*sin_θm 0.0 sin_θm cos_θm 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
    ]
    state_ = rot_mat2 * state_emrot
    
    state_EarthIne = state_ + [μ2*cos(θm), μ2*sin(θm), 0, 0, 0, 0]

    
    return state_EarthIne
end


function transform_EearthIne_to_sb1(state_earthIne::Vector, θm::Real, ωm::Real, μ2::Real, as::Real)
    # EarthIne -> EMrot
    # move the orign from earth to B1
    state = state_earthIne - [μ2*cos(θm), μ2*sin(θm), 0, 0, 0, 0]
        
    # construct transformation matrix
    cos_θm = cos(θm)
    sin_θm = sin(θm)
    rotmat = inv(
        [
            cos_θm -sin_θm 0.0 0.0 0.0 0.0
            sin_θm cos_θm 0.0 0.0 0.0 0.0
            0.0 0.0 1.0 0.0 0.0 0.0
            -ωm*sin_θm -ωm*cos_θm 0.0 cos_θm -sin_θm 0.0
            ωm*cos_θm -ωm*sin_θm 0.0 sin_θm cos_θm 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ],
    )
    # transform
    state_em = rotmat * state

    # EMrot -> SB1
    C = [
        cos_θm      -sin_θm   0 0       0      0
        sin_θm      cos_θm    0 0       0      0
        0           0         1 0       0      0
        -ωm*sin_θm -ωm*cos_θm 0 cos_θm -sin_θm 0 
         ωm*cos_θm -ωm*sin_θm 0 sin_θm  cos_θm 0
         0          0         0 0       0      1
    ]
    state_ = C * state_em
    state_sb1 = state_ + [as, 0,0,0,0,0]
    
    return state_sb1
end


"""
convert the initial parameters (rp, ra, alpha) into the initial state in the SB1 frame 
"""
function paramIni_to_sb1(rp::Real, α::Real, ra::Real, θm::Real, ωm ::Real, μ2::Real, as::Real)
    # build an initial state in Earth inertial frame
    sma = (rp + ra) / 2
    v = sqrt((1-μ2) * (2/rp - 1/sma))
    
    state_i_EIne = [
        rp*cos(α), 
        rp*sin(α), 
        0,
        -v*sin(α),
        v*cos(α),
        0
    ] 
    println("state_EIne: ", state_i_EIne)
    
    # EarthIne -> Sunb1
    state_i_sb1 = SailorMoon.transform_EearthIne_to_sb1(state_i_EIne, θm, ωm, μ2, as)
    
    return state_i_sb1
end



"""
transformation from Cartesian frame to cylindrical frame, changing only positions but not velocity 
"""
function cart2cylind_only_pos(state::Vector)
    x = state[1]
    y = state[2]
    xdot = state[4]
    ydot = state[5]
    r = sqrt(x^2 + y^2)
    θ = atan(y,x)

    # might want to use it later...?
    # rdot = (x*xdot + y*ydot) / sqrt(x^2 + y^2)
    # θdot = (x*ydot - xdot*y) / (x^2+y^2)

    return vcat([r, θ], state[3:6])
end


"""
transformation from Cartesian frame to cylindrical frame, changing only positions but not velocity 
"""
function cylind2cart_only_pos(state::Vector)
    r = state[1]
    θ = state[2]
    rdot = state[4]
    θdot = state[5]

    x = r * cos(θ)
    y = r * sin(θ)

    # might want to use it later...?
    # xdot = rdot * cos(θ) - r * θdot * sin(θ)
    # ydot = rdot * sin(θ) + r * θdot * cos(θ)

    return vcat([x, y], state[3:6])
end
