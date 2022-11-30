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