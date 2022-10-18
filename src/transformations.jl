"""
Transformations
"""


"""
Transform state from Earth-Moon rotating frame to Sun-B1 rotating frame.
Careful with sign of ωs!! (should be negative)
"""
function transform_EMrot_to_SunB1(
    state,
    θs::Real,
    ωs::Real,
)
    ωm = -ωs
    θm = π - θs
    cos_θm = cos(θm)
    sin_θm = sin(θm)
    C = [
        cos_θm -sin_θm 0
        sin_θm  cos_θm 0
        0       0      1
    ]
    Cdot = [
        -ωm*sin_θm -ωm*cos_θm 0
         ωm*cos_θm -ωm*sin_θm 0
         0          0         1
    ]
    state_conv = vcat(
        C*state[1:3],
        C*state[4:6] + Cdot*state[1:3],
    )
    return state_conv
end



"""
Transform state from Earth-Moon rotating frame to Sun-B1 rotating frame.
Careful with sign of ωm!! (should be positive)
"""
function transform_SunB1_to_EMrot(state, )
end
