"""
    matchpoint_residual(c_mp::Vector, SFset::SimsFlanaganSettings, sv_fwd::Vector, sv_bck::Vector)

Compute residual at match-point
"""
function matchpoint_residual(
    sv_fwd::Vector,
    sv_bck::Vector,
)
    return sv_bck[:] - sv_fwd[:]
end
