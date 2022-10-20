"""
Dynamics functions
"""

abstract type AbstractParameterType end

Base.@kwdef struct dynamics_params <: AbstractParameterType
    mu1::Real
    mu2::Real
    mus::Real
    lstar::Real
    tstar::Real
    as::Real
    oms::Real
    oml::Real
end

"""
Get dynamics parameters for Earth-Moon-Sun system
"""
function dyanmics_parameters()
    mu2    = 1.215058560962404E-2     # Moon
    mu1    = 1 - mu2                  # Earth
    gm_em  = 4.0350323550225981E+05   # GM of Earth-Moon system
    gm_sun = 1.3271244004193938E+11   # GM of Sun

    t_sidereal = 27.3217  *86400        # sec
    t_synodic  = 29.530588*86400      # sec

    tstar = t_sidereal / (2π)         # sec
    lstar = (tstar^2 * gm_em)^(1/3)   # km

    mus   = gm_sun/gm_em
    as    = 1.000003 * 1.495978707e8 / lstar

    oms   = -2π/(t_synodic/tstar)     # rad/[canonical time]
    oml   = 2π/(t_sidereal/tstar)     # rad/[canonical time]
    return dynamics_params(
        mu1, mu2, mus, lstar, tstar, as, oms, oml
    )
end



"""
    lagrange_points(μ::Float64)

Function computes lagrange points from CR3BP paraneter μ.
Returns a 5 x 6 matrix, each row corresopnding to state of lagrange point.
"""
function lagrange_points(μ::Float64)
    # place-holder
    l = 1 - μ;
    # collinear points
    fl1(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)+μ-l)*x^2 + (μ^2*l^2+2*(l^2+μ^2))*x + (μ^3-l^3);
    xl1 = find_zero(fl1, (0,2), Bisection());

    fl2(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)-μ-l)*x^2 + (μ^2*l^2+2*(l^2-μ^2))*x - (μ^3+l^3);
    xl2 = find_zero(fl2, (0,2), Bisection());

    fl3(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)+μ+l)*x^2 + (μ^2*l^2+2*(μ^2-l^2))*x + (μ^3+l^3);
    xl3 = find_zero(fl3, (-2,0), Bisection());

    LP = zeros(5, 6);
    LP[1,1] = xl1;
    LP[2,1] = xl2;
    LP[3,1] = xl3;
    # equilateral points
    LP[4,1] = cos(pi/3)-μ;
    LP[4,2] = sin(pi/3);
    LP[5,1] = cos(pi/3)-μ;
    LP[5,2] = -sin(pi/3);
    return LP
end


"""
    rhs_cr3bp!(du,u,p,t)

Right-hand side expression for state-vector in CR3BP

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_cr3bp!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 + z^2 );
    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # derivatives of velocities
    du[4] =  2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[5] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    du[6] =            -((1-p[1])/r1^3)*z - (p[1]/r2^3)*z;
end



"""

Right-hand side expression for state-vector in BCR4BP

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p = [μ2, μS, as, θ0, ωM, τ, γ, β]
        m* : mE + mL (6.0455 * 10^22)
        μ2 : mL / m* (0.012150585609624)
        μS : mS / m* (0.00000303951)
        as : (sun-B1 distance) / l* (388.709677419)
            l* : Earth-moon distance (384,400 km)
        ωM : Earth-Moon line's angular velocity around E-M barycenter
        τ  : thrust magnitude (0~1)
        γ  : thrust angle 1
        β  : thrust angle 2
        mdot : mass-flow rate
    - `t`: time
"""
function rhs_bcr4bp_with_mass!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]
    τ, γ, β, mdot = p[6], p[7], p[8], p[9]

    θ = θ0 + ωM * t

    # create Thrust term
    # To Do: we want to change this to the other cooridnate frame (directing towards Sun's tidal force etc.)
    #        in the future
    T = dv_inertial_angles(μ, [x,y,z], [τ,γ,β])
    Tx, Ty, Tz = T[1], T[2], T[3]

    # compute distances
    r30 = sqrt((as - x)^2 + y^2 + z^2)
    r31 = sqrt((-μ2*cos(θ) - x)^2 + (-μ2*sin(θ) - y)^2 + z^2)
    r32 = sqrt(((1-μ2)*cos(θ) - x)^2 + ((1-μ2)*sin(θ) - y)^2 + z^2)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # forces applied (Sun, Earth, Moon, and thrust term)
    Fx = μS / r30^3 * (-as - x) + (1-μ2) / r31^3 * (-μ2*cos(θ) - x)    + μ2 / r32^3 * ((1-μ2)*cos(θ) - x) + Tx
    Fy = μS / r30^3 * (-y)      + (1-μ2) / r31^3 * (-μ2*sin(θ) - y)    + μ2 / r32^3 * ((1-μ2)*sin(θ) - y) + Ty
    Fz = μS / r30^3 * (-z)      + (1-μ2) / r31^3 * (-z)                + μ2 / r32^3 * (-z)                + Tz

    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =             Fz

    du[7] = -mdot * τ
end

"""

Right-hand side expression for state-vector in BCR4BP, thrusting predefined direction

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p = [μ2, μS, as, θ0, ωM, τ, γ, β]
        m* : mE + mL (6.0455 * 10^22)
        μ2 : mL / m* (0.012150585609624)
        μS : mS / m* (0.00000303951)
        as : (sun-B1 distance) / l* (388.709677419)
            l* : Earth-moon distance (384,400 km)
        ωM : Earth-Moon line's angular velocity around E-M barycenter
        τ  : thrust magnitude (0~1)
        mdot : mass-flow rate
    - `t`: time
"""
function rhs_bcr4bp_predifined_dir!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]
    τ, mdot = p[6], p[7]

    θ = θ0 + ωM * t

    # create Thrust term
    # To Do: we want to change this to the other cooridnate frame (directing towards Sun's tidal force etc.)
    #        in the future
    T = dv_sun_dir_angles(μ, as, [x,y,z], τ)
    Tx, Ty, Tz = T[1], T[2], T[3]

    # compute distances
    r30 = sqrt((as - x)^2 + y^2 + z^2)
    r31 = sqrt((-μ2*cos(θ) - x)^2 + (-μ2*sin(θ) - y)^2 + z^2)
    r32 = sqrt(((1-μ2)*cos(θ) - x)^2 + ((1-μ2)*sin(θ) - y)^2 + z^2)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # forces applied (Sun, Earth, Moon, and thrust term)
    Fx = μS / r30^3 * (-as - x) + (1-μ2) / r31^3 * (-μ2*cos(θ) - x)    + μ2 / r32^3 * ((1-μ2)*cos(θ) - x) + Tx
    Fy = μS / r30^3 * (-y)      + (1-μ2) / r31^3 * (-μ2*sin(θ) - y)    + μ2 / r32^3 * ((1-μ2)*sin(θ) - y) + Ty
    Fz = μS / r30^3 * (-z)      + (1-μ2) / r31^3 * (-z)                + μ2 / r32^3 * (-z)                + Tz

    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =             Fz

    du[7] = -mdot * τ
end
