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
    r_park::Real
    v_park::Real
end


"""
Get dynamics parameters for Earth-Moon-Sun system
"""
function dyanmics_parameters(use_sun::Bool=true)
    mu2    = 1.215058560962404E-2     # Moon
    mu1    = 1 - mu2                  # Earth
    gm_em  = 4.0350323550225981E+05   # GM of Earth-Moon system
    gm_sun = 1.3271244004193938E+11   # GM of Sun

    t_sidereal = 27.3217  *86400        # sec
    t_synodic  = 29.530588*86400      # sec

    tstar = t_sidereal / (2π)         # sec
    lstar = (tstar^2 * gm_em)^(1/3)   # km

    if use_sun == true
        mus   = gm_sun/gm_em
        as    = 1.000003 * 1.495978707e8 / lstar
    else
        mus = 0.0
        as  = 1.0
    end

    oms   = -2π/(t_synodic/tstar)     # rad/[canonical time]
    oml   =  2π/(t_synodic/tstar)     # rad/[canonical time]
    # oml  = abs(oms)/(1-abs(oms))      # this might be correct? (Boudad 2022 PhD thesis Eq. C.19)

    # parking radius
    r_park_km = 6378 + 200            # parking radius, km
    r_park = r_park_km/lstar          # parking radius, canonical
    v_park = sqrt(mu1/r_park)         # parking velocity, canonical
    return dynamics_params(
        mu1, mu2, mus, lstar, tstar, as, oms, oml, r_park, v_park
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
    - `t`: time
"""
function rhs_bcr4bp!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM, τ, γ, β = p

    θ = θ0 + ωM * t

    # compute distances
    r30 = sqrt((as - x)^2 + y^2 + z^2)
    r31 = sqrt((-μ2*cos(θ) - x)^2 + (-μ2*sin(θ) - y)^2 + z^2)
    r32 = sqrt(((1-μ2)*cos(θ) - x)^2 + ((1-μ2)*sin(θ) - y)^2 + z^2)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # forces applied (Sun, Earth, Moon, and thrust term)
    Fx = μS / r30^3 * (-as - x) + (1-μ2) / r31^3 * (-μ2*cos(θ) - x)    + μ2 / r32^3 * ((1-μ2)*cos(θ) - x)
    Fy = μS / r30^3 * (-y)      + (1-μ2) / r31^3 * (-μ2*sin(θ) - y)    + μ2 / r32^3 * ((1-μ2)*sin(θ) - y)
    Fz = μS / r30^3 * (-z)      + (1-μ2) / r31^3 * (-z)                + μ2 / r32^3 * (-z)

    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =             Fz
end


"""

Right-hand side expression for state-vector in BCR4BP

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p = [μ2, μS, as, θ0, ωM, τ, γ, β, mdot, tmax]
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
        tmax : max thrust
    - `t`: time
"""
function rhs_bcr4bp_with_mass!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]
    τ, γ, β, mdot, tmax = p[6], p[7], p[8], p[9], p[10]

    θ = θ0 + ωM * t  # moon angle

    # create Thrust term
    T = dv_inertial_angles([τ,γ,β])
    Tx, Ty, Tz = T * tmax / u[7]  # mdot

    # compute distances
    r30 = sqrt((as + x)^2 + y^2 + z^2)
    r31 = sqrt((-μ2*cos(θ) - x)^2 + (-μ2*sin(θ) - y)^2 + z^2)
    r32 = sqrt(((1-μ2)*cos(θ) - x)^2 + ((1-μ2)*sin(θ) - y)^2 + z^2)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # forces applied (Sun, Earth, Moon, and thrust term)
    Fx = μS / r30^3 * (-as-x) + μS/as^3 * as  + (1-μ2) / r31^3 * (-μ2*cos(θ) - x)    + μ2 / r32^3 * ((1-μ2)*cos(θ) - x) + Tx
    Fy = μS / r30^3 * (-y)      + (1-μ2) / r31^3 * (-μ2*sin(θ) - y)    + μ2 / r32^3 * ((1-μ2)*sin(θ) - y) + Ty
    Fz = μS / r30^3 * (-z)      + (1-μ2) / r31^3 * (-z)                + μ2 / r32^3 * (-z)                + Tz

    # println(μS / r30^3)
    println(μS / r30^3 * (-as-x)+ μS/as^3 * as , " ", (1-μ2) / r31^3 * (-μ2*cos(θ) - x), " ",  μ2 / r32^3 * ((1-μ2)*cos(θ) - x))

    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =           + Fz

    du[7] = -mdot * τ
end

"""

Right-hand side expression for state-vector in BCR4BP

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p = [μ2, μS, as, θ0, ωM, τ, γ, β, mdot, tmax]
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
        tmax : max thrust
    - `t`: time
"""
function rhs_bcr4bp_sb1frame!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]
    τ, γ, β, mdot, tmax = p[6], p[7], p[8], p[9], p[10]

    θ = θ0 + ωM * t  # moon angle

    # create Thrust term
    T = dv_inertial_angles([τ,γ,β])
    Tx, Ty, Tz = T * tmax / u[7]  # mdot


    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # earth and moon location
    xe = μS/(μS+1) - μ2/as *cos(θ)
    ye = μS/(μS+1) - μ2/as *sin(θ)
    ze = 0

    xm = μS/(μS+1) + (1-μ2)/as *cos(θ)
    ym = μS/(μS+1) + (1-μ2)/as *sin(θ)
    zm = 0

    # compute distances
    r30 = sqrt((-1/(μS+1)- x)^2 + y^2 + z^2)
    r31 = sqrt((xe - x)^2 + (ye - y)^2 + (ze-z)^2)
    r32 = sqrt((xm - x)^2 + (ym - y)^2 + (zm-z)^2)
    println(r30, " ", r31, " ", r32)


    Fx = -(μS/(μS+1))*(x-1/(μS+1))/r30^3 - (1-μ2)/(μS+1)*(x-xe)/r31^3 - μ2/(μS+1)*(x-xm)/r32^3 + Tx
    Fy = -(μS/(μS+1))*(y-1/(μS+1))/r30^3 - (1-μ2)/(μS+1)*(y-ye)/r31^3 - μ2/(μS+1)*(y-ym)/r32^3 + Ty
    Fz = -(μS/(μS+1))*(z-1/(μS+1))/r30^3 - (1-μ2)/(μS+1)*(z-ze)/r31^3 - μ2/(μS+1)*(z-zm)/r32^3 + Tz
    # println(-(μS/(μS+1))*(x-1/(μS+1))/r30^3 , " ", - (1-μ2)/(μS+1)*(x-xe)/r31^3, " ",  - μ2/(μS+1)*(x-xm)/r32^3)


    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =           + Fz

    du[7] = -mdot * τ
end


"""

Right-hand side expression for state-vector in BCR4BP, thrusting predefined direction

# Arguments
    - `du`: cache array of duative of state-vector, mutated
    - `u`: state-vector
    - `p`: parameters, where p = [μ2, μS, as, θ0, ωM, τ, mdot, tmax]
        m* : mE + mL (6.0455 * 10^22)
        μ2 : mL / m* (0.012150585609624)
        μS : mS / m* (0.00000303951)
        as : (sun-B1 distance) / l* (388.709677419)
            l* : Earth-moon distance (384,400 km)
        ωM : Earth-Moon line's angular velocity around E-M barycenter
        τ  : thrust magnitude (0~1)
        mdot : mass-flow rate
        tmax : max thrust
    - `t`: time
"""
function rhs_bcr4bp_predifined_dir!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]
    τ, mdot, tmax = p[6], p[7], p[8]

    θ = θ0 + ωM * t

    # create Thrust term
    T = dv_sun_dir_angles(μ, as, [x,y,z], τ)  # dimensionless vector
    Tx, Ty, Tz = T * tmax / u[7]  # make into [m/s^2]

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

    du[4] =  2*vy + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =           + Fz

    du[7] = -mdot * τ
end



"""
    rhs_bcr4bp_thrust!(du,u,p,t)

BCR4BP equation of motion, centered at Earth-Moon barycenter.

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector, including mass
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s,
                             p[6] = τ, p[7] = θ, p[8] = β, p[9] = mdot, p[10] = tmax
    - `t`: time
"""
function rhs_bcr4bp_emframe_thrust!(du,u,p,t)
    # unpack arguments
    μ, μ_3, t0, a_s, ω_s, τ, γ, β, mdot, tmax, dv_fun = p
    # decompose state
    x, y, z, vx, vy, vz, mass = u

    # calculate radii
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 )
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 )

    # sun position
    θ = ω_s*t + t0
    xs = a_s*cos(θ)
    ys = a_s*sin(θ)
    zs = 0.0
    r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

    # compute delta-V vector
    dir_v = dv_fun(μ_3, a_s, θ, u[1:6], [τ, γ, β])
    ts = tmax*dir_v

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-μ)/r1^3)*(μ+x)  + (μ/r2^3)*(1-μ-x)  + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs ) + ts[1]/mass
    du[5] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys ) + ts[2]/mass
    du[6] =           - ((1-μ)/r1^3)*z      - (μ/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs ) + ts[3]/mass

    # mass derivative
    du[7] = -mdot*τ
end
