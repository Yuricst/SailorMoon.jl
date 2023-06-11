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
    omb::Real
    r_park::Real
    v_park::Real
end


"""
Get dynamics parameters for Earth-Moon-Sun system
"""
function dynamics_parameters(use_sun::Bool=true)
    mu2    = 1.215058560962404E-2     # Moon
    mu1    = 1 - mu2                  # Earth
    gm_em  = 4.0350323550225981E+05   # GM of Earth-Moon system
    gm_sun = 1.3271244004193938E+11   # GM of Sun

    t_sidereal = 27.3217  *86400      # sec
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

    oms   = -2π/(t_synodic/tstar)     # rad/[canonical time]  # Sun-B1 rotating velocity in EM rot frame
    oml   =  2π/(t_synodic/tstar)     # rad/[canonical time]  # E-M rotating velocity in S-B1 rot frame
    omb   = 2*π * (t_synodic-t_sidereal) / (t_sidereal*t_synodic) * tstar   # Sun-B1 rotating velocity in S-B1 rot frame

    # parking radius
    r_park_km = 6378 + 200            # parking radius, km
    r_park = r_park_km/lstar          # parking radius, canonical
    v_park = sqrt(mu1/r_park)         # parking velocity, canonical
    return dynamics_params(
        mu1, mu2, mus, lstar, tstar, as, oms, oml, omb, r_park, v_park
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

Right-hand side expression for state-vector in BCR4BP at Earth-Moon barycenter

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

Right-hand side expression for state-vector in BCR4BP in Sun-B1 frame
(origin = B2 = Sun-B1 barycenter)

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
            (= 2π/(t_synodic/t_star))
        ωb : rotating velocity of S-B1 frame
            (= 2π*(t_synodic - t_sidereal)/(t_sidereal*t_synodic)*t_star)
        τ  : thrust magnitude (0~1)
        γ  : thrust angle 1
        β  : thrust angle 2
        mdot : mass-flow rate
        tmax : max thrust
    - `t`: time

*If you normalize l* = E-M distance, m* = mE + mM, and t* = t_sidereal/2pi,
then ωb = 0.074800, ωM = 0.9252. These values should work in this project.
"""
function rhs_bcr4bp_sb1frame!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM, ωb = p[1], p[2], p[3], p[4], p[5], p[6]
    τ, γ, β, mdot, tmax = p[7], p[8], p[9], p[10], p[11]

    θ = θ0 + ωM * t  # moon angle

    # create Thrust term
    #  T = dv_inertial_angles([τ,γ,β])
    T = [0.0, 0.0, 0.0]
    Tx, Ty, Tz = T * tmax / u[7]  # mdot


    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # earth and moon location
    xe = - μ2 *cos(θ) + as * μS/(μS+1)
    ye = - μ2 *sin(θ)
    ze = 0

    xm = (1-μ2) *cos(θ) + as * μS/(μS+1)
    ym = (1-μ2) *sin(θ)
    zm = 0
    # println(xe, " ", ye, " ", ze)

    # compute distances
    r30 = sqrt((-as/(μS+1) - x)^2 + y^2 + z^2)      # SC-Sun
    r31 = sqrt((xe - x)^2 + (ye - y)^2 + (ze-z)^2)  # SC-earth
    r32 = sqrt((xm - x)^2 + (ym - y)^2 + (zm-z)^2)  # SC-Moon
    # println(r30, " ", r31, " ", r32)

    Fx = -(μS)*(x+as/(μS+1))/r30^3 - (1-μ2)*(x-xe)/r31^3 - μ2*(x-xm)/r32^3 + Tx
    Fy = -(μS)* y /r30^3           - (1-μ2)*(y-ye)/r31^3 - μ2*(y-ym)/r32^3 + Ty
    Fz = -(μS)* z /r30^3           - (1-μ2)*(z-ze)/r31^3 - μ2*(z-zm)/r32^3 + Tz
    # println(-μS*(x-1/(μS+1))/r30^3 , " ", - (1-μ2)*(x-xe)/r31^3, " ",  - μ2*(x-xm)/r32^3)

    du[4] =  2*ωb*vy + ωb^2*x + Fx
    du[5] = -2*ωb*vx + ωb^2*y + Fy
    du[6] =                   + Fz

    du[7] = -mdot * τ
end


"""

Right-hand side expression for state-vector in "reduced" BCR4BP in Sun-B1 frame
(origin = Sun). This EoM is exactly the same to the rhs_bcr4bp_emframe_thrust! after the adequate
coorinate transformation. Note that E-M rotating frame's BCR4BP and the original S-B1 rotating frame's BCR4BP
are different dynamics, strictly speaking.

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
            (= 2π/(t_synodic/t_star))
        ωb : rotating velocity of S-B1 frame
            (= 2π*(t_synodic - t_sidereal)/(t_sidereal*t_synodic)*t_star)
        τ  : thrust magnitude (0~1)
        γ  : thrust angle 1
        β  : thrust angle 2
        mdot : mass-flow rate
        tmax : max thrust
    - `t`: time

*If you normalize l* = E-M distance, m* = mE + mM, and t* = t_sidereal/2pi,
then ωb = 0.074800, ωM = 0.9252. These values should work in this project.
"""
function rhs_bcr4bp_sb1frame2!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM, ωb = p[1], p[2], p[3], p[4], p[5], p[6]
    τ, mdot, tmax = p[7], p[8], p[9], p[10], p[11]

    θ = θ0 + ωM * t  # moon angle

    # create Thrust term
#     T = dv_inertial_angles([τ,γ,β])
    T = [0.0, 0.0, 0.0]
    Tx, Ty, Tz = T * tmax / u[7]  # mdot


    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # earth and moon location
    xe = - μ2 *cos(θ) + as
    ye = - μ2 *sin(θ)
    ze = 0

    xm = (1-μ2) *cos(θ) + as
    ym = (1-μ2) *sin(θ)
    zm = 0
#     println(xe, " ", ye, " ", ze)

    # compute distances
    r30 = sqrt(x^2 + y^2 + z^2)      # SC-Sun
    r31 = sqrt((xe - x)^2 + (ye - y)^2 + (ze-z)^2)  # SC-earth
    r32 = sqrt((xm - x)^2 + (ym - y)^2 + (zm-z)^2)  # SC-Moon
#     println(r30, " ", r31, " ", r32)

    Fx = -(μS)*(x)/r30^3 - (1-μ2)*(x-xe)/r31^3 - μ2*(x-xm)/r32^3 + Tx
    Fy = -(μS)* y /r30^3           - (1-μ2)*(y-ye)/r31^3 - μ2*(y-ym)/r32^3 + Ty
    Fz = -(μS)* z /r30^3           - (1-μ2)*(z-ze)/r31^3 - μ2*(z-zm)/r32^3 + Tz
#     println(-μS*(x-1/(μS+1))/r30^3 , " ", - (1-μ2)*(x-xe)/r31^3, " ",  - μ2*(x-xm)/r32^3)

    du[4] =  2*ωb*vy + ωb^2*x + Fx
    du[5] = -2*ωb*vx + ωb^2*y + Fy
    du[6] =                   + Fz

    du[7] = -mdot * τ
end

"""

Right-hand side expression for state-vector in "reduced" BCR4BP in Sun-B1 frame
(origin = Sun). This EoM is exactly the same to the rhs_bcr4bp_emframe_thrust! after the adequate
coorinate transformation. Note that E-M rotating frame's BCR4BP and the original S-B1 rotating frame's BCR4BP
are different dynamics, strictly speaking.

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
            (= 2π/(t_synodic/t_star))
        ωb : rotating velocity of S-B1 frame
            (= 2π*(t_synodic - t_sidereal)/(t_sidereal*t_synodic)*t_star)
        τ  : thrust magnitude (0~1)
        γ  : thrust angle 1
        β  : thrust angle 2
        mdot : mass-flow rate
        tmax : max thrust
    - `t`: time

*If you normalize l* = E-M distance, m* = mE + mM, and t* = t_sidereal/2pi,
then ωb = 0.074800, ωM = 0.9252. These values should work in this project.
"""
function rhs_bcr4bp_sb1frame2_thrust!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM, ωb = p[1], p[2], p[3], p[4], p[5], p[6]
    τ, γ, β, mdot, tmax = p[7], p[8], p[9], p[10], p[11]
    dv_fun = p[12]

    θ = θ0 + ωM * t  # moon angle

    # create Thrust term
    dir_v = dv_fun(μS, as, θ, ωM, u[1:6], [τ, γ, β]) 
    Tx, Ty, Tz = dir_v * tmax / u[7]  # mdot
    # println(T)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # earth and moon location
    xe = - μ2 *cos(θ) + as
    ye = - μ2 *sin(θ)
    ze = 0

    xm = (1-μ2) *cos(θ) + as
    ym = (1-μ2) *sin(θ)
    zm = 0
#     println(xe, " ", ye, " ", ze)

    # compute distances
    r30 = sqrt(x^2 + y^2 + z^2)      # SC-Sun
    r31 = sqrt((xe - x)^2 + (ye - y)^2 + (ze-z)^2)  # SC-earth
    r32 = sqrt((xm - x)^2 + (ym - y)^2 + (zm-z)^2)  # SC-Moon
#     println(r30, " ", r31, " ", r32)

    Fx = -(μS)*(x)/r30^3 - (1-μ2)*(x-xe)/r31^3 - μ2*(x-xm)/r32^3 + Tx
    Fy = -(μS)* y /r30^3           - (1-μ2)*(y-ye)/r31^3 - μ2*(y-ym)/r32^3 + Ty
    Fz = -(μS)* z /r30^3           - (1-μ2)*(z-ze)/r31^3 - μ2*(z-zm)/r32^3 + Tz
#     println(-μS*(x-1/(μS+1))/r30^3 , " ", - (1-μ2)*(x-xe)/r31^3, " ",  - μ2*(x-xm)/r32^3)

    du[4] =  2*ωb*vy + ωb^2*x + Fx
    du[5] = -2*ωb*vx + ωb^2*y + Fy
    du[6] =                   + Fz

    du[7] = -mdot * τ
end


"""
    adding the coasting term... 
"""
function rhs_bcr4bp_sb1frame2_thrust_bal!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ2, μS, as, θ0, ωM, ωb = p[1], p[2], p[3], p[4], p[5], p[6]
    τ, γ, β, mdot, tmax = p[7], p[8], p[9], p[10], p[11]
    dv_fun = p[12]
    tstar = p[13]

    θ = θ0 + ωM * t  # moon angle

    # 10 days 
    t_ballistic = 10 * 24*60*60 / tstar
    # create Thrust term
    if abs(t) < t_ballistic
        dir_v = [0.0, 0.0, 0.0]
        mdot = 0.0

        Tx = 0.0
        Ty = 0.0
        Tz = 0.0
    else 
        dir_v = dv_fun(μS, as, θ, ωM, u[1:6], [τ, γ, β]) 
        Tx, Ty, Tz = dir_v * tmax / u[7]  # mdot
    end

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]

    # earth and moon location
    xe = - μ2 *cos(θ) + as
    ye = - μ2 *sin(θ)
    ze = 0

    xm = (1-μ2) *cos(θ) + as
    ym = (1-μ2) *sin(θ)
    zm = 0
#     println(xe, " ", ye, " ", ze)

    # compute distances
    r30 = sqrt(x^2 + y^2 + z^2)      # SC-Sun
    r31 = sqrt((xe - x)^2 + (ye - y)^2 + (ze-z)^2)  # SC-earth
    r32 = sqrt((xm - x)^2 + (ym - y)^2 + (zm-z)^2)  # SC-Moon
#     println(r30, " ", r31, " ", r32)

    Fx = -(μS)*(x)/r30^3 - (1-μ2)*(x-xe)/r31^3 - μ2*(x-xm)/r32^3 + Tx
    Fy = -(μS)* y /r30^3           - (1-μ2)*(y-ye)/r31^3 - μ2*(y-ym)/r32^3 + Ty
    Fz = -(μS)* z /r30^3           - (1-μ2)*(z-ze)/r31^3 - μ2*(z-zm)/r32^3 + Tz
#     println(-μS*(x-1/(μS+1))/r30^3 , " ", - (1-μ2)*(x-xe)/r31^3, " ",  - μ2*(x-xm)/r32^3)

    du[4] =  2*ωb*vy + ωb^2*x + Fx
    du[5] = -2*ωb*vx + ωb^2*y + Fy
    du[6] =                   + Fz

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
    dir_v = dv_fun(μS, as, θ, ωM, u[1:6], [τ, γ, β]) 
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


