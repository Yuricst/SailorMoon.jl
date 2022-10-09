"""
Dynamics functions
"""


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
    - `p`: parameters, where p = [μ1, μ2, as, ωM, τ, β, γ]
        m* : mE + mL (6.0455 * 10^22)
        μ1 : mL / m* (0.012150585609624)
        μ2 : mS / m* (0.00000303951)
        as : (sun-B1 distance) / l* (388.709677419)
            l* : Earth-moon distance (384,400 km)
        ωM : Earth-Moon line's angular velocity around E-M barycenter
        τ  : thrust magnitude (0~1)
        β  : thrust angle 1
        γ  : thrust angle 2
    - `t`: time
"""
function rhs_bcr4bp!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    μ1, μ2, as, θ0, ωM = p[1], p[2], p[3], p[4], p[5]

    θ = θ0 + ωM * t

    # compute distances
    r30 = sqrt((as - x)^2 + y^2 + z^2)
    r31 = sqrt((-μ1*cos(θ) - x)^2 + (-μ1*sin(θ) - y)^2 + z^2)
    r32 = sqrt(((1-μ1)*cos(θ) - x)^2 + ((1-μ1)*sin(θ) - y)^2 + z^2)

    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    
    # forces applied (Sun, Earth, Moon, and thrust term)
    Fx = μ2 / r30^3 * (-as - x) + (1-μ1) / r31^3 * (-μ1*cos(θ) - x)    + μ1 / r32^3 * ((1-μ1)*cos(θ) - x) + Tx
    Fy = μ2 / r30^3 * (-y)      + (1-μ1) / r31^3 * (-μ1*sin(θ) - y)    + μ1 / r32^3 * ((1-μ1)*sin(θ) - y) + Ty
    Fz = μ2 / r30^3 * (-z)      + (1-μ1) / r31^3 * (-z)                 + μ1 / r32^3 * (-z)               + Tz

    du[4] = 2*vy  + x + Fx
    du[5] = -2*vx + y + Fy
    du[6] =             Fz
end