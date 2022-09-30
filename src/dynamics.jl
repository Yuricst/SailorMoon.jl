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
    - `p`: parameters, where p = [μ, θ0, ωM]
    - `t`: time
"""
function rhs_cr3bp!(du,u,p,t)
end