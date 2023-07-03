"""
    Setting up the initial guess and ub/lb of the variables
"""


"""
    updated version along with multishoot_trajectory2
"""
function make_ig_bounds2(row, τ_ig, n_arc::Int64, scale::Bool=false)
    # change the coordinates into cylindrical (only position)
    # svm_mid_cyl = vcat(cart2cylind_only_pos(sv_mid_cart), row.m_ra)


    tof_leo2mid = row.dt2
    tof_mid2lpo = row.dt1
    m_rp  = row.m_rp
    m_lpo = row.m_lpo
    tof   = row.tof

    θsf = row.thetasf
    ϕ0  = row.phi0

    x_lr    = row.x_lr
    y_lr    = row.y_lr
    z_lr    = row.z_lr
    xdot_lr = row.xdot_lr
    ydot_lr = row.ydot_lr
    zdot_lr = row.zdot_lr
    m_lr    = row.m_lr
    t_lr    = row.t_lr

    if scale
        factor  = 10
        factort = 10
        # scaling problem: special treatment for x_lr and x_ra
        x_lr = x_lr - param3b.as
        svm_mid_cart = [
            row.x_ra-param3b.as, row.y_ra, 0.0, 
            row.xdot_ra, row.ydot_ra, 0.0, row.m_ra
            ] ./ factor

    else 
        factor  = 1
        factort = 1
        svm_mid_cart = [row.x_ra, row.y_ra, 0.0, row.xdot_ra, row.ydot_ra, 0.0, row.m_ra] 
    end

    rE = 6375 # km

    # x_lr = [x,y,z,vx,vy,vz, m, tof_back, tof_fwd, controls...] (9 + 6*n_arc)
    ig_x_lr = vcat(
        x_lr/factor, y_lr/factor, 0.0, xdot_lr/factor, ydot_lr/factor, 0.0, m_lr,
        (tof - t_lr)/factort, (tof_leo2mid/2 + t_lr - tof)/factort,
        vcat([[τ_ig,0,0] for i = 1:2*n_arc]...)
    )

    # x_mid = [r,theta,z, vx,vy,vz, m, tof_back, tof_forward, controls...] (9 + 6*n_arc)
    ig_x_mid = vcat(
        svm_mid_cart, 
        tof_leo2mid/2/factort, tof_mid2lpo/2/factort, 
        vcat([[τ_ig,0,0] for i = 1:2n_arc]...)
    )

    # x_LPO = [θf, ϕ, mf, tof, controls...] (4 + 3*n_arc)
    ig_x_LPO = vcat(
        [θsf, ϕ0, m_lpo, tof_mid2lpo/2/factort],
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    x0 = vcat(ig_x_lr, ig_x_mid, ig_x_LPO)

    ### lb, ub of variables 
    # for the lyapunov orbit, everything is 2D, so just set β=0 in [τ, γ, β]

    dstate = 0.3 / factor 

    lx_lr = vcat(
        ig_x_lr[1] - dstate, ig_x_lr[2] - dstate, 0.0, -1.5/factor, -1.5/factor, -0.0, 1.0,
        0.8*(tof - t_lr)/factort, 0.8*(tof_leo2mid/2 + t_lr - tof)/factort,
        vcat([[0.0, 0.0, 0.0] for i = 1:n_arc]..., [[0.0, -pi, 0.0] for i = 1:n_arc]...)
    )

    ux_lr = vcat(
        ig_x_lr[1] + dstate, ig_x_lr[2] + dstate, 0, 1.5/factor, 1.5/factor, 0.0, 1.3,
        1.2*(tof - t_lr)/factort, 1.2*(tof_leo2mid/2 + t_lr - tof)/factort,
        vcat([[0.0, 0.0, 0.0] for i = 1:n_arc]..., [[1.0, pi, 0.0] for i = 1:n_arc]...)
    )

    lx_mid = vcat(
        ig_x_mid[1] - dstate, ig_x_mid[2] - dstate, 0.0, -1.0/factor, -1.0/factor, 0.0, 1.0, 
        0.7*tof_leo2mid/2/factort, 0.7*tof_mid2lpo/2/factort, 
        vcat([[0.0, -pi, 0.0] for i = 1:2n_arc]...)
    )
    ux_mid = vcat(
        ig_x_mid[1] + dstate, ig_x_mid[2] + dstate, 0.0, 1.0/factor, 1.0/factor, 0.0, 1.2,
        1.3*tof_leo2mid/2/factort, 1.3*tof_mid2lpo/2/factort, 
        vcat([[1.0, pi, 0.0] for i = 1:2n_arc]...)
    )

    lx_lpo = vcat(
        [θsf-pi/4, ϕ0-pi/4, 1.000, 0.7*tof_mid2lpo/2/factort],
        vcat([[0.0, -pi, 0.0] for i = 1:n_arc]...)
    )
    ux_lpo = vcat(
        [θsf+pi/4, ϕ0+pi/4, 1.000, 1.3*tof_mid2lpo/2/factort],
        vcat([[1.0, pi, 0.0] for i = 1:n_arc]...)
    )

    lx = vcat(lx_lr, lx_mid, lx_lpo)
    ux = vcat(ux_lr, ux_mid, ux_lpo)

    return x0, lx, ux
end


"""
    for multishoot_trajectory2. 
    The input row[3:end] is the raw variable of x0,x1,...xN.
    used in opt_from_output.jl
"""
function make_ig_bounds2_raw(row, τ_ig, n_arc::Int64, scale::Bool=false)

    x0 = collect(values(row[4:end])) 

    # svm_mid_cyl = x0[10+6*n_arc:16+6*n_arc]

    tof_leo2mid0 = x0[8]   # leo -> lr direction tof 
    tof_leo2mid1 = x0[9]   # lr -> mid direction tof 
    tof_leo2mid2 = x0[17+6*n_arc]   # mid -> leo direction tof 
    tof_mid2lpo1 = x0[18+6*n_arc]   # mid -> lpo direction tof
    tof_mid2lpo2 = x0[22+12*n_arc]   # lpo => mid direction tof 
    tof = tof_leo2mid0 + tof_leo2mid1 + tof_leo2mid2 + tof_mid2lpo1 + tof_mid2lpo2

    θsf  = x0[19+12*n_arc]
    ϕ0   = x0[20+12*n_arc]
    m_rp = x0[21+12*n_arc]

    x_lr    = x0[1]
    y_lr    = x0[2]
    z_lr    = x0[3]
    xdot_lr = x0[4]
    ydot_lr = x0[5]
    zdot_lr = x0[6]
    m_lr    = x0[7]
    # t_lr    = row.t_lr

    x_ra    = x0[10+6*n_arc]
    y_ra    = x0[11+6*n_arc]

    rE = 6375 # km

    if scale
        factor  = 10
        factort = 10 

        # state
        x_lr = x_lr - param3b.as
        x_ra = x_ra - param3b.as
        x0[1]          = x_lr  
        x0[10+6*n_arc] = x_ra 
        x0[1:6] = x0[1:6] / factor
        x0[10+6n_arc:15+6n_arc] = x0[10+6n_arc:15+6n_arc] / factor

        # time
        x0[8:9] = x0[8:9]./factort   
        x0[17+6*n_arc:18+6*n_arc] = x0[17+6*n_arc:18+6*n_arc] ./ factort  
        x0[22+12*n_arc] = x0[22+12*n_arc] ./ factort
    else
        factor  = 1
        factort = 1
    end

    dstate = 0.3/factor

    ### lb, ub of variables 
    lx_lr = vcat(
        x0[1] - dstate, x0[2] - dstate, 0, -1.5/factor, -1.5/factor, 0.0, 1.0,
        0.9*tof_leo2mid0/factort, 0.8*tof_leo2mid1/factort,
        vcat([[0.0,0.0,0.0] for i = 1:2n_arc]...)
    )
    ux_lr = vcat(
        x0[1] + dstate, x0[2] + dstate, 0, 1.5/factor, 1.5/factor, 0.0, 1.3,
        1.1*tof_leo2mid0/factort, 1.2*tof_leo2mid1/factort,
        vcat([[0.0,0.0,0.0] for i = 1:2n_arc]...)
    )

    lx_mid = vcat(
        x0[10+6*n_arc] - dstate,  x0[11+6*n_arc]-dstate, 0.0, -1.0/factor, -1.0/factor, 0.0, 1.0, 
        0.8*tof_leo2mid2/factort, 0.8*tof_mid2lpo1/factort, 
        vcat([[0.0,-pi,0.0] for i = 1:2n_arc]...)
    )
    ux_mid = vcat(
        x0[10+6*n_arc] + dstate,  x0[11+6*n_arc] + dstate, 0.0, 1.0/factor, 1.0/factor, 0.0, 1.2, 
        1.2*tof_leo2mid2/factort, 1.2*tof_mid2lpo1/factort,
        vcat([[1.0,pi,0.0] for i = 1:2n_arc]...)
    )

    lx_lpo = vcat(
        [θsf-pi/4, ϕ0-pi/4, 1.000, 0.7*tof_mid2lpo2/factort],
        vcat([[0.0,-pi,0.0] for i = 1:n_arc]...)
    )
    ux_lpo = vcat(
        [θsf+pi/4, ϕ0+pi/4, 1.000, 1.3*tof_mid2lpo2/factort],
        vcat([[1.0,pi,0.0] for i = 1:n_arc]...)
    )

    lx = vcat(lx_lr, lx_mid, lx_lpo)
    ux = vcat(ux_lr, ux_mid, ux_lpo)


    return x0, lx, ux
end




# ======== Not used right now =================================================================================

function make_ig_bounds(row, τ_ig, n_arc::Int64)
    sv_mid_cart = [row.x_ra, row.y_ra, row.z_ra, row.xdot_ra, row.ydot_ra, row.zdot_ra]
    # change the coordinates into cylindrical (only position)
    svm_mid_cyl = vcat(cart2cylind_only_pos(sv_mid_cart), row.m_ra)

    tof_leo2mid = row.dt2
    tof_mid2lpo = row.dt1
    rp = row.rp_kep
    ra = row.ra_kep
    α = row.alpha
    m_rp = row.m_rp

    θsf = row.thetasf
    ϕ0  = row.phi0

    rE = 6375 # km

    # x_LEO = [rp, ra, α, m0, tof, controls...]
    ig_x_LEO = vcat(
        [rp, ra,  α, 1.0, tof_leo2mid/2],
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    # x_mid = [r,theta,z, vx,vy,vz, m, tof_back, tof_forward, controls...]
    ig_x_mid = vcat(
        svm_mid_cyl, tof_leo2mid/2, tof_mid2lpo/2, 
        vcat([[τ_ig,0,0] for i = 1:2n_arc]...)
    )

    # x_LPO = [θf, ϕ, mf, tof, controls...]
    ig_x_LPO = vcat(
        [θsf, ϕ0, m_rp, tof_mid2lpo/2],
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    x0 = vcat(ig_x_LEO, ig_x_mid, ig_x_LPO)


    ### lb, ub of variables 
    lx_leo = vcat(
        [(rE+200)/param3b.lstar, 0.0,  -pi, 0.7*m_rp, 0.7*tof_leo2mid/2],
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )

    ux_leo = vcat(
        [(rE+2000)/param3b.lstar, 10.0,  pi, 1.2*m_rp, 1.3*tof_leo2mid/2],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx_mid = vcat(
        0.8*svm_mid_cyl[1], (svm_mid_cyl[2]-pi/12), -1.0, -2.0, -2.0, -0.0, 1.0, 
        0.7*tof_leo2mid/2, 0.7*tof_mid2lpo/2, 
        vcat([[0.0,-pi,-pi] for i = 1:2n_arc]...)
    )
    ux_mid = vcat(
        1.2*svm_mid_cyl[1], (svm_mid_cyl[2]+pi/12), 1.0, 1.0, 1.0, 0.0, 1.2*m_rp,
        1.3*tof_leo2mid/2, 1.3*tof_mid2lpo/2, 
        vcat([[1.0,0,0] for i = 1:2n_arc]...)
    )

    # for psi and theta, initial guess is generated in range of 0 ~ 2pi -> -pi ~ 3pi ub/lb
    lx_lpo = vcat(
        [-pi, -pi, 0.0, 0.7*tof_mid2lpo/2],
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )
    ux_lpo = vcat(
        [3*pi, 3*pi, 1.5, 1.3*tof_mid2lpo/2],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx = vcat(lx_leo, lx_mid, lx_lpo)
    ux = vcat(ux_leo, ux_mid, ux_lpo)

    return x0, lx, ux
end


"""
    updated version along with multishoot_trajectory3
"""
function make_ig_bounds3(row, τ_ig, n_arc::Int64)
    sv_mid_cart = [row.x_ra, row.y_ra, row.z_ra, row.xdot_ra, row.ydot_ra, row.zdot_ra]
    # change the coordinates into cylindrical (only position)
    svm_mid_cyl = vcat(cart2cylind_only_pos(sv_mid_cart), row.m_ra)

    tof_leo2mid = row.dt2
    tof_mid2lpo = row.dt1
    rp    = row.rp_kep
    ra    = row.ra_kep
    α     = row.alpha
    m_rp  = row.m_rp
    tof   = row.tof

    θsf = row.thetasf
    ϕ0  = row.phi0

    x_lr    = row.x_lr
    y_lr    = row.y_lr
    z_lr    = row.z_lr
    xdot_lr = row.xdot_lr
    ydot_lr = row.ydot_lr
    zdot_lr = row.zdot_lr
    m_lr    = row.m_lr
    t_lr    = row.t_lr

    rE = 6375 # km

    # x_LEO = [rp, ra, α, m0, tof_fwd, controls...]  (5 + 3*n_arc)
    ig_x_LEO = vcat(
        rp, ra,  α, 1.0, tof - t_lr,
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    # x_lr = [x,y,z,vx,vy,vz, m, tof_fwd, controls...] (8 + 3*n_arc)
    ig_x_lr = vcat(
        x_lr, y_lr, z_lr, xdot_lr, ydot_lr, zdot_lr, 1.0,
        tof_leo2mid/2 + t_lr - tof,
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    # x_mid = [r,theta,z, vx,vy,vz, m, tof_back, tof_fwd, controls...] (9 + 6*n_arc)
    ig_x_mid = vcat(
        svm_mid_cyl, tof_leo2mid/2, tof_mid2lpo/2, 
        vcat([[τ_ig,0,0] for i = 1:2*n_arc]...)
    )

    # x_LPO = [θf, ϕ, mf, tof_back, controls...] (4 + 3*n_arc)
    ig_x_LPO = vcat(
        [θsf, ϕ0, m_rp, tof_mid2lpo/2],
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    x0 = vcat(ig_x_LEO, ig_x_lr, ig_x_mid, ig_x_LPO)

    ### lb, ub of variables 
    lx_leo = vcat(
        [(rE+200)/param3b.lstar, 0.0,  -pi, 0.7*m_rp, 0.8*(tof - t_lr)],
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )

    ux_leo = vcat(
        [(rE+2000)/param3b.lstar, 10.0,  pi, 1.2*m_rp, 1.2*(tof_leo2mid/2)],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx_lr = vcat(
        x_lr - 0.3, y_lr - 0.3, z_lr - 0.2 , -Inf, -Inf, -Inf, 1.0,
        0.8*(tof_leo2mid/2 + t_lr - tof),
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )

    ux_lr = vcat(
        x_lr + 0.3, y_lr + 0.3, z_lr + 0.2 , Inf, Inf, Inf, 1.2*m_rp,
        1.2*(tof_leo2mid/2 + t_lr - tof),
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx_mid = vcat(
        0.8*svm_mid_cyl[1], (svm_mid_cyl[2]-pi/12), -1.0, -Inf, -Inf, -Inf, 1.0, 
        0.7*tof_leo2mid/2, 0.7*tof_mid2lpo/2, 
        vcat([[0.0,-pi,-pi] for i = 1:2n_arc]...)
    )
    ux_mid = vcat(
        1.2*svm_mid_cyl[1], (svm_mid_cyl[2]+pi/12), 1.0, Inf, Inf, Inf, 1.2*m_rp,
        1.3*tof_leo2mid/2, 1.3*tof_mid2lpo/2, 
        vcat([[1.0,0,0] for i = 1:2n_arc]...)
    )

    lx_lpo = vcat(
        [-pi, -pi, 0.0, 0.7*tof_mid2lpo/2],
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )
    ux_lpo = vcat(
        [3*pi, 3*pi, 1.5, 1.3*tof_mid2lpo/2],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx = vcat(lx_leo, lx_lr, lx_mid, lx_lpo)
    ux = vcat(ux_leo, ux_lr, ux_mid, ux_lpo)

    return x0, lx, ux
end



"""
    updated version along with multishoot_trajectory2. 
    including tof as an explicit variable (did not go well. 5/19/2023)
"""
function make_ig_bounds2plus(row, τ_ig, n_arc::Int64)
    sv_mid_cart = [row.x_ra, row.y_ra, row.z_ra, row.xdot_ra, row.ydot_ra, row.zdot_ra]
    # change the coordinates into cylindrical (only position)
    svm_mid_cyl = vcat(cart2cylind_only_pos(sv_mid_cart), row.m_ra)

    tof_leo2mid = row.dt2
    tof_mid2lpo = row.dt1
    m_rp  = row.m_rp
    m_lpo = row.m_lpo
    tof   = row.tof

    θsf = row.thetasf
    ϕ0  = row.phi0

    x_lr    = row.x_lr
    y_lr    = row.y_lr
    z_lr    = row.z_lr
    xdot_lr = row.xdot_lr
    ydot_lr = row.ydot_lr
    zdot_lr = row.zdot_lr
    m_lr    = row.m_lr
    t_lr    = row.t_lr

    rE = 6375 # km

    # x_lr = [x,y,z,vx,vy,vz, m, tof_back, tof_fwd, controls...] (9 + 6*n_arc)
    ig_x_lr = vcat(
        x_lr, y_lr, z_lr, xdot_lr, ydot_lr, zdot_lr, m_lr,
        tof - t_lr, tof_leo2mid/2 + t_lr - tof,
        vcat([[τ_ig,0,0] for i = 1:2*n_arc]...)
    )

    # x_mid = [r,theta,z, vx,vy,vz, m, tof_back, tof_forward, controls...] (9 + 6*n_arc)
    ig_x_mid = vcat(
        svm_mid_cyl, tof, tof_mid2lpo/2, 
        vcat([[τ_ig,0,0] for i = 1:2n_arc]...)
    )

    # x_LPO = [θf, ϕ, mf, tof, controls...] (4 + 3*n_arc)
    ig_x_LPO = vcat(
        [θsf, ϕ0, m_lpo, tof_mid2lpo/2],
        vcat([[τ_ig,0,0] for i = 1:n_arc]...)
    )

    x0 = vcat(ig_x_lr, ig_x_mid, ig_x_LPO)

    ### lb, ub of variables 
    lx_lr = vcat(
        x_lr - 0.3, y_lr - 0.3, z_lr - 0.2 , -Inf, -Inf, -Inf, 1.0,
        0.8*(tof - t_lr), 0.8*(tof_leo2mid/2 + t_lr - tof),
        vcat([[0.0,-pi,-pi] for i = 1:2n_arc]...)
    )

    ux_lr = vcat(
        x_lr + 0.3, y_lr + 0.3, z_lr + 0.2 , Inf, Inf, Inf, 1.5*m_lr,
        1.2*(tof - t_lr), 1.2*(tof_leo2mid/2 + t_lr - tof),
        vcat([[1.0,pi,pi] for i = 1:2n_arc]...)
    )

    lx_mid = vcat(
        0.8*svm_mid_cyl[1], (svm_mid_cyl[2]-pi/12), -1.0, -Inf, -Inf, -Inf, 1.0, 
        0.8*tof, 0.7*tof_mid2lpo/2, 
        vcat([[0.0,-pi,-pi] for i = 1:2n_arc]...)
    )
    ux_mid = vcat(
        1.2*svm_mid_cyl[1], (svm_mid_cyl[2]+pi/12), 1.0, Inf, Inf, Inf, 1.2*m_rp,
        1.2*tof, 1.3*tof_mid2lpo/2, 
        vcat([[1.0,pi,pi] for i = 1:2n_arc]...)
    )

    lx_lpo = vcat(
        [-pi, -pi, 1.0, 0.7*tof_mid2lpo/2],
        vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
    )
    ux_lpo = vcat(
        [3*pi, 3*pi, 1.0, 1.3*tof_mid2lpo/2],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx = vcat(lx_lr, lx_mid, lx_lpo)
    ux = vcat(ux_lr, ux_mid, ux_lpo)

    return x0, lx, ux
end
