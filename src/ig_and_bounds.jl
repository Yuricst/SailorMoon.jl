

function make_ig_bounds(row, τ_ig, n_arc::Int64)
    sv_mid_cart = [row.x_ra, row.y_ra, row.z_ra, row.xdot_ra, row.ydot_ra, row.zdot_ra]
    # change the coordinates into cylindrical (only position)
    svm_mid_cyl = vcat(SailorMoon.cart2cylind_only_pos(sv_mid_cart), row.m_ra)

    tof_leo2mid = row.dt2
    tof_mid2lpo = row.dt1
    rp = row.rp_kep
    ra = row.ra_kep
    α = row.alpha
    m_rp = row.m_rp

    θsf = row.thetaf
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
        [(rE+600)/param3b.lstar, 10.0,  pi, 1.2*m_rp, 1.3*tof_leo2mid/2],
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
        [pi, pi, 1.5, 1.3*tof_mid2lpo/2],
        vcat([[1.0,pi,pi] for i = 1:n_arc]...)
    )

    lx = vcat(lx_leo, lx_mid, lx_lpo)
    ux = vcat(ux_leo, ux_mid, ux_lpo)

    return x0, lx, ux
end
