"""
 set up the basic parameters for the multiple shooting
"""

abstract type AbstractParameterType end

Base.@kwdef struct multishoot_params <: AbstractParameterType
    dt::Real
    ballistic_time::Real
    ballistic_time_back::Real
    tmax::Real
    mdot::Real
    n_arc::Int
end


function multi_shoot_parameters(param3b::dynamics_params)
    tmax_si = 400e-3   # N
    isp_si  = 2500   # sec
    mdot_si = tmax_si / (isp_si * 9.81)
    mstar = 2500  # kg
    dt    = 0.005
    n_arc = 5

    ballistic_time      = 1 * 86400 / param3b.tstar
    ballistic_time_back = 10 * 86400 / param3b.tstar

    tmax = AstrodynamicsBase.dimensional2canonical_thrust(
        tmax_si, mstar, param3b.lstar, param3b.tstar
        )
    mdot = AstrodynamicsBase.dimensional2canonical_mdot(
            mdot_si, mstar, param3b.tstar
        )

    return multishoot_params(dt, ballistic_time, ballistic_time_back, tmax, mdot, n_arc)

end

lpo_filename = joinpath(dirname(@__FILE__), "data/lpo_L2_1200km_S.bson")
LPOArrival = load_lpo(lpo_filename, 2)

# dummy params for initialization
params = [
    param3b.mu2, param3b.mus, param3b.as, 0.0, param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 0.0, 0.0,
    dv_no_thrust
]
_prob_base = ODEProblem(rhs_bcr4bp_sb1frame2_thrust!, zeros(7,1), [0, -10.0], params);



"""
    unpacking variables along with multishoot_trajectory2
"""
function unpack_x2(x::AbstractVector{T}, n_arc::Int, verbose::Bool=false, scale::Bool=false) where T
    
    factor  = 10
    factort = 10

    # unpack
    nx = length(x)
    x_lr  = x[1 : 9+6*n_arc]
    x_mid = x[10+6*n_arc : 18+12*n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[19+12*n_arc : 22+15*n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    if scale
        # state
        x_lr[1:6]  = x_lr[1:6]  * factor
        x_mid[1:6] = x_mid[1:6] * factor
        x_lr[1]    = x_lr[1]  + param3b.as
        x_mid[1]   = x_mid[1] + param3b.as

        # time 
        x_lr[8:9]  = x_lr[8:9] * factor
        x_mid[8:9] = x_mid[8:9] * factor
        x_LPO[4]   = x_LPO[4] *factort   # lpo => mid direction tof

        # println("variables are scaled!")
    end

    # get time of flights
    tofs = [x_lr[8], x_lr[9], x_mid[8], x_mid[9], x_LPO[4]]
    θf = x_LPO[1]
    θs = [
        θf - param3b.oms*sum(broadcast(abs, tofs[end-3:end])),
        θf - param3b.oms*sum(broadcast(abs, tofs[end-1:end])),
        θf
    ]
    # print message
    if verbose
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
        # println("θf: ", θf)
    end
    return x_lr, x_mid, x_LPO, tofs, θs

end

"""
    as we scaled x0s, we would like to convert from old x to new x
"""
function oldx2newx(x, n_arc)
    factor = 10

    # unpack
    nx = length(x)
    x_lr  = x[1 : 9+6*n_arc]
    x_mid = x[10+6*n_arc : 18+12*n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[19+12*n_arc : 22+15*n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    # length 
    x_lr[1]    = x_lr[1]  - param3b.as
    x_mid[1]   = x_mid[1] - param3b.as
    x_lr[1:6]  = x_lr[1:6]  / factor
    x_mid[1:6] = x_mid[1:6] / factor

    # time
    x_lr[8:9]  = x_lr[8:9] / factor
    x_mid[8:9] = x_mid[8:9] / factor
    
    xnew = vcat(x_lr, x_mid, x_LPO)


end



function get_LEO_state(x_LEO, θs, param_multi::multishoot_params, verbose::Bool=false)
    rp_input, ra_input, α = x_LEO[1:3]
    sv0_sunb1 = vcat(
        paramIni_to_sb1(rp_input, α, ra_input, pi-θs[1], param3b.oml, param3b.mu2, param3b.as),
        x_LEO[4]
        )

    if verbose
        println("SMA [km]: ", (rp_input + ra_input)*param3b.lstar)
        println("ra  [km]: ", x_LEO[1]*param3b.lstar)
        println("Energy in inertial frame: ", -(param3b.mu2)/(rp_input + ra_input))
        println("sv0_sb1: ", (sv0_sunb1[1:3]-[param3b.as,0,0])*param3b.lstar, sv0_sunb1[4:6]*param3b.lstar/param3b.tstar)
    end


    # ballistic propagation with small time-steps
    params = [
        param3b.mu2, param3b.mus, param3b.as, pi-θs[1], 
        param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 
        param_multi.mdot, param_multi.tmax,
        dv_no_thrust
    ]

    _prob = remake(_prob_base; tspan=[0, param_multi.ballistic_time], u0 = sv0_sunb1, p = params)
    sol_ballistic_fwd = integrate_rk4(_prob, 0.001);
    θ_iter = θs[1] + param3b.oms*sol_ballistic_fwd.t[end]


    return sol_ballistic_fwd.u[end], θ_iter, sol_ballistic_fwd
end

function get_LPO_state(x_LPO, θs, param_multi::multishoot_params, verbose::Bool=false)
    xf = set_terminal_state2(x_LPO[2], θs[3], param3b, LPOArrival)

    # transform to Sun-B1 frame
    svf_sunb1 = vcat(
        transform_EMrot_to_SunB1(xf, pi - θs[3], param3b.oml, param3b.as),
        x_LPO[3]   # mf
    )

    # ballistic propagation with small time-steps
    params = [
        param3b.mu2, param3b.mus, param3b.as, pi - θs[3], 
        param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 
        param_multi.mdot, param_multi.tmax,
        dv_no_thrust
    ]
    # println("params: ", params)
    _prob = remake(
        _prob_base; tspan=[0,-param_multi.ballistic_time_back],
        u0 = svf_sunb1, p = params
    )
    # sol_ballistic_bck = integrate_rk4(_prob, 0.001);
    sol_ballistic_bck = solve(_prob, Tsit5(); reltol=1e-12, abstol=1e-12);

    θ_iter = θs[3] - param3b.oms * param_multi.ballistic_time_back

    if verbose
        println("svf_sunb1: \n", svf_sunb1)
    end
    return sol_ballistic_bck.u[end], θ_iter, sol_ballistic_bck
end

function propagate_arc!(sv0, θ0, tspan, x_control, dir_func, param_multi::multishoot_params, get_sols::Bool, sol_param_list, name::String)
    sv_iter = [el for el in sv0]
    θ_iter = θ0
    for i = 1:param_multi.n_arc
        τ, γ, β = x_control[1+3*(i-1) : 3*i]

        params = [
            param3b.mu2, param3b.mus, param3b.as, pi - θ_iter, 
            param3b.oml, param3b.omb, τ, γ, β, 
            param_multi.mdot, param_multi.tmax,
            dir_func
        ]

        # Tsit5()
        _prob = remake(_prob_base; tspan=tspan, u0=sv_iter, p=params)
        sol = solve(_prob, AutoTsit5(Rosenbrock23()); reltol=1e-12, abstol=1e-12)

        # rk4 in DifferentialEquations
        # _prob = remake(_prob_base; tspan=tspan, u0=sv_iter, p=params)
        # sol = solve(_prob, RK4(); dt=0.005, adaptive=false)

        # rk4() - original
        # _prob = remake(_prob_base; tspan=tspan, u0 = sv_iter, p = params)
        # sol = integrate_rk4(_prob, param_multi.dt);



        if get_sols
            push!(sol_param_list, [sol, params, name])
        end
        # update θ0 and sv0
        θ_iter += param3b.oms*sol.t[end]
        sv_iter = sol.u[end]
    end
    return sv_iter
end

## ---- mulitple shooting -----------------------------------------------------------

"""
    Updated version of multishoot_trajectory.
    arcs are: LEO <- x_lr -> <- apogee -> <- LPO
"""
function multishoot_trajectory2(
    x::AbstractVector{T},
    dir_func, 
    param_multi::multishoot_params,
    get_sols::Bool=false, 
    verbose::Bool=false,
    scale::Bool=false
    ) where T

    # unpack decision vector
    x_lr, x_mid, x_LPO, tofs, θs = unpack_x2(x, param_multi.n_arc, false, scale) 

    # initialize storage
    sol_param_list = []

    sv_lr = x_lr[1:6]  # state-vector at midpoint (position is cylindrical frame)
    svm_lr = vcat(sv_lr, x_lr[7])

    # propagate midpoint backward (-> LEO)
    # FIXME: we probably want to add "coasting" for the final TBD sec. (1 day?)
    svf_lr_bck = propagate_arc!(
        svm_lr, θs[1], [0, -tofs[1]/param_multi.n_arc], x_lr[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "xlr_bck_arc"
    )

    # propaagte midpoint forward
    svf_lr_fwd = propagate_arc!(
        svm_lr, θs[1], [0, tofs[2]/param_multi.n_arc], x_lr[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "xlr_fwd_arc"
    )

    # propagate midpoint backward
    sv0 = x_mid[1:6]  # state-vector at midpoint (position is cylindrical frame)
    # sv0_cart = cylind2cart_only_pos(sv0_cyl)  # convert to Cartesian coordinate
    svm0 = vcat(sv0, x_mid[7])

    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[3]/param_multi.n_arc], x_mid[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[4]/param_multi.n_arc], x_mid[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, param_multi, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, -(tofs[5] - param_multi.ballistic_time_back)/param_multi.n_arc], x_LPO[5 : end], 
        dir_func, param_multi,
        get_sols, sol_param_list, "lpo_arc"
    )

    # periapsis 
    rp_target = (6375 + 500) / param3b.lstar
    θs0 = θs[end] - param3b.oms * sum(tofs)
    θe0 = 2*pi - θs0    # earth angle at LEO

    # Earth -> SC vector
    sc_earth = [
        svf_lr_bck[1] - (param3b.as + param3b.mu2*cos(θe0)),
        svf_lr_bck[2] - param3b.mu2*sin(θe0),
        svf_lr_bck[3]
    ]
    peri_cond = norm(sc_earth) - rp_target
    # println("|sc_earth| = ", norm(sc_earth)*param3b.lstar)
    # println("alt: ", peri_cond)

    # SC needs to be launched tangentially to the Earth: dot(r_{E,SC}, v_SC) = 0 
    dep_LEO = sc_earth[1]*svf_lr_bck[4] + sc_earth[2]*svf_lr_bck[5] + sc_earth[3]*svf_lr_bck[6]
    # dep_LEO = 0.0

    # residuals  
    res = vcat(
        svf_mid_bck[1:2] - svf_lr_fwd[1:2],
        svf_mid_bck[4:5] - svf_lr_fwd[4:5],
        svf_mid_bck[7] - svf_lr_fwd[7],
        svf_lpo[1:2] - svf_mid_fwd[1:2],
        svf_lpo[4:5] - svf_mid_fwd[4:5],
        svf_lpo[7] - svf_mid_fwd[7],
        peri_cond, dep_LEO)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_bck], tofs
    end
end



"""
    Updated version of multishoot_trajectory2.
    fixed TOF
    arcs are: LEO <- x_lr -> <- apogee -> <- LPO
"""
function multishoot_trajectory4(
    x::AbstractVector{T},
    dir_func, 
    param_multi::multishoot_params,
    tof_target::Real,
    get_sols::Bool=false, 
    verbose::Bool=false,
    scale::Bool=false
    ) where T

    # unpack decision vector
    x_lr, x_mid, x_LPO, tofs, θs = unpack_x2(x, param_multi.n_arc, false, scale)
    
    # println("x_lr: ",  [round(el,digits=3) for el in x_lr])
    # println("x_mind: ", [round(el,digits=3) for el in x_mid])
    # println("x_LPO: ",  [round(el,digits=3) for el in x_LPO])
    # println("tofs: ",  [round(el,digits=3) for el in tofs])

    # initialize storage
    sol_param_list = []

    sv_lr = x_lr[1:6]  # state-vector at midpoint (position is cylindrical frame)
    svm_lr = vcat(sv_lr, x_lr[7])

    # propagate midpoint backward (-> LEO)
    # FIXME: we probably want to add "coasting" for the final TBD sec. (1 day?)
    svf_lr_bck = propagate_arc!(
        svm_lr, θs[1], [0, -tofs[1]/param_multi.n_arc], x_lr[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "xlr_bck_arc"
    )

    # propaagte midpoint forward
    svf_lr_fwd = propagate_arc!(
        svm_lr, θs[1], [0, tofs[2]/param_multi.n_arc], x_lr[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "xlr_fwd_arc"
    )

    # propagate midpoint backward
    sv0 = x_mid[1:6]  # state-vector at midpoint (position is cylindrical frame)
    # sv0_cart = cylind2cart_only_pos(sv0_cyl)  # convert to Cartesian coordinate
    svm0 = vcat(sv0, x_mid[7])

    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[3]/param_multi.n_arc], x_mid[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[4]/param_multi.n_arc], x_mid[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, param_multi, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, -(tofs[5] - param_multi.ballistic_time_back)/param_multi.n_arc], x_LPO[5 : end], 
        dir_func, param_multi,
        get_sols, sol_param_list, "lpo_arc"
    )

    # periapsis 
    alt = (6375 + 500) / param3b.lstar
    θs0 = θs[end] - param3b.oms * sum(tofs)
    θe0 = 2*pi - θs0    # earth angle at LEO
    sc_earth = [
        svf_lr_bck[1] - (param3b.as + param3b.mu2*cos(θe0)),
        svf_lr_bck[2] - param3b.mu2*sin(θe0),
        svf_lr_bck[3]
    ]
    peri_cond = norm(sc_earth) - alt
    # println("alt: ", peri_cond)

    # println(tofs)
    tof_cond = tof_target - sum(tofs)

    # LEO tangential departure condition
    dep_LEO = sc_earth[1]*svf_lr_bck[4] + sc_earth[2]*svf_lr_bck[5] + sc_earth[3]*svf_lr_bck[6]
    # dep_LEO = 0.0

    # residuals    
    res = vcat(
        svf_mid_bck[1:2] - svf_lr_fwd[1:2],
        svf_mid_bck[4:5] - svf_lr_fwd[4:5],
        svf_mid_bck[7] - svf_lr_fwd[7],
        svf_lpo[1:2] - svf_mid_fwd[1:2],
        svf_lpo[4:5] - svf_mid_fwd[4:5],
        svf_lpo[7] - svf_mid_fwd[7],
        peri_cond, tof_cond, dep_LEO)[:]

    # println(res)
    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_bck], tofs
    end
end


"""
    Updated version of multishoot_trajectory2.
    fixed m_LEO
    arcs are: LEO <- x_lr -> <- apogee -> <- LPO
"""
function multishoot_trajectory5(
    x::AbstractVector{T},
    dir_func, 
    param_multi::multishoot_params,
    mleo_target::Real,
    get_sols::Bool=false, 
    verbose::Bool=false,
    scale::Bool=false
    ) where T

    # unpack decision vector
    x_lr, x_mid, x_LPO, tofs, θs = unpack_x2(x, param_multi.n_arc, false, scale)

    println("x_lr: ", x_lr)
    println("x_mind: ", x_mid)
    println("x_LPO: ", x_LPO)
    println("tofs: ", tofs)

    # initialize storage
    sol_param_list = []

    sv_lr = x_lr[1:6]  # state-vector at midpoint (position is cylindrical frame)
    svm_lr = vcat(sv_lr, x_lr[7])

    # propagate midpoint backward (-> LEO)
    # FIXME: we probably want to add "coasting" for the final TBD sec. (1 day?)
    svf_lr_bck = propagate_arc!(
        svm_lr, θs[1], [0, -tofs[1]/param_multi.n_arc], x_lr[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "xlr_bck_arc"
    )

    # propaagte midpoint forward
    svf_lr_fwd = propagate_arc!(
        svm_lr, θs[1], [0, tofs[2]/param_multi.n_arc], x_lr[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "xlr_fwd_arc"
    )

    # propagate midpoint backward
    sv0 = x_mid[1:6]  # state-vector at midpoint (position is cylindrical frame)
    # sv0_cart = cylind2cart_only_pos(sv0_cyl)  # convert to Cartesian coordinate
    svm0 = vcat(sv0, x_mid[7])

    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[3]/param_multi.n_arc], x_mid[10 : 9+3*param_multi.n_arc], dir_func, param_multi, 
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[4]/param_multi.n_arc], x_mid[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, param_multi, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, -(tofs[5] - param_multi.ballistic_time_back)/param_multi.n_arc], x_LPO[5 : end], 
        dir_func, param_multi,
        get_sols, sol_param_list, "lpo_arc"
    )

    # periapsis 
    alt = (6375 + 500) / param3b.lstar
    θs0 = θs[end] - param3b.oms * sum(tofs)
    θe0 = 2*pi - θs0    # earth angle at LEO
    sc_earth = [
        svf_lr_bck[1] - (param3b.as + param3b.mu2*cos(θe0)),
        svf_lr_bck[2] - param3b.mu2*sin(θe0),
        svf_lr_bck[3]
    ]
    peri_cond = norm(sc_earth) - alt
    # println("alt: ", peri_cond)

    m_leo_cond = mleo_target - svf_lr_bck[7]

    # LEO tangential departure condition
    dep_LEO = sc_earth[1]*svf_lr_bck[4] + sc_earth[2]*svf_lr_bck[5] + sc_earth[3]*svf_lr_bck[6]
    # dep_LEO = 0.0

    # residuals    
    res = vcat(
        svf_mid_bck[1:2] - svf_lr_fwd[1:2],
        svf_mid_bck[4:5] - svf_lr_fwd[4:5],
        svf_mid_bck[7] - svf_lr_fwd[7],
        svf_lpo[1:2] - svf_mid_fwd[1:2],
        svf_lpo[4:5] - svf_mid_fwd[4:5],
        svf_lpo[7] - svf_mid_fwd[7],
        peri_cond, m_leo_cond, dep_LEO)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_bck], tofs
    end
end

# function to retreive the time series data (t, u(t)) from x0 
function x2time_series(
    x::Vector, 
    dir_func, 
    param_multi::multishoot_params,
    scale::Bool=false,
    )

    θm_lpo = x[19+12*param_multi.n_arc]  
    u  = []
    t  = []
    th = []

    x_lr, x_mid, x_LPO, tofs, θs = unpack_x2(x, param_multi.n_arc, false, scale)

    res, sol_param_list, sols_ballistic, tofs = multishoot_trajectory2(x, dir_func, param_multi, true, false, scale)
        
    # ballistic legs
    for sol_ballistic in sols_ballistic
        u = sol_ballistic.u[:]
        t = sol_ballistic.t[:]
        th = [[0,0,0] for i in collect(1:size(u,1))]

    end
    
    # nonbalistic legs: lr_bck -> lr_fwd -> mid_bck -> mid_fwd -> lpo_bck
    for j = 1:Int(floor(length(sol_param_list)/param_multi.n_arc))
        
        for k = 1:param_multi.n_arc
            
            if mod(j,2) == 1
                # backward propagation
                sol, _, name = sol_param_list[length(sol_param_list) - j*param_multi.n_arc + k]
                t_append = sol.t[:] .+ t[end]
                u_append = sol.u[:]

                # somehow, get a thrust parameter (1x3 double)
                if j == 5      # lr_bck
                    params = x_lr[10 : 9+3*param_multi.n_arc]
                elseif j == 3  # mid_bck
                    params = x_mid[10 : 9+3*param_multi.n_arc]                    
                elseif j == 1  # lpo_bck
                    params = x_LPO[5 : end]
                end 

                # println("test: ", params[13:15])
                thrust_angle = params[3*(param_multi.n_arc-k+1)-2 : 3*(param_multi.n_arc-k+1)]
                
            else
                # forward propagation
                sol, _, name = sol_param_list[length(sol_param_list)-(j-1)*param_multi.n_arc - k + 1]
                u_append = sol.u[end:-1:1, :]
                t_append = sol.t[end:-1:1] .- sol.t[end] .+ t[end]

                # somehow, get a thrust parameter (1x3 double)
                if j == 4      # lr_fwd
                    params = x_lr[10+3*param_multi.n_arc : end]
                elseif j == 2  # mid_fwd
                    params = x_mid[10+3*param_multi.n_arc : end]
                end 
                
                thrust_angle = params[3*k-2 : 3*k]

            end

            th_append = get_thrust(t_append, u_append, thrust_angle, θm_lpo, dir_func)
            u  = vcat(u, u_append)
            t  = vcat(t, t_append)
            th = vcat(th, th_append)

        end
        
    end
    
    u = Base.Array(u)
    u = hcat(u...)[:,:]
    th = hcat(th...)

    return t, u, th
end

# thrust_param: [τ, γ, β]
function get_thrust(tlist, ulist, thrust_param, θm_lpo, dir_func)
    thrusts = []
    for (idx, tnow) in enumerate(tlist)
        θm    = θm_lpo - param3b.oml*tnow
        unow  = ulist[idx]
        thvec = dir_func(param3b.mus, param3b.as, θm, param3b.oml, unow, thrust_param)
        # thrusts = vcat(thrusts, transpose(thvec))
        push!(thrusts, thvec)

    end 

    return thrusts
end





## ----- Not used as of now --------------------------------------------

function unpack_x(x::AbstractVector{T}, n_arc::Int, verbose::Bool=false) where T
    # unpack
    nx = length(x)
    x_LEO = x[1:5+3n_arc]
    x_mid = x[6+3n_arc:14+9n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[15+9n_arc:18+12n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    # get time of flights
    tofs = [x_LEO[5], x_mid[8], x_mid[9], x_LPO[4]]
    θf = x_LPO[1]
    θs = [
        θf - param3b.oms*sum(broadcast(abs, tofs)),
        θf - param3b.oms*sum(broadcast(abs, tofs[3:4])),
        θf
    ]
    # print message
    if verbose
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
        println("θf: ", θf)
    end
    return x_LEO, x_mid, x_LPO, tofs, θs
end

function unpack_x2plus(x::AbstractVector{T}, n_arc::Int, verbose::Bool=false) where T
    # unpack
    nx = length(x)
    x_lr  = x[1:9+6*n_arc]
    x_mid = x[10+6*n_arc:18+12*n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[19+12*n_arc:22+15*n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    # get time of flights
    tofs = [x_lr[8], x_lr[9], x_mid[8] - (x_lr[8]+x_lr[9]+x_mid[9]+x_LPO[4]), x_mid[9], x_LPO[4]]
    θf = x_LPO[1]
    θs = [
        θf - param3b.oms*sum(broadcast(abs, tofs[end-3:end])),
        θf - param3b.oms*sum(broadcast(abs, tofs[end-1:end])),
        θf
    ]
    # print message
    if verbose
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
        # println("θf: ", θf)
    end
    return x_lr, x_mid, x_LPO, tofs, θs
end

function unpack_x3(x::AbstractVector{T}, n_arc::Int, verbose::Bool=false) where T
    # unpack
    nx = length(x)
    x_LEO = x[1 : 5+3*n_arc]
    x_lr  = x[6+3*n_arc : 13+6*n_arc]
    x_mid = x[14+6*n_arc : 22+12*n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[23+12*n_arc : 26+15*n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    # get time of flights
    tofs = [x_LEO[5], x_lr[8], x_mid[8], x_mid[9], x_LPO[4]]
    θf = x_LPO[1]
    θs = [
        θf - param3b.oms*sum(broadcast(abs, tofs[end-3:end])),
        θf - param3b.oms*sum(broadcast(abs, tofs[end-1:end])),
        θf
    ]
    # print message
    if verbose
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
        # println("θf: ", θf)
    end
    return x_LEO, x_lr, x_mid, x_LPO, tofs, θs
end


"""
    arc: LEO -> <- apogee -> <- LPO
"""

function multishoot_trajectory(
    x::AbstractVector{T},
    dir_func, 
    param_multi::multishoot_params,
    get_sols::Bool=false, 
    verbose::Bool=false
    ) where T
    
    # unpack decision vector
    x_LEO, x_mid, x_LPO, tofs, θs = unpack_x(x, param_multi.n_arc)

    # initialize storage
    sol_param_list = []

    # propagate from LEO forward
    sv0_LEO, θ0_leo, sol_ballistic_fwd = get_LEO_state(x_LEO, θs, ballistic_time, verbose)
    svf_LEO = propagate_arc!(
        sv0_LEO, θ0_leo, [0, (tofs[1] - param_multi.ballistic_time)/param_multi.n_arc], x_LEO[5 : end], dir_func, param_multi,
        get_sols, sol_param_list, "leo_arc"
    )

    # propagate midpoint backward
    sv0_cyl = x_mid[1:6]  # state-vector at midpoint (position is cylindrical frame)
    sv0_cart = cylind2cart_only_pos(sv0_cyl)  # convert to Cartesian coordinate
    svm0 = vcat(sv0_cart, x_mid[7])

    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[2]/param_multi.n_arc], x_mid[10 : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[3]/param_multi.n_arc], x_mid[10+3n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, param_multi, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, (-tofs[4] + param_multi.ballistic_time_back)/param_multi.n_arc], 
        x_LPO[5 : end], dir_func, param_multi,
        get_sols, sol_param_list, "lpo_arc"
    )

    # residuals
    res = vcat(svf_mid_bck - svf_LEO, svf_lpo - svf_mid_fwd)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_fwd, sol_ballistic_bck], tofs
    end
end


"""
    arcs are: LEO -> x_lr -> <- apogee -> <- LPO
"""
function multishoot_trajectory3(
    x::AbstractVector{T},
    dir_func, 
    param_multi::multishoot_params,
    get_sols::Bool=false, 
    verbose::Bool=false
    ) where T

    # unpack decision vector
    x_LEO, x_lr, x_mid, x_LPO, tofs, θs = unpack_x3(x, param_multi.n_arc)
    # initialize storage
    sol_param_list = []

    # propagate from LEO forward
    sv0_LEO, θ0_leo, sol_ballistic_fwd = get_LEO_state(x_LEO, θs, param_multi, verbose)
    svf_LEO = propagate_arc!(
        sv0_LEO, θ0_leo, [0, (tofs[1] - param_multi.ballistic_time)/param_multi.n_arc], x_LEO[6 : end], dir_func, param_multi,
        get_sols, sol_param_list, "leo_arc"
    )

    sv_lr = x_lr[1:6]  # state-vector at midpoint (position is cylindrical frame)
    svm_lr = vcat(sv_lr, x_lr[7])

    # propaagte midpox_lrint forward
    svf_lr_fwd = propagate_arc!(
        svm_lr, θs[1], [0, tofs[2]/param_multi.n_arc], x_lr[9 : end], dir_func, param_multi,
        get_sols, sol_param_list, "xlr_fwd_arc"
    )

    # propagate midpoint backward
    sv0_cyl = x_mid[1:6]  # state-vector at midpoint (position is cylindrical frame)
    sv0_cart = cylind2cart_only_pos(sv0_cyl)  # convert to Cartesian coordinate
    svm0 = vcat(sv0_cart, x_mid[7])

    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[3]/param_multi.n_arc], x_mid[10 : 9+3*param_multi.n_arc], dir_func, param_multi,
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[4]/param_multi.n_arc], x_mid[10+3*param_multi.n_arc : end], dir_func, param_multi,
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, param_multi, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, (-tofs[4] + param_multi.ballistic_time_back)/param_multi.n_arc], x_LPO[5 : end], dir_func, param_multi,
        get_sols, sol_param_list, "lpo_arc"
    )

    # residuals
    res = vcat(svf_LEO - svm_lr, svf_mid_bck - svf_lr_fwd, svf_lpo - svf_mid_fwd)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_bck], tofs
    end
end


