using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
# import joptimise
using Printf
using JSON

include("../../julia-r3bp/R3BP/src/R3BP.jl")

# set up the basic parameters 
param3b = dyanmics_parameters()
tmax_si = 280e-3 * 4   # N
isp_si = 1800   # sec
mdot_si = tmax_si / (isp_si * 9.81)
mstar = 2500  # kg

tmax = AstrodynamicsBase.dimensional2canonical_thrust(
    tmax_si, mstar, param3b.lstar, param3b.tstar
)
mdot = AstrodynamicsBase.dimensional2canonical_mdot(
    mdot_si, mstar, param3b.tstar
)


function unpack_x(x::AbstractVector{T}, verbose::Bool=true) where T
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

function get_LEO_state(x_LEO, θs, verbose::Bool=false) 
    rp_input, ra_input, α = x_LEO[1:3]
    sv0_sunb1 = vcat(
        paramIni_to_sb1(rp, α, ra, θm0, param3b.oml, param3b.mu2, param3b.as), 
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
        param3b.mu2, param3b.mus, param3b.as, pi-θs[1], param3b.oml, param3b.omb, 0.0, 0.0, 0.0, mdot, tmax,
        dv_sun_dir_angles
    ]
    
    _prob = remake(_prob_base; tspan=[0,ballistic_time], u0 = sv0_sunb1, p = params)
    sol_ballistic_fwd = integrate_rk4(_prob, 0.001);
    θ_iter = θs[1] + param3b.oms*sol_ballistic_fwd.t[end]

    
    return sol_ballistic_fwd.u[end], θ_iter, sol_ballistic_fwd
end

function get_LPO_state(x_LPO, θs, verbose::Bool=false)
    xf = set_terminal_state2(x_LPO[2], θs[3], param3b, LPOArrival)

    # transform to Sun-B1 frame
    svf_sunb1 = vcat(
        transform_EMrot_to_SunB1(xf, θs[3], param3b.oms, param3b.as),
        x_LPO[3]   # mf
    )    
        
    # ballistic propagation with small time-steps
    params = [
        param3b.mu2, param3b.mus, param3b.as, pi - θs[3], param3b.oml, param3b.omb, 0.0, 0.0, 0.0, mdot, tmax,
        dv_no_thrust
    ]
    println("params: ", params)
    _prob = remake(
        _prob_base; tspan=[0,-ballistic_time_back], 
        u0 = svf_sunb1, p = params
    )
    sol_ballistic_bck = integrate_rk4(_prob, 0.001);
    θ_iter = θs[3] - param3b.oms*ballistic_time_back
    
    if verbose
        println("svf_sunb1: \n", svf_sunb1)
    end
    return sol_ballistic_bck.u[end], θ_iter, sol_ballistic_bck
end

function propagate_arc!(sv0, θ0, tspan, dt, x_control, get_sols::Bool, sol_param_list, name::String)
    sv_iter = [el for el in sv0]
    θ_iter = 1*θ0
    for i = 1:n_arc
        τ, γ, β = x_control[1+3*(i-1) : 3*i]
        
        params = [
            param3b.mu2, param3b.mus, param3b.as, pi - θ_iter, param3b.oml, param3b.omb, τ, γ, β, mdot, tmax,
            dv_sun_dir_angles
        ]
        _prob = remake(_prob_base; tspan=tspan, u0 = sv_iter, p = params)
        sol = integrate_rk4(_prob, dt);
        #sol = DifferentialEquations.solve(_prob, RK4(), reltol=1e-10, abstol=1e-10)
        if get_sols
            push!(sol_param_list, [sol, params, name])
        end
        # update θ0 and sv0
        θ_iter += param3b.oms*sol.t[end]
        sv_iter = sol.u[end]
    end
    return sv_iter
end

multishoot_trajectory = function (x::AbstractVector{T}, get_sols::Bool=false, verbose::Bool=false) where T
    # unpack decision vector
    x_LEO, x_mid, x_LPO, tofs, θs = unpack_x(x)
    
    # initialize storage
    sol_param_list = []
    
    # propagate from LEO forward
    sv0_LEO, θ0_leo, sol_ballistic_fwd = get_LEO_state(x_LEO, θs, verbose)
    svf_LEO = propagate_arc!(
        sv0_LEO, θ0_leo, [0, tofs[1]/n_arc], dt, x_LEO[5 : end],
        get_sols, sol_param_list, "leo_arc"
    )
    
    # propagate midpoint backward
    svm0 = x_mid[1:7]  # state-vector at midpoint, [r,v,mass]
    svf_mid_bck = propagate_arc!(
        svm0, θs[2], [0, -tofs[2]/n_arc], dt, x_mid[10 : end],
        get_sols, sol_param_list, "mid_bck_arc"
    )
    
    # propaagte midpoint forward
    svf_mid_fwd = propagate_arc!(
        svm0, θs[2], [0, tofs[3]/n_arc], dt, x_mid[10+3n_arc : end],
        get_sols, sol_param_list, "mid_fwd_arc"
    )
    
    # propagate from LPO backward
    sv0_LPO, θ0_lpo, sol_ballistic_bck = get_LPO_state(x_LPO, θs, verbose)
    svf_lpo = propagate_arc!(
        sol_ballistic_bck.u[end], θ0_lpo, [0, -tofs[4]/n_arc], dt, x_LPO[5 : end],
        get_sols, sol_param_list, "lpo_arc"
    )
    
    # residuals
    res = vcat(svf_mid_bck - svf_LEO, svf_lpo - svf_mid_fwd)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_fwd,sol_ballistic_bck], tofs
    end
end
























dt = 0.001

param3b = dyanmics_parameters()
lps = lagrange_points(param3b.mu2)

lp = 2
Az_km = 1200.0
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)#, saveat=LinRange(0, period, n+1))
monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
ys0 = R3BP.get_eigenvector(monodromy, true, 1);

# arrival LPO object
LPOArrival = SailorMoon.CR3BPLPO(
    res.x0, res.period, ys0, prob_cr3bp_stm, 1e-6, Tsit5(), 1e-12, 1e-12, 0.005
);

