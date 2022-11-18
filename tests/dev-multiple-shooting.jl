using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
import joptimise
using Printf
using JSON
using BSON

plotly()

include("../../julia-r3bp/R3BP/src/R3BP.jl")
include("../src/SailorMoon.jl")   # relative path to main file of module

function cart2spherical(sv_cartesian::Array{<:Real,1})
    # unpack state-vector
    x,y,z,vx,vy,vz = sv_cartesian
    r = sqrt(x^2 + y^2 + z^2)
    sv_spherical = [
        r,
        atan(y,x),
        asin(z/r),
        (x*vx+y*vy+z*vz)/r,
        (vx*y - x*vy)/(x^2+y^2),
        (z*(x*vx+y*vy) - (x^2+y^2)*vz)/((x^2+y^2+z^2)*sqrt(x^2+y^2))
    ]
    return sv_spherical
end

# solver settings within fitness function
# https://diffeq.sciml.ai/stable/solvers/dynamical_solve/#Symplectic-Integrators
#method = RK4()  # CalvoSanz4()
#reltol = 1e-10
#abstol = 1e-10
dt = 0.001

param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

# lp = 2
# Az_km = 1200.0
# println("Halo guess Az_km: $Az_km")
# northsouth = 3   # 1 or 3
# guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
# res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")
#
# x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
#prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
# sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)#, saveat=LinRange(0, period, n+1))
# monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
# ys0 = R3BP.get_eigenvector(monodromy, true, 1);
#
# # arrival LPO object
# LPOArrival = SailorMoon.CR3BPLPO(
#     res.x0, res.period, ys0, prob_cr3bp_stm, 1e-6, Tsit5(), 1e-12, 1e-12, 0.005
# );

# load LPO
_load_dict = BSON.load("lpo.bson")
x0_stm = vcat(_load_dict[:x0], reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
LPOArrival = SailorMoon.CR3BPLPO(
    [el for el in _load_dict[:x0]],
    _load_dict[:period],
    [el for el in _load_dict[:ys0]],
    prob_cr3bp_stm, 1e-6, Tsit5(), 1e-12, 1e-12, 0.005
);

tmax_si = 0.3  # N
isp_si = 3500  # sec
mdot_si = tmax_si / (isp_si * 9.81)
mstar = 2000  # kg
rp_parking = (6378+200)/param3b.lstar   # parking orbit radius

tmax = AstrodynamicsBase.dimensional2canonical_thrust(
    tmax_si, mstar, param3b.lstar, param3b.tstar
)
mdot = AstrodynamicsBase.dimensional2canonical_mdot(
    mdot_si, mstar, param3b.tstar
)

params = [
    param3b.mu2, param3b.mus, 0.0, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0,
    SailorMoon.dv_sun_dir_angles
]
#_prob_base = ODEProblem(R3BP.rhs_bcr4bp_thrust!, rand(7), [0,1], params);
_prob_base = ODEProblem(SailorMoon.rhs_bcr4bp_emframe_thrust!, rand(7), [0,1], params);


n_arc = 5

function unpack_x(x::AbstractVector{T}, verbose::Bool=false) where T
    # unpack
    nx = length(x)
    x_LEO = x[1:4+3n_arc]
    x_mid = x[5+3n_arc:13+9n_arc]    # x[5+3n_arc:4+3n_arc+9+6n_arc]
    x_LPO = x[14+9n_arc:17+12n_arc]  # x[14+9n_arc:13+9n_arc+4+3n_arc]

    # get time of flights
    tofs = [x_LEO[4], x_mid[8], x_mid[9], x_LPO[4]]
    θf = x_LPO[1]
    θs = [
        θf - param3b.oms*sum(broadcast(abs, tofs)),
        θf - param3b.oms*sum(broadcast(abs, tofs[3:4])),
        θf
    ]
    # print message
    if verbose
        @printf("Initial mass : %2.4f\n", x_LEO[3])
        @printf("Final mass   : %2.4f\n", x_LPO[3])
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
    end
    return x_LEO, x_mid, x_LPO, tofs, θs
end

dt = 0.01

# ballistic time right after launch
ballistic_time = 1*86400 / param3b.tstar
ballistic_time_back = 1*86400 / param3b.tstar

function propagate_arc!(sv0, θ0, tspan, x_control, get_sols::Bool, sol_param_list, name::String)
    sv_iter = [el for el in sv0]
    θ_iter = 1*θ0
    for i = 1:n_arc
        τ, γ, β = x_control[1+3*(i-1) : 3*i]
        params = [
            param3b.mu2, param3b.mus, θ_iter, param3b.as, param3b.oms, τ, γ, β, mdot, tmax,
            SailorMoon.dv_sun_dir_angles
        ]
        _prob = remake(_prob_base; tspan=tspan, u0 = sv_iter, p = params)
        sol = SailorMoon.integrate_rk4(_prob, dt);
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

multishoot_trajectory = function (x::AbstractVector{T}, get_sols::Bool=false) where T
    # unpack
    x_LEO, x_mid, x_LPO, tofs, θs = unpack_x(x)

    # construct initial state
    sma = (rp_parking + x_LEO[1])/2
    ecc = (x_LEO[1] - rp_parking)/(x_LEO[1] + rp_parking)
    sv0_kep = [sma, ecc, 0.0, x_LEO[2], 0.0, 0.0]
    sv0_i = AstrodynamicsBase.kep2cart(sv0_kep, param3b.mu1)
    sv0 = vcat(inertial2rotating(sv0_i, θs[1], 1.0) + [-param3b.mu2,0,0,0,0,0], x_LEO[3])
    # ballistic propagation with small time-steps
    params = [
        param3b.mu2, param3b.mus, θs[1]-param3b.oms*ballistic_time, param3b.as, param3b.oms,
        0.0, 0.0, 0.0, mdot, tmax, SailorMoon.dv_sun_dir_angles
    ]
    _prob = remake(_prob_base; tspan=[0,ballistic_time], u0 = sv0, p = params)
    sol_ballistic_fwd = SailorMoon.integrate_rk4(_prob, 0.001);

    # construct midpoint state
    svm = x_mid[1:7]  # state-vector at midpoint, [r,v,mass]

    # construct final state
    svf = vcat(
        SailorMoon.set_terminal_state(x_LPO[2], param3b, LPOArrival, true),
        x_LPO[3]
    )
    # ballistic propagation with small time-steps
    params = [
        param3b.mu2, param3b.mus, θs[3]+param3b.oms*ballistic_time_back, param3b.as, param3b.oms,
        0.0, 0.0, 0.0, mdot, tmax, SailorMoon.dv_sun_dir_angles
    ]
    _prob = remake(_prob_base; tspan=[0,-ballistic_time_back], u0 = svf, p = params)
    sol_ballistic_bck = SailorMoon.integrate_rk4(_prob, 0.001);

    # initialize storage
    sol_param_list = []

    # forward propagation
    sv_leo_mp = propagate_arc!(
        sol_ballistic_fwd.u[end], θs[1], [0, tofs[1]/n_arc], x_LEO[5 : end],
        get_sols, sol_param_list, "leo_arc"
    )

    # middle point propagation backward
    sv_mid_bck_mp = propagate_arc!(
        svm, θs[2], [0, -tofs[2]/n_arc], x_mid[10 : end],
        get_sols, sol_param_list, "mid_bck_arc"
    )

    # middle point propagation forward
    sv_mid_fwd_mp = propagate_arc!(
        svm, θs[2], [0, tofs[3]/n_arc], x_mid[10+3n_arc : end],
        get_sols, sol_param_list, "mid_fwd_arc"
    )

    # back propagation
    sv_lpo_mp = propagate_arc!(
        sol_ballistic_bck.u[end], θs[3], [0, -tofs[4]/n_arc], x_LPO[5 : end],
        get_sols, sol_param_list, "lpo_arc"
    )

    # residuals
    res = vcat(sv_mid_bck_mp - sv_leo_mp, sv_lpo_mp - sv_mid_fwd_mp)[:]

    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list, [sol_ballistic_fwd,sol_ballistic_bck], tofs
    end
end


sv_mid = [
    -2.58,
    -3.39,
    0.0,
    -3.330371473282953,
    2.4064127167314444,
    0.0,
    1.0,
]  # mid point state-vector

# create test decision vector
τ_ig = 0.0
ig_x_LEO = vcat([4.541281, 2.76442069, 1.0, 5.0], vcat([[τ_ig,0,0] for i = 1:n_arc]...));
ig_x_mid = vcat(sv_mid, 3.0, 5.6, vcat([[τ_ig,0,0] for i = 1:2n_arc]...));
ig_x_LPO = vcat([3.1523571, -0.0069668, 1.0, 10.0], vcat([[τ_ig,0,0] for i = 1:n_arc]...));
ig_x = vcat(ig_x_LEO, ig_x_mid, ig_x_LPO);

lx_LEO = vcat(
    [4.2, -4π, 1.0, 4.0], vcat([[0,-π,-π] for i = 1:n_arc]...)
);
ux_LEO = vcat(
    [4.8, 4π, 3.0, 8.0], vcat([[0, π, π] for i = 1:n_arc]...)
);

lx_mid = vcat(
    [-5, -5, 0, -4, -4, 0, 1.0, 3.0, 3.0], vcat([[0,-π,-π] for i = 1:2n_arc]...)
);
ux_mid = vcat(
    [ 5,  5, 0,  4,  4, 0, 3.0, 10.0, 10.0], vcat([[0, π, π] for i = 1:2n_arc]...)
);

lx_LPO = vcat(
    [-2π, -π, 1.0, 4.0], vcat([[0,-π,-π] for i = 1:n_arc]...)
);
ux_LPO = vcat(
    [ 4π, π, 1.0, 12.0], vcat([[0, π, π] for i = 1:n_arc]...)
);

lx = vcat(lx_LEO, lx_mid, lx_LPO);
ux = vcat(ux_LEO, ux_mid, ux_LPO);


x_LEO, x_mid, x_LPO, tofs, θs = unpack_x(ig_x, true);

@time hmp_ig, sol_param_list_ig, sols_ballistic_ig, _ = multishoot_trajectory(ig_x, true);

arcs_color = Dict(
    "leo_arc" => :navy,
    "mid_bck_arc" => :red1,
    "mid_fwd_arc" => :darkorange,
    "lpo_arc" => :gold
)

function plot_trajectory!(pl, all_grey::Bool, sols_ballistic, sol_param_list)
    # ballistic leg
    if all_grey
        c = :grey
    else
        c = :dodgerblue
    end
    for sol_ballistic in sols_ballistic
        plot!(pl, hcat(sol_ballistic.u...)[1,:], hcat(sol_ballistic.u...)[2,:], c=c)
    end
    # trajectory
    for i = 1:length(sol_param_list)
        if all_grey
            c = :grey
        else
            c = arcs_color[name]
        end
        sol, _, name = sol_param_list[i]
        plot!(pl, SailorMoon.Array(sol)[1,:], SailorMoon.Array(sol)[2,:],
            linewidth=1.5, label="$name", c=c)
    end
    return
end


# problem settings
ng = 14
lg = [0.0 for idx=1:ng];
ug = [0.0 for idx=1:ng];

fitness! = function (g, x)
    hmp, _, _, tofs = multishoot_trajectory(x, true)
    g[:] = hmp[:]
    f = x[3]   # minimize initial mass
    return sum(tofs)
end

gfoo  = zeros(ng)
fitness!(gfoo, ig_x)
gfoo

n_arc

ip_options = Dict(
    "max_iter" => 5,   # approx 100
    "print_level" => 5,
    "acceptable_tol" => 1e-6,
    "constr_viol_tol" => 1e-6,
)
x0 = [el for el in ig_x]
xopt, fopt, info = joptimise.minimize(
    fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug,
    #derivatives=joptimise.ForwardAD(),
    options=ip_options,
);


hmp, sol_param_list, sols_ballistic, _ = multishoot_trajectory(xopt, true);
x_LEO_opt, x_mid_opt, _, _, _ = unpack_x(xopt, true);

arcs_color = Dict(
    "leo_arc" => :navy,
    "mid_bck_arc" => :red1,
    "mid_fwd_arc" => :darkorange,
    "lpo_arc" => :gold,
)

pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)
scatter!(pcart, lps[:,1], lps[:,2], marker=:diamond, color=:red, label="LPs")

# initial guess
plot_trajectory!(pcart, true, sols_ballistic_ig, sol_param_list_ig)

# ballistic legs
for sol_ballistic in sols_ballistic
    plot!(pcart, hcat(sol_ballistic.u...)[1,:], hcat(sol_ballistic.u...)[2,:], c=:dodgerblue, label="Ballistic")
end
# trajectory
for i = 1:length(sol_param_list)
    sol, _, name = sol_param_list[i]
    plot!(pcart, SailorMoon.Array(sol)[1,:], SailorMoon.Array(sol)[2,:],
        linewidth=1.5, label="$name", c=arcs_color[name])
end
# control node
scatter!(pcart, [x_mid[1]], [x_mid[2]], marker=:circle, color=:black, label="MP-init")
scatter!(pcart, [x_mid_opt[1]], [x_mid_opt[2]], marker=:diamond, color=:deeppink, label="MP-opt")
pcart
