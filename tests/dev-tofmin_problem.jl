using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
import joptimise
using Printf
using JSON

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
method = RK4()  # CalvoSanz4()
reltol = 1e-10
abstol = 1e-10
#dt = 0.005

param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

lp = 2
Az_km = 1200.0
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
sol = solve(prob_cr3bp_stm, method, reltol=reltol, abstol=abstol)#, saveat=LinRange(0, period, n+1))
monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
ys0 = R3BP.get_eigenvector(monodromy, true, 1);

# arrival LPO object
LPOArrival = SailorMoon.CR3BPLPO(
    res.x0, res.period, ys0, prob_cr3bp_stm, 1e-6, Tsit5(), 1e-12, 1e-12
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

params = [param3b.mu2, param3b.mus, 0.0, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0]
_prob_base = ODEProblem(R3BP.rhs_bcr4bp_thrust!, rand(7), [0,1], params)

propagate_trajectory = function (x::AbstractVector{T}, get_sols::Bool=false) where T
    # unpack
    nx = length(x)
    θf, tof, eta, r_apogee, raan, ϕ, m0, mf = x[1:8]  # θf: Sun angle at final time
    tau1     = x[9 : floor(Int, (nx-8)/2 + 8)]  # discretization numbers are the same for the first & second arc
    tau2     = x[floor(Int, (nx-8)/2 + 9) : end]
    tof_fwd = tof * eta
    tof_bck = tof * (1 - eta)

    # construct initial state
    sma = (rp_parking + r_apogee)/2
    ecc = (r_apogee - rp_parking)/(r_apogee + rp_parking)
    sv0_kep = [sma, ecc, 0.0, raan, 0.0, 0.0]
    θ0 = θf - param3b.oms*(tof_fwd + tof_bck)   # initial Sun angle
    sv0_i = AstrodynamicsBase.kep2cart(sv0_kep, param3b.mu1)
    sv0 = vcat(inertial2rotating(sv0_i, θ0, 1.0) + [-param3b.mu2,0,0,0,0,0], m0)

    # construct final state
    #x0_stm = vcat(LPOArrival.x0, reshape(I(6), (36,)))
    #tspan = [0, ϕ*LPOArrival.period]
    #prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, tspan, [param3b.mu2,])
    #sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)
    svf = vcat(SailorMoon.set_terminal_state(ϕ, param3b, LPOArrival), mf)

    # initialize storage
    sols_fwd, sols_bck = [], []
    params_fwd, params_bck = [], []

    # forward propagation
    nhalf = Int(n/2)
    tspan_fwd = [0, tof_fwd/nhalf]
    for i = 1:nhalf
        τ, γ, β = tau1[3*i-2 : 3*i]
        params = [param3b.mu2, param3b.mus, θ0, param3b.as, param3b.oms, τ, γ, β, mdot, tmax]
        _prob = remake(
            _prob_base; tspan=tspan_fwd,
            u0 = sv0,
            p = params,
        )
        sol = DifferentialEquations.solve(_prob, method, reltol=reltol, abstol=abstol)
        if get_sols
            push!(sols_fwd, sol)
            push!(params_fwd, params)
        end
        # update θ0 and sv0
        θ0 += param3b.oms*sol.t[end]
        sv0 = sol.u[end]
    end

    # back propagation
    tspan_bck = [0, -tof_bck/(n-nhalf)]
    for i = 1:n-nhalf
        τ, γ, β = tau2[3*i-2 : 3*i]
        params = [param3b.mu2, param3b.mus, θf, param3b.as, param3b.oms, τ, γ, β, mdot, tmax]
        _prob = remake(
            _prob_base; tspan=tspan_bck,
            u0 = svf,
            p = params,
        )
        sol = DifferentialEquations.solve(_prob, method, reltol=reltol, abstol=abstol)
        if get_sols
            push!(sols_bck, sol)
            push!(params_bck, params)
        end
        # update θf (note sol.t[end] is negative so += is correct) and svf
        θf += param3b.oms*sol.t[end]
        svf = sol.u[end]
    end

    # residual
    if get_sols == false
        return svf - sv0
    else
        return sols_fwd, sols_bck, params_fwd, params_bck
    end
end

function xprint(x)
    θf, tof, eta, r_apogee, raan, ϕ, m0, mf = x[1:8]
    @printf("Launch RA   : %1.4f\n", r_apogee)
    @printf("Launch RAAN : %3.4f\n", rad2deg(raan))
    @printf("TOF [day]   : %3.4f\n", tof*param3b.tstar/86400)
    @printf("TOF [TU]    : %3.4f\n", tof)
    @printf("m0          : %2.4f\n", m0)
    @printf("mf          : %2.4f\n", mf)
end

function get_controls(x)
    nx = length(x)
    tau1     = x[9 : floor(Int, (nx-8)/2 + 8)]  # discretization numbers are the same for the first & second arc
    tau2     = x[floor(Int, (nx-8)/2 + 9) : end]
    tau1_list = [tau1[3*i-2 : 3*i] for i = 1:Int(n/2)]
    tau2_list = [tau2[3*i-2 : 3*i] for i = 1:n-Int(n/2)]
    tau1 = hcat(tau1_list...)
    tau2 = hcat(tau2_list...)
    return tau1, tau2
end


## Solve problem
n = 20

bounds_tau = [0,1]
bounds_γ   = [-π, π]
bounds_β   = [-π, π]
# θf, tof, eta, sma, raan, ϕ, m0, mf
xtest = [
    3.171156742785936,
    24.678179489698685,
    0.4188256425390488,
    4.966847320088643,
    10.446876808409478,
    -0.0124160745030169,
    5.031998701118364,
    1.0,
]
for i = 1:n
    global xtest = vcat(xtest, [0,0,0])
end

# θf, tof, eta, sma, raan, ϕ, m0, mf
lx = [
    2.6, 18, 0.4, 4.0, 0.0, -1.0, 1.0, 1.0
]
ux = [
    3.2, 27, 0.42, 5.4, 2π, 1.0, 10.0, 1.0
]
for i = 1:n
    global lx = vcat(lx, [bounds_tau[1],bounds_γ[1],bounds_β[1]])
    global ux = vcat(ux, [bounds_tau[2],bounds_γ[2],bounds_β[2]])
end


sols_fwd, sols_bck, params_fwd, params_bck = propagate_trajectory(xtest, true);
propagate_trajectory(xtest, false)

cs_fwd = palette([:darkred, :navyblue], max(Int(n/2),2))
cs_bck = palette([:darkgreen, :goldenrod], max(n-Int(n/2),2))

pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)
scatter!(pcart, lps[:,1], lps[:,2], marker=:diamond, color=:red, label="LPs")
# trajectory
for (ifwd,sol_fwd) in enumerate(sols_fwd)
    plot!(pcart, Array(sol_fwd)[1,:], Array(sol_fwd)[2,:], color=cs_fwd[ifwd],
        linewidth=1.5, label="fwd $ifwd", linestyle=:dashdot)
end
for (ibck,sol_bck) in enumerate(sols_bck)
    plot!(pcart, Array(sol_bck)[1,:], Array(sol_bck)[2,:], color=cs_bck[ibck],
        linewidth=1.5, label="bck $ibck", linestyle=:solid)
end
plot!(pcart; title="Initial guess")

# problem settings
ng = 7
lg = [0.0 for idx=1:ng];
ug = [0.0 for idx=1:ng];

fitness! = function (g, x)
    # evaluate objective & objective gradient (trivial)
    #f = 1       # whichever x corresponds to e.g. mass at LEO
    #g[:] = propagate_trajectory(x, false)

    sols_fwd, sols_bck = propagate_trajectory(x, true)
    sol_fwd = sols_fwd[end]
    sol_bck = sols_bck[end]

    # convert to spherical coordinates
    svf_spherical = vcat(cart2spherical(sol_fwd.u[end][1:6]), sol_fwd.u[end][7])
    svb_spherical = vcat(cart2spherical(sol_bck.u[end][1:6]), sol_bck.u[end][7])
    g[:] = svb_spherical - svf_spherical #sol_bck.u[end] - sol_fwd.u[end]

    # minimize initial mass + tof
    f = sols_fwd[1].u[1][7] + 0.05 * x[2]
    return f
end

ip_options = Dict(
    "max_iter" => 200,   # approx 100
    "print_level" => 5,
    "acceptable_tol" => 1e-5,
    "constr_viol_tol" => 1e-5,
)
x0 = [el for el in xtest]
xopt, fopt, info = joptimise.minimize(
    fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug,
    #derivatives=joptimise.ForwardAD(),
    options=ip_options,
);


sols_fwd, sols_bck, params_fwd, params_bck = propagate_trajectory(xopt, true);
propagate_trajectory(xopt, false)

dr = sols_bck[end].u[end][1:3] - sols_fwd[end].u[end][1:3]
dv = sols_bck[end].u[end][4:6] - sols_fwd[end].u[end][4:6]
gfoo  = zeros(ng)
ffoo = fitness!(gfoo, xopt)
println("ffoo: $ffoo")
println("gfoo: $gfoo")
@printf("Position offset = %1.4f [DU]\n", norm(dr))
@printf("Velocity offset = %1.4f [DU/TU]\n", norm(dv))

xprint(xopt)


## Plot
cs_fwd = palette([:darkred, :navyblue], max(Int(n/2),2))
cs_bck = palette([:darkgreen, :goldenrod], max(n-Int(n/2),2))
pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)
scatter!(pcart, lps[:,1], lps[:,2], marker=:diamond, color=:red, label="LPs")
# trajectory
for (ifwd,sol_fwd) in enumerate(sols_fwd)
    plot!(pcart, Array(sol_fwd)[1,:], Array(sol_fwd)[2,:], color=cs_fwd[ifwd],
        linewidth=1.5, label="fwd $ifwd", linestyle=:dashdot)
end
for (ibck,sol_bck) in enumerate(sols_bck)
    plot!(pcart, Array(sol_bck)[1,:], Array(sol_bck)[2,:], color=cs_bck[ibck],
        linewidth=1.5, label="bck $ibck", linestyle=:solid)
end
plot!(pcart; title="Solved")
display(pcart)

# save solution to file
solution_dict = Dict(
    "n" => n,
    "tmax" => tmax,
    "mdot" => mdot,
    "lpo_x0" => res.x0,
    "lpo_period" => res.period,
    "xopt" => xopt,
)
save_filename = joinpath("..", "results", "test_save.json")
open(save_filename, "w") do f
    JSON.print(f, solution_dict, 4)
end