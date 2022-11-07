using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
import joptimise
using Printf

plotly()

include("../../julia-r3bp/R3BP/src/R3BP.jl")
include("../src/SailorMoon.jl")   # relative path to main file of module

param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

lp = 2
Az_km = 1200.0
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")
res.flag

x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)#, saveat=LinRange(0, period, n+1))
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
    θf, tof, eta, sma, ecc, raan, ϕ = x[1:7]  # θf: Sun angle at final time
    tau1     = x[8 : floor(Int, (nx-7)/2 + 7)]  # discretization numbers are the same for the first & second arc
    tau2     = x[floor(Int, (nx-7)/2 + 8) : end]
    tof_fwd = tof * eta
    tof_bck = tof * (1 - eta)
    
    # construct initial state
    sv0_kep = [sma, ecc, 0.0, raan, 0.0, 0.0]
    θ0 = θf - param3b.oms*(tof_fwd + tof_bck)   # initial Sun angle
    sv0_i = AstrodynamicsBase.kep2cart(sv0_kep, param3b.mu1)
    sv0 = vcat(inertial2rotating(sv0_i, θ0, 1.0) + [-param3b.mu2,0,0,0,0,0], 1.0)
    
    # construct final state
    #x0_stm = vcat(LPOArrival.x0, reshape(I(6), (36,)))
    #tspan = [0, ϕ*LPOArrival.period]
    #prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, tspan, [param3b.mu2,])
    #sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)
    svf = vcat(SailorMoon.set_terminal_state(ϕ, param3b, LPOArrival), 1.0)
    
    # initialize storage
    sols_fwd, sols_bck = [], []
    
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
        sol = DifferentialEquations.solve(_prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        if get_sols
            push!(sols_fwd, sol)
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
        sol = DifferentialEquations.solve(_prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        if get_sols
            push!(sols_bck, sol)
        end
        # update θf (note sol.t[end] is negative so += is correct) and svf
        θf += param3b.oms*sol.t[end]
        svf = sol.u[end]
    end
    
    # residual
    if get_sols == false
        return svf - sv0
    else
        return sols_fwd, sols_bck
    end
end

n = 20

bounds_tau = [0,1]
bounds_γ   = [-π, π]
bounds_β   = [-π, π]
# θf, tof, eta, sma, ecc, raan = x
xtest = [
    2.4, 22, 0.4, 1.8, 0.87, 4.5, 0.02
]
for i = 1:n
    global xtest = vcat(xtest, [0.1,0,0])
end

lx = [
    2.6, 18, 0.3, 1.9, 0.7, 0.0, -1.0
]
ux = [
    3.2, 27, 0.7, 2.4, 0.995, 2π, 1.0
]
for i = 1:n
    global lx = vcat(lx, [0,bounds_γ[1],bounds_β[1]])
    global ux = vcat(ux, [0,bounds_γ[1],bounds_β[2]])
end

sols_fwd, sols_bck = propagate_trajectory(xtest, true);
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
pcart

# problem settings
ng = 3
lg = [0.0 for idx=1:ng];
ug = [0.0 for idx=1:ng];

fitness! = function (g, x)
    # evaluate objective & objective gradient (trivial)
    #f = 1       # whichever x corresponds to e.g. mass at LEO
    #g[:] = propagate_trajectory(x, false)
    
    sols_fwd, sols_bck = propagate_trajectory(x, true)
    sol_fwd = sols_fwd[end]
    sol_bck = sols_bck[end]
    g[:] = sol_bck.u[end][1:3] - sol_fwd.u[end][1:3]
    f = norm(sol_bck.u[end][4:6] - sol_fwd.u[end][4:6])
    return f
end

gfoo  = zeros(ng)
fitness!(gfoo, xtest)
gfoo

ip_options = Dict(
    "max_iter" => 100,   # approx 100
    "print_level" => 5,
    "acceptable_tol" => 1e-6,
    "constr_viol_tol" => 1e-6,
)
x0 = [el for el in xtest]
xopt, fopt, info = joptimise.minimize(
    fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug,
    #derivatives=joptimise.ForwardAD(),
    options=ip_options,
);

info




