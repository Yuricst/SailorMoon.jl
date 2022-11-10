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
_prob_base = ODEProblem(R3BP.rhs_bcr4bp_thrust!, rand(7), [0,1], params);

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
        @printf("ToF per arc  : %3.3f, %3.3f, %3.3f, %3.3f\n", tofs...)
        @printf("Phase angles : %3.3f, %3.3f, %3.3f\n", θs...)
    end
    return x_LEO, x_mid, x_LPO, tofs, θs
end

function propagate_arc!(sv0, θ0, tspan, x_control, get_sols::Bool, sol_param_list, name::String)
    sv_iter = [el for el in sv0]
    θ_iter = 1*θ0
    for i = 1:n_arc
        τ, γ, β = x_control[1+3*(i-1) : 3*i]
        params = [param3b.mu2, param3b.mus, θ_iter, param3b.as, param3b.oms, τ, γ, β, mdot, tmax]
        _prob = remake(_prob_base; tspan=tspan, u0 = sv_iter, p = params)
        sol = DifferentialEquations.solve(_prob, method, reltol=reltol, abstol=abstol)
        if get_sols
            push!(sol_param_list, [sol,params, name])
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
    
    # construct midpoint state
    svm = x_mid[1:7]  # state-vector at midpoint, [r,v,mass]
    
    # construct final state
    svf = vcat(SailorMoon.set_terminal_state(x_LPO[2], param3b, LPOArrival), x_LPO[3])

    # initialize storage
    sol_param_list = []

    # forward propagation
    sv_leo_mp = propagate_arc!(
        sv0, θs[1], [0, tofs[1]/n_arc], x_LEO[5 : end], 
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
        svf, θs[3], [0, -tofs[4]/n_arc], x_LPO[5 : end], 
        get_sols, sol_param_list, "lpo_arc"
    )

    # residuals
    res = vcat(sv_mid_bck_mp - sv_leo_mp, sv_lpo_mp - sv_mid_fwd_mp)[:]
    
    # output
    if get_sols == false
        return res
    else
        return res, sol_param_list
    end
end

sv_mid = [
    -2.592332271045137,
    -3.4055462335079487,
    0.0,
    -3.330371473282953,
    2.4064127167314444,
    0.0,
    1.0,
]  # mid point state-vector

# create test decision vector
ig_x_LEO = vcat([4.541281, 2.76442069, 1.0, 5.0], vcat([[0,0,0] for i = 1:n_arc]...));
ig_x_mid = vcat(sv_mid, 3.055608, 5.87818, vcat([[0,0,0] for i = 1:2n_arc]...));
ig_x_LPO = vcat([3.1523571, -0.0069668, 1.0, 10.0], vcat([[0,0,0] for i = 1:n_arc]...));
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

check = 0
for i = 1:length(lx)
    if lx[i] <= ig_x[i]  <= ux[i]
        check += 1
    else
        println("idx $i")
    end
end
check == length(lx)

ig_x[27]

n_arc

x_LEO, x_mid, x_LPO, tofs, θs = unpack_x(ig_x, true);

hmp, sol_param_list = multishoot_trajectory(ig_x, true);
hmp

arcs_color = Dict(
    "leo_arc" => :navy, 
    "mid_bck_arc" => :red1, 
    "mid_fwd_arc" => :darkorange, 
    "lpo_arc" => :gold
)

pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)
scatter!(pcart, lps[:,1], lps[:,2], marker=:diamond, color=:red, label="LPs")
# trajectory
for i = 1:length(sol_param_list)
    sol, _, name = sol_param_list[i]
    plot!(pcart, Array(sol)[1,:], Array(sol)[2,:], 
        linewidth=1.5, label="$name", c=arcs_color[name])
end
# control node
scatter!(pcart, [x_mid[1]], [x_mid[2]], marker=:circle, color=:black, label="MP")
pcart

# problem settings
ng = 14
lg = [0.0 for idx=1:ng];
ug = [0.0 for idx=1:ng];

fitness! = function (g, x)
    hmp = multishoot_trajectory(x, false)
    g[:] = hmp[:]
    f = x[3]   # minimize initial mass
    return f
end

gfoo  = zeros(ng)
fitness!(gfoo, ig_x)
gfoo

ip_options = Dict(
    "max_iter" => 5,   # approx 100
    "print_level" => 4,
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



hmp, sol_param_list = multishoot_trajectory(xopt, true);
hmp

arcs_color = Dict(
    "leo_arc" => :navy, 
    "mid_bck_arc" => :red1, 
    "mid_fwd_arc" => :darkorange, 
    "lpo_arc" => :gold
)

pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)
scatter!(pcart, lps[:,1], lps[:,2], marker=:diamond, color=:red, label="LPs")
# trajectory
for i = 1:length(sol_param_list)
    sol, _, name = sol_param_list[i]
    plot!(pcart, Array(sol)[1,:], Array(sol)[2,:], 
        linewidth=1.5, label="$name", c=arcs_color[name])
end
# control node
scatter!(pcart, [x_mid[1]], [x_mid[2]], marker=:circle, color=:black, label="MP")
pcart


