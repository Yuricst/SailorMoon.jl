"""
Given initial guess, earth θf and ϕ0 that matches the initial condition at the LEO.
x = [ϕ0, θf]
"""

#push!(LOAD_PATH,"../src/")
using joptimise

include("../../julia-r3bp/R3BP/src/R3BP.jl")
include("../src/SailorMoon.jl")   # relative path to main file of module

param3b = SailorMoon.dyanmics_parameters()
lps = SailorMoon.lagrange_points(param3b.mu2)

## set up of initial condition (Lyapunov orbit)
lp = 2
Az_km = 1200.0
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
# for Halo propagation, keep the tol as tight as possible 
sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12) #, saveat=LinRange(0, period, n+1))
monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
ys0 = R3BP.get_eigenvector(monodromy, true, 1) # monodromy eigenvector

ϵr = 1e-6
ϵv = 1e-6

LPOArrival = SailorMoon.CR3BPLPO2(
    res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12
);

function objective(x)
    # using x, generate the trajectory 
    ϕ0, θf = x[1], x[2]

    xf = vcat(SailorMoon.set_terminal_state2(ϕ0, θf, param3b, LPOArrival), 1.0)

    




    rp_target = 200 + 6375  # altitude of LEO 
    rp = 
    f = rp - rp_target 
    return f
end


function obj_con(g, x)
    # compute objective
    f = objective(x)
    g[1] = 0.0
    return f
end

# initial guess
x0 = [4.0; 4.0]
# bounds on variables
lx = [-5.0; -5.0]
ux = [5.0; 5.0]

# bounds on constriants
lg = [0.0]
ug = [0.0]
# number of constraints
ng = 1


## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "tol" => 1e-6
)

xopt, fopt, info = minimize(obj_con, x0, ng; lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt", options=ip_options, derivatives=ForwardAD());

println("Done with IPOPT!")
println(info)
println(xopt)