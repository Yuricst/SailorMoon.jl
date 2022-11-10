using DifferentialEquations
using Plots
using LinearAlgebra
import ForwardDiff
import DiffResults
using AstrodynamicsBase
# import joptimise
using Printf
using DataFrames
using JSON
using CSV
plotly()


#### callback functions ###################
function terminate_condition(u,t,int)
    # Hit earth
    cond1 = sqrt((u[1] - param3b.mu2) ^2 + u[2]^2 + u[3]^2) - (6600 / param3b.lstar)
    # Hit moon
    cond2 = sqrt((u[1] - (1-param3b.mu2)) ^2 + u[2]^2 + u[3]^2) - (1770 / param3b.lstar)
    return - cond1 * cond2
end
# terminate condition: hit earth or moon
function terminate_affect!(int)
    terminate!(int)
end

# store the apoapsis value
function apoapsis_cond(u,t,int)
    # condition 1: SC is sufficiently far from the moon
    cond1 = sqrt((u[1] - (1-param3b.mu2))^2 + u[2]^2 + u[3]^2) - 5000 / param3b.lstar 
    
    # condition 2: dot product of velocity and position is zero
    cond2 = dot(u[1:3], u[4:6]) == 0
    
    return cond1 * cond2
end
function apoapsis_affect!(int)
    return NaN
end

# store the periapsis value and terminate
function periapsis_cond(u,t,int)
    # 0.95 * moon sma < SC distantce from B1 < 1.05 * moon sma
    cond1 = sqrt(u[1]^2 + u[2]^2 + u[3]^2) - (1 - param3b.mu2)*0.95 
    cond2 = - sqrt(u[1]^2 + u[2]^2 + u[3]^2) + (1 - param3b.mu2)*1.05
    return cond1 * cond2
end

function periapsis_affect!(int)
    terminate!(int)
end
###############################################



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
sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)#, saveat=LinRange(0, period, n+1))
monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
ys0 = R3BP.get_eigenvector(monodromy, true, 1) # monodromy eigenvector

## Grid search parameters: CHANGE HERE
ϕ_vec    = [0.0]    #range(0, stop=2*pi, length=20)
epsr_vec = [1e-6]   #10 .^(-12:-6)
epsv_vec = [1e-6]   #10 .^(-12:-6)
θ_vec    = [0.0]    #range(0, stop=2*pi, length=20)
tof_bck  = 30

## make initial conditions 
grids = []
for ϕ0 in ϕ_vec
    for ϵr in epsr_vec
        for ϵv in epsv_vec
            for θf in θ_vec
                
                βf = π - θf                
                # arrival LPO object
                LPOArrival = SailorMoon.CR3BPLPO2(
                    res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12
                );
                
                xf = vcat(SailorMoon.set_terminal_state2(ϕ0, θf, param3b, LPOArrival), 1.0)
                push!(grids, [ϕ0, ϵr, ϵv, θf, xf])
                
            end
        end
    end
end

# println("grids: $grids")

## initialize problem 

# include callback functions 
terminate_cb = ContinuousCallback(terminate_condition,terminate_affect!)
apoapsis_cb  = ContinuousCallback(apoapsis_cond, apoapsis_affect!; rootfind=true, save_positions=(true,false))
periapsis_cb = ContinuousCallback(apoapsis_cond, apoapsis_affect!)
cbs = CallbackSet(terminate_cb, apoapsis_cb, periapsis_cb)

svf_ = zeros(Float64, 1, 7)
params = [param3b.mu2, param3b.mus, 0.0, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0]
tspan = [0, -tof_bck]

prob = ODEProblem(R3BP.rhs_bcr4bp_thrust!, svf_, tspan, params)
# sol_bck = solve(prob_bck, Tsit5(), reltol=1e-12, abstol=1e-12);

## make ensemble problems
function prob_func(prob, i, repeat)
    print("\rproblem # $i")
    remake(prob, u0=grids[i][5], p=[param3b.mu2, param3b.mus, grids[i][4], param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0])
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=length(grids),
            callback=cbs, reltol=1e-12, abstol=1e-12,
            save_everystep=true);


## data extraction and make csv
# make dataframe
entries = ["phi0", "epsr", "epsv", "thetaf", "rp", "ra", "tof"]
df = DataFrame([ name =>[] for name in entries])


pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)

# extract the ensemble simulation
for (idx,sol) in enumerate(sim)
    
    ϕ0 = grids[idx][1]
    ϵr = grids[idx][2]
    ϵv = grids[idx][3]
    θf = grids[idx][4]

    plot!(pcart, Array(sol)[1,:], Array(sol)[2,:], color=:blue, linewidth=1.5, label="sol", linestyle=:dashdot)
    display(pcart)


    if sol.retcode == :Terminated
        if periapsis_cond(sol.u[end], 0.0, 0.0) == true
            # got to lunar orbit!
            rp = sqrt(sol.u[end][1]^2 + sol.u[end][2]^2 + sol[end][3]^2)
            ra = sqrt(sol.u[end-1][2]^2 + sol.u[end-1][2]^2 + sol[end-1][3]^2)
            tof = sol.t[end]
            push!(df, [ϕ0, ϵr, ϵv, θf, rp, ra, tof])
            println("idx $idx is a success!")
        end
    end
end

print(df)
CSV.write("grid_search.csv", df)


