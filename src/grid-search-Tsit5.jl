using Distributed
@everywhere begin
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
end

plotly()

# addprocs(10)
@show procs()

@everywhere  begin

    include("../../julia-r3bp/R3BP/src/R3BP.jl")
    include("../src/SailorMoon.jl")   # relative path to main file of module

    # integrator parameters
    alg = Tsit5()
    rtol = 1e-12
    atol = 1e-12

    function terminate_condition(u,t,int)
        # Hit earth
        cond1 = sqrt((u[1] + param3b.mu2) ^2 + u[2]^2 + u[3]^2) - (1000 / param3b.lstar)
        # Hit moon
        cond2 = sqrt((u[1] - (1-param3b.mu2)) ^2 + u[2]^2 + u[3]^2) - (870 / param3b.lstar)
        return - cond1 * cond2
    end
    
    # store the apoapsis value
    function apoapsis_cond(u,t,int)
        # condition 1: SC is sufficiently far from the moon
        r = sqrt((u[1] - (1-param3b.mu2))^2 + u[2]^2 + u[3]^2) # SC-Moon distance
        moon_soi = 5000 / param3b.lstar # define "sphere of influence"
        
        if sqrt(u[1]^2 + u[2]^2 + u[3]^2) > 1.5 
            # condition 2: dot product of velocity and position is zero
            return dot((u[1:3] + [param3b.mu2, 0.0, 0.0]), u[4:6])
        else
            return NaN
        end
    end
    
    # store the periapsis value and terminate
    function periapsis_cond(u,t,int)
        r = sqrt((u[1] + param3b.mu2)^2 + u[2]^2 + u[3]^2)  # SC-earth distance
        earth_leo_ub = (6357 + 1500) / param3b.lstar
        earth_leo_lb = 3000 / param3b.lstar
        
        if earth_leo_lb < r < earth_leo_ub
            return dot((u[1:3] + [param3b.mu2, 0.0, 0.0]), u[4:6])
        else 
            return NaN
        end
    end
    
    # generalized aps condition 
    function aps_cond(u,t,int)
        pos = u[1:3] + [param3b.mu2, 0.0, 0.0]
        return dot(pos, u[4:6]) 
    end
    
    # perilune
    function perilune_cond(u,t,int)
        r = sqrt((u[1] - (1 - param3b.mu2))^2 + u[2]^2 + u[3]^2)  # SC-Moon distance
        moon_soi = 66000 / param3b.lstar
        
        if r < moon_soi && -t > (10*86400 / param3b.tstar)
            return dot((u[1:3] - [1 - param3b.mu2, 0.0, 0.0]), u[4:6])
        else 
            return NaN
        end
    end
    
    
    # affect!
    terminate_affect!(int) = terminate!(int)
    no_affect!(int) = NaN

    function plot_circle(radius, x, y, n=50)
        circle = zeros(2,n)
        thetas = LinRange(0.0, 2π, n)
        for i = 1:n
            circle[1,i] = radius*cos(thetas[i]) + x
            circle[2,i] = radius*sin(thetas[i]) + y
        end
        return circle
    end
end
##########################################

@everywhere begin
    t0 = time()
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

    ## Grid search parameters: CHANGE HERE
    n = 100
    ϕ_vec    = LinRange(0, 2*pi, n+1)[1:n] #
    epsr_vec = 10 .^(-9:-6)
    epsv_vec = 10 .^(-9:-6)
    θ_vec    = LinRange(0, 2*pi, n+1)[1:n]  # [2.890265, 3.015929]     #
    tof_bck  = 150 * 86400 / param3b.tstar


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
end

# println("grids: $grids")

## initialize problem 

# include callback functions 
@everywhere begin
    terminate_cb = ContinuousCallback(terminate_condition, terminate_affect!; rootfind=false)
    apoapsis_cb  = ContinuousCallback(apoapsis_cond, no_affect!; rootfind=false, save_positions=(false,true))
    periapsis_cb = ContinuousCallback(periapsis_cond, terminate_affect!)
    # aps_cb       = ContinuousCallback(aps_cond, no_affect!; rootfind=false, save_positions=(false,true))
    perilune_cb  = ContinuousCallback(perilune_cond, no_affect!; rootfind=false, save_positions=(false,true))
    
    cbs = CallbackSet(terminate_cb, apoapsis_cb, periapsis_cb, perilune_cb)
    # print(cbs)


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
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=length(grids),
            callback=cbs, reltol=1e-12, abstol=1e-12,
            save_everystep=true);


## data extraction and make csv
# make dataframe
entries = ["id", "phi0", "epsr", "epsv", "thetaf", "ra", "rp", "dt1", "dt2", "x_ra", "x_rp", "tof", "lfb"]
df = DataFrame([ name =>[] for name in entries])

pcart = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.4)

id = 1
# extract the ensemble simulation
for (idx,sol) in enumerate(sim)
    global lfb_count = 0
    
    ϕ0 = grids[idx][1]
    ϵr = grids[idx][2]
    ϵv = grids[idx][3]
    θf = grids[idx][4]
    
    
    if sol.retcode == :Terminated
#         println(periapsis_cond(sol.u[end], 0.0, 0.0))
        if ~isnan(periapsis_cond(sol.u[end], 0.0, 0.0))  # == true
            t_vec = -sol.t[:]
            rp = sqrt(sol.u[end][1]^2 + sol.u[end][2]^2 + sol[end][3]^2)
            x_rp = sol.u[end]
            
            # find apoapsis
            r_vec = sqrt.((Array(sol)[1,:] .- [param3b.mu2]).^2+Array(sol)[2,:].^2+Array(sol)[3,:].^2)
            ra, id_ra = findmax(r_vec)
            x_ra = sol.u[id_ra]
            dt_ra = - sol.t[id_ra]
            dt_rp = -sol.t[end] - dt_ra 
            
            # flag: lunar flyby? 
            rm_vec = sqrt.((Array(sol)[1,:] .- [1-param3b.mu2]).^2+Array(sol)[2,:].^2+Array(sol)[3,:].^2)
            id_lfb = findall(rm_vec .< (66000 / param3b.lstar) .&& t_vec .> (10*86400/param3b.tstar))
            
            if ~isempty(id_lfb)
#                 println("there might be a perilune")
                for k in id_lfb
                    if ~isnan(perilune_cond(sol.u[k], sol.t[k], 0.0))
                        global lfb_count += 1
                    end
                end    
            end

            tof = -sol.t[end]
            push!(df, [id, ϕ0, ϵr, ϵv, θf, rp, ra, dt_ra, dt_rp, x_ra, x_rp, tof, lfb_count])
            println("idx $idx is a success!")
            plot!(pcart, Array(sol)[1,:], Array(sol)[2,:], color=:blue, linewidth=1.0, label="sol $idx", linestyle=:solid)
#             println(sol.u, sol.t)
            global id += 1
        end
    end
end

# print(df)
CSV.write("grid_search.csv", df)

moon = plot_circle(1738/param3b.lstar, 1-param3b.mu2, 0.0)
earth = plot_circle(6375/param3b.lstar, -param3b.mu2, 0.0)
moon_soi = plot_circle(66000/param3b.lstar, 1-param3b.mu2, 0.0)
leo_lb = plot_circle((3000)/param3b.lstar, -param3b.mu2, 0.0)
leo_ub = plot_circle((6375+1500)/param3b.lstar, -param3b.mu2, 0.0)
plot!(pcart, leo_lb[1,:], leo_lb[2,:], c=:black, lw=1.0, label="LEO lb")
plot!(pcart, leo_ub[1,:], leo_ub[2,:], c=:black, lw=1.0, label="LEO ub")
plot!(pcart, earth[1,:], earth[2,:], c=:green, lw=1.0, label="earth")
plot!(pcart, moon[1,:], moon[2,:], c=:orange, lw=1.0, label="moon")
plot!(pcart, moon_soi[1,:], moon_soi[2,:], c=:black, lw=1.0, label="moon soi")


display(pcart)

