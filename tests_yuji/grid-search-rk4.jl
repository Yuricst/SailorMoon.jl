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
    plotly()


    include("../../julia-r3bp/R3BP/src/R3BP.jl")
    include("../src/SailorMoon.jl")   # relative path to main file of module

    #### CALLBACK FUNCTIONS #############################################
    # store the apoapsis value
    function apoapsis_cond(u,t,int)
        θm0 = prob_base.p[4]
        θm = θm0 + param3b.oml * t
        
        # for convenience, instead of taking the distance of SC-moon, assume r_a > 2.0
        if sqrt((u[1]-param3b.as)^2 + u[2]^2 + u[3]^2) > 2.0
            # condition 2: dot product of velocity and position is zero
            return dot([u[1] - param3b.as - (-param3b.mu2 * cos(θm)), u[2] - (-param3b.mu2 * sin(θm)), u[3]], u[4:6])
        else
            return NaN
        end
    end

    # store the periapsis value and terminate
    function periapsis_cond(u,t,int)
        θm0 = prob_base.p[4]
        θm = θm0 + param3b.oml * t
        r = sqrt((u[1] - param3b.as - (-param3b.mu2 * cos(θm)))^2 + (u[2] - (-param3b.mu2 * sin(θm))) ^2 + u[3]^2)  # SC-earth distance
    #     r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        earth_leo_ub = 3 * param3b.mu2   # 10000 / param3b.lstar  # km
        earth_leo_lb = 200  / param3b.lstar  # km
        
        if earth_leo_lb < r < earth_leo_ub
            return dot([u[1] - param3b.as - (-param3b.mu2 * cos(θm)), u[2] - (-param3b.mu2 * sin(θm)), u[3]], u[4:6])
        else 
            return NaN
        end
    end

    function perilune_cond(u,t,int)
        θm0 = prob_base.p[4]
        θm  = θm0 + param3b.oml * t
        
        # moon-SC distance
        r = sqrt((u[1] - param3b.as - (1-param3b.mu2)*cos(θm))^2 + (u[2] - (1-param3b.mu2)*sin(θm)) ^2 + u[3]^2)
        moon_soi = 66100 / param3b.lstar
        v_moon   = (1-param3b.mu2)*param3b.oml * [-sin(θm), cos(θm), 0] 
        
        if r < moon_soi #&& -t > (10*86400 / param3b.tstar)
            return dot([u[1] - param3b.as - (1-param3b.mu2) * cos(θm), u[2] - (1-param3b.mu2) * sin(θm), u[3]], u[4:6]-v_moon)
        else 
            return NaN
        end
        
    end
        
    terminate_affect!() = true
    no_affect!() = false


    function rhs_bcr4bp_thrust!(du,u,p,t)
        # unpack arguments
        μ, μ_3, t0, a_s, ω_s, τ, θ, β, mdot, tmax = p
        # decompose state
        x, y, z, vx, vy, vz, mass = u

        # calculate radii
        r1 = sqrt( (x+μ)^2 + y^2 + z^2 )
        r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 )

        # sun position
        xs = a_s*cos(ω_s*t + t0)
        ys = a_s*sin(ω_s*t + t0)
        zs = 0.0
        r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

        # compute delta-V vector
    #     dir_v = dv_lvlh2inertial(u[1:6], [τ, θ, β])
        dir_v = [0.0, 0.0, 0.0]
        ts = tmax*dir_v

        # position-state derivative
        du[1] = vx
        du[2] = vy
        du[3] = vz

        # velocity derivatives: Earth, moon, sun
        du[4] =  2*vy + x - ((1-μ)/r1^3)*(μ+x)  + (μ/r2^3)*(1-μ-x)  + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs ) + ts[1]/mass
        du[5] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys ) + ts[2]/mass
        du[6] =           - ((1-μ)/r1^3)*z      - (μ/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs ) + ts[3]/mass
        
    #     println(r1, " ", r2, " ", r3)
    #     println(((1-μ)/r1^3)*(μ+x), " ", (μ/r2^3)*(1-μ-x), " ",  -(μ_3/r3^3)*(x-xs) )

        # mass derivative
        du[7] = -mdot*τ
    end

    function plot_circle(radius, x, y, n=50)
        circle = zeros(2,n)
        thetas = LinRange(0.0, 2π, n)
        for i = 1:n
            circle[1,i] = radius*cos(thetas[i]) + x
            circle[2,i] = radius*sin(thetas[i]) + y
        end
        return circle
    end


    param3b = SailorMoon.dynamics_parameters()
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
    n = 50
    θs_vec   = [3.76991118430775] #LinRange(0, 2*pi, n+1)[1:n]  #    #
    ϕ_vec    =  [0.628318530717958] #LinRange(0, 2*pi, n+1)[1:n] #
    epsr_vec = [1e-6]   #10 .^(-12:-6)
    epsv_vec = [1e-6]   #10 .^(-12:-6)
    tof_bck  = 120 * 86400 / param3b.tstar

    # include callback functions 
    apoapsis_cb  = ContinuousCallback(apoapsis_cond, no_affect!)
    periapsis_cb = ContinuousCallback(periapsis_cond, terminate_affect!)
    perilune_cb  = ContinuousCallback(perilune_cond, no_affect!)
    cbs = [apoapsis_cb, periapsis_cb, perilune_cb];



    ## make initial conditions 
    grids = []
    for ϕ0 in ϕ_vec
        for ϵr in epsr_vec
            for ϵv in epsv_vec
                for θsf in θs_vec
                    
                    θmf = π - θsf                
                    # arrival LPO object
                    LPOArrival = SailorMoon.CR3BPLPO2(
                        res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12, 0.005
                    );
                    
                    xf = SailorMoon.set_terminal_state2(ϕ0, θmf, param3b, LPOArrival)
                    # in Sun-B1 frame
                    xf_sb1 = vcat(SailorMoon.transform_EMrot_to_SunB1(xf, θsf, param3b.oms, param3b.as), 1.0)
                    
                    push!(grids, [ϕ0, ϵr, ϵv, θsf, xf_sb1])
    #                 println(xf)
                    
                end
            end
        end
    end

    svf_ = grids[1][5]
    params = [param3b.mu2, param3b.mus, param3b.as, pi - grids[1][4], param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 0.0, 0.0]
    tspan = [0, -tof_bck]

    prob_base = ODEProblem(SailorMoon.rhs_bcr4bp_sb1frame2!, svf_, tspan, params)
    # sol_bck = SailorMoon.integrate_rk4(prob_base, 0.001, cbs)

    # sol = solve(prob_base, Tsit5(), reltol=1e-12, abstol=1e-12);

    ptraj = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.2)

    ## data extraction and make csv
    # make dataframe
    entries = ["id", "phi0", "epsr", "epsv", "thetaf", "ra", "rp", "dt1", "dt2", "x_ra", "x_rp", "tof", "lfb"]
    df = DataFrame([ name =>[] for name in entries])
    id = 1

    for i in 1:size(grids)[1]
        global lfb_count = 0
        global prob_base = remake(prob_base, u0=grids[i][5], p=[param3b.mu2, param3b.mus, param3b.as, pi - grids[i][4], param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 0.0, 0.0])
        global sol = SailorMoon.integrate_rk4(prob_base, 0.001, cbs, false)
        # global sol2 = SailorMoon.integrate_rk4(prob_base, 0.001, cbs, true)

        if sol.t[end] > prob_base.tspan[2]
            println("$i is terminated!")
            t_vec = -sol.t
            
            rp = sqrt((sol.u[end][1]-param3b.as)^2 + sol.u[end][2]^2 + sol.u[end][3]^2)
            x_rp = sol.u[end]

            # find apoapsis      
            r_vec = sqrt.((hcat(sol.u...)[1,:] .+ param3b.mu2.*cos.(prob_base.p[4] .+ param3b.oml .* sol.t) .- [param3b.as]).^2
                        .+ (hcat(sol.u...)[2,:] .+ param3b.mu2.*sin.(prob_base.p[4] .+ param3b.oml .* sol.t)).^2
                        .+  hcat(sol.u...)[3,:].^2)

            ra, id_ra = findmax(r_vec)
            x_ra  = sol.u[id_ra]
            dt_ra = - sol.t[id_ra]
            dt_rp = - sol.t[end] - dt_ra 

            # flag: lunar flyby? 
            rm_vec = sqrt.((hcat(sol.u...)[1,:] .- (1-param3b.mu2).*cos.(prob_base.p[4].+param3b.oml.*sol.t) .- [param3b.as]).^2
                        + (hcat(sol.u...)[2,:] .- (1-param3b.mu2).*sin.(prob_base.p[4].+param3b.oml.*sol.t)).^2
                        +  hcat(sol.u...)[3,:].^2)
            id_lfb = findall(rm_vec .< (66100 / param3b.lstar) .&& t_vec .> (10*86400/param3b.tstar))
            
    #         println(sol.event_states)
        
            if ~isempty(id_lfb)
                for k in id_lfb
    #                 print("\r$k   ")
    #                 println(perilune_cond(sol.u[k], sol.t[k], 0.0))
                    if ~isnan(perilune_cond(sol.u[k], sol.t[k], 0.0))
                        global lfb_count += 1
                    end
                end    
                println("  -> lfb_count: $lfb_count")
            end

            tof = -sol.t[end]
            
            if ra > 2.0
                # scatter!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], color=:blue, shape=:circle, markersize=2.0, label="flyby?")
                # plot!(ptraj, hcat(sol2.u...)[1,:], hcat(sol2.u...)[2,:], label="sol $i")
                    
                ϕ0  = grids[i][1]
                ϵr  = grids[i][2]
                ϵv  = grids[i][3]
                θsf = grids[i][4]
                push!(df, [id, ϕ0, ϵr, ϵv, θsf, rp, ra, dt_ra, dt_rp, x_ra, x_rp, tof, lfb_count])
                global id += 1
            
            end
        end
    end

end # end @everywhere




# println(df)
CSV.write("grid_search0130.csv", df)
    
# moon = plot_circle(1-param3b.mu2, param3b.as , 0.0)
# earth = plot_circle(param3b.mu2, param3b.as, 0.0)
# LEO_ub = plot_circle(3*param3b.mu2, param3b.as, 0.0)
# moon_soi_ub = plot_circle(1-param3b.mu2+66000/param3b.lstar, param3b.as, 0.0)
# moon_soi_lb = plot_circle(1-param3b.mu2-66000/param3b.lstar, param3b.as, 0.0)

# plot!(ptraj, earth[1,:], earth[2,:], c=:green, lw=1.0, label="earth")
# plot!(ptraj, moon[1,:], moon[2,:], c=:orange, lw=1.0, label="moon")
# plot!(ptraj, moon_soi_ub[1,:], moon_soi_ub[2,:], c=:grey, lw=1.0, label="moon_soi_ub")
# plot!(ptraj, moon_soi_lb[1,:], moon_soi_lb[2,:], c=:grey, lw=1.0, label="moon_soi_lb")

# gui(ptraj)
# savefig(ptraj, "plot_test.png")