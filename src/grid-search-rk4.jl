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

    param3b = SailorMoon.dyanmics_parameters()
    lps = SailorMoon.lagrange_points(param3b.mu2)


    # some inputs needed for the thrust profile
    tmax_si = 280e-3 * 4  # N
    isp_si = 1800 # sec
    mdot_si = tmax_si / (isp_si * 9.81)
    mstar = 2500  # kg

    tmax = AstrodynamicsBase.dimensional2canonical_thrust(
        tmax_si, mstar, param3b.lstar, param3b.tstar
    )
    mdot = AstrodynamicsBase.dimensional2canonical_mdot(
        mdot_si, mstar, param3b.tstar
    )

    dv_fun = SailorMoon.dv_no_thrust

    global earth_leo_ub = 9000 / param3b.lstar   # 10000 / param3b.lstar  # km
    global earth_leo_lb = 4500 / param3b.lstar  # km

    #### CALLBACK FUNCTIONS #################
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
        ub = earth_leo_ub 
        lb = earth_leo_lb 
        
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

    # terminate if the SC is sufficiently close to the earth 
    function proximity_earth_cond(u,t,int)
        θm0 = prob_base.p[4]
        θm = θm0 + param3b.oml * t
        r = sqrt((u[1] - param3b.as - (-param3b.mu2 * cos(θm)))^2 + (u[2] - (-param3b.mu2 * sin(θm))) ^2 + u[3]^2)  # SC-earth distance
        return r - 12000 / param3b.lstar
    end
    

        
    terminate_affect!() = true
    no_affect!() = false

    function plot_circle(radius, x, y, n=50)
        circle = zeros(2,n)
        thetas = LinRange(0.0, 2π, n)
        for i = 1:n
            circle[1,i] = radius*cos(thetas[i]) + x
            circle[2,i] = radius*sin(thetas[i]) + y
        end
        return circle
    end



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
    n = 120
    θs_vec   = [180/180*pi] #LinRange(0, 2*pi, n+1)[1:n]  # [3.76991118430775]   #[180/180*pi]  #
    ϕ_vec    = [0.0] # LinRange(0, 2*pi, n+1)[1:n] # [0.628318530717958]  [0.0]    #
    epsr_vec = 10.0 .^(-6)
    epsv_vec = 10.0 .^(-6)
    tof_bck  = 120 * 86400 / param3b.tstar

    # include callback functions 
    apoapsis_cb  = ContinuousCallback(apoapsis_cond, no_affect!)
    periapsis_cb = ContinuousCallback(periapsis_cond, terminate_affect!)
    perilune_cb  = ContinuousCallback(perilune_cond, no_affect!)
    prixmity_earth_cb = ContinuousCallback(proximity_earth_cond, terminate_affect!)
    cbs = [apoapsis_cb, prixmity_earth_cb, perilune_cb];

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
    tspan = [0, -tof_bck]

    params = [param3b.mu2, param3b.mus, param3b.as, pi - grids[1][4], param3b.oml, param3b.omb, 1.0, 0.0, 0.0, 0.01, 0.01, SailorMoon.dv_sun_dir_angles2]
    prob_base = ODEProblem(SailorMoon.rhs_bcr4bp_sb1frame2_thrust!, svf_, tspan, params)

    ptraj = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.2)

    ## data extraction and make csv
    # make dataframe
    entries = ["id", "phi0", "epsr", "epsv", "thetaf",
              "rp_kep", "ra_kep", "alpha",  # important guys
              "ra", "dt1", "dt2",
              "x_ini", "y_ini", "z_ini", "xdot_ini", "ydot_ini", "zdot_ini", "m_ini",
              "x_ra", "y_ra", "z_ra", "xdot_ra", "ydot_ra", "zdot_ra", "m_ra",
              "x_rp", "y_rp", "z_rp", "xdot_rp", "ydot_rp", "zdot_rp", "m_rp",
              "tof", "m0", "lfb"]
    df = DataFrame([ name =>[] for name in entries])
    id = 1

    for i in 1:size(grids)[1]
        global lfb_count = 0
        global prob_base  = remake(prob_base, u0=grids[i][5], p=[param3b.mu2, param3b.mus, param3b.as, pi - grids[i][4], param3b.oml, param3b.omb, 1.0, 0.0, 0.0, mdot, tmax, dv_fun])
        global sol  = SailorMoon.integrate_rk4(prob_base, 0.001, cbs, false)
        
        # did SC hit proximity of the earth? 
        if sol.t[end] > prob_base.tspan[2]
            println("$i is terminated!")
            
            # Using Keplar 2 body problem, find the rp analytically
            θsf = grids[i][4]
            r_entry = sol.u[end][1:6]
            θmf = pi - θsf
            θm0 = θmf + param3b.oml * sol.t[end]

            r_entry_EIne = SailorMoon.transform_sb1_to_EearthIne(r_entry, θm0, param3b.oml, param3b.mu2, param3b.as)
            h_entry = cross(r_entry_EIne[1:3], r_entry_EIne[4:6])
            println(θm0)
            println(r_entry_EIne[1:3], r_entry_EIne[4:6])
            println(sol.u[end])

            # we want SC to leave from Earth in CCW direction 
            if h_entry[3] > 0.0
                coe_entry = cart2kep(r_entry_EIne, param3b.mu1)
                sma, ecc, inc, OMEGA, omega, nu = coe_entry

                rp_kep = sma * (1-ecc)
                ra_kep = sma * (1+ecc)
                println("rp_kep", rp_kep)
                println("ra_kep", ra_kep)

                if ra_kep > 2.0
                    # generate state @ periapsis
                    state_rp = kep2cart([sma, ecc, inc, OMEGA, omega, 0.0], param3b.mu1)
                    state_rp = SailorMoon.transform_EearthIne_to_sb1(state_rp, θm0, param3b.oml, param3b.mu2, param3b.as)
                    
                    # obtian the eccentric anomaly & mean anomaly at entrance
                    cosE = (ecc + cos(nu)) / (1 + ecc*cos(nu))
                    sinE = sqrt(1-ecc^2) * sin(nu) / (1 + cos(nu))
                    E = atan(sinE, cosE)
                    M = E - ecc*sin(E)
                    n = sqrt(param3b.mu1 / sma^3) 
                    tof_finale = abs(M / n)

                    tof_tot = -sol.t[end] + tof_finale

                    t_vec = -sol.t
                    # rp = sqrt((sol.u[end][1]-param3b.as)^2 + sol.u[end][2]^2 + sol.u[end][3]^2)
                    # x_rp = sol.u[end]

                    # find apoapsis      
                    r_vec = sqrt.((hcat(sol.u...)[1,:] .+ param3b.mu2.*cos.(prob_base.p[4] .+ param3b.oml .* sol.t) .- [param3b.as]).^2
                                .+ (hcat(sol.u...)[2,:] .+ param3b.mu2.*sin.(prob_base.p[4] .+ param3b.oml .* sol.t)).^2
                                .+  hcat(sol.u...)[3,:].^2)

                    ra, id_ra = findmax(r_vec)
                    dt_ra = - sol.t[id_ra]
                    dt_rp = tof_tot - dt_ra 

                    x_ra = sol.u[id_ra][1]
                    y_ra = sol.u[id_ra][2]
                    z_ra = sol.u[id_ra][3]
                    xdot_ra = sol.u[id_ra][4]
                    ydot_ra = sol.u[id_ra][5]
                    zdot_ra = sol.u[id_ra][6]
                    m_ra = sol.u[id_ra][7]                   

                    # flag: lunar flyby? 
                    rm_vec = sqrt.((hcat(sol.u...)[1,:] .- (1-param3b.mu2).*cos.(prob_base.p[4].+param3b.oml.*sol.t) .- [param3b.as]).^2
                                +  (hcat(sol.u...)[2,:] .- (1-param3b.mu2).*sin.(prob_base.p[4].+param3b.oml.*sol.t)).^2
                                +   hcat(sol.u...)[3,:].^2)
                    id_lfb = findall(rm_vec .< (66100 / param3b.lstar) .* t_vec .> (10*86400/param3b.tstar))
                            
                    if ~isempty(id_lfb)
                        for k in id_lfb
                            if ~isnan(perilune_cond(sol.u[k], sol.t[k], 0.0))
                                global lfb_count += 1
                            end
                        end    
                        # print(" ->> lfb_count; $lfb_count")
                    end

                    # obtain α
                    θs0 = θsf - param3b.oms * tof_tot
                    θm0 = π - θs0
                    rE = [
                        param3b.as - param3b.mu2 * cos(θm0),
                        -param3b.mu2 * sin(θm0),
                        0.0
                    ]

                    # r_sc - r_E
                    vec = state_rp[1:3] - rE 
                    x_unit = [1.0, 0.0, 0.0]
                    α = acos(dot(vec, x_unit) / norm(vec))

                    if cross(x_unit, vec)[3] <= 0
                        α = -α
                    end

                    m0 = sol.u[end][end]
            
                    ϕ0  = grids[i][1]
                    ϵr  = grids[i][2]
                    ϵv  = grids[i][3]
                    x_ini = sol.u[1][1]
                    y_ini = sol.u[1][2]
                    z_ini = sol.u[1][3]
                    xdot_ini = sol.u[1][4]
                    ydot_ini = sol.u[1][5]
                    zdot_ini = sol.u[1][6]
                    m_ini = sol.u[1][7]

                    x_rp = state_rp[1]
                    y_rp = state_rp[2]
                    z_rp = state_rp[3]
                    xdot_rp = state_rp[4]
                    ydot_rp = state_rp[5]
                    zdot_rp = state_rp[6]
                    m_rp = sol.u[end][7]

                    # scatter!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], color=:blue, shape=:circle, markersize=2.0, label="event?")
                    # plot!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], label="no thrust")
                
                    push!(df, [id, ϕ0, ϵr, ϵv, θsf, 
                               rp_kep, ra_kep, α, 
                               ra, dt_ra, dt_rp, 
                               x_ini, y_ini, z_ini, xdot_ini, ydot_ini, zdot_ini, m_ini,
                               x_ra, y_ra, z_ra, xdot_ra, ydot_ra, zdot_ra, m_ra,
                               x_rp, y_rp, z_rp, xdot_rp, ydot_rp, zdot_rp, m_rp,
                               tof_tot,  m0, lfb_count])
                    global id += 1
                end
            end
        end
    end

    # println(df)
    CSV.write("grid_search1129.csv", df)
        
    # moon = plot_circle(1-param3b.mu2, param3b.as , 0.0)
    # earth = plot_circle(param3b.mu2, param3b.as, 0.0)
    # LEO_ub = plot_circle(3*param3b.mu2, param3b.as, 0.0)
    # moon_soi_ub = plot_circle(1-param3b.mu2+66000/param3b.lstar, param3b.as, 0.0)
    # moon_soi_lb = plot_circle(1-param3b.mu2-66000/param3b.lstar, param3b.as, 0.0)


    # plot!(ptraj, hcat(sol_bck.u...)[1,:], hcat(sol_bck.u...)[2,:], color=:deeppink, label="RK4")
    # plot!(ptraj, Array(sol)[1,:], Array(sol)[2,:], color=:blue, linewidth=1.0, label="sol", linestyle=:solid)

    # plot!(ptraj, earth[1,:], earth[2,:], c=:green, lw=1.0, label="earth")
    # plot!(ptraj, moon[1,:], moon[2,:], c=:orange, lw=1.0, label="moon")
    # plot!(ptraj, moon_soi_ub[1,:], moon_soi_ub[2,:], c=:grey, lw=1.0, label="moon_soi_ub")
    # plot!(ptraj, moon_soi_lb[1,:], moon_soi_lb[2,:], c=:grey, lw=1.0, label="moon_soi_lb")

    # display(ptraj)

end

