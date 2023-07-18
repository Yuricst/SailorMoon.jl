using Distributed

# addprocs()
# @show procs()

@everywhere  begin
    using DifferentialEquations
    using Plots
    using LinearAlgebra
    import ForwardDiff
    import DiffResults
    using AstrodynamicsBase
    using Printf
    using DataFrames
    using JSON
    using CSV
    using ODEInterfaceDiffEq

    include("../src/SailorMoon.jl")
    include("../../julia-r3bp/R3BP/src/R3BP.jl")


    ## ====== something you want to change =====================
    
    out_fname = "data/grid_search_Tsit5_0717_tidal_Thrust.csv"
    dv_fun = SailorMoon.dv_tidal_dir_sb1frame

    ## =========================================================

    param3b = SailorMoon.dynamics_parameters()

    # some inputs needed for the thrust profile
    tmax_si = 400e-3   # N
    isp_si = 2500 # sec
    mdot_si = tmax_si / (isp_si * 9.81)
    mstar = 2500  # kg

    global earth_leo_ub = 30000 / param3b.lstar  # km
    global earth_leo_lb = 3000 / param3b.lstar  # km

    tmax = AstrodynamicsBase.dimensional2canonical_thrust(
        tmax_si, mstar, param3b.lstar, param3b.tstar
    )
    mdot = AstrodynamicsBase.dimensional2canonical_mdot(
        mdot_si, mstar, param3b.tstar
    )

    if dv_fun == SailorMoon.dv_no_thrust 
        tmax = 0.0
        mdot = 0.0 
    end 

    #### CALLBACK FUNCTIONS #################
    # store the apoapsis value
    function apoapsis_cond(u,t,int)
        θmLPO = int.p[4]
        θm = θmLPO + param3b.oml * t
        
        # for convenience, instead of taking the distance of SC-moon, assume r_a > 2.0
        if sqrt((u[1]-param3b.as)^2 + u[2]^2 + u[3]^2) > 2.0
            # condition 2: dot product of velocity and position is zero
            return dot([u[1] - param3b.as - (-param3b.mu2 * cos(θm)), u[2] - (-param3b.mu2 * sin(θm)), u[3]], u[4:6])
        else
            return NaN
        end
    end

    # store the periapsis around Earth value and terminate
    function periapsis_LEO_cond(u,t,int)

        ub = earth_leo_ub 
        r = sqrt((u[1] - param3b.as)^2 + u[2] ^2 + u[3]^2)  
        val = NaN 

        if r < ub
            θmLPO = int.p[4]
            θm = θmLPO + param3b.oml * t

            r_sc_earth = [u[1] - param3b.as - (-param3b.mu2 * cos(θm)), u[2] - (-param3b.mu2 * sin(θm)), u[3]]
            v_sc = u[4:6]

            val = dot(r_sc_earth, v_sc) 
        end

        return val 
    end

    # general periapsis condition w.r.t. Earth, inside the Earth 
    function periapsis_cond(u,t,int)

        val = NaN 
        r = sqrt((u[1] - param3b.as)^2 + u[2] ^2 + u[3]^2)  

        if abs(t) > 3  && r < param3b.mu1  # second condiiton enforces the SC to be inside the lunar radius

            θmLPO = int.p[4]
            θm = θmLPO + param3b.oml * t

            r_sc_earth = [u[1] - param3b.as - (-param3b.mu2 * cos(θm)), u[2] - (-param3b.mu2 * sin(θm)), u[3]]
            v_sc = u[4:6]

            val = dot(r_sc_earth, v_sc) 
        end

        return val 
    end

    # count the lunar radius crossing but not flyby
    function lunar_radius_cond(u,t,int)

        if isa(int, AbstractFloat)
            θmLPO = int
        else
            θmLPO = int.p[4]
        end
        θm  = θmLPO + param3b.oml * t
        
        # moon-SC distance
        r = sqrt((u[1] - param3b.as - (1-param3b.mu2)*cos(θm))^2 + (u[2] - (1-param3b.mu2)*sin(θm)) ^2 + u[3]^2)
        moon_soi = 66100 / param3b.lstar

        # usually LPO ~ apogee is about Dt = 12~14
        if abs(t) > 8 && r - moon_soi > 0 
            return sqrt((u[1] - param3b.as)^2 + u[2]^2 + u[3]^2) - param3b.mu1 
        else
            return NaN
        end

    end

    # count the lunar radius crossing but not flyby
    function lunar_radius_cond2(u,t,int)
        
        # pick up if SC get closer to the lunar SMA in S-B1 frame (consider the SOI of the moon 
        # so that we can avoid the integration "skipping" the r < lunar raidus domain)
        if abs(t) < 5  
            r = norm(u[1:3] - [param3b.as, 0, 0])
            return r - param3b.mu1 - (66100 / 4 / param3b.lstar)
        else
            return NaN
        end

    end


    function perilune_cond(u,t,int)
        # if abs(t) > 1
            if isa(int, AbstractFloat)
                θmLPO = int
            else
                θmLPO = int.p[4]
            end

            θm  = θmLPO + param3b.oml * t
            
            # moon-SC distance
            r = sqrt((u[1] - param3b.as - (1-param3b.mu2)*cos(θm))^2 + (u[2] - (1-param3b.mu2)*sin(θm)) ^2 + u[3]^2)
            moon_soi = 66100 / param3b.lstar
            v_moon   = (1-param3b.mu2)*param3b.oml * [-sin(θm), cos(θm), 0] 
            
            # strict perilune condition
            # if r < moon_soi #&& -t > (10*86400 / param3b.tstar)
            #     return dot([u[1] - param3b.as - (1-param3b.mu2) * cos(θm), u[2] - (1-param3b.mu2) * sin(θm), u[3]], u[4:6]-v_moon)
            # else 
            #     return NaN
            # end

            # check when the SC passes the boundary of the lunar SOI
            return moon_soi - r
        # else 
            # return NaN
        # end
        
    end

    function switch2ballistic_cond(u,t,int)
        if abs(t) > 8 
            r = norm(u[1:3] - [param3b.as, 0, 0])
            return r - param3b.mu1
        else 
            return NaN 
        end
    end


    # switch on the engine after escaping from the "invariant manifold"
    function wbs_boundary_cond(u,t,int)
        coast_time = 10 * (24*60*60) / param3b.tstar
        return abs(t) - coast_time
    end 

    
    function collision_cond(u,t,int)
        val = NaN 
        r = sqrt((u[1] - param3b.as)^2 + u[2] ^2 + u[3]^2)  # earth-SC distance 
        return (r - earth_leo_lb*0.8) * (r - 8)  # terminate if it's too close to Earth, and remove if it's too far from Earth
    end

    ### affect!
    terminate_affect!(int) = terminate!(int)
    no_affect!(int) = NaN

    function turn_off_engine_affect!(int)
        int.p[7]  = 0.0    # τ
        int.p[10] = 0.0    # mdot
        int.p[11] = 0.0    # tmax
        # println(" engine turned off at t=", int.t)
    end

    function turn_on_engine_affect!(int)
        int.p[7]  = 1.0     # τ
        int.p[10] = mdot    # mdot
        int.p[11] = tmax    # tmax
        # println(" engine turned on at t=", int.t)
    end

end

@everywhere begin
    t0 = time()
    lps = SailorMoon.lagrange_points(param3b.mu2)

    ## set up of initial condition (Lyapunov orbit)
    lp = 2
    Az_km = 0.0
    println("Halo guess Az_km: $Az_km")
    northsouth = 3   # 1 or 3
    guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
    res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

    x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
    prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
    # for Halo propagation, keep the tol as tight as possible 
    sol = solve(prob_cr3bp_stm, Tsit5(); reltol=1e-12, abstol=1e-12) #, saveat=LinRange(0, period, n+1))

    monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
    ys0 = R3BP.get_eigenvector(monodromy, true, 1) # monodromy eigenvector

    ### Grid search parameters
    ### OPTION 1; grid generations
    n = 60
    m = 300
    ϕ_vec    = LinRange(0, 2*pi, m+1)[1:m]  # [0.335103216] [0.0]    
    θs_vec   = LinRange(0, 2*pi, n+1)[1:n]  # [0.104719755] [180/180*pi]   # this should be 1 but not 2pi!!!!!!!
    epsr_vec = 10.0 .^(-5)
    epsv_vec = 10.0 .^(-5)
    tof_bck  = 120 * 86400 / param3b.tstar

    grids = []
    for ϕ0 in ϕ_vec
        for ϵr in epsr_vec
            for ϵv in epsv_vec
                for θsf in θs_vec
                    
                    # arrival LPO object
                    LPOArrival = SailorMoon.CR3BPLPO2(
                        res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12, 0.005
                    );
                    
                    xf = SailorMoon.set_terminal_state2(ϕ0, pi-θsf, param3b, LPOArrival)
                    # in Sun-B1 frame
                    xf_sb1 = vcat(SailorMoon.transform_EMrot_to_SunB1(xf, pi-θsf, param3b.oml, param3b.as), 1.0)
                    
                    # println("xf_sb1: ", xf_sb1)
                    push!(grids, [ϕ0, ϵr, ϵv, θsf, xf_sb1])
                    
                end
            end
        end
    end

    ### OPTION 2;  a bit of cheating... re-filtering the existing file 
    # filename = "data/grid_search_Tsit5_0525_EMrotThrust.csv"
    # df = DataFrame(CSV.File(filename))
    # ϕ_vec    = df.phi0
    # θs_vec   = df.thetasf
    # epsr_vec = 10.0 .^(-5)
    # epsv_vec = 10.0 .^(-5)
    # tof_bck  = 120 * 86400 / param3b.tstar
    # grids = []

    # for i in collect(1:1:length(ϕ_vec))
    #     ϕ0  = ϕ_vec[i]
    #     θsf = θs_vec[i]
    #     # arrival LPO object
    #     LPOArrival = SailorMoon.CR3BPLPO2(
    #         res.x0, res.period, ys0, prob_cr3bp_stm, epsr_vec, epsv_vec, Tsit5(), 1e-12, 1e-12, 0.005
    #     );
        
    #     xf = SailorMoon.set_terminal_state2(ϕ0, pi-θsf, param3b, LPOArrival)
    #     xf_sb1 = vcat(SailorMoon.transform_EMrot_to_SunB1(xf, pi-θsf, param3b.oml, param3b.as), 1.0)  # in Sun-B1 frame
    #     push!(grids, [ϕ0, epsr_vec, epsv_vec, θsf, xf_sb1])
    
    # end

end


# include callback functions 
@everywhere begin
    ## affect!: upcrossing (neg->pos), affect_neg: downcrossing (pos->neg)

    apoapsis_cb    = ContinuousCallback(apoapsis_cond, no_affect!; rootfind=false, save_positions=(false,true))
    # periapsis_cb = ContinuousCallback(periapsis_LEO_cond, terminate_affect!)
    periapsis_cb   = ContinuousCallback(periapsis_cond, terminate_affect!)
    # aps_cb       = ContinuousCallback(aps_cond, no_affect!; rootfind=false, save_positions=(false,true))
    perilune_cb    = ContinuousCallback(perilune_cond, nothing, no_affect!; rootfind=false, save_positions=(true,false))
    # lunar_rad_cb = ContinuousCallback(lunar_radius_cond, nothing, no_affect!; rootfind=false, save_positions=(false, true))
    lunar_rad_cb   = ContinuousCallback(lunar_radius_cond2, nothing, no_affect!; rootfind=false, save_positions=(true, false))
    ballistic_cb   = ContinuousCallback(switch2ballistic_cond, nothing, turn_off_engine_affect!; rootfind=false, save_positions=(true, false))
    # wbs_cb         = ContinuousCallback(wbs_boundary_cond, turn_on_engine_affect!; rootfind=false, save_positions=(false,false))
    collision_cb   = ContinuousCallback(collision_cond, terminate_affect!; rootfind=false, save_positions=(false, true))

    cbs = CallbackSet(apoapsis_cb, periapsis_cb, perilune_cb, lunar_rad_cb, ballistic_cb, collision_cb)  # wbs_cb

    svf_ = zeros(Float64, 1, 7)
    tspan = [0, -tof_bck]

    # EMrot EOM
    # params = [param3b.mu2, param3b.mus, 0.0, param3b.as, param3b.oms, 0.0, 0.0, 0.0, 0.0, 0.0]
    # prob = ODEProblem(R3BP.rhs_bcr4bp_thrust!, svf_, tspan, params)

    # SB1rot EOM
    params = [param3b.mu2, param3b.mus, param3b.as, pi, param3b.oml, param3b.omb, 1.0, 0.0, 0.0, 0.0, 0.0, dv_fun, param3b.tstar]
    prob = ODEProblem(SailorMoon.rhs_bcr4bp_sb1frame2_thrust_bal!, svf_, tspan, params)

    ## make ensemble problems
    function prob_func(prob, i, repeat)
        print("\rproblem # $i")
        remake(prob, u0=grids[i][5], p=[param3b.mu2, param3b.mus, param3b.as, pi - grids[i][4], param3b.oml, param3b.omb, 1.0, 0.0, 0.0, mdot, tmax, dv_fun, param3b.tstar])
        # remake(prob, u0=grids[i][5], p=[param3b.mu2, param3b.mus, param3b.as, pi - grids[i][4], param3b.oml, param3b.omb, 0.0, 0.0, 0.0, 0.0, 0.0, dv_fun, param3b.tstar])
    end
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

sim = solve(ensemble_prob, AutoTsit5(Rosenbrock23()), EnsembleThreads(); trajectories=length(grids),
            callback=cbs, reltol=1e-12, abstol=1e-12,
            save_everystep=true);
tofs = [sol.t[end] for sol in sim]

# sim = solve(ensemble_prob, RK4(),  dt=0.005, adaptive=false, EnsembleThreads(), trajectories=length(grids),
#             callback=cbs,
#             save_everystep=true);

ptraj = plot(size=(700,500), frame_style=:box, aspect_ratio=:equal, grid=0.2, xlim=[385,392], ylim=[-3,3])
# plot!(ptraj, xlims=(385,395), ylims=(-8,8))
cmap = :viridis


## data extraction and make csv
# make dataframe
entries = [
    "id", "phi0", "epsr", "epsv", "thetasf",
    "rp_kep", "ra_kep", "alpha", 
    "ra", "dt1", "dt2",
    "x_lpo", "y_lpo", "z_lpo", "xdot_lpo", "ydot_lpo", "zdot_lpo", "m_lpo",
    "x_ra", "y_ra", "z_ra", "xdot_ra", "ydot_ra", "zdot_ra", "m_ra",
    "x_rp", "y_rp", "z_rp", "xdot_rp", "ydot_rp", "zdot_rp", "m_rp",
    "x_lr", "y_lr", "z_lr", "xdot_lr", "ydot_lr", "zdot_lr", "m_lr", "t_lr",
    "tof", "lfb"
]
df = DataFrame([ name =>[] for name in entries])

global id = 1
color_start = "orange"
color_end = "blue"
color_gradation = cgrad([color_start, color_end], tof_bck)
# extract the ensemble simulation
for (i,sol) in enumerate(sim)
    global lfb_count = 0
    global lradius_cross = 0
    global flag = false  

    # println(sol.t[end] > -tof_bck)
    
    if sol.t[end] > -tof_bck
        # println("sol $i : ", sol.retcode)

        # the terminal point (i.e., first (time backward) periapsis w.r.t. Earth) is in the range of LEO? 
        r_entry = sol.u[end][1:6]
        θsf = grids[i][4]
        θmf = pi - θsf
        θm0 = θmf + param3b.oml * sol.t[end]  # moon angle at the periapsis
        r_sc_earth = norm([r_entry[1] - param3b.as - (-param3b.mu2 * cos(θm0)), r_entry[2] - (-param3b.mu2 * sin(θm0)), r_entry[3]])

        if earth_leo_lb < r_sc_earth && r_sc_earth < earth_leo_ub

            # Using Keplar 2 body problem, find the rp analytically
            r_entry_EIne = SailorMoon.transform_sb1_to_EearthIne(r_entry, θm0, param3b.oml, param3b.mu2, param3b.as)
            h_entry = cross(r_entry_EIne[1:3], r_entry_EIne[4:6])

            # we want SC to leave from Earth in CCW direction 
            if h_entry[3] > 0.0

                coe_entry = cart2kep(r_entry_EIne, param3b.mu1)
                sma, ecc, inc, OMEGA, omega, nu = coe_entry
                
                # 06/15/2023: I don't like this if statement, but this is enevitable... 
                if ecc < 1 + 1e-6
                    rp_kep = sma * (1-ecc)
                    ra_kep = sma * (1+ecc)
                    
                    # generate state @ periapsis
                    state_rp = kep2cart([sma, ecc, inc, OMEGA, omega, 0.0], param3b.mu1)
                    state_rp = SailorMoon.transform_EearthIne_to_sb1(state_rp, θm0, param3b.oml, param3b.mu2, param3b.as)
                    # println("state_rp: ", state_rp)

                    # obtian the eccentric anomaly & mean anomaly at entrance
                    cosE = (ecc + cos(nu)) / (1 + ecc*cos(nu))
                    sinE = sqrt(1-ecc^2) * sin(nu) / (1 + cos(nu))
                    E = atan(sinE, cosE)
                    M = E - ecc*sin(E)
                    n_ = sqrt(param3b.mu1 / sma^3) 
                    tof_finale = abs(M / n_)

                    tof_tot = -sol.t[end] + tof_finale

                    t_vec = -sol.t
                    # rp = sqrt((sol.u[end][1]-param3b.as)^2 + sol.u[end][2]^2 + sol.u[end][3]^2)
                    # x_rp = sol.u[end]

                    # find apoapsis      
                    # relative to Earth
                    # r_vec = sqrt.((hcat(sol.u...)[1,:] .+ param3b.mu2.*cos.(θmf .+ param3b.oml .* sol.t) .- [param3b.as]).^2
                    #             .+ (hcat(sol.u...)[2,:] .+ param3b.mu2.*sin.(θmf .+ param3b.oml .* sol.t)).^2
                    #             .+  hcat(sol.u...)[3,:].^2)
                    
                    # relative to B1 
                    r_vec = sqrt.((hcat(sol.u...)[1,:] .- [param3b.as]).^2
                                .+ hcat(sol.u...)[2,:].^2
                                .+ hcat(sol.u...)[3,:].^2)

                    ra, id_ra = findmax(r_vec)
                    dt_ra = - sol.t[id_ra]
                    dt_rp = tof_tot - dt_ra 

                    x_ra    = sol.u[id_ra][1]
                    y_ra    = sol.u[id_ra][2]
                    z_ra    = sol.u[id_ra][3]
                    xdot_ra = sol.u[id_ra][4]
                    ydot_ra = sol.u[id_ra][5]
                    zdot_ra = sol.u[id_ra][6]
                    m_ra    = sol.u[id_ra][7]                   

                    # # flag: lunar flyby? (cf. moon SOI = 66100 km)
                    # rm_vec = sqrt.((hcat(sol.u...)[1,:] .- (1-param3b.mu2).*cos.(θmf .+ param3b.oml.*sol.t) .- [param3b.as]).^2
                    #             +  (hcat(sol.u...)[2,:] .- (1-param3b.mu2).*sin.(θmf .+ param3b.oml.*sol.t)).^2
                    #             +   hcat(sol.u...)[3,:].^2)
                    # id_lfb = findall(rm_vec .< (66100 / param3b.lstar) .* t_vec .> (10*86400/param3b.tstar))
                            
                    # if ~isempty(id_lfb)
                    #     for k in id_lfb
                    #         if ~isnan(perilune_cond(sol.u[k], sol.t[k], pi - θsf))
                    #             global lfb_count += 1
                    #         end
                    #     end    
                    # end


                    for (i, x) in enumerate(sol.u)
                        t = sol.t[i]
                        c1 = perilune_cond(x, t, pi-θsf)
                        c2 = lunar_radius_cond2(x, t, pi-θsf)
                        
                        # println("prilune cond: ", c1, " lunar rad cond: ", c2, " , t = ", t)
                        if  abs(c1) < 5e-2 
                            global lfb_count += 1
                        end

                        if abs(c2) < 3e-2
                            # println("watch out, lunar flyby is happening too early... ")
                            global flag = true 
                        end 

                    end

                    if flag == false 

                        r_vec[1:id_ra] = 100 * ones(Float64, (1,id_ra)) # dummy variables so that the id_lunar_rad occurs after the apoapsis
                        id_lunar_rad   = findmin(abs.(r_vec .- param3b.mu1))
                        id_lunar_rad   = id_lunar_rad[2]
                        x_l    = sol.u[id_lunar_rad][1]
                        y_l    = sol.u[id_lunar_rad][2]
                        z_l    = sol.u[id_lunar_rad][3]
                        xdot_l = sol.u[id_lunar_rad][4]
                        ydot_l = sol.u[id_lunar_rad][5]
                        zdot_l = sol.u[id_lunar_rad][6]
                        m_l    = sol.u[id_lunar_rad][7]  
                        t_lrad = -sol.t[id_lunar_rad]

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
                        α = acos(dot(vec, x_unit) / norm(vec))  # velocity vector angle at rp 

                        if cross(x_unit, vec)[3] <= 0
                            α = -α
                        end

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
                        # plot!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:])

                        push!(df, [id, ϕ0, ϵr, ϵv, θsf, 
                                rp_kep, ra_kep, α, 
                                ra, dt_ra, dt_rp, 
                                x_ini, y_ini, z_ini, xdot_ini, ydot_ini, zdot_ini, m_ini,
                                x_ra, y_ra, z_ra, xdot_ra, ydot_ra, zdot_ra, m_ra,
                                x_rp, y_rp, z_rp, xdot_rp, ydot_rp, zdot_rp, m_rp,
                                x_l, y_l, z_l, xdot_l, ydot_l, zdot_l, m_l, t_lrad,
                                tof_tot, lfb_count])
                        println("idx $i is a success!")

                        color = color_gradation[round(-sol.t[end])]
                        # plot!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], color=color, label="", linewidth=0.8)

                        global id += 1

                        plot!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], label="", linewidth=0.8)

                    end
                end
            end
        end
    end

    plot!(ptraj, hcat(sol.u...)[1,:], hcat(sol.u...)[2,:], label="", linewidth=0.8)

end

scatter!(ptraj, [param3b.as], [0.0], label="")  # roughly earth
circle = SailorMoon.plot_circle(1-param3b.mu2, param3b.as, 0.0)  # moon
plot!(ptraj, circle[1,:], circle[2,:], label="")
display(ptraj)

# print(df)
CSV.write(out_fname, df)




