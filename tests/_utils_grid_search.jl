"""
Utilities for grid-search function
"""

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

function earth_escape(u,t,int)
    r = sqrt((u[1] - param3b.as)^2 + u[2] ^2 + u[3]^2)  
    return r - 1.2 * 145*6378/param3b.lstar
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


apoapsis_cb    = ContinuousCallback(apoapsis_cond, no_affect!; rootfind=false, save_positions=(false,true))
# periapsis_cb = ContinuousCallback(periapsis_LEO_cond, terminate_affect!)
periapsis_cb   = ContinuousCallback(periapsis_cond, terminate_affect!)
# aps_cb       = ContinuousCallback(aps_cond, no_affect!; rootfind=false, save_positions=(false,true))
perilune_cb    = ContinuousCallback(perilune_cond, nothing, no_affect!; rootfind=false, save_positions=(true,false))
# lunar_rad_cb = ContinuousCallback(lunar_radius_cond, nothing, no_affect!; rootfind=false, save_positions=(false, true))
lunar_rad_cb   = ContinuousCallback(lunar_radius_cond2, nothing, no_affect!; rootfind=false, save_positions=(true, false))
ballistic_cb   = ContinuousCallback(switch2ballistic_cond, nothing, turn_off_engine_affect!; rootfind=false, save_positions=(true, false))
# wbs_cb         = ContinuousCallback(wbs_boundary_cond, turn_on_engine_affect!; rootfind=false, save_positions=(false,false))
earth_escape_cb = ContinuousCallback(earth_escape, terminate_affect!; rootfind=false, save_positions=(false,false))


cbs = CallbackSet(apoapsis_cb, periapsis_cb, perilune_cb, lunar_rad_cb, ballistic_cb, earth_escape_cb)  # wbs_cb
#cbs = CallbackSet(apoapsis_cb, perilune_cb, lunar_rad_cb, ballistic_cb, earth_escape_cb)  # wbs_cb
