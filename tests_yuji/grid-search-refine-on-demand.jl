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
using DataFrames


using ProgressMeter
using Printf
using Roots

include("../src/SailorMoon.jl")
include("../../julia-R3BP/R3BP/src/R3BP.jl")

param3b = SailorMoon.dynamics_parameters()
out_fname = "data/grid_search_Tsit5_0614_velThrust.csv"


# --------------------------------------------------------------------------- #
# some inputs needed for the thrust profile
tmax_si = 400e-3   # N
isp_si = 2500 # sec
mdot_si = tmax_si / (isp_si * 9.81)
mstar = 2500  # kg

global earth_leo_ub = 0.4   # 15000 / param3b.lstar  # km
global earth_leo_lb = 3000 / param3b.lstar  # km

tmax = AstrodynamicsBase.dimensional2canonical_thrust(
    tmax_si, mstar, param3b.lstar, param3b.tstar
)
mdot = AstrodynamicsBase.dimensional2canonical_mdot(
    mdot_si, mstar, param3b.tstar
)

# Choose DV transcription
# dv_no_thrust
# dv_tidal_dir_sb1frame
# dv_vel_dir_sb1frame (*)
# dv_EMrotdir_sb1frame (*)
dv_fun = SailorMoon.dv_vel_dir_sb1frame

if dv_fun == SailorMoon.dv_no_thrust
    tmax = 0.0
    mdot = 0.0 
end 

# define callbacks 
include("../tests/_utils_grid_search.jl")

# --------------------------------------------------------------------------- #

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
global id = 0



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
n = 10 # 60
m = 10 #300
#ϕ_vec    = LinRange(0, 2*pi, m+1)[1:m]  # [0.335103216] [0.0]    
θs_vec   = LinRange(0, 2*pi, n+1)[1:n]  # [0.104719755] [180/180*pi] 
ϵr, ϵv = 1e-5, 1e-5
tof_bck  = 120 * 86400 / param3b.tstar
tspan = [0, -tof_bck]


function find_perigee_SunB1(sol)
    r1min = 1e8
    for u in sol.u
        r_Earth = sqrt((u[1] - param3b.as)^2 + u[2] ^2 + u[3]^2)  
        if r_Earth < r1min
            r1min = r_Earth
        end
    end
    return r1min
end


"""
For given combination of θs and ϕ0, return perigee and ODEProblem's solution
"""
function get_perigee_sol(θsf::Real, ϕ0::Real=0.0)
    LPOArrival = SailorMoon.CR3BPLPO2(
        res.x0, res.period, ys0, prob_cr3bp_stm, ϵr, ϵv, Tsit5(), 1e-12, 1e-12, 0.005
    );

    # LOI state
    xf_sb1 = vcat(
        SailorMoon.transform_EMrot_to_SunB1(
            SailorMoon.set_terminal_state2(ϕ0, pi-θsf, param3b, LPOArrival),
            pi-θsf,
            param3b.oml,
            param3b.as
        ),
        1.0
    )
    params = [
        param3b.mu2, param3b.mus, param3b.as,
        pi - θsf,
        param3b.oml, param3b.omb, 1.0, 0.0, 0.0, mdot, tmax, dv_fun, param3b.tstar
    ]


    # ODEProblem
    prob = ODEProblem(
        SailorMoon.rhs_bcr4bp_sb1frame2_thrust_bal!,
        xf_sb1,
        tspan,
        params,
    )
    _sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, callback=cbs)
    _perigee = find_perigee_SunB1(_sol)
    return _sol, _perigee 
end

"""
Run `get_perigee_sol()` over a grid of θsfs, for fixed θ0
"""
function grid_perigees_sols(θsfs::Vector, ϕ0::Real)
    sols = []
    perigees = []
    @showprogress for θsf in θsfs
        _sol, _perigee = get_perigee_sol(θsf, ϕ0)
        push!(sols, _sol)
        push!(perigees, _perigee)
    end
    return sols, perigees
end

# 1. run to get perigee ----------------------------------------------------- #
target_perigee = (6378 + 200)/param3b.lstar

ϕ0 = deg2rad(145)  # fixed 

n_θsf = 180
θsfs = LinRange(0, 2π, n_θsf+1)[1:n_θsf]
θsfs = [el for el in θsfs]
#θsfs = LinRange(1, 2, n_θsf+1)[1:n_θsf]
#θsfs = LinRange(0, deg2rad(30), n_θsf+1)[1:n_θsf]
sols, _perigees = grid_perigees_sols(θsfs[:], ϕ0)
println("Minimum perigee: ", minimum(_perigees)*param3b.lstar, " km")

# 2. zoom into local minimum of perigee ------------------------------------ #
θsf_interest = Real[]
rp_threshold = 30000 / param3b.lstar
global Δdir = 0.0  # initialize 
for idx = 1:length(_perigees)-1
    if (_perigees[idx+1] - _perigees[idx])*Δdir < 0.0 && _perigees[idx] < rp_threshold
        # detected change of sign
        push!(θsf_interest, θsfs[idx])
    end
    global Δdir = _perigees[idx+1] - _perigees[idx]
end

# 3. Re-search by zooming into regions of θsf of interest -------------------- #
@printf("Number of Sun-angle candidates: %1.0f\n", length(θsf_interest))
n_refined = 30
refined_sols_perigees = []
for idx = 1:length(θsf_interest)
    dθs = deg2rad(10)
    θsf = θsf_interest[idx]
    θsfs_refined = LinRange(θsf - dθs, θsf + dθs, n_refined)
    θsfs_refined = [el for el in θsfs_refined]
    _sols_refined, _perigees_refined = grid_perigees_sols(θsfs_refined, ϕ0)
    # storage
    push!(refined_sols_perigees, [θsfs_refined, _sols_refined, _perigees_refined])
end
@printf("refined_sols_perigees = %1.0f\n", length(refined_sols_perigees))

# 4. Run regula-falsi ------------------------------------------------------ #
"""Residual function for root-solving at perigee"""
res_func = function (θsf, get_sol::Bool=false)
    _sol, rp = get_perigee_sol(θsf, ϕ0)
    if get_sol 
        return rp - target_perigee, _sol
    else
        return rp - target_perigee
    end
end

println("Refining solution via Regula-Falsi...")
solved_sols_perigees = []
@showprogress for ref in refined_sols_perigees
    θsfs_refined, _, perigees_refined = ref     # unpack
    # check through θsf 
    previous_sign = 0.0
    for (idx,rp) in enumerate(perigees_refined[1:end-1])
        rp_next = perigees_refined[idx+1]
        current_sign = (rp_next - target_perigee) * (rp - target_perigee)
        if current_sign < 0.0
            # detected change of sign, solve root-solving
            println("Run root-solving problem! θsfs_refined[idx] = ", θsfs_refined[idx]*180/π)
            println("Currnet rp - target_radius = ", rp - target_perigee)
            println("Next rp - target_radius = ", rp_next - target_perigee)

            θsf_solved = find_zero(
                res_func,
                (θsfs_refined[idx], θsfs_refined[idx+1]),
                Bisection()
            )
            _residual, solved_sol = res_func(θsf_solved, true)
            push!(
                solved_sols_perigees,
                [θsf_solved, solved_sol, target_perigee + _residual]
            )

            # store the results
            # the terminal point (i.e., first (time backward) periapsis w.r.t. Earth) is in the range of LEO? 
            state_rp = solved_sol.u[end][1:6]
            θsf      = θsf_solved
            # θmf = pi - θsf_solved
            # θm0 = θmf + param3b.oml * sol.t[end]  # moon angle at the periapsis
            # r_sc_earth = norm([r_entry[1] - param3b.as - (-param3b.mu2 * cos(θm0)), r_entry[2] - (-param3b.mu2 * sin(θm0)), r_entry[3]])

            # # Using Keplar 2 body problem, find the rp analytically
            # r_entry_EIne = SailorMoon.transform_sb1_to_EearthIne(r_entry, θm0, param3b.oml, param3b.mu2, param3b.as)
            # h_entry = cross(r_entry_EIne[1:3], r_entry_EIne[4:6])

            # # we want SC to leave from Earth in CCW direction 
            # # we may or may not need this 
            # # if h_entry[3] > 0.0

            #     coe_entry = cart2kep(r_entry_EIne, param3b.mu1)
            #     println("KOE: ", coe_entry)
            #     sma, ecc, inc, OMEGA, omega, nu = coe_entry

            #     rp_kep = sma * (1-ecc)
            #     ra_kep = sma * (1+ecc)
                
            #     # generate state @ periapsis
            #     state_rp = kep2cart([sma, ecc, inc, OMEGA, omega, 0.0], param3b.mu1)
            #     state_rp = SailorMoon.transform_EearthIne_to_sb1(state_rp, θm0, param3b.oml, param3b.mu2, param3b.as)
            #     # println("state_rp: ", state_rp)

            #     # obtian the eccentric anomaly & mean anomaly at entrance
            #     cosE = (ecc + cos(nu)) / (1 + ecc*cos(nu))
            #     sinE = sqrt(1-ecc^2) * sin(nu) / (1 + cos(nu))
            #     E = atan(sinE, cosE)
            #     M = E - ecc*sin(E)
            #     n_ = sqrt(param3b.mu1 / sma^3) 
            #     tof_finale = abs(M / n_)

            tof_tot = -solved_sol.t[end]  # + tof_finale
            println("tof_tot = ", tof_tot )

            # find apoapsis relative to B1 
            r_vec = sqrt.((hcat(solved_sol.u...)[1,:] .- [param3b.as]).^2 .+ hcat(solved_sol.u...)[2,:].^2 .+ hcat(solved_sol.u...)[3,:].^2)

            ra, id_ra = findmax(r_vec)
            dt_ra = - solved_sol.t[id_ra]
            dt_rp = tof_tot - dt_ra 

            x_ra    = solved_sol.u[id_ra][1]
            y_ra    = solved_sol.u[id_ra][2]
            z_ra    = solved_sol.u[id_ra][3]
            xdot_ra = solved_sol.u[id_ra][4]
            ydot_ra = solved_sol.u[id_ra][5]
            zdot_ra = solved_sol.u[id_ra][6]
            m_ra    = solved_sol.u[id_ra][7]                   

            r_vec[1:id_ra] = 100 * ones(Float64, (1,id_ra)) # dummy variables so that the id_lunar_rad occurs after the apoapsis
            id_lunar_rad   = findmin(abs.(r_vec .- param3b.mu1))
            id_lunar_rad   = id_lunar_rad[2]
            x_l    = solved_sol.u[id_lunar_rad][1]
            y_l    = solved_sol.u[id_lunar_rad][2]
            z_l    = solved_sol.u[id_lunar_rad][3]
            xdot_l = solved_sol.u[id_lunar_rad][4]
            ydot_l = solved_sol.u[id_lunar_rad][5]
            zdot_l = solved_sol.u[id_lunar_rad][6]
            m_l    = solved_sol.u[id_lunar_rad][7]  
            t_lrad = -solved_sol.t[id_lunar_rad]

            # obtain α
            θs0 = θsf_solved - param3b.oms * tof_tot
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

            # ϕ0  = grids[i][1]
            x_ini = solved_sol.u[1][1]
            y_ini = solved_sol.u[1][2]
            z_ini = solved_sol.u[1][3]
            xdot_ini = solved_sol.u[1][4]
            ydot_ini = solved_sol.u[1][5]
            zdot_ini = solved_sol.u[1][6]
            m_ini = solved_sol.u[1][7]

            x_rp = state_rp[1]
            y_rp = state_rp[2]
            z_rp = state_rp[3]
            xdot_rp = state_rp[4]
            ydot_rp = state_rp[5]
            zdot_rp = state_rp[6]
            m_rp = solved_sol.u[end][7]

            global id += 1 
            
            rp_kep, ra_kep = NaN, NaN
            lfb_count = NaN

            # println(sol.u)
            # println(sol.t)

            push!(df, [id, ϕ0, ϵr, ϵv, θsf_solved, 
                    rp_kep, ra_kep, α, 
                    ra, dt_ra, dt_rp, 
                    x_ini, y_ini, z_ini, xdot_ini, ydot_ini, zdot_ini, m_ini,
                    x_ra, y_ra, z_ra, xdot_ra, ydot_ra, zdot_ra, m_ra,
                    x_rp, y_rp, z_rp, xdot_rp, ydot_rp, zdot_rp, m_rp,
                    x_l, y_l, z_l, xdot_l, ydot_l, zdot_l, m_l, t_lrad,
                    tof_tot, lfb_count])
            # println("idx $i is a success!")
                
            # end
        end
    end
end

CSV.write(out_fname, df)
println("Database successfully written, ", id, " solutions found")

# 5. plot results --------------------------------------------------------- #
# A. plot trajectory
ptraj = plot(size=(600,600), frame_style=:box, aspect_ratio=:equal, grid_alpha=0.5, legend=false,
    xlabel="x, LU", ylabel="y, LU")

for (idx,ref) in enumerate(refined_sols_perigees)
    θsfs_refined, sols_refined, _ = ref  # unpack
    for (idx,sol) in enumerate(sols_refined)
        scatter!(ptraj, [Array(sol)[1,1]], [Array(sol)[2,1]], marker=:circle, color=:blue, label="LOI")
        plot!(ptraj, Array(sol)[1,:], Array(sol)[2,:], label="Traj", lw=0.7, color=:blue)
    end
end

# solved trajectory
for solved in solved_sols_perigees
    _, solved_sol, _ = solved     # unpack
    scatter!(ptraj, [Array(solved_sol)[1,1]], [Array(solved_sol)[2,1]], marker=:circle, color=:blue, label="LOI")
    plot!(ptraj, Array(solved_sol)[1,:], Array(solved_sol)[2,:], label="Traj", lw=0.7, color="lime")
end


# # B. plot perigee distribution
# ptrend = plot(size=(800,600), frame_style=:box, grid_alpha=0.5, legend=:topright,
#     xlabel="θsf", ylabel="Perigee, km")
# plot!(ptrend, θsfs*180/π, _perigees*param3b.lstar, marker=:circle, color=:blue, label="ELET Perigee")

# for ref in refined_sols_perigees
#     θsfs_refined, _, perigees_refined = ref     # unpack
#     plot!(ptrend, θsfs_refined*180/π, perigees_refined*param3b.lstar,
#         marker=:circle, color=:deeppink, label="Refined")
# end

# for solved in solved_sols_perigees
#     θsfs_solved, _, perigees_solved = solved     # unpack
#     plot!(ptrend, [θsfs_solved*180/π], [perigees_solved*param3b.lstar],
#         marker=:circle, color=:lime, label="Solved")
# end

# hline!(ptrend, [target_perigee*param3b.lstar,], color=:red, label="Target Perigee")

# # display plots 
# display(plot(ptraj, ptrend; size=(1400,600)))


# 6. Process results ------------------------------------------------------- #
# The solved solutions are in the `solved_sols_perigees` list; the second entry
# in each entry of the list is the ODESolution struct.
