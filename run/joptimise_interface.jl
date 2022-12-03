using ForwardDiff
using Suppressor
using CSV
using DataFrames

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

include("../src/SailorMoon.jl")

### PARAMETERS ###################################
# csv file to load the initial solution
filename = "grid_search1129.csv"
# dv_dir function corresponding to the csv file 
dir_func = SailorMoon.dv_no_thrust 

n_arc = 5

##################################################

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess
df = DataFrame(CSV.File(filename))

# maybe want to use "for row in eachrow(df)" to automate the process...? 
row = df[1,:]

sv_mid_cart = [row.x_ra, row.y_ra, row.z_ra, row.xdot_ra, row.ydot_ra, row.zdot_ra]
# change the coordinates into cylindrical (only position)
svm_mid = vcat(SailorMoon.cart2cylind_only_pos(sv_mid_cart), row.m_ra)

tof_leo2mid = row.dt2
tof_mid2lpo = row.dt1
rp = row.rp_kep
ra = row.ra_kep
α = row.alpha
m_rp = row.m_rp

θsf = row.thetaf
ϕ0  = row.phi0

# x_LEO = [ra, rp, α, m0, tof, controls...]
ig_x_LEO = vcat(
    [rp, ra,  α, 1.0, tof_leo2mid/2],
    vcat([[τ_ig,0,0] for i = 1:n_arc]...)
)

ig_x_mid = vcat(
    sv_mid, tof_leo2mid/2, tof_mid2lpo/2, 
    vcat([[τ_ig,0,0] for i = 1:2n_arc]...)
)

# x_LPO = [θf, ϕ, mf, tof, controls...]
ig_x_LPO = vcat(
    [θsf, ϕ0, m_rp, tof_mid2lpo/2],
    vcat([[τ_ig,0,0] for i = 1:n_arc]...)
)

x0 = vcat(ig_x_LEO, ig_x_mid, ig_x_LPO)
###########################################################
# lb, ub of variables 
lx_leo = vcat(
    [200/para.m3b.lstar, 0.8*ra,  0.7*α, 1.2*m_rp, 0.7*tof_leo2mid/2],
    vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
)

ux_leo = vcat(
    [250/param3b.lstar, 1.2*ra,  1.3*α, 0.7*m_rp, 1.3*tof_leo2mid/2],
    vcat([[1.0,pi,pi] for i = 1:n_arc]...)
)


lx_mid = vcat(
    0.9*svm_mid[1], (svm_mid[2]-pi/12), -inf, -inf, -inf, -inf, 1.0, 
    0.7*tof_leo2mid/2, 0.7*tof_mid2lpo/2, 
    vcat([[0.0,-pi,-pi] for i = 1:2n_arc]...)
)
lx_mid = vcat(
    1.1*svm_mid[1], (svm_mid[2]-pi/12), inf, inf, inf, inf, 1.2*m_rp,
    1.3*tof_leo2mid/2, 1.3*tof_mid2lpo/2, 
    vcat([[1.0,0,0] for i = 1:2n_arc]...)
)

lx_lpo = vcat(
    [-pi, -pi, 1.0, 0.7*tof_mid2lpo/2],
    vcat([[0.0,-pi,-pi] for i = 1:n_arc]...)
)
ux_lpo = vcat(
    [pi, pi, 1.0, 1.3*tof_mid2lpo/2],
    vcat([[1.0,pi,pi] for i = 1:n_arc]...)
)

lx = vcat(lx_leo, lx_mid, lx_lpo)
ux = vcat(ux_leo, ux_mid, ux_lpo)

lg = []
ug = []
# number of constraints
ng = 1


## run minimizer with IPOPT
# ip_options = Dict(
#     "max_iter" => 2500,   # 1500 ~ 2500
#     "tol" => 1e-6,
#     "output_file" => "test_hogehoge_ipopt.out"
# )

# xopt, fopt, info = joptimise.minimize(rosenbrock!, x0, ng;
#     lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
#     options=ip_options, outputfile=true);

# println("Done with IPOPT!")
# println(info)
# println(xopt)


res = SailorMoon.multishoot_trajectory(ig_x, dir_func, n_arc, false, false) 

println("res: ", res)
