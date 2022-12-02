using ForwardDiff
using Suppressor

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

# csv file to load the initial solution
filename = "hogehoge.csv"
# dv_dir function corresponding to the csv file 
dir_func = dv_no_thrust   

# load initial guess
df = DataFrame(CSV.File(filename))

# maybe want to use "for row in eachrow(df)" to automate the process...? 
row = df[1,:]

sv_mid_cart = row.x_ra
# change the coordinates into cylindrical (only position)
sv_mid = cart2cylind_only_pos(sv_mid_car)

tof_leo2mid = row.dt2
tof_mid2lpo = row.dt1
rp = row.rp_kep
ra = row.ra_kep

θsf = row.thetaf
ϕ0  = row.phi0


# create test decision vector
τ_ig = 0.0
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
    [θsf, ϕ0, 1.0, tof_mid2lpo/2],
    vcat([[τ_ig,0,0] for i = 1:n_arc]...)
)

ig_x = vcat(ig_x_LEO, ig_x_mid, ig_x_LPO)

res = multishoot_trajectory(ig_x, dir_func, false, false) 

