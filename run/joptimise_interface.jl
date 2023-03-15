using ForwardDiff
using Suppressor
using CSV
using DataFrames
using LinearAlgebra

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

include("../src/SailorMoon.jl")

### PARAMETERS ###################################
# csv file to load the initial solution
filename = "grid_search.csv"
# dv_dir function corresponding to the csv file 
dir_func = SailorMoon.dv_no_thrust 

# 3body parameter
param3b = SailorMoon.dynamics_parameters()

n_arc = 5

# run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "tol" => 1e-6,
    "output_file" => "ELET_ipopt.out"
)
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

x0, lx, ux = SailorMoon.make_ig_bounds2(row, τ_ig, n_arc)

# res = SailorMoon.multishoot_trajectory(x0, dir_func, n_arc)
# print("size of res: ", size(res))

# use get_fitness or get_fitness2
fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2(n_arc, dir_func, x0)

# checking if the initial guess is good enough
res = eval_sft(x0)
println("ub - x0: ", ux - x0)
println("x0 - lb: ", x0 - lx)
println("x0: ", x0)
println("residual (needs to be 0): ", res)

# make sure the initial guess is inbetween ub & lb
vec = vcat(ux - x0, x0 - lx)
if any(vec .< 0.0)
    error("Error: (At least one element of) initial guess is infinging the defined ub/lb.") 
end

xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
    options=ip_options, outputfile=true, 
)  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...

println(Info)

# println(xopt)
# res = SailorMoon.multishoot_trajectory(ig_x, dir_func, n_arc, false, false) 

# println("res: ", res)
