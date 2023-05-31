using ForwardDiff
using Suppressor
using CSV
using DataFrames
using LinearAlgebra

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

include("../src/SailorMoon.jl")

### INPUTS ###################################
# csv file to load the initial solution
filename = "data/output_0519.csv"
# dv_dir function corresponding to the csv file 
dir_func = SailorMoon.dv_EMrotdir_sb1frame 

# 3body parameter
param3b = SailorMoon.dynamics_parameters()

# multiple shooting parameter
paramMulti = SailorMoon.multi_shoot_parameters(param3b)

# run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "tol" => 1e-2,
    "constr_viol_tol" => 1e-4,
    "dual_inf_tol" => 1e-1,
    "output_file" => "ELET_ipopt.out",
    "mu_strategy" => "adaptive",
    "acceptable_constr_viol_tol" => 1e-4
)

sn_options = Dict(
    "Major feasibility tolerance" => 1.e-6,
    "Major optimality tolerance"  => 1.e-6,
    "Minor feasibility tolerance" => 1.e-6,
    "Major iterations limit" => 1000,
    "Major print level" => 1,
    # "printfile" => "snopt_print.out",
)

# arc design (1 or 2 or 3)
arc_design = 2

output_fname = "data/output_opt_0519.csv"

########################################## 

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess
df = DataFrame(CSV.File(filename))

# maybe want to use "for row in eachrow(df)" to automate the process...? 
row = df[1,:]

x0, lx, ux = SailorMoon.make_ig_bounds2_raw(row, τ_ig, paramMulti.n_arc)

# obtain m_LEO 
_, sol_list, _, _ = SailorMoon.multishoot_trajectory2(x0, dir_func, paramMulti, true, false)
sol, _, _  = sol_list[1]
m_leo = sol[end, end]
println("m_leo; ", m_leo)

fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2_fixmf_mintof(dir_func, paramMulti, x0, m_leo)

# checking if the initial guess is good enough
res = eval_sft(x0)
# println("ub - x0: ", ux - x0)
# println("x0 - lb: ", x0 - lx)
# println("ub - lb; ", ux-lx)
# println("x0: ", x0)
# println("residual (needs to be 0): ", res)


# inital guess
xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug, solver="snopt",
    options=sn_options, outputfile=true, 
)  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...

fixed_tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]
vec = vcat(fixed_tof, xopt[7], xopt)  # tof, m_leo
CSV.write(output_fname,  Tables.table(transpose(vec)), writeheader=false, append=true)


println(Info)
println("Now, using the initial guess, we reoptimize...")

# while true
#     global xopt 

#     # transfer the TOF (fixed value)
#     fixed_tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]

#     # change the value a little bit... 
#     fixed_tof = fixed_tof - 0.05

#     # redefine the equality constriant
#     fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2_fixToF(dir_func, paramMulti, x0, fixed_tof)

#     xopt, fopt, Info = joptimise.minimize(fitness!, xopt, ng;
#         lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
#         options=ip_options, outputfile=false, 
#     ) 

#     if Info == :Solve_Succeeded
#         vec = vcat(fixed_tof, xopt[7], xopt)
#         CSV.write(output_fname,  Tables.table(transpose(vec)), writeheader=false, append=true)
#     else
#         println("Optimization couldn't succeed. Terminated... ")
#         break
#     end
# end

