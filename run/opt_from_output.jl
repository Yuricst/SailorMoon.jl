using ForwardDiff
using Suppressor
using CSV
using DataFrames
using LinearAlgebra
using Printf

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

include("../src/SailorMoon.jl")


## === INPUTS ==================================================
# csv file to load the initial solution
filename = "data/diffcorr_0619_EMrotThrust2.csv"
dir_func = SailorMoon.dv_EMrotdir_sb1frame
output_fname = "data/opt_0619_EMrotThrust.csv"
optim_solver = "snopt"
## =============================================================

# 3body parameter
param3b = SailorMoon.dynamics_parameters()

# multiple shooting parameter
paramMulti = SailorMoon.multi_shoot_parameters(param3b)

# run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 50,   # 1500 ~ 2500
    "tol" => 1e-2,
    "constr_viol_tol" => 1e-4,
    "dual_inf_tol" => 1e-1,
    "output_file" => "ELET_ipopt.out",
    "mu_strategy" => "adaptive",
    "acceptable_constr_viol_tol" => 1e-4
)

sn_options = Dict(
    "Major feasibility tolerance" => 1.e-6,
    "Major optimality tolerance"  => 1.e-4,
    "Minor feasibility tolerance" => 1.e-6,
    "Major iterations limit" => 400,
    "Major print level" => 1,
    "Major step limit" => 0.002,   # 0.1 - 0.01? # default; 2
    "central difference interval" => 1e-6
    # "printfile" => "snopt_print.out",
)

# arc design (1 or 2 or 3)
arc_design = 2

########################################## 

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess
df = CSV.read(filename, DataFrame; header=0);

# maybe want to use "for row in eachrow(df)" to automate the process...? 
row = df[1,:]

x0, lx, ux = SailorMoon.make_ig_bounds2_raw(row, τ_ig, paramMulti.n_arc, true)
println("x_lr: ",  [round(el,digits=3) for el in x0[1 : 9]])
println("x_mid: ", [round(el,digits=3) for el in x0[10+6*paramMulti.n_arc : 18+6*paramMulti.n_arc]])
println("x_lpo: ", [round(el,digits=3) for el in x0[19+12*paramMulti.n_arc : 22+12*paramMulti.n_arc]])

# obtain m_LEO 
_, sol_list, _, tofs = SailorMoon.multishoot_trajectory2(x0, dir_func, paramMulti, true, false, true)
sol, _, _  = sol_list[1]
m_leo = sol[end, end]
println("m_leo; ", round(m_leo, digits=4))

# fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness5_minToF_fixmleo(dir_func, paramMulti, x0, m_leo)
fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness4_minmleo_fixToF(dir_func, paramMulti, x0, sum(tofs), true)

# checking if the initial guess is good enough
res = eval_sft(x0)
# println("ub - x0: ", ux - x0)
# println("x0 - lb: ", x0 - lx)
# println("ub - lb; ", ux-lx)
println("residual (needs to be 0): ", [round(el,digits=3) for el in res])

if optim_solver == "ipopt"
    xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
    options=ip_options, outputfile=true,
    )  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...
elseif optim_solver == "snopt"
    xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
        lx=lx, ux=ux, lg=lg, ug=ug, solver="snopt",
        options=sn_options, outputfile=true, lencw=100000, iSumm=6
    )  #  derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...

else 
    error("optim_solver needs to be ipopt or snopt")
end


fixed_tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]
vec = vcat(1, fixed_tof, xopt[7], xopt)  # tof, m_leo
CSV.write(output_fname,  Tables.table(transpose(vec)), writeheader=false, append=true)

println(Info)
println(x0)
# println("Now, using the initial guess, we reoptimize...")

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

