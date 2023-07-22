# using Distributed

# @everywhere begin

using ForwardDiff
using Suppressor
using CSV
using DataFrames
using LinearAlgebra
using Printf
using Base.Threads

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise

include("../src/SailorMoon.jl")

## === INPUTS ==================================================
# csv file to load the initial solution
filename = "data/opt_0717_EMrotThrust.csv"
dir_func = SailorMoon.dv_EMrotdir_sb1frame
output_fname = "data/opt_0717_EMrotThrust_retry.csv"
optim_solver = "snopt"
## =============================================================

# 3body parameter
param3b = SailorMoon.dynamics_parameters()

# multiple shooting parameter
paramMulti = SailorMoon.multi_shoot_parameters(param3b)


# run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 100,   # 1500 ~ 2500
    "tol" => 1e-4,
    "constr_viol_tol" => 1e-6,
    "dual_inf_tol" => 1e-1,
    "output_file" => "ELET_ipopt.out",
    "mu_strategy" => "adaptive",
    "acceptable_constr_viol_tol" => 1e-4,
    "print_level" => 0,
)

sn_options = Dict(
    "Major feasibility tolerance" => 1.e-6,
    "Major optimality tolerance"  => 1.e-1,
    "Minor feasibility tolerance" => 1.e-6,
    "Major iterations limit" => 50,
    "Major print level" => 1,
    # "Major step limit" => 0.01,   # 0.1 - 0.01? # default; 2,  0.001 ;looks working in general 
    # "linesearch tolerance" => 0.95, 
    "printfile" => "snopt_opt.out",
)


# ===============================================================

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess ( "grid_serach_XXX.csv" )
df = DataFrame(CSV.File(filename))
# df = df[394:395, :]
height = size(df,1)
height = 2
# end

# for (m, row) in enumerate( eachrow( df ) ) 
# @distributed for m in collect(1:height)
# for m in collect(1:height)
row = df[1,:]

# perform the differential correction only if there is no flyby

x0, lx, ux = SailorMoon.make_ig_bounds2_raw(row, τ_ig, paramMulti.n_arc, false, paramMulti.n_arc)
fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2(dir_func, paramMulti, x0)
# fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2_fixToF(dir_func, paramMulti, x0, row.tof)

# checking if the initial guess is good enough
res = eval_sft(x0)
# println("altitude difference: ", res[end-1]*param3b.lstar, " km")
# println("ub - x0: ", ux - x0)
# println("x0 - lb: ", x0 - lx)
# println("ub - lb; ", ux-lx)
# println("x0: ", x0)
println("residual (needs to be 0): ", res)

# _, sol_list, sol_bal, _ = SailorMoon.multishoot_trajectory2(x0, dir_func, paramMulti, true, false)
# for j = 1:length(sol_list)
#     sol, _, name = sol_list[j]
#     mf = sol[end, :]
#     println("m's; ", sol[end, 1], "  ", sol[end, end])
# end

# inital guess

# xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
# lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
# options=ip_options, outputfile=true,
# )  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...

xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug, solver="snopt",
    options=sn_options, outputfile=true, lencw=50000, iSumm=6,
)  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...


if Info == "Finished successfully: optimality conditions satisfied"  || Info == :Solve_Succeeded
    # println("solved! row added for ig #", m, "/", height)

    fixed_tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]
    sol_vec = vcat(row[1], fixed_tof, xopt[7], xopt)  # original id, tof, mass@LEO, xopt
    sol_vec = transpose(sol_vec)

    CSV.write(output_fname,  Tables.table(sol_vec), writeheader=false, append=true)

end 


    # println("optimization #", m)

    #     println("First differential correction is done. Now, using the initial guess, we reoptimize...")
    #     while true
    #         global xopt 

    #         # transfer the TOF (fixed value)
    #         \
            
    #         fixed_tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]

    #         # change the value a little bit... 
    #         fixed_tof = fixed_tof - 0.05

    #         # redefine the equality constriant
    #         fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2_fixToF(dir_func, paramMulti, x0, fixed_tof)

    #         if optim_solver == "ipopt"
    #             xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    #             lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
    #             options=ip_options, outputfile=true,
    #             )  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...
    #         elseif optim_solver == "snopt"
    #             xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
    #                 lx=lx, ux=ux, lg=lg, ug=ug, solver="snopt",
    #                 options=sn_options, outputfile=true, lencw=5000, iSumm=6,
    #             )  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...
    #         else 
    #             error("optim_solver needs to be ipopt or snopt")
    #         end

    #         if Info == :Solve_Succeeded
    #             sol_vec = vcat(fixed_tof, xopt[7], xopt)
    #             CSV.write(output_fname,  Tables.table(transpose(sol_vec)), writeheader=false, append=true)
    #         else
    #             println("Optimization couldn't succeed. Terminated... ")
    #             break
    #         end
    #     end

# end

