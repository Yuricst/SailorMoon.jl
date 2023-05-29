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
filename = "data/grid_search_Tsit5_0525_EMrotThrust.csv"
# dv_dir function corresponding to the csv file 
dir_func = SailorMoon.dv_EMrotdir_sb1frame 

# 3body parameter
param3b = SailorMoon.dynamics_parameters()

# multiple shooting parameter
paramMulti = SailorMoon.multi_shoot_parameters(param3b)

# run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 100,   # 1500 ~ 2500
    "tol" => 1e-2,
    "constr_viol_tol" => 1e-4,
    "dual_inf_tol" => 1e-1,
    "output_file" => "ELET_ipopt.out",
    "mu_strategy" => "adaptive",
    "acceptable_constr_viol_tol" => 1e-4
)

# arc design (1 or 2 or 3)
arc_design = 2

output_fname = "output_opt_from_gs_0529.csv"

### PARAMETERS #################################

if dir_func == SailorMoon.dv_no_thrust
    τ_ig = 0.0
else 
    τ_ig = 1.0
end

# load initial guess
df = DataFrame(CSV.File(filename))

# for (m, row) in enumerate( eachrow( df ) ) 

    # maybe want to use "for row in eachrow(df)" to automate the process...? 
    row = df[10,:]

    x0, lx, ux = SailorMoon.make_ig_bounds2(row, τ_ig, paramMulti.n_arc)

    fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2(dir_func, paramMulti, x0)

    # println("initial tof: ", [x0[8], x0[9],  x0[17+6n_arc], x0[18+6n_arc], x0[22+12n_arc]])

    # checking if the initial guess is good enough
    res = eval_sft(x0)
    # println("ub - x0: ", ux - x0)
    # println("x0 - lb: ", x0 - lx)
    # println("ub - lb; ", ux-lx)
    # println("x0: ", x0)
    println("residual (needs to be 0): ", res)

    # inital guess
    xopt, fopt, Info = joptimise.minimize(fitness!, x0, ng;
        lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
        options=ip_options, outputfile=true, 
    )  # derivatives=joptimise.UserDeriv());  # to use AD, need this additional parameter...

    tof = xopt[8] + xopt[9] + xopt[17+6*paramMulti.n_arc] + xopt[18+6*paramMulti.n_arc] + xopt[22+12*paramMulti.n_arc]
    vec = vcat(tof, xopt[7], xopt)
    CSV.write(output_fname,  Tables.table(transpose(vec)), writeheader=false, append=true)

    println(Info)
    println("Now, using the initial guess, we reoptimize...")


    # first, obtain the terminal mass
    _, sol_param_list, _, _ = SailorMoon.multishoot_trajectory2(xopt, dir_func, paramMulti, true, false) 
    svf_lr_bck = sol_param_list[1]  # lr -> LEO arc 
    sol, _, _ = svf_lr_bck
    fixed_mf = sol.u[end][7]
    fitness!, ng, lg, ug, eval_sft = SailorMoon.get_fitness2_fixmf_mintof(dir_func, paramMulti, x0, fixed_mf)

    xopt, fopt, Info = joptimise.minimize(fitness!, xopt, ng;
        lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
        options=ip_options, outputfile=false, 
    ) 

    if Info == :Solve_Succeeded
        vec = vcat(fixed_tof, xopt[7], xopt)
        CSV.write(output_fname,  Tables.table(transpose(vec)), writeheader=false, append=true)
    else
        println("Optimization couldn't succeed. Terminated... ")
        break
    end


# end 