"""
Test minimize function
"""

using ForwardDiff
using Suppressor

push!(LOAD_PATH,"../../joptimise/src/")
using joptimise


function rosenbrock(x)
    f = (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
    return f
end


function rosenbrock!(g, x)
    # compute objective
    f = rosenbrock(x)
    # constraint
    g[1] = x[1]^2 + x[2]^2 - 1.0
    return f
end

# initial guess
x0 = [4.0; 4.0]
# bounds on variables
lx = [-5.0; -5.0]
ux = [5.0; 5.0]
# bounds on constriants
lg = [0.0]
ug = [0.0]
# number of constraints
ng = 1


## run minimizer with IPOPT
ip_options = Dict(
    "max_iter" => 2500,   # 1500 ~ 2500
    "tol" => 1e-6,
    "output_file" => "test_hogehoge_ipopt.out"
)

xopt, fopt, info = joptimise.minimize(rosenbrock!, x0, ng;
    lx=lx, ux=ux, lg=lg, ug=ug, solver="ipopt",
    options=ip_options, outputfile=true);

println("Done with IPOPT!")
println(info)
println(xopt)


## run minimizer with SNOPT
# sn_options = Dict(
#     "Major feasibility tolerance" => 1.e-6,
#     "Major optimality tolerance"  => 1.e-6,
#     "Minor feasibility tolerance" => 1.e-6,
#     "Major iterations limit" => 1000,
#     "Major print level" => 1,
#     "printfile" => "snopt_print.out",
# )

# # initial guess
# output = @capture_out begin
# x0 = [4.0; 4.0]
# xopt, fopt, info = joptimise.minimize(rosenbrock!, x0, ng; lx=lx, ux=ux, lg=lg, ug=ug, solver="snopt", options=sn_options);
# end

# println("Done with SNOPT!")
# println(info)
# println(xopt)
