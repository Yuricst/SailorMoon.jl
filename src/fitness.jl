"""
Generate fitness function

# Arguments
    - `dir_func`: reference direction of thrust 
    - `param_multi`: parameters for multiple shooting 
"""

"""
    Differential correction (Objective is constant). 
    free ToF and m_LEO
"""
function get_fitness2(
    dir_func,
    param_multi::multishoot_params,
    x
)
    # number of constraints: 7 states (pos,vel,mass) * 2 + 2 (rp & tangential departure)
    ng = 16

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory2(x, dir_func, param_multi, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 22 + 15*param_multi.n_arc  # number of decision variables 

    # storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    # df_onehot = zeros(nx)
    # df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = 1.0
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end

"""
    minimization of ToF.
    free ToF and m_LEO
"""
function get_fitness2_minToF(
    dir_func,
    param_multi::multishoot_params,
    x
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    ng = 16

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory2(x, dir_func, param_multi, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 22 + 15*param_multi.n_arc  # number of decision variables 

    # storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    # df_onehot = zeros(nx)
    # df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = x[8] + x[9] + x[17+6*param_multi.n_arc] + x[18+6*param_multi.n_arc] + x[22+12*param_multi.n_arc]    # tof
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end

"""
    fix tof, differential correction
"""
function get_fitness2_fixToF(
    dir_func,
    param_multi::multishoot_params,
    x, 
    tof
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    ng = 17

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory4(x, dir_func, param_multi, tof, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 22 + 15*param_multi.n_arc  # number of decision variables 
    storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    df_onehot = zeros(nx)
    df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = 1.0  #  x[8] + x[9] + x[17+6*param_multi.n_arc] + x[18+6*param_multi.n_arc] + x[22+12*param_multi.n_arc]    # tof
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end

"""
    fix m_LEO, minimize ToF
"""
function get_fitness2_fixmf_mintof(
    dir_func,
    param_multi::multishoot_params,
    x, 
    mf
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    ng = 17

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory5(x, dir_func, param_multi, mf, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 22 + 15*param_multi.n_arc  # number of decision variables 
    storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    df_onehot = zeros(nx)
    df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = x[8] + x[9] + x[17+6*param_multi.n_arc] + x[18+6*param_multi.n_arc] + x[22+12*param_multi.n_arc]  
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end




## ---- not used as of now --------------------------

function get_fitness(
    dir_func, 
    param_multi::multishoot_params,
    x
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    ng = 14

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory(x, dir_func, param_multi, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 18 + 12*param_multi.n_arc  # number of decision variables 
    storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    df_onehot = zeros(nx)
    df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = 1 #- x[17 + 9*n]    # currently thinking of maximization of mf, given m0 = 1.0
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end


function get_fitness3(
    dir_func,
    param_multi::multishoot_params,
    x
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    ng = 21

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory3(x, dir_func, param_multi, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 26 + 15*param_multi.n_arc   # number of decision variables 
    storage_ad = DiffResults.JacobianResult(x)  # initialize storage for AD
    df_onehot = zeros(nx)
    df_onehot[2] = 1.0   # insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # # create objective function
    # fitness! = function (g, df, dg, x::AbstractVector{T}) where T
    #     # evaluate objective & objective gradient (trivial)
    #     f = x[2]       # whichever x corresponds to e.g. mass at LEO
    #     df[1:nx] = df_onehot[:]
    #     # evalue constraint & constraint gradient
    #     ForwardDiff.jacobian!(storage_ad, eval_sft, x)
    #     g[:] = storage_ad.value
    #     # constraints gradients jacobian
    #     dg[:,:] = DiffResults.jacobian(storage_ad)
    #     return f
    # end

    # create objective function
    fitness! = function (g, x::AbstractVector{T}) where T
        # evaluate objective & objective gradient (trivial)
        f = 1  # - x[17 + 9*n]    # currently thinking of maximization of mf, given m0 = 1.0
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        # println("g(res): ", round.(g[:], digits=3))
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug, eval_sft
end