"""
Generate fitness function

# Arguments
    - `n` : number of sf transcription (discretization per arc)
    - `dir_func`: reference direction of thrust 
"""
function get_fitness(
    n::Int,
    dir_func,
    x
)
    # number of constraints: 7 states (pos,vel,mass) * 2 
    # FIXME ... check if this is correct!
    ng = 14

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res = multishoot_trajectory(x, dir_func, n, false, false) 
        
        # compute constraints
        # residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize (for AD)
        residuals = zeros(ng)
        residuals[:] = res[:]
        return residuals
    end

    nx = 18 + 12*n  # number of decision variables 
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
        f = - x[17 + 9*n]    # currently thinking of maximization of mf, given m0 = 1.0
        g[:] = eval_sft(x)       # constraints (that need to be zero)
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints

    return fitness!, ng, lg, ug
end
