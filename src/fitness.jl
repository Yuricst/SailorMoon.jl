"""
Generate fitness function

# Arguments
    - `n` : number of sf transcription (discretization per arc)
    - `SCParam` : [mdot, tmax]
    - `Propagator` : 
    - `param3b` : 
    - `LPOArrival` : 
"""
function get_fitness(
    n::Int,
    SCparam::Vector{Float64},
    Propagator::ODEPropagator,
    param3b::AbstractParameterType,
    LPOArrival::AbstractTerminalType
)
    # number of constraints FIXME ... check if this is correct!
    ng = 7

    # function that computes constraints of SFT
    eval_sft = function (x::AbstractVector{T}) where T
        # unpack decision vector & residual
        res, _, _, _, _ = sf_propagate(x,n,SCparam,Propagator,param3b,LPOArrival)

        # compute constraints
        residuals = ForwardDiff.Dual[0 for i = 1:ng]   # initialize
        residuals[:] = res[:]
        return residuals
    end

    nx = 8 + 6*n  # number of decision variables (tof(1),c_launch(4),c_arr(3),tau1(3n),tau2(3n))
    storage_ad = DiffResults.JacobianResult(x0)  # initialize storage for AD
    df_onehot = zeros(nx)
    df_onehot[3] = 1.0   # FIXME ... insert 1 to whichever index of x corresponding to e.g. mass at LEO

    # create objective function
    fitness! = function (g, df, dg, x)
        # evaluate objective & objective gradient (trivial)
        f = x[3]       # FIXME ... whichever x corresponds to e.g. mass at LEO
        df[1:nx] = df_onehot[:]
        # evalue constraint & constraint gradient
        ForwardDiff.jacobian!(storage_ad, eval_sft, x)
        g[:] = storage_ad.value
        # constraints gradients jacobian
        dg[:,:] = DiffResults.jacobian(storage_ad)
        return f
    end

    # problem bounds
    lg = [0.0 for idx=1:ng]   # lower bounds on constraints
    ug = [0.0 for idx=1:ng]   # upper bounds on constraints
    lx = [
        1.0,          # lower bounds on decision variables
    ]
    ux = [
        1.0,
    ]                 # upper bounds on decision variables
    return fitness!, ng, lx, ux, lg, ug
end
