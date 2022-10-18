"""
SFT objective, setting, propagator struct
"""


## Struct for objective of optimization
abstract type AbstractObjectiveType end

"""
Final mass mazimization problem
"""
struct MassOptimal <: AbstractObjectiveType
    objective_str::String
end

"""
Total time of flight minimization problem
"""
struct TofOptimal <: AbstractObjectiveType
    objective_str::String
end

"""
Lateness minimization problem (w.r.t. reference arrival epoch)
"""
struct LatenessOptimal <: AbstractObjectiveType
    objective_str::String
    nominal_arrival::Float64
end

"""
Initial mass minimization problem
"""
struct WetMassOptimal <: AbstractObjectiveType
    objective_str::String
end

"""
Get objective struct from string
"""
function construct_objective_type(objective_str::String; kwargs...)
    if cmp(objective_str, "mass") == 0 || cmp(objective_str, "mf") == 0
        Objective = MassOptimal(objective_str)

    elseif cmp(objective_str, "tof") == 0
        Objective = TofOptimal(objective_str)

    elseif cmp(objective_str, "lateness") == 0
        Objective = LatenessOptimal(objective_str, Dict(kwargs)[:nominal_arrival])

    elseif cmp(objective_str, "wetmass") == 0 || cmp(objective_str, "m0") == 0
        Objective = WetMassOptimal(objective_str)

    end
    return Objective
end


"""
Serialize objective struct
"""
function serialize_obejctive(Objective::AbstractObjectiveType)
    return struct2dict(Objective)  # convert struct to dict
end


"""
Reconstruct objective object from dictionary
"""
function construct_objective_type(Objective_dict::Dict)
    if cmp(Objective_dict["objective_str"], "mass") == 0 || cmp(Objective_dict["objective_str"], "mf") == 0
        Objective = MassOptimal(Objective_dict["objective_str"])

    elseif cmp(Objective_dict["objective_str"], "tof") == 0
        Objective = TofOptimal(Objective_dict["objective_str"])

    elseif cmp(Objective_dict["objective_str"], "lateness") == 0
        Objective = LatenessOptimal(
            Objective_dict["objective_str"],
            Objective_dict["nominal_arrival"],
        )

    elseif cmp(Objective_dict["objective_str"], "wetmass") == 0 || cmp(Objective_dict["objective_str"], "m0") == 0
        Objective = WetMassOptimal(Objective_dict["objective_str"])

    end
    return Objective
end



## Struct to hold main options of Sims-Flanagan problem
abstract type SimsFlanaganSettings end


"""
Structure holds Sims-Flanagan transcription problem settings with constant number of segments
"""
struct SimsFlanaganConstantSegments <: SimsFlanaganSettings
    mu::Float64
    n::Int
    np::Int
    phase_bodies::Array
    etmin::Float64
    etmax::Float64
    tof_total_max::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    process_dv_fun::Function
    mmin::Float64
    et_scale::Float64
    et0_interp::Float64
    et1_interp::Float64
    dt_interp::Float64
    control_node_type::String
    verbosity::Int
    zero_soi::Bool
    center_body::Int
    cmp_type::String
    cmp_weights::Vector{Float64}
    h_safe_factor::Float64
    use_launch_vehicle::Bool
end


"""
Structure holds Sims-Flanagan transcription problem settings with variable number of segments
"""
struct SimsFlanaganVariableSegments <: SimsFlanaganSettings
    mu::Float64
    n::Vector{Int}
    np::Int
    phase_bodies::Array
    etmin::Float64
    etmax::Float64
    tof_total_max::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    process_dv_fun::Function
    mmin::Float64
    et_scale::Float64
    et0_interp::Float64
    et1_interp::Float64
    dt_interp::Float64
    control_node_type::String
    verbosity::Int
    zero_soi::Bool
    center_body::Int
    cmp_type::String
    cmp_weights::Vector{Float64}
    h_safe_factor::Float64
    use_launch_vehicle::Bool
end


"""
Structure holds Chebyshev-CBCBC formulation settings
"""
struct SimsFlanaganChebyCBCBC <: SimsFlanaganSettings
    mu::Float64
    np::Int
    phase_bodies::Array
    etmin::Float64
    etmax::Float64
    tof_total_max::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    process_dv_fun::Function
    mmin::Float64
    et_scale::Float64
    et0_interp::Float64
    et1_interp::Float64
    dt_interp::Float64
    control_node_type::String
    verbosity::Int
    zero_soi::Bool
    center_body::Int
    cmp_type::String
    cmp_weights::Vector{Float64}
    h_safe_factor::Float64
    use_launch_vehicle::Bool
    # NEW THINGS!
    t_spacing::Float64
    n_coefficients::Int
    max_forward_bias::Float64
end


"""
Structure holds Chebyshev-CBarbitrary formulation settings (CBCB...BC).
`nb` holds the number of bang arcs, which must be an odd number.
"""
struct SimsFlanaganChebyOddBang <: SimsFlanaganSettings
    mu::Float64
    np::Int
    phase_bodies::Array
    etmin::Float64
    etmax::Float64
    tof_total_max::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    process_dv_fun::Function
    mmin::Float64
    et_scale::Float64
    et0_interp::Float64
    et1_interp::Float64
    dt_interp::Float64
    control_node_type::String
    verbosity::Int
    zero_soi::Bool
    center_body::Int
    cmp_type::String
    cmp_weights::Vector{Float64}
    h_safe_factor::Float64
    use_launch_vehicle::Bool
    # NEW THINGS!
    t_spacing::Float64
    n_coefficients::Int
    nb::Int  # number of bang-arcs, odd number
end


"""
Structure holds Chebyshev-CBarbitrary formulation settings (CBCB...BC).
`nb` holds the number of bang arcs, which must be an even number.
"""
struct SimsFlanaganChebyEvenBang <: SimsFlanaganSettings
    mu::Float64
    np::Int
    phase_bodies::Array
    etmin::Float64
    etmax::Float64
    tof_total_max::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    process_dv_fun::Function
    mmin::Float64
    et_scale::Float64
    et0_interp::Float64
    et1_interp::Float64
    dt_interp::Float64
    control_node_type::String
    verbosity::Int
    zero_soi::Bool
    center_body::Int
    cmp_type::String
    cmp_weights::Vector{Float64}
    h_safe_factor::Float64
    use_launch_vehicle::Bool
    # NEW THINGS!
    t_spacing::Float64
    n_coefficients::Int
    nb::Int  # number of bang-arcs, even number
end

# union of chebyshev formulation
SimsFlanaganCheby = Union{
    SimsFlanaganChebyCBCBC,
    SimsFlanaganChebyOddBang,
    SimsFlanaganChebyEvenBang
}


"""
    serialize_SimsFlanaganSettings(SFset::SimsFlanaganSettings)

Serialize Sims-Flanagan settings to dictionary
"""
function serialize_SimsFlanaganSettings(SFset::SimsFlanaganSettings)
    SFset_dict = struct2dict(SFset)  # convert struct to dict
    SFset_dict["process_dv_fun"] = string(SFset_dict["process_dv_fun"])
    return SFset_dict
end


"""
    construct_SimsFlanaganSettings(SFset_dict::Dict)

Construct Sims-Flanagan settings object
"""
function construct_SimsFlanaganSettings(SFset_dict::Dict)
    # fixes for backward compatibility
    if SFset_dict["process_dv_fun"] == "dv_lvlh2inertial"
        SFset_dict["process_dv_fun"] = dv_lvlh2inertial
    end
    if haskey(SFset_dict, "cmp_weights") == false # backward compatibility
        SFset_dict["cmp_weights"] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    end
    if haskey(SFset_dict, "h_safe_factor") == false # backward compatibility
        SFset_dict["h_safe_factor"] = 0.1
    end
    if haskey(SFset_dict, "use_launch_vehicle") == false # backward compatibility
        SFset_dict["use_launch_vehicle"] = false
    end
    if haskey(SFset_dict, "n") == true
        if typeof(SFset_dict["n"]) == Int64
            return tostruct(SimsFlanaganConstantSegments, SFset_dict)
        else
            return tostruct(SimsFlanaganVariableSegments, SFset_dict)
        end
    else
        # backward compatibility
        if haskey(SFset_dict, "n") == false
            SFset_dict["max_forward_bias"] = 0.7
        end
        return tostruct(SimsFlanaganChebyCBCBC, SFset_dict)
    end
end




## Struct for propagator settings
# helper function to interpret solver name from DifferentialEquations
"""
    interpret_solver_name(solver_str)

Helper function to interpret solver name when using ODE-based SFT legs
"""
function interpret_solver_name(solver_str::String)
    if cmp(solver_str, "Tsit5") == 0 || startswith(solver_str, "OrdinaryDiffEq.Tsit5") == true
        solver = Tsit5()
    elseif cmp(solver_str, "lsoda") == 0 || startswith(solver_str, "OrdinaryDiffEq.lsoda") == true
        solver = lsoda()
    elseif cmp(solver_str, "BS3") == 0 ||
        startswith(solver_str, "OrdinaryDiffEq.BS3") == true ||
        startswith(solver_str, "BS3") == true
        solver = BS3()
    elseif cmp(solver_str, "DP5") == 0 || startswith(solver_str, "OrdinaryDiffEq.DP5") == true
        solver = DP5()
    elseif cmp(solver_str, "DP8") == 0 || startswith(solver_str, "OrdinaryDiffEq.DP8") == true
        solver = DP8()
    elseif cmp(solver_str, "Vern7") == 0 || startswith(solver_str, "OrdinaryDiffEq.Vern7") == true
        solver = Vern7()
    elseif cmp(solver_str, "Rodas4") == 0 || startswith(solver_str, "OrdinaryDiffEq.Rodas4") == true
        solver = Rodas4()
    elseif cmp(solver_str, "Trapezoid") == 0 ||
           startswith(solver_str, "OrdinaryDiffEq.Trapezoid") == true
        solver = Trapezoid()
    else
        println("Solver name $solver_str is not compatible!")
    end
    return solver
end

# struct for propagator option
abstract type AbstractPropagatorType end

"""
Propagator class for Kepler's equation method
"""
struct KeplerPropagator <: AbstractPropagatorType
    tol::Float64
    maxiter::Int
end

"""
Propagator class for ODE integration method
"""
struct ODEPropagator <: AbstractPropagatorType
    method::Any
    reltol::Float64
    abstol::Float64
    dt::Any
    prob_base::ODEProblem
    tol::Float64
    maxiter::Int
end


"""
Create Propagator object

# Args:
    - `ArrivalCondition::AbstractArrivalCondition`: arrival condition
    - `propagator::String`: propagator type, "kepler" or "ode"
    - `propagator_opts::Dict`: propagator options dictionary

# Returns
    - `AbstractPropagatorType`: Kepler or ODE Propagator object
"""
function construct_propagator(
    ArrivalCondition::AbstractArrivalCondition,
    propagator::String,
    propagator_opts = nothing,
    chebyshev_problem::Bool=false,
)
    # set default options
    if isnothing(propagator_opts) == true
        propagator_opts = Dict(
            "type" => propagator,
            "tol" => 1.e-14,
            "maxiter" => 20,   # for Kepler's method
            "reltol" => 1.e-12,
            "abstol" => 1.e-12,
            "method" => "Tsit5",    # need to fix lsoda
            "dt" => nothing,        # nothing or a float
        )
    end

    # construct Propagator object
    if cmp(propagator, "kepler") == 0
        tol = propagator_opts["tol"]
        maxiter = propagator_opts["maxiter"]
        Propagator = KeplerPropagator(tol, maxiter)
    else
        # construct prob_base
        if chebyshev_problem == false
            if typeof(ArrivalCondition) == ACHyperbolic
                prob_base = ODEProblem(
                    twobody_fixthrust!,
                    [1.0, 0.0, 0.01, 0.0, 1.0, 0.0],
                    (0.0, 1.0),
                    (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                )
                # elseif cmp(bctype, "l1") == 0 || cmp(bctype, "l2") == 0 || cmp(bctype, "lpo1") == 0 || cmp(bctype, "lpo2") == 0
            else
                prob_base = ODEProblem(
                    twobody_fixthrust_thirdbody!,
                    [1.0, 0.0, 0.01, 0.0, 1.0, 0.0],
                    (0.0, 1.0),
                    (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 432000.0, [], 0.0, 5.0),
                )
                # CR3BP parameters
                #ω_cr3bp = 1.0
                #CR3BP_param = get_cr3bp_param(10, SFset.phase_bodies[end])
            end

        else   # chebyshev problem ODE
            prob_base = ODEProblem(
                twobody_callable_thrust!,
                [1.0, 0.0, 0.01, 0.0, 1.0, 0.0],
                (0.0, 1.0),
                (1.0, 0.0, 0.0, ones(3), ones(3), 1.0),
            )
        end

        method = interpret_solver_name(propagator_opts["method"])
        reltol = propagator_opts["reltol"]
        abstol = propagator_opts["abstol"]
        dt = propagator_opts["dt"]
        Propagator = ODEPropagator(method, reltol, abstol, dt, prob_base,
            1.e-14, 20)
    end
    return Propagator
end


"""
    serialize_propagator(Propagator::AbstractPropagatorType)

Serialize propagator into dictionary
"""
function serialize_propagator(Propagator::AbstractPropagatorType)
    Propagator_dict = struct2dict(Propagator)
    if typeof(Propagator) == KeplerPropagator
        Propagator_dict["type"] = "kepler"
    else
        Propagator_dict["type"] = "ode"
        Propagator_dict["method"] = string(Propagator.method)[1:end-2]
        Propagator_dict["prob_base"] = 0
    end
    return Propagator_dict
end


"""
    construct_propagator(
        Propagator_dict::Dict,
        ArrivalCondition::AbstractArrivalCondition,
        SFset::Union{SimsFlanaganConstantSegments,SimsFlanaganVariableSegments},
    )

Construct propagator from dictionary
"""
function construct_propagator(
    Propagator_dict::Dict,
    ArrivalCondition::AbstractArrivalCondition,
    SFset::Union{SimsFlanaganConstantSegments,SimsFlanaganVariableSegments},
)
    if cmp(Propagator_dict["type"], "kepler") == 0
        Propagator = KeplerPropagator(Propagator_dict["tol"], Propagator_dict["maxiter"])

    elseif cmp(Propagator_dict["type"], "ode") == 0
        # construct prob_base
        if typeof(ArrivalCondition) == ACHyperbolic
            prob_base = ODEProblem(
                twobody_fixthrust!,
                [1.0, 0.0, 0.01, 0.0, 1.0, 0.0],
                (0.0, 1.0),
                (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            )
            # elseif cmp(bctype, "l1") == 0 || cmp(bctype, "l2") == 0 || cmp(bctype, "lpo1") == 0 || cmp(bctype, "lpo2") == 0
        else
            prob_base = ODEProblem(
                twobody_fixthrust_thirdbody!,
                [1.0, 0.0, 0.01, 0.0, 1.0, 0.0],
                (0.0, 1.0),
                (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 432000.0, [], 0.0, 5.0),
            )
            # CR3BP parameters
            #ω_cr3bp = 1.0
            #CR3BP_param = get_cr3bp_param(10, SFset.phase_bodies[end])
        end

        # construct propagator object
        Propagator = ODEPropagator(
            interpret_solver_name(Propagator_dict["method"]),
            Propagator_dict["reltol"],
            Propagator_dict["abstol"],
            Propagator_dict["dt"],
            prob_base,
            1.e-14, 20
        )
    end
    return Propagator
end



"""
    construct_propagator(Propagator_dict::Dict, ArrivalCondition::AbstractArrivalCondition,
        SFset::Union{SimsFlanaganConstantSegments,SimsFlanaganVariableSegments},
    )

Construct propagator from dictionary
"""
function construct_propagator(
    Propagator_dict::Dict,
    ArrivalCondition::AbstractArrivalCondition,
    SFset::SimsFlanaganChebyCBCBC,
)
    if cmp(Propagator_dict["type"], "kepler") == 0
        Propagator = KeplerPropagator(Propagator_dict["tol"], Propagator_dict["maxiter"])

    elseif cmp(Propagator_dict["type"], "ode") == 0
        # construct prob_base
        prob_base = ODEProblem(
            twobody_chebyshev_extension!,
            [1.0, 0.0, 0.01, 0.0, 1.0, 0.0, 0.0, 0.0],  # [r,v,m,θ,β]
            (0.0, 1.0),  # [t0,t1]
            (1.0, 0.0, 0.0, ones(SFset.n_coefficients), ones(SFset.n_coefficients), 1.0),
        )
        if typeof(Propagator_dict["method"]) == String
            method = interpret_solver_name(Propagator_dict["method"])
        else
            method = Propagator_dict["method"]
        end

        # construct propagator object
        Propagator = ODEPropagator(
            method,
            Propagator_dict["reltol"],
            Propagator_dict["abstol"],
            Propagator_dict["dt"],
            prob_base,
            1.e-14, 20
        )
    end
    return Propagator
end
