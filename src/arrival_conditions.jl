"""
Arrival conditions handling
"""


abstract type AbstractArrivalCondition end

"""
Regular SFT arrival condition to planet position with >=0 v-infinity
"""
struct ACHyperbolic <: AbstractArrivalCondition
    bctype::String   # FIXME -- remove later!
    vinf_arrival_max::Float64
end


"""
Arrival to R3BP Libration point
"""
struct ACL1 <: AbstractArrivalCondition
    bctype::String
    lp::Int   # 1 ~ 5
    target_r3bp::Array
    CR3BP_param::Any
end

struct ACL2 <: AbstractArrivalCondition
    bctype::String
    lp::Int   # 1 ~ 5
    target_r3bp::Array
    CR3BP_param::Any
end

struct ACL3 <: AbstractArrivalCondition
    bctype::String
    lp::Int   # 1 ~ 5
    target_r3bp::Array
    CR3BP_param::Any
end

struct ACL4 <: AbstractArrivalCondition
    bctype::String
    lp::Int   # 1 ~ 5
    target_r3bp::Array
    CR3BP_param::Any
end

struct ACL5 <: AbstractArrivalCondition
    bctype::String
    lp::Int   # 1 ~ 5
    target_r3bp::Array
    CR3BP_param::Any
end

ACLagrangePoint = Union{ACL1,ACL2,ACL3,ACL4,ACL5}


"""
Arriving to single LPO's manifold's poincare section
"""
struct AClpo1 <: AbstractArrivalCondition
    bctype::String
    manifold_interp_func::Any  # callable
    t_manifold::Float64
    target_lpo::Dict
    CR3BP_param::Any
end


"""
Arriving to single LPO at varying location along its manifold
"""
struct AClpo2 <: AbstractArrivalCondition
    bctype::String
    sim_manifold::EnsembleSolution
    t_manifold_min::Float64
    t_manifold_max::Float64
    target_lpo::Dict
    CR3BP_param::Any
end


# """
# Arriving to single LPO, with varying manifold propagation time
# """
# struct AClpo3 <: AbstractArrivalCondition
#     bctype::String
#     prob_lpo_stm::ODEProblem  # ODE problem with STM and initial condition of LPO
#     prob_branch::ODEProblem   # ODE Problem without STM, for branch
#     y0::Array   # stable eigenvector
#     ϵ::Float64  # perturbation scalar (with sign)
#     period::Float64
#     t_manifold_min::Float64
#     t_manifold_max::Float64
# end


ThreeBodyCapture = Union{ACL1,ACL2,ACL3,ACL4,ACL5,AClpo1,AClpo2}


"""
    construct_arrival_condition(bctype::String, arrival_ID::Int, vinf_arrival_max::Float64; kwargs...)

Construct arrival condition based on given parameters
"""
function construct_arrival_condition(
    bctype::String,
    arrival_ID::Int,
    vinf_arrival_max::Float64;
    kwargs...,
)
    if bctype == "hyperbolic"
        ArrivalCondition = ACHyperbolic(bctype, vinf_arrival_max)

    elseif bctype == "l1"
        CR3BP_param = get_cr3bp_param(10, arrival_ID)
        lps = lagrangePoints(CR3BP_param.mu)
        target_r3bp = lps[1, :]
        ArrivalCondition = ACL1(bctype, 1, target_r3bp, CR3BP_param)

    elseif bctype == "l2"
        CR3BP_param = get_cr3bp_param(10, arrival_ID)
        lps = lagrangePoints(CR3BP_param.mu)
        target_r3bp = lps[2, :]
        ArrivalCondition = ACL2(bctype, 2, target_r3bp, CR3BP_param)

    elseif bctype == "lpo1"
        CR3BP_param = get_cr3bp_param(10, arrival_ID)
        ArrivalCondition = AClpo1(
            bctype,
            Dict(kwargs)[:manifold_interp_func],
            Dict(kwargs)[:t_manifold],
            Dict(kwargs)[:target_lpo],
            CR3BP_param,
        )

    elseif bctype == "lpo2"
        CR3BP_param = get_cr3bp_param(10, arrival_ID)
        ArrivalCondition = AClpo2(
            bctype,
            Dict(kwargs)[:sim_manifold],
            Dict(kwargs)[:t_manifold_min],
            Dict(kwargs)[:t_manifold_max],
            Dict(kwargs)[:target_lpo],
            CR3BP_param,
        )

    end
    return ArrivalCondition
end


"""
    construct_arrival_condition(bctype::String, arrival_ID::CartesianReferenceNode, vinf_arrival_max::Float64; kwargs...)

Construct arrival condition based on given parameters
"""
function construct_arrival_condition(
    bctype::String,
    arrival_ID::CartesianReferenceNode,
    vinf_arrival_max::Float64;
    kwargs...,
)
    if bctype == "hyperbolic"
        ArrivalCondition = ACHyperbolic(bctype, vinf_arrival_max)
    else
        println("Error: if phase_bodies are CartesianReferenceNode, bctype must be hyperbolic!")
    end
    return ArrivalCondition
end


"""
    serialize_arrival_condition(ArrivalCondition::AbstractArrivalCondition)

Serialize arrival condition from object to dictionary
"""
function serialize_arrival_condition(ArrivalCondition::AbstractArrivalCondition)
    ArrivalCondition_Dict = struct2dict(ArrivalCondition)  # convert struct to dict
    delete!(ArrivalCondition_Dict, "sim_manifold")
    # convert CR3BP-parameter into dictionary
    #ArrivalCondition_Dict["CR3BP_param"] = struct2dict(ArrivalCondition_Dict["CR3BP_param"])
    return ArrivalCondition_Dict
end


"""
    construct_arrival_condition(ACDict::Dict, arrival_ID::Int, center_body_ID::Int=10)

Construct arrival condition from dictionary

# Arguments
    - `ACDict::Dict`: Arrival condition dictionary
    - `arrival_ID::Int`: arrival body ID, i.e. final control node body
    - `center_body_ID::Int`: center body ID, default is 10 (Solar System Barycenter)

# Returns
    - `AbstractArrivalCondition`: arrival condition object
"""
function construct_arrival_condition(
    ACDict::Dict,
    arrival_ID::Int,
    center_body_ID::Int = 10,
)
    if cmp(ACDict["bctype"], "hyperbolic") == 0
        ArrivalCondition = ACHyperbolic(ACDict["bctype"], ACDict["vinf_arrival_max"])

    elseif cmp(ACDict["bctype"], "l1") == 0
        # Get CR3BP parameters
        CR3BP_param = get_cr3bp_param(center_body_ID, arrival_ID)
        # lps = lagrangePoints(CR3BP_param.mu)
        # target_r3bp = lps[2,:]
        # construct arrival condition object
        ArrivalCondition = ACL1(ACDict["bctype"], 1, ACDict["target_r3bp"], CR3BP_param)

    elseif cmp(ACDict["bctype"], "l2") == 0
        # Get CR3BP parameters
        CR3BP_param = get_cr3bp_param(center_body_ID, arrival_ID)
        # construct arrival condition object
        ArrivalCondition = ACL2(ACDict["bctype"], 2, ACDict["target_r3bp"], CR3BP_param)

    elseif cmp(ACDict["bctype"], "lpo1") == 0
        # Get CR3BP parameters
        CR3BP_param = get_cr3bp_param(center_body_ID, arrival_ID)
        # get manifold
        target_lpo = ACDict["target_lpo"]
        target_lpo["x0"] = [el for el in target_lpo["x0"]]  # overwrite to list of float

        _, _, manifold_interp_func = GALT.get_manifold_interp_func(
            arrival_ID,
            target_lpo["x0"],
            target_lpo["period"],
            ACDict["t_manifold"],
            ϵ = target_lpo["ϵ"],
            n = target_lpo["num_branch"],
        )

        # construct arrival condition object
        ArrivalCondition = AClpo1(
            ACDict["bctype"],
            manifold_interp_func,
            ACDict["t_manifold"],
            ACDict["target_lpo"],
            CR3BP_param,
        )

    elseif cmp(ACDict["bctype"], "lpo2") == 0
        # Get CR3BP parameters
        CR3BP_param = get_cr3bp_param(center_body_ID, arrival_ID)
        # get manifold
        target_lpo = ACDict["target_lpo"]
        target_lpo["x0"] = [el for el in target_lpo["x0"]]  # overwrite to list of float
        sim_manifold, _ = R3BP.get_manifold(
            target_lpo["mu"],
            target_lpo["x0"],
            target_lpo["period"],
            target_lpo["tprop_cr3bp"],
            true,
            n = target_lpo["num_branch"],
            ϵ = target_lpo["ϵ"],
            xdir = target_lpo["xdir"],
            reltol = target_lpo["reltol"],
            abstol = target_lpo["abstol"],
            method = Tsit5(),
        )
        # construct arrival condition object
        ArrivalCondition = AClpo2(
            ACDict["bctype"],
            sim_manifold,
            ACDict["t_manifold_min"],
            ACDict["t_manifold_max"],
            ACDict["target_lpo"],
            CR3BP_param,
        )
    end
    return ArrivalCondition
end


"""
    construct_arrival_condition(ACDict::Dict, arrival_ID::Int, center_body_ID::Int=10)

Construct arrival condition from dictionary

# Arguments
    - `ACDict::Dict`: Arrival condition dictionary
    - `arrival_ID::Int`: arrival body ID, i.e. final control node body
    - `center_body_ID::Int`: center body ID, default is 10 (Solar System Barycenter)

# Returns
    - `AbstractArrivalCondition`: arrival condition object
"""
function construct_arrival_condition(
    ACDict::Dict,
    arrival_ID::CartesianReferenceNode,
    center_body_ID::Int = 10,
)
    ArrivalCondition = ACHyperbolic(ACDict["bctype"], ACDict["vinf_arrival_max"])
    return ArrivalCondition
end
