"""
Control node functions
"""


"""
    get_controlnode_range(
        et0::Float64,
        et1::Float64,
        bodyID::Int,
        frame::String,
        lstar::Float64,
        star::Float64,
        steps::Int
    )

Wrapper for get_controlnode over range of epoch
"""
function get_controlnode_range(
    et0::Float64,
    et1::Float64,
    bodyID::Int,
    frame::String,
    lstar::Float64,
    vstar::Float64,
    steps::Int,
    center_body::Int = 10,
)
    ets = LinRange(et0, et1, steps)
    svs = zeros(6, steps)
    for i = 1:steps
        sv = spice_spkssb(ets[i], bodyID, frame, lstar, vstar, center_body)
        svs[:, i] = sv
    end
    return svs
end


"""
Get control-node positions for `CartesianReferenceNode`
"""
function get_controlnode_range(
    et0::Float64,
    et1::Float64,
    node::CartesianReferenceNode,
    steps::Int,
    mu::Float64=1.0,
)
    ets = LinRange(et0, et1, steps)
    svs = zeros(6, steps)
    for (i,et) in enumerate(ets)
        t = et - node.et_ref
        svs[:, i] = keplerder_nostm(mu, node.sv_ref, 0.0, t, 1.e-13, 20)
    end
    return svs
end


# """
#     get_controlnode(et::Float64, ctuples)
#
# Get state-vector of body based on Spline1D objects obtained from interp_spkssb()
#
# Args:
#     `et::Float64`: epoch to extract state
#     `ctuple`: tuple of Spline1D objects for all 6 states
#
# Returns:
#     (Array): state evaluated via interpolation at epoch et
# """
# function get_controlnode(et::Float64, ctuples)
#     # construct array
#     sv = [
#         evaluate(ctuples[1], et),
#         evaluate(ctuples[2], et),
#         evaluate(ctuples[3], et),
#         evaluate(ctuples[4], et),
#         evaluate(ctuples[5], et),
#         evaluate(ctuples[6], et),
#     ]
#     return sv
# end


"""
    get_control_nodes(SFset::SimsFlanaganSettings, i::Int, et_fwd::Float64, et_bck::Float64)

Get control nodes directly from SPICE. Provide epochs in canoncical terms.
"""
function get_control_nodes(
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
)
    node_fwd = spice_spkssb(
        et_fwd * SFset.tstar,
        SFset.phase_bodies[i],
        "ECLIPJ2000",
        SFset.lstar,
        SFset.vstar,
        SFset.center_body,
    )
    node_bck = spice_spkssb(
        et_bck * SFset.tstar,
        SFset.phase_bodies[i+1],
        "ECLIPJ2000",
        SFset.lstar,
        SFset.vstar,
        SFset.center_body,
    )
    return node_fwd, node_bck
end


"""
    get_control_nodes(
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        node::Vector{SPICEControlNode},
    )

Get control node for SPICE-based control node
"""
function get_control_nodes(
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    node::Vector{SPICEControlNode},
)
    if length(phase_bodies[i]) == 1
        node_fwd = spice_spkssb(
            et_fwd * SFset.et_scale,
            SFset.phase_bodies[i],
            "ECLIPJ2000",
            SFset.lstar,
            SFset.vstar,
        )
    else
        node_fwd = phase_bodies[i]
    end
    if length(phase_bodies[i+1]) == 1
        node_bck = spice_spkssb(
            et_bck * SFset.et_scale,
            SFset.phase_bodies[i+1],
            "ECLIPJ2000",
            SFset.lstar,
            SFset.vstar,
        )
    else
        node_bck = phase_bodies[i+1]
    end
    return node_fwd, node_bck
end


"""
    get_control_nodes(
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        nodes::Vector{SPICEInterpControlNode},
    )

Get control node for SPICE-interpolated control node.
"""
function get_control_nodes(
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    nodes::Vector{SPICEInterpControlNode},
)
    if length(SFset.phase_bodies[i]) == 1
        node_fwd = [
            evaluate(nodes[i].ctuple[1], et_fwd),
            evaluate(nodes[i].ctuple[2], et_fwd),
            evaluate(nodes[i].ctuple[3], et_fwd),
            evaluate(nodes[i].ctuple[4], et_fwd),
            evaluate(nodes[i].ctuple[5], et_fwd),
            evaluate(nodes[i].ctuple[6], et_fwd),
        ]
    else
        node_fwd = phase_bodies[i]
    end
    if length(SFset.phase_bodies[i+1]) == 1
        node_bck = [
            evaluate(nodes[i+1].ctuple[1], et_bck),
            evaluate(nodes[i+1].ctuple[2], et_bck),
            evaluate(nodes[i+1].ctuple[3], et_bck),
            evaluate(nodes[i+1].ctuple[4], et_bck),
            evaluate(nodes[i+1].ctuple[5], et_bck),
            evaluate(nodes[i+1].ctuple[6], et_bck),
        ]
    else
        node_bck = phase_bodies[i+1]
    end
    return node_fwd, node_bck
end


"""
    get_control_nodes(
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        nodes::Vector{CartesianReferenceNode},
    )

Get control node for SPICE-interpolated control node.
"""
function get_control_nodes(
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    nodes::Vector{CartesianReferenceNode},
)
    # forward node
    t_fwd = et_fwd - nodes[i].et_ref
    node_fwd = keplerder_nostm(SFset.mu, nodes[i].sv_ref, 0.0, t_fwd, 1.e-13, 20)
    # backward node
    t_bck = et_bck - nodes[i+1].et_ref
    node_bck = keplerder_nostm(SFset.mu, nodes[i+1].sv_ref, 0.0, t_bck, 1.e-13, 20)
    return node_fwd, node_bck
end


"""
    function process_controlnodes(
        ArrivalCondition::ACHyperbolic,
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        vinf::Vector{Float64},
        α::Vector{Float64},
        δ::Vector{Float64},
        nodes::Union{Vector{SPICEInterpControlNode},Vector{SPICEControlNode}}
    )

Process control-nodes with v-infinity's
"""
function process_controlnodes(
    ArrivalCondition::ACHyperbolic,
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    vinf::Vector{Float64},
    α::Vector{Float64},
    δ::Vector{Float64},
    nodes::Union{Vector{SPICEInterpControlNode},Vector{SPICEControlNode},Vector{CartesianReferenceNode}}
)
    # control nodes
    node_fwd, node_bck = get_control_nodes(SFset, i, et_fwd, et_bck, nodes)
    # if cmp(SFset.control_node_type, "interp") == 0
    #     node_fwd, node_bck =
    #         get_control_nodes(SFset, i, et_fwd, et_bck, ctuples)
    # elseif cmp(SFset.control_node_type, "spice") == 0
    #     node_fwd, node_bck = get_control_nodes(SFset, i, et_fwd, et_bck)
    # end

    # compute & append Delta-V   # FIXME RA/DEC implementation!
    vinf_fwd =
        SFset.process_dv_fun(SFset.mu, node_fwd, [vinf[i], α[2i-1], δ[2i-1]])
    vinf_bck =
        SFset.process_dv_fun(SFset.mu, node_bck, [vinf[i+1], α[2i], δ[2i]])
    if SFset.zero_soi == true
        control_fwd = node_fwd + vcat(zeros(3, 1), vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
        control_bck = node_bck + vcat(zeros(3, 1), vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
    else
        rSOI_vinf_fwd =
            get_body_soi(SFset.phase_bodies[i], SFset.center_body) /
            SFset.lstar * vinf_fwd / norm(vinf_fwd)
        control_fwd = node_fwd + vcat(rSOI_vinf_fwd, vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
        rSOI_vinf_bck =
            get_body_soi(SFset.phase_bodies[i+1], SFset.center_body) /
            SFset.lstar * vinf_bck / norm(vinf_bck)
        control_bck = node_bck + vcat(rSOI_vinf_bck, vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
    end
    # return state
    return node_fwd, node_bck, control_fwd, control_bck
end


"""
    function process_controlnodes(
        ArrivalCondition::Union{ACL1,ACL2,ACL3,ACL4,ACL5},
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        vinf::Vector{Float64},
        α::Vector{Float64},
        δ::Vector{Float64},
        ctuples::Union{Nothing,Vector},
    )

Process control-nodes with v-infinity's
"""
function process_controlnodes(
    ArrivalCondition::Union{ACL1,ACL2,ACL3,ACL4,ACL5},
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    vinf::Vector{Float64},
    α::Vector{Float64},
    δ::Vector{Float64},
    ctuples::Union{Nothing,Vector},
    nodes::Vector{AbstractControlNodeType},
)
    # control nodes
    node_fwd, node_bck = get_control_nodes(SFset, i, et_fwd, et_bck, ctuples)
    # compute & append Delta-V   # FIXME RA/DEC implementation!
    vinf_fwd =
        SFset.process_dv_fun(SFset.mu, node_fwd, [vinf[i], α[2i-1], δ[2i-1]])
    if SFset.zero_soi == true
        control_fwd = node_fwd + vcat(zeros(3, 1), vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    else
        rSOI_vinf_fwd =
            get_body_soi(SFset.phase_bodies[i], SFset.center_body) /
            SFset.lstar * vinf_fwd / norm(vinf_fwd)
        control_fwd = node_fwd + vcat(rSOI_vinf_fwd, vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    end

    # check if arriving to final body
    if i != SFset.np
        vinf_bck =
            SFset.process_dv_fun(SFset.mu, node_bck, [vinf[i+1], α[2i], δ[2i]])
        if SFset.zero_soi == true
            control_bck =
                node_bck + vcat(zeros(3, 1), vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        else
            rSOI_vinf_bck =
                get_body_soi(SFset.phase_bodies[i+1], SFset.center_body) /
                SFset.lstar * vinf_bck / norm(vinf_bck)
            control_bck =
                node_bck + vcat(rSOI_vinf_bck, vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        end
    else   #if arriving to final body, overwrite control_bck
        control_bck = get_lagrange_point_eclipj2000(
            SFset,
            ArrivalCondition,
            et_bck,
            SFset.phase_bodies[i+1],
            1.0,  # FIXME - ω_cr3bp,
        )
    end
    # return state
    return node_fwd, node_bck, control_fwd, control_bck
end



"""
    function process_controlnodes(
        ArrivalCondition::AClpo1,
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        vinf::Vector{Float64},
        α::Vector{Float64},
        δ::Vector{Float64},
        ctuples::Union{Nothing,Vector},
    )

Process control-nodes with v-infinity's
"""
function process_controlnodes(
    ArrivalCondition::AClpo1,
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    vinf::Vector{Float64},
    α::Vector{Float64},
    δ::Vector{Float64},
    ctuples::Union{Nothing,Vector},
    nodes::Vector{AbstractControlNodeType},
)
    # control nodes
    node_fwd, node_bck = get_control_nodes(SFset, i, et_fwd, et_bck, ctuples)
    # compute & append Delta-V   # FIXME RA/DEC implementation!
    vinf_fwd =
        SFset.process_dv_fun(SFset.mu, node_fwd, [vinf[i], α[2i-1], δ[2i-1]])
    if SFset.zero_soi == true
        control_fwd = node_fwd + vcat(zeros(3, 1), vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    else
        rSOI_vinf_fwd =
            get_body_soi(SFset.phase_bodies[i], SFset.center_body) /
            SFset.lstar * vinf_fwd / norm(vinf_fwd)
        control_fwd = node_fwd + vcat(rSOI_vinf_fwd, vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    end

    # check if arriving to final body
    if i != SFset.np
        vinf_bck =
            SFset.process_dv_fun(SFset.mu, node_bck, [vinf[i+1], α[2i], δ[2i]])
        if SFset.zero_soi == true
            control_bck =
                node_bck + vcat(zeros(3, 1), vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        else
            rSOI_vinf_bck =
                get_body_soi(SFset.phase_bodies[i+1], SFset.center_body) /
                SFset.lstar * vinf_bck / norm(vinf_bck)
            control_bck =
                node_bck + vcat(rSOI_vinf_bck, vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        end
    else   #if arriving to final body, overwrite control_bck
        control_bck = get_lpo_patchpoint_eclipj2000(
            SFset,
            ArrivalCondition,
            et_bck,
            SFset.phase_bodies[i+1],
            α[end],
        ) #, ArrivalCondition.CR3BP_param, ω_cr3bp)
    end
    # return state
    return node_fwd, node_bck, control_fwd, control_bck
end


"""
    function process_controlnodes(
        ArrivalCondition::AClpo2,
        SFset::SimsFlanaganSettings,
        i::Int,
        et_fwd::Float64,
        et_bck::Float64,
        vinf::Vector{Float64},
        α::Vector{Float64},
        δ::Vector{Float64},
        ctuples::Union{Nothing,Vector},
    )

Process control-nodes with v-infinity's
"""
function process_controlnodes(
    ArrivalCondition::AClpo2,
    SFset::SimsFlanaganSettings,
    i::Int,
    et_fwd::Float64,
    et_bck::Float64,
    vinf::Vector{Float64},
    α::Vector{Float64},
    δ::Vector{Float64},
    ctuples::Union{Nothing,Vector},
    nodes::Vector{AbstractControlNodeType},
)
    # control nodes
    node_fwd, node_bck = get_control_nodes(SFset, i, et_fwd, et_bck, ctuples)
    # compute & append Delta-V   # FIXME RA/DEC implementation!
    vinf_fwd =
        SFset.process_dv_fun(SFset.mu, node_fwd, [vinf[i], α[2i-1], δ[2i-1]])
    if SFset.zero_soi == true
        control_fwd = node_fwd + vcat(zeros(3, 1), vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    else
        rSOI_vinf_fwd =
            get_body_soi(SFset.phase_bodies[i], SFset.center_body) /
            SFset.lstar * vinf_fwd / norm(vinf_fwd)
        control_fwd = node_fwd + vcat(rSOI_vinf_fwd, vinf_fwd ./ SFset.vstar)[:]   # convert km/s -> canonical
    end

    # check if arriving to final body
    if i != SFset.np
        vinf_bck =
            SFset.process_dv_fun(SFset.mu, node_bck, [vinf[i+1], α[2i], δ[2i]])
        if SFset.zero_soi == true
            control_bck =
                node_bck + vcat(zeros(3, 1), vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        else
            rSOI_vinf_bck =
                get_body_soi(SFset.phase_bodies[i+1], SFset.center_body) /
                SFset.lstar * vinf_bck / norm(vinf_bck)
            control_bck =
                node_bck + vcat(rSOI_vinf_bck, vinf_bck ./ SFset.vstar)[:]   # convert km/s -> canonical
        end

    else   #if arriving to final body, overwrite control_bck
        # construct manifold interpolation
        manifold_interp_func =
            get_manifold_interp_func(ArrivalCondition.sim_manifold, δ[end])
        # construct control node
        control_bck = get_lpo_patchpoint_eclipj2000(
            SFset,
            ArrivalCondition,
            manifold_interp_func,
            et_bck,
            SFset.phase_bodies[i+1],
            α[end],
        )  #, ArrivalCondition.CR3BP_param, ω_cr3bp)
    end

    # return state
    return node_fwd, node_bck, control_fwd, control_bck
end
