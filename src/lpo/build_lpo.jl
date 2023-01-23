"""
build a LPO object and save it as a BSON file 
"""

include("../../../julia-r3bp/R3BP/src/R3BP.jl")


function build_lpo(lp::Int=2, Az_km::Real=1200.0, northsouth::Integer=3, dt::Real=0.005)

    param3b = dynamics_parameters()
    lps = lagrange_points(param3b.mu2)

    println("Halo guess Az_km: $Az_km")
    guess0 = R3BP.halo_analytical_construct(param3b.mu2, lp, Az_km, param3b.lstar, northsouth)
    res = R3BP.ssdc_periodic_xzplane([param3b.mu2,], guess0.x0, guess0.period, fix="period")

    x0_stm = vcat(res.x0, reshape(I(6), (6^2,)))[:]
    prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, res.period, (param3b.mu2))
    sol = solve(prob_cr3bp_stm, Tsit5(), reltol=1e-12, abstol=1e-12)#, saveat=LinRange(0, period, n+1))
    monodromy = R3BP.get_stm(sol, 6)   # get monodromy matrix
    ys0 = R3BP.get_eigenvector(monodromy, true, 1);

    # build a dictionary: res.x0, res.period, ys0
    if northsouth == 1
        ns = "N"
    elseif northsouth == 3
        ns = "S"
    end 

    filename = "lpo_L" * string(lp) * "_" * string(Int(Az_km)) * "km_" * ns * ".bson"
    bson((filename), Dict(:x0 => res.x0, :period => res.period, :ys0 => ys0))

end



function load_lpo(filename, type::Int=1, epsr::Real=1e-6, epsv::Real=1e-6, abstol::Real=1e-12, reltol::Real=1e-12, dt::Real=0.005)
    _load_dict = BSON.load(filename)
    x0_stm = vcat(_load_dict[:x0], reshape(I(6), (6^2,)))[:]
    prob_cr3bp_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, _load_dict[:period], (param3b.mu2))

    if type == 1 
        # CR3BPLPO object 
        LPOArrival = SailorMoon.CR3BPLPO(
            [el for el in _load_dict[:x0]],
            _load_dict[:period],
            [el for el in _load_dict[:ys0]],
            prob_cr3bp_stm, epsr, Tsit5(), abstol, reltol, dt
        )

    elseif type == 2
        # CR3BPLPO2 object
        LPOArrival = SailorMoon.CR3BPLPO2(
            [el for el in _load_dict[:x0]],
            _load_dict[:period],
            [el for el in _load_dict[:ys0]],
            prob_cr3bp_stm, epsr, epsv, Tsit5(), abstol, reltol, dt
        )
    end

end
