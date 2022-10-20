
μ1, μ2, μS, lstar, tstar, as, oms, oml = dyanmics_parameters()

# elliptical GTO
# initial state generation
period = 10.5 * 60 * 60  # s 
e = 0.73 
a = (μ1 / (4*pi^2) * period^2)^(1/3)

rp = a*(1-e)
vp = sqrt(μ*(1/rp - 1/(2*a)))

# make a state vector
state_p = [rp, 0, 0, 0, vp, 0]

state_p_rot = transform_earthIne_to_EMrot(state_p, 0, ωm)



