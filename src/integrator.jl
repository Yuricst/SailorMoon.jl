"""
Direct implementation of fixed-step RK4 integrator
"""

Base.@kwdef struct RK4Solution
    u
    t
    prob::ODEProblem
    retcode
    alg="RK4"
    event_states
    event_times
end


"""
    integrate_rk4(
        prob::ODEProblem,
        dt::Real,
        cbs::Union{Nothing,Vector}=nothing,
        save_all::Bool=true,
        exact_tf::Bool=true,
    )

Fixed-step RK4 integration.
"""
function integrate_rk4(
    prob::ODEProblem,
    dt::Real,
    cbs::Union{Nothing,Vector}=nothing,
    save_all::Bool=true,
    exact_tf::Bool=true,
)
    # get number of steps
    nsteps = Int(ceil(abs((prob.tspan[2] - prob.tspan[1])/dt)))
    direction = sign(prob.tspan[2] - prob.tspan[1])
    dt = direction * dt
    # initialize
    u_iter = prob.u0
    t_iter = prob.tspan[1]
    us = [u_iter,]
    ts = Real[t_iter,]
    nx = length(prob.u0)
    k1, k2, k3, k4 = zeros(nx), zeros(nx), zeros(nx), zeros(nx)
    # for checking callback change of sign
    if isnothing(cbs) == false
        n_cb = length(cbs)
        _check_sign_cbs = [0.0 for i = 1:n_cb]
        _event_states = [[] for i = 1:n_cb]
        _event_times  = [[] for i = 1:n_cb]
    else
        _event_states = nothing
        _event_times  = nothing
    end
    retcode = :Initialized
    # iterate
    for i = 1:nsteps
        # adjust time-step
        if exact_tf && abs(t_iter + dt) > abs(prob.tspan[2])
            dt = prob.tspan[2] - t_iter
        end

        prob.f(k1, u_iter, prob.p, t_iter)
        prob.f(k2, u_iter+0.5dt*k1, prob.p, t_iter+0.5dt)
        prob.f(k3, u_iter+0.5dt*k2, prob.p, t_iter+0.5dt)
        prob.f(k4, u_iter+dt*k3, prob.p, t_iter+dt)

        # update state and time
        u_iter = u_iter + 1/6 * dt * (k1 + 2k2 + 2k3 + k4)
        t_iter += dt

        if save_all == true || i == nsteps
            push!(us, u_iter)
            push!(ts, t_iter)
        end

        # check callback
        if isnothing(cbs) == false
            for (icb,cb) in enumerate(cbs)
                update = cb.condition(u_iter, t_iter, 0)  # 0 is a placeholder for integrator
                if update * _check_sign_cbs[icb] < 0
                    # store event states
                    push!(_event_states[icb], u_iter)
                    push!(_event_times[icb],  t_iter)
                    # check whether to terminate
                    if cb.affect!() == true
                        # terminate
                        retcode = :Terminated
                        if save_all == false  # if not saving all the steps, make sure it is saved
                            push!(us, u_iter)
                            push!(ts, t_iter)
                        end
                        break
                    end
                else
                    if isnan(update)
                        _check_sign_cbs[icb] = 0
                    else
                        _check_sign_cbs[icb] = update  # update value of callback
                    end
                end
            end
        end
        # break if necessary
        if retcode == :Terminated
            break
        end
    end
    # create solution object
    if abs(ts[end]) >= abs(prob.tspan[2])
        retcode = :Success
    else
        retcode = :PrematureEnd
    end
    res = RK4Solution(
        us, ts, prob, retcode, "RK4", _event_states, _event_times
    )
    return res
end
