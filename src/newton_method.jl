"""
	differential correction via Newton method 
	currently Jacobian is computed via finite difference
"""

function fitness(x)
	# evaluate equality constraint violations
	return c_eq
end

function get_jacobian(func::Function, x, nc::Int, h::Float64=1e-8)
	# for example, using central differencing (self-implemented)
	DF = zeros(nc,size(x,1))
	for i = 1:size(x,1)
		x_ptb_pls = [el for el in x]   # not sure if `x_ptb_plus = x` is ok. (in python we would need copy.copy())
		x_ptb_min = [el for el in x]
		x_ptb_pls[i] += h
		x_ptb_min[i] -= h
		DF[:,i] = (func(x_ptb_pls) - func(x_ptb_min))/(2h)
	end

	# of course, we could also use things like
	# FiniteDifferences.jl or ForwardDiff.jl etc.

	return DF
end


"""
	Considering only [x_lr, y_lr, xdot_lr, ydot_lr] for the variables
	corresponding to multishoot_trajectory2
	xfix: pass the original x0 array (including [x_lr, y_lr, xdot_lr, ydot_lr])
	xvar: [x_lr, y_lr, xdot_lr, ydot_lr]
"""
function get_jacobian2(func::Function, xfix::Vector, xvar::Vector, nc::Int, h::Float64=1e-7)
	# for example, using central differencing (self-implemented)
	DF = zeros(nc,size(xvar,1))
	x_ptb_pls = [el for el in xfix]   # not sure if `x_ptb_plus = x` is ok. (in python we would need copy.copy())
	x_ptb_min = [el for el in xfix]

	for i = 1:size(xvar,1)
		xvar[i] += h
		if i == 1
			x_ptb_pls[1] = xvar[i]
		elseif i == 2
			x_ptb_pls[2] = xvar[i]
		elseif i == 3
			x_ptb_pls[4] = xvar[i]
		elseif i == 4
			x_ptb_pls[5] = xvar[i]
		end
		
		xvar[i] -= 2h
		if i == 1
			x_ptb_min[1] = xvar[i]
		elseif i == 2
			x_ptb_min[2] = xvar[i]
		elseif i == 3
			x_ptb_min[4] = xvar[i]
		elseif i == 4
			x_ptb_min[5] = xvar[i]
		end

		DF[:,i] = (func(x_ptb_pls) - func(x_ptb_min))/(2h)
	end

	# of course, we could also use things like
	# FiniteDifferences.jl or ForwardDiff.jl etc.

	return DF
end


function newton_update(x, f, DF)
	nc, nx = size(DF)
	if nc == nx
		# Newton-Raphson
		x_new = x - inv(DF)*f
	elseif nc < nx
		# minimum-norm
		x_new = x - transpose(DF) * inv(DF*transpose(DF)) * f
	elseif nc > nx
		# least square
		x_new = x - inv(transpose(DF)*DF) * transpose(DF) * f
	end
	return x_new
end



function differential_correction(x0, fitness::Function, max_iter::Int=20, ctol::Float64=1e-6)
	# initialize
	x_iter = x0
	xs, fs = Vector[], Vector[]
	converged = false

	for i = 1:max_iter
		# evaluate fitness
		f_iter = fitness(x_iter)
		# store
		push!(xs, x_iter)
		push!(fs, f_iter)
		
		# check convergence
		if norm(f_iter) < ctol
			break
		end

		#  evaluate gradient
		nc = size(f_iter,1)   # number of equality constraints
		# println("nc:", nc)
		DF = get_jacobian(fitness, x_iter, nc)

		# update
		x_iter = newton_update(x_iter, f_iter, DF)

		println("iter #", i)
		println(f_iter)
		diff = x_iter - xs[end]
		println("update: ", [round(el, sigdigits=3) for el in diff])

	end
	return xs, fs, converged
end



function differential_correction2(x0, fitness::Function, max_iter::Int=20, ctol::Float64=1e-6)
	# initialize
	x_full = x0
	x_var = vcat(x0[1:2], x0[4:5])
	xs, fs = Vector[], Vector[]
	converged = false

	for i = 1:max_iter
		# evaluate fitness
		f_iter = fitness(x_full)
		# store
		push!(xs, x_full)
		push!(fs, f_iter)
		
		# check convergence
		if norm(f_iter) < ctol
			break
		end

		#  evaluate gradient
		nc = size(f_iter,1)   # number of equality constraints
		# println("nc:", nc)
		DF = get_jacobian2(fitness, x_full, x_var, nc)

		# update
		x_var = newton_update(x_var, f_iter, DF)

		println("iter #", i)
		println(f_iter)
		diff = x_var - xs[end]
		println("update: ", [round(el, sigdigits=3) for el in diff])

		x_full[1] = x_var[1]
		x_full[2] = x_var[2]
		x_full[4] = x_var[3]
		x_full[5] = x_var[4]

	end
	return xs, fs, converged
end