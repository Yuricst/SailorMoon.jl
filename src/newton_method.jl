"""
	differential correction via Newton method 
	currently Jacobian is computed via finite difference
"""

function fitness(x)
	# evaluate equality constraint violations
	return c_eq
end

function get_jacobian(func::Function, x, nc::Int, h::Float64=1e-6)
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

	end
	return xs, fs, converged
end