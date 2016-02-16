__precompile__()
module LoewnerJohnEllipsoids

using Convex
export inner_ellipsoid,
	   outer_ellipsoid,
	   box,
	   barrier_value

type LinearConstraint{T<:Real}
	# Represents the constraint a * x <= b
	a::Vector{T}
	b::T
end

type Ellipsoid{T}
	# Represents an ellipsoid as the set 
	# { x | x = C*u + d, ||u||_2 <= 1 }
	# or the set
	# { x | ||(x - d)' * P' * P * (x - d) <= 1 }
	# where P = inv(C)' * inv(C)
	C::Array{T, 2}
	d::Array{T, 1}
	P::Array{T, 2}
end

Ellipsoid{T}(C::Array{T, 2}, d::Array{T, 1}) = begin
	C_inv = inv(C)
	Ellipsoid(C, d, C_inv' * C_inv)
end

Ellipsoid{T}(C::Array{T,2}, d::Array{T,2}) = Ellipsoid(C, vec(d))

function barrier_value{T}(ell::Ellipsoid{T}, x)
	# Evaluates the ellipsoid's barrier function, which is defined to be < 0
	# inside the ellipsoid, 0 on the boundary, and > 0 outside
	return ((x - ell.d)' * ell.P * (x - ell.d))[1] - 1
end

type Polyhedron{T}
	constraints::Vector{LinearConstraint{T}}
end

function box{T}(lower::Vector{T}, upper::Vector{T})
	@assert length(lower) == length(upper)
	dimension = length(lower)
	constraints = Array(LinearConstraint{T}, 2 * dimension)
	for i = 1:dimension
		a = zeros(T, dimension)
		a[i] = 1
		constraints[i] = LinearConstraint(a, upper[i])
	end
	for i = 1:dimension
		a = zeros(T, dimension)
		a[i] = -1
		constraints[i + dimension] = LinearConstraint(a, -lower[i])
	end
	Polyhedron(constraints)
end

function inner_ellipsoid{T<:Real}(polyhedron::Polyhedron{T})
	@assert length(polyhedron.constraints) > 2
	dimension = length(polyhedron.constraints[1].a)
	C = Variable(dimension, dimension)
	d = Variable(dimension)
	problem = maximize(logdet(C))
	for constraint in polyhedron.constraints
		problem.constraints += norm(C*constraint.a) <= constraint.b - dot(constraint.a, d)
	end
	solve!(problem)
	return Ellipsoid(C.value, d.value)
end

function outer_ellipsoid{T}(points::Vector{T})
	@assert length(points) > 1
	dimension = length(points[1])
	P = Variable(dimension, dimension)
	q = Variable(dimension)
	problem = maximize(logdet(P))
	for point in points
		problem.constraints += norm(P*point - q) <= 1
	end
	solve!(problem)
	P, q = P.value, q.value
	C = inv(P)
	return Ellipsoid(C, C * q)
end

end # module
