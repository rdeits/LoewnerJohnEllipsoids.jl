using LoewnerJohnEllipsoids
using Base.Test

# write your own tests here
function test_unit_circle()
	polyhedron = box([-1.,-1], [1., 1])
	e = inner_ellipsoid(polyhedron)
	@show e
	@test isapprox(e.d, [0,0], atol=1e-4)
	@test isapprox(e.C, [1 0; 0 1], atol=1e-4)
end
test_unit_circle()

function test_shifted_scaled_box()
	polyhedron = box([1,1], [2,2])
	e = inner_ellipsoid(polyhedron)
	@show e
	@test isapprox(e.d, [1.5, 1.5], atol=1e-4)
	@test isapprox(e.C, [0.5 0; 0 0.5], atol=1e-4)
end
test_shifted_scaled_box()

function test_outer_unit_box()
	points = Any[[-1,-1], [-1,1], [1,-1], [1,1]]
	e = outer_ellipsoid(points)
	@show e
	@test isapprox(e.d, [0, 0], atol=1e-4)
	@test isapprox(e.C, [sqrt(2) 0; 0 sqrt(2)], atol=1e-4)
end
test_outer_unit_box()

function test_barrier()
	polyhedron = box([1,1], [2,2])
	e = inner_ellipsoid(polyhedron)
	@test isapprox(barrier_value(e, [1,1.5]), 0, atol=1e-4)
	@test isapprox(barrier_value(e, [2,1.5]), 0, atol=1e-4)
	@test barrier_value(e, [1.5, 1.5]) < 0
	@test barrier_value(e, [2.5, 2.5]) > 0
end
test_barrier()

