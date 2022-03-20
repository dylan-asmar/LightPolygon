
using Test
using LightPolygon
using PolygonOps
import PolygonOps: inpolygon

"""
Use `inpolygon` function from PolygonOps.jl to check if a point is in the polygon.
"""
function inpolygon(pt::Tuple{Real,Real}, poly::polygon)
    vec_of_nodes = Vector{Vector{Real}}(undef, length(poly.nodes) + 1)
    for (i, n) in enumerate(poly.nodes)
        vec_of_nodes[i] = [n.x, n.y]
    end
    vec_of_nodes[end] = [poly.nodes[1].x, poly.nodes[1].y]
    return inpolygon(pt, vec_of_nodes)
end

@testset verbose=true "All Tests" begin
    @testset "Areas" begin 
        pentagon1 = regular_poly((2., 3.), 5, 5, 0.)
        pentagon2 = regular_poly((-3., 2.), 5, 5, 2π*(rand()-0.5))
        @test isapprox(area(pentagon1), area(pentagon2), atol=1e-7)
        @test isapprox(area(pentagon1), 59.4410323, atol=1e-7)

        rect = rectangle(Tuple(rand(2)), 3.0, 6.0, 2π*(rand()-0.5))
        @test isapprox(area(rect), 18.0, atol=1e-7)

        square4 = square(Tuple(rand(2)), 2.0, 2π*(rand()-0.5))
        @test isapprox(area(square4), 4.0, atol=1e-7)

        reg99 = regular_poly(Tuple(rand(2)), 1000, 3.5, 2π*(rand()-0.5))
        @test isapprox(area(reg99), π*3.5^2, atol=1e-3)

        center = Tuple(rand(2) .- 0.5)
        rect_border = rectangle(center, 5., 10., 0.)
        hole1 = rectangle(center .+ (-0.5, 0.5), 0.5, 1., 0.)
        hole2 = regular_poly(center .+ (1., 3.), 3, 1.5, 0.)
        polygon_with_holes = polygon_w_holes(rect_border, [hole1, hole2])
        calculated_area = 50. - 0.5 - 2.92283574
        @test isapprox(area(polygon_with_holes), calculated_area, atol=1e-7)
    end

    @testset "Within Polygon Compare to PolygonOps" begin
        center = Tuple(rand(2))
        rect_border = rectangle(center, 5., 10., 0.)
        pts = []
        for _ = 1:10
            push!(pts, Tuple(rand(2)).*10 .- 5. .+ center)
        end
        push!(pts, center .+ (2.5, 0.))
        push!(pts, center .+ (-2.5, 0.))
        push!(pts, center .+ (0., 5.))
        push!(pts, center .+ (0., -5.))
        push!(pts, center .+ (2.5, 5.))
        push!(pts, center .+ (-2.5, -5.))
        for pt ∈ pts
            poly_ops = inpolygon(pt, rect_border)
            light_poly = in_polygon(pt, rect_border)
            @test poly_ops == light_poly
        end
    end

    @testset "Within Polygon with Holes" begin
        center = Tuple(rand(2) .- 0.5)
        rect_border = rectangle(center, 5., 10., 0.)
        hole1 = rectangle(center .+ (-1., -1.), 0.5, 1., 0.)
        hole2 = regular_poly(center .+ (0.75, 3.), 3, 1.5, 0.)
        polygon_with_holes = polygon_w_holes(rect_border, [hole1, hole2])
        
        pt = center .+ (-2.4, -4.9)
        @test in_polygon(pt, polygon_with_holes) == 1
        pt = center .+ (-1., -1.)
        @test in_polygon(pt, polygon_with_holes) == 0
        pt = center .+ (0.75, 3.)
        @test in_polygon(pt, polygon_with_holes) == 0

        pt = center .+ (2.5, 5.)
        @test in_polygon(pt, polygon_with_holes) == -1
        pt = center .+ (-1., -1.) .- (0.25, 0.)
        @test in_polygon(pt, polygon_with_holes) == -1
        pt = center .+ (0.75, 3.) .+ (1.5, 0.)
        @test in_polygon(pt, polygon_with_holes) == -1
    end

    @testset "Vis polygons" begin
        square_border = square((-4.5, 13.), 20., 0.)
        hole1 = square((-4.5, 18.), 2., 0.)
        hole2 = square((0.5, 13.), 2., 0.)
        polygon_with_holes = polygon_w_holes(square_border, [hole1, hole2])
        
        vis_poly = visibility_polygon((-4.5, 13.), square_border)
        @test area(square_border) == area(vis_poly)
        
        vis_poly = visibility_polygon((-4.5, 13.), polygon_with_holes)
        calculated_area = 400. - (25. - 4.) - (25. - 4.)
        @test calculated_area == area(vis_poly)
    end
end