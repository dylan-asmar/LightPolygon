
using Test
using PolygonOps
using LightPolygon

"""
Use `inpolygon` function from PolygonOps.jl to check if a point is in the polygon.
"""
function inpolygon(p::Tuple{Real,Real}, poly::polygon)
    vec_of_nodes = Vector{Vector{Real}}(undef, length(poly.nodes) + 1)
    for (i, n) in enumerate(poly.nodes)
        vec_of_nodes[i] = [n.x, n.y]
    end
    vec_of_nodes[end] = [poly.nodes[1].x, poly.nodes[1].y]
    return inpolygon(p, vec_of_nodes)
end

@test 1 == 1