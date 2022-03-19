module LightPolygon

using Reexport
using LinearAlgebra
using Plots
@reexport import Plots.plot
@reexport import Plots.plot!

export 
    node,
    node_p,
    edge,
    edge_p,
    polygon,
    polygon_w_holes,
    in_polygon,
    regular_poly,
    rectangle,
    square,
    area,
    visibility_polygon


include("poly_basics.jl")
include("poly_functions.jl")

end