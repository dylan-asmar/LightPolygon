
using Plots
using LightPolygon

function run_vp_poly_test(view_pt=nothing)
    # Jagged polygon for border
    border_poly = polygon([
        (2, 2), (0.5, 1.7), (0.3, 1.65), (0.5, 2.1), (0.5, 2.3), 
        (-2, 2), (-2.3, 1.8), (-2.3, 1.2), (-2.1, 1.2), (-2.1, 1.7), (-1.9, 1.6),
        (-2, -2), (-1, -2), (-1, -2.3), (0, -2),
        (2, -2), (2, -1), (2.2, 0), (2.2, 1), (2.3, 1.7), (2.0, 1.2)]
    )

    # Vector of holes
    ph = Vector{polygon}()
    push!(ph, polygon([(-1.0, 0.2), (-1.05, 0.5), (-0.9, 0.51), (-1.15, 0.75),]))
    push!(ph, regular_poly((-1, 1), 8, 0.25))
    push!(ph, regular_poly((-0.79, -1.0), 8, 0.15))
    push!(ph, regular_poly((-1.5, 0.0), 8, 0.15))
    push!(ph, regular_poly((-.6, 0.9), 8, 0.15))
    push!(ph, regular_poly((-0.85, -0.5), 5, 0.3))
    push!(ph, square((1.05, -0.5), 0.5, pi/8))
    push!(ph, square((0.5, 1), 0.5, -pi/3))
    push!(ph, square((-1.5, -0.5), 0.5, pi/3.9))
    push!(ph, square((1.25, -2.0), 1.0, pi/7))
    push!(ph, square((0.0, -1.5), 1.0))
    push!(ph, square((0.75, 0.0), 0.5)) 
    push!(ph, rectangle((-0.1, -0.5), 0.75, 0.25, pi/5))

    # Creating polygon with holes
    poly_h = polygon_w_holes(border_poly, ph)

    # Select a random point if one isn't provided
    if isnothing(view_pt)
        view_pt = Tuple(rand(1,2)*5 .- 2.5)
        while in_polygon(view_pt, poly_h) == 0
            view_pt = Tuple(rand(1,2)*5 .- 2.5)
        end
    end

    # Get the visibility polygon
    vp = visibility_polygon(view_pt, poly_h)

    # Plot
    plt = plot(border_poly, aspect_ratio=:equal; xlim=[-2.5,2.5], ylim=[-2.5, 2.5], 
                legend=false, color=:black)
    for ph_i in ph
        plot!(ph_i, color=:white)
    end
    plt = plot!(vp, opacity=0.7, color=palette(:default)[5])
    plt = scatter!(view_pt, marker=:xcross, color=:red, markersize=5)
    display(plt)

    println("Point: $(view_pt)")
    println("Area of VP: $(area(vp))")
end

run_vp_poly_test()


