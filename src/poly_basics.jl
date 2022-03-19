

"""
Node/vertex to use in a polygon
"""
struct node
    x::Real
    y::Real
end

"""
Polar node/vertex to use in a polygon
"""
struct node_p
    r::Real
    θ::Real
end

"""
Edge of polygon.
    s - sart node
    e - end node
"""
struct edge
    s::node
    e::node
end

"""
Edge of polygon with polar nodes

    s - sart node
    e - end node
"""
struct edge_p
    s::node_p
    e::node_p
end

"""
Basic polygon

    nodes - vector of nodes
    edges - vector of edges
"""
struct polygon
    nodes::Vector{node}
    edges::Vector{edge}
end

"""
Basic polygon with holes

    border - basic polygon
    holes - Vector of polygons defining the holes
    edges - Vector of all edges. Composed of `border` in CCW order and edges/nodes of `holes` in CW order
"""
struct polygon_w_holes
    border::polygon
    holes::Vector{polygon}
    edges::Vector{edge}
end

"""
Create polygon from a vector of edges
"""
function polygon(edges::Vector{edge})
    poly_edges = Vector{edge}()
    for ii in 1:length(edges)
        if ii > 1 && edges[ii].s != edges[ii-1].e
            push!(poly_edges, edge(edges[ii-1].e, edges[ii].s))
        end
        push!(poly_edges, edges[ii])
    end
    if edges[1].s != edges[end].e
        push!(poly_edges, edge(edges[end].e, edges[1].s))
    end

    nodes = Vector{node}(undef, length(poly_edges))
    for (ii, ed) in enumerate(poly_edges)
        nodes[ii] = ed.s
    end
    check_edges!(poly_edges)
    return polygon(nodes, poly_edges)
end

"""
Create polygon from vector of (x,y) tuples
"""
function polygon(node_vec::Vector{Tuple{T,A}}) where {T,A}
    nodes = Vector{node}(undef, length(node_vec))
    edges = Vector{edge}(undef, length(node_vec))
    for (i, n) in enumerate(node_vec)
        nodes[i] = node(n[1], n[2])
        if i > 1
            edges[i-1] = edge(nodes[i-1], nodes[i])
        end
    end
    edges[end] = edge(nodes[end], nodes[1])
    check_edges!(edges)
    return polygon(nodes, edges)
end

"""
Create a polygon with holes

    - border: basic polygon
    - holes: vector of basic polygons

The holes vector must be basic polygons with nodes in CCW order
"""
function polygon_w_holes(border::polygon, holes::Vector{polygon})
    holes′ = condition_polys(holes)
    mod_holes = Vector{polygon}()
    for hole_poly in holes′
        mod_edges = non_overlap_poly_border(border, hole_poly)
        if !isempty(mod_edges)
            push!(mod_holes, polygon(mod_edges))
        end
    end
    num_edges = length(border.edges) + sum([length(p.edges) for p in mod_holes])
    edges = Vector{edge}(undef, num_edges)
    edges[1:length(border.edges)] = border.edges
    cnt = length(border.edges)
    for hole_poly in mod_holes
        for ii in length(hole_poly.edges):-1:1
            cnt += 1
            edges[cnt] = edge(hole_poly.edges[ii].e, hole_poly.edges[ii].s)
        end
    end
    return polygon_w_holes(border, holes, edges)
end

"""
Create regular polygon
"""
function regular_poly(center::Tuple{Real,Real}, num_sides::Int, radius::Real, orientation::Real=π / 2)
    num_sides > 2 || error("Regular polygon does not exist for $num_sides sides")
    cart_pts = Vector{Tuple{Real,Real}}(undef, num_sides)
    ϕ = 2π / num_sides
    for ii in 0:num_sides-1
        cart_pts[ii+1] = pol2cart(radius, orientation + ϕ * ii) .+ center
    end
    return polygon(cart_pts)
end

"""
Create rectangle from center, width, height, and rotation
"""
function rectangle(center::Tuple{Real,Real}, width::Real, height::Real, orientation::Real=0.0)
    pts = Vector{Tuple{Real,Real}}(undef, 4)
    pts[1] = (width / 2, height / 2)
    pts[2] = (-width / 2, height / 2)
    pts[3] = (-width / 2, -height / 2)
    pts[4] = (width / 2, -height / 2)

    for ii in 1:4
        pts[ii] = Tuple(rotate_vector(collect(pts[ii]), orientation)) .+ center
    end
    return polygon(pts)
end

"""
Create a square given a center, width, and orientation
"""
function square(center::Tuple{Real,Real}, width::Real, orientation::Real=0.0)
    return rectangle(center, width, width, orientation)
end

"""
Function to rotate a 2D vector 
"""
function rotate_vector(vector_in::Vector{T}, θ::A) where {T,A}
    length(vector_in) == 2 || error("Only implemented for 2D vectors")
    rot_matrix = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    return rot_matrix * vector_in
end

"""
Helper function to remove edges with the same start and end nodes.
"""
function check_edges!(edges::Vector{edge})
    idxs = Vector{Int}()
    for (ii, ed) in enumerate(edges)
        if is_equal(ed.s, ed.e)
            push!(idxs, ii)
        end
    end
    deleteat!(edges, idxs)
end

"""
Test whether two nodes are equal
"""
function is_equal(n1::node, n2::node)
    return (n1.x == n2.x) && (n1.y == n2.y)
end

"""
Get modified edges of p2 based on overlaps with p1
"""
function non_overlap_poly_border(p1::polygon, p2::polygon)
    mod_edges = Vector{edge}()
    for ed in p2.edges
        s_in = in_polygon((ed.s.x, ed.s.y), p1)
        e_in = in_polygon((ed.e.x, ed.e.y), p1)
        if s_in != -1 && e_in != -1
            if s_in && e_in
                push!(mod_edges, ed)
            elseif s_in || e_in
                px, py = get_edge_intx(ed, p1)
                mod_node = node(px, py)
                mod_ed = s_in ? edge(ed.s, mod_node) : edge(mod_node, ed.e)
                push!(mod_edges, mod_ed)
            elseif (s_in == 0 && e_in == 0)

            end
        elseif (s_in == 1 && e_in == -1) || (e_in == 1 && s_in == -1)
            push!(mod_edges, ed)
        end
    end
    return mod_edges
end

"""
Calculate intercept of an edge with another polygon
"""
function get_edge_intx(ed::edge, poly::polygon)
    for p_ed in poly.edges
        px, py, t = intersection_from_pts(ed.s.x, ed.s.y, ed.e.x, ed.e.y,
            p_ed.s.x, p_ed.s.y, p_ed.e.x, p_ed.e.y)
        if t >= 0 && t <= 1.0
            return (px, py)
        end
    end
    error("Did not find a valid intersection")
end

"""
For a vector of polygons, condition so no polygons overlap
"""
function condition_polys(polys::Vector{polygon})
    new_polys = Vector{polygon}()
    n = length(polys)
    for ii = 1:n
        changed = false
        mod_edges = polys[ii].edges
        mod_poly = polygon(mod_edges)
        for jj = 1:n
            if jj == ii
                continue
            end
            for n in mod_poly.nodes
                if in_polygon((n.x, n.y), polys[jj]) == 1
                    changed = true
                    mod_edges = non_overlap_poly(polys[jj], polygon(mod_edges))
                    break
                end
            end
            if changed
                mod_poly = polygon(mod_edges)
            end
        end
        push!(new_polys, mod_poly)
    end
    return new_polys
end

function non_overlap_poly(p1::polygon, p2::polygon)
    mod_edges = Vector{edge}()
    for ed in p2.edges
        s_in = in_polygon((ed.s.x, ed.s.y), p1)
        e_in = in_polygon((ed.e.x, ed.e.y), p1)
        if s_in != -1 && e_in != -1
            if !s_in && !e_in
                push!(mod_edges, ed)
            elseif s_in || e_in
                px, py = get_edge_intx(ed, p1)
                mod_node = node(px, py)
                mod_ed = s_in ? edge(mod_node, ed.e) : edge(ed.s, mod_node)
                push!(mod_edges, mod_ed)
            end
        elseif e_in < 1 || s_in < 1
            push!(mod_edges, ed)
        end
    end
    return mod_edges
end

function plot(poly::polygon; kwargs...)
    x_coords = [n.x for n ∈ poly.nodes]
    y_coords = [n.y for n ∈ poly.nodes]
    plt = plot(Shape(x_coords, y_coords); kwargs...)
    return plt
end

function plot_border(poly::polygon; kwargs...)
    x_coords = [n.x for n ∈ poly.nodes]
    x_coords = [x_coords; poly.nodes[1].x]
    y_coords = [n.y for n ∈ poly.nodes]
    y_coords = [y_coords; poly.nodes[1].y]
    plot(x_coords, y_coords; kwargs...)
end

function plot_border!(poly::polygon; kwargs...)
    x_coords = [n.x for n ∈ poly.nodes]
    x_coords = [x_coords; poly.nodes[1].x]
    y_coords = [n.y for n ∈ poly.nodes]
    y_coords = [y_coords; poly.nodes[1].y]
    plot!(x_coords, y_coords; kwargs...)
end

function plot!(poly::polygon; kwargs...)
    x_coords = [n.x for n ∈ poly.nodes]
    y_coords = [n.y for n ∈ poly.nodes]
    plot!(Shape(x_coords, y_coords); kwargs...)
end