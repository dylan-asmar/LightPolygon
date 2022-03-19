
"""
Calculate the area of basic polygon
"""
function area(poly::polygon)
    A = 0.
    n = length(poly.nodes)
    for ii = 1:n-1
        A += 1/2*(poly.nodes[ii].y + poly.nodes[ii+1].y)*(poly.nodes[ii].x - poly.nodes[ii+1].x)
    end
    A += 1/2*(poly.nodes[end].y + poly.nodes[1].y)*(poly.nodes[end].x - poly.nodes[1].x)
    return A
end

"""
Calculate the area of a polygon with holes. Assumes the holes to not overlap.
"""
function area(poly_h::polygon_w_holes)
    A_border = area(poly_h.border)
    
    # Assumes the holes do not overlap
    A_holes = 0.
    for hole in poly_h.holes
        A_holes += area(hole)
    end
    return A_border - A_holes
end

"""
Determine if point is in a polygon. Works by using a verticle line and horizontal line and counting
    the number of intersections with the polygon. Edge cases of when edges and/or nodes are aligned 
    with the vertical and horizontal lines when counting intercepts.

    Returns 
        - 1 if inside the polygon
        - 0 if outside the polygon
        - -1 if on the border of the polygon
"""
function in_polygon(p::Tuple{Real,Real}, poly::polygon)
    # Horizontal line
    cnt = 0
    cnt_edges = 0
    on_edge = false
    first_edge = false
    for e in poly.edges
        sx = e.s.x - p[1]
        ex = e.e.x - p[1]
        sy = e.s.y - p[2]
        ey = e.e.y - p[2]
        if (sx >= 0 || ex >= 0) && ((sy >= 0 && ey <= 0) || (sy <= 0 && ey >= 0))
            x_axis_intercept = (sx - ex) * sy / (ey - sy) + sx
            if isnan(x_axis_intercept)
                on_edge = true
            elseif isapprox(x_axis_intercept, 0, atol=10e-10)
                on_edge = true
            elseif x_axis_intercept > 0
                cnt += 1
                itx_end_v = (isapprox(x_axis_intercept + p[1], e.e.x) && isapprox(p[2], e.e.y))
                itx_str_v = (isapprox(x_axis_intercept + p[1], e.s.x) && isapprox(p[2], e.s.y))
                if (itx_end_v || itx_str_v) && !check_colinear([[e.s.x, e.s.y], [e.e.x, e.e.y], collect(p)])
                    if !first_edge
                        first_edge = true
                    elseif iseven(cnt)
                        cnt -= 1
                        first_edge = false
                    else
                        first_edge = false
                    end
                end
            end
        end
    end
    if on_edge
        return -1
    end
    horiz_line = isodd(cnt)

    # Vertical line
    cnt = 0
    cnt_edges = 0
    on_edge = false
    first_edge = false
    for e in poly.edges
        sx = e.s.x - p[1]
        ex = e.e.x - p[1]
        sy = e.s.y - p[2]
        ey = e.e.y - p[2]
        if (sy >= 0 || ey >= 0) && (sx >= 0 && ex <= 0) || (sx <= 0 && ex >= 0)
            y_axis_intercept = (sy - ey) * sx / (ex - sx) + sy
            if isnan(y_axis_intercept)
                if sign(ey) * sign(sy) < 1
                    on_edge = true
                elseif ey > 0 && sy > 0
                    cnt += 1
                end
            elseif isapprox(y_axis_intercept, 0, atol=10e-10)
                on_edge = true
            elseif y_axis_intercept > 0
                cnt += 1
                itx_end_v = (isapprox(y_axis_intercept + p[2], e.e.y) && isapprox(p[1], e.e.x))
                itx_str_v = (isapprox(y_axis_intercept + p[2], e.s.y) && isapprox(p[1], e.s.x))
                if (itx_end_v || itx_str_v) &&
                   !check_colinear([[e.s.x, e.s.y], [e.e.x, e.e.y], collect(p)])
                    if !first_edge
                        first_edge = true
                    elseif iseven(cnt)
                        cnt -= 1
                        first_edge = false
                    else
                        first_edge = false
                    end
                end
            end
        end
    end
    if on_edge
        return -1
    end
    vert_line = isodd(cnt)

    return horiz_line && vert_line
end

# TODO: Extend this function to return the proper values if on edges vs just true or false.
"""
Determine if point is in a polygon with holes. Must be on the border or inside the border
    and not inside a hole. If a point is on the border of a hole, it is still inside counted
    as inside. 
"""
function in_polygon(p::Tuple{Real,Real}, poly::polygon_w_holes)
    if in_polygon(p, poly.border) == 0
        return false
    end
    for ph in poly.holes
        if in_polygon(p, ph) == 1
            return false
        end
    end
    return true
end

"""
Calculate the visibility polygon from a point inside a polygon
"""
function visibility_polygon(pt::Tuple{Real,Real}, poly::polygon)
    return visibility_polygon(pt, poly.edges, poly)
end
function visibility_polygon(pt::Tuple{Real,Real}, poly::polygon_w_holes)
    return visibility_polygon(pt, poly.edges, poly.border)
end

function visibility_polygon(pt::Tuple{Real,Real}, edges::Vector{edge}, border::polygon)
    edges_polar = Vector{edge_p}(undef, length(edges) * 2)
    edges′ = Vector{edge}(undef, length(edges) * 2)
    cnt = 0
    for (ii, ed) in enumerate(edges)
        split_edge = false
        xs′ = ed.s.x - pt[1]
        ys′ = ed.s.y - pt[2]
        xe′ = ed.e.x - pt[1]
        ye′ = ed.e.y - pt[2]
        rm, θm1, θm2 = 0.0, 0.0, 2π

        # If we are crossing angle 0 on this edge, break up the edge at the 0 angle
        if (xs′ > 0 || xe′ > 0) && (ys′ > 0 && ye′ < 0) || (ys′ < 0 && ye′ > 0)
            x_intercept = (xs′ - xe′) * ys′ / (ye′ - ys′) + xs′
            if x_intercept > 0
                split_edge = true
                rm = x_intercept
                if ye′ > 0
                    θm1, θm2 = 2π, 0.0
                end
            end
        end

        (rs, θs) = cart2pol(xs′, ys′)
        (re, θe) = cart2pol(xe′, ye′)

        if (isapprox(ys′, 0) && isapprox(xs′, 0)) || (isapprox(ye′, 0) && isapprox(xe′, 0))
            if isapprox(ys′, 0)
                rs = xs′
                θs = ye′ > 0 ? 0.0 : 2π
            end
            if isapprox(ye′, 0)
                re = xe′
                θe = ys′ > 0 ? 0.0 : 2π
            end
        end

        θs_near_0 = isapprox(θs, 0) || isapprox(θs, 2π)
        θe_near_0 = isapprox(θe, 0) || isapprox(θe, 2π)
        if θs_near_0
            θs = ye′ > 0 ? 0.0 : 2π
        end
        if θe_near_0
            θe = ys′ > 0 ? 0.0 : 2π
        end

        if isapprox(xs′, 0) && isapprox(ys′, 0)
            θs = θe
        elseif isapprox(xe′, 0) && isapprox(ye′, 0)
            θe = θs
        end

        if split_edge
            if θs <= θm1
                cnt += 1
                edges_polar[cnt] = edge_p(node_p(rs, θs), node_p(rm, θm1))
                edges′[cnt] = edge(node(xs′, ys′), node(rm, 0.0))
            end
            if θm2 <= θe
                cnt += 1
                edges_polar[cnt] = edge_p(node_p(rm, θm2), node_p(re, θe))
                edges′[cnt] = edge(node(rm, 0.0), node(xe′, ye′))
            end
        else
            if θs <= θe
                cnt += 1
                edges_polar[cnt] = edge_p(node_p(rs, θs), node_p(re, θe))
                edges′[cnt] = edge(node(xs′, ys′), node(xe′, ye′))
            end
        end
    end
    resize!(edges_polar, cnt)
    resize!(edges′, cnt)

    Q_mat = Matrix{Float64}(undef, cnt * 2, 6)
    cnt = 0
    for ii in 1:length(edges_polar)
        cnt += 1
        Q_mat[cnt, 1] = edges_polar[ii].s.r
        Q_mat[cnt, 2] = edges_polar[ii].s.θ
        Q_mat[cnt, 3] = 0.0
        Q_mat[cnt, 4] = ii
        Q_mat[cnt, 5] = edges′[ii].s.x
        Q_mat[cnt, 6] = edges′[ii].s.y
        cnt += 1
        Q_mat[cnt, 1] = edges_polar[ii].e.r
        Q_mat[cnt, 2] = edges_polar[ii].e.θ
        Q_mat[cnt, 3] = 1.0
        Q_mat[cnt, 4] = ii
        Q_mat[cnt, 5] = edges′[ii].e.x
        Q_mat[cnt, 6] = edges′[ii].e.y
    end

    for ii in [3, 1, 2]
        Q_mat = sort_by_column(Q_mat, ii)
    end

    SL = Vector{Real}()
    SLr = Vector{Real}()
    PV = Vector{Vector{Real}}()
    active_edge = Vector{Float64}(undef, 1)
    active_edge_un = Vector{Float64}()

    # current_vertex = Q_mat[1, :]
    active_edge[1] = Q_mat[1, 4]
    push!(active_edge_un, active_edge[1])
    pushfirst!(PV, Q_mat[1, [1, 2, 5, 6]])

    for ii in 2:size(Q_mat, 1)
        current_vertex = Q_mat[ii, :]
        if current_vertex[3] == 1
            deleteat!(active_edge_un, findall(x -> x == current_vertex[4], active_edge_un))
        end

        handle_event_point(current_vertex, active_edge, PV, SL, SLr, edges_polar, edges′, active_edge_un, border, pt)
        if length(PV) >= 2
            while length(PV) >= 2 && isapprox(PV[1][1], PV[2][1]) && isapprox(PV[1][2], PV[2][2])
                deleteat!(PV, 1)
            end
            # If reinstating, need to change the call to colinear
            # while length(PV) >=3 && check_colinear([pv[3:4] for pv in PV[1:3]])
            #     deleteat!(PV, 2)
            # end
        end
        if length(PV) >= 4
            angles = [PV[ii][2] for ii in length(PV):-1:2]
            angles_traversed = mod(sum(diff(angles)), 2π)
            full_circle = isapprox(angles_traversed, 0.0, atol=10e-10)
            full_circle = full_circle || isapprox(angles_traversed, 2π, atol=10e-10)
            if full_circle
                popfirst!(PV)
                break
            end
        end
    end
    PV_tup = [Tuple(pv[3:4]) .+ pt for pv in PV]
    return polygon(PV_tup)
end

function handle_event_point(current_vertex::Vector{T}, active_edge::Vector{T},
    PV::Vector{Vector{Real}},
    SL::Vector{Real}, SLr::Vector{Real},
    edges_polar::Vector{edge_p}, edges::Vector{edge},
    active_edge_un::Vector{T},
    border::polygon, pt::Tuple{Real,Real}) where {T}

    # If current vertex is an end vertex and current vertex is part of the active edge
    if current_vertex[3] == 1 && current_vertex[4] == active_edge[1]
        second_check = true
        pushfirst!(PV, current_vertex[[1, 2, 5, 6]])

        while !isempty(SL)
            angs, Rs = Vector{Real}(), Vector{Real}()
            intx_cs = Vector{Tuple{Real,Real,Real}}()
            intx_ps = Vector{Tuple{Real,Real}}()
            remove_idxs = Vector{Int}()

            for (ii, ed_idx) in enumerate(SL)
                if edges_polar[Int(ed_idx)].e.θ > current_vertex[2]
                    intx_c, intx_p = intersection_pt(edges[Int(ed_idx)], current_vertex[5:6])
                    push!(angs, intx_p[2])
                    push!(Rs, intx_p[1])
                    push!(intx_cs, intx_c)
                    push!(intx_ps, intx_p)
                else
                    push!(remove_idxs, ii)
                end
            end

            deleteat!(SL, remove_idxs)
            deleteat!(SLr, remove_idxs)

            if !isempty(angs)
                idx = 0
                bool_vec = minimum(angs) .≈ angs
                angs_idx = collect(1:length(angs))[bool_vec]
                if length(angs_idx) > 1
                    Rs_idx = argmin(Rs[angs_idx])
                    idx = angs_idx[Rs_idx]
                else
                    idx = angs_idx[1]
                end

                second_check = false
                intx_edge = popat!(SL, idx)
                intx_r = popat!(SLr, idx)

                intx_c, intx_p = intx_cs[idx], intx_ps[idx]

                pushfirst!(PV, [collect(intx_p); collect(intx_c[1:2])])
                active_edge[1] = intx_edge
                push!(active_edge_un, active_edge[1])
                break
            end
        end

        if !isempty(active_edge_un) && second_check
            intx_c, intx_p = intersection_pt(edges[Int(active_edge_un[end])], current_vertex[5:6])
            if ((intx_p[1] > current_vertex[1] || isapprox(intx_p[1], current_vertex[1]))
                && intx_c[3] >= 0 && intx_c[3] <= 1)

                pushfirst!(PV, [collect(intx_p); collect(intx_c[1:2])])
                active_edge[1] = active_edge_un[end]
            end
        end
    end

    # If current vertex is a start vertex
    if current_vertex[3] == 0
        intx_c, intx_p = intersection_pt(edges[Int(active_edge[1])], current_vertex[5:6])
        if (intx_p[1] < current_vertex[1] - 10e-10) && !isapprox(intx_p[1], 0)
            # idx = searchsortedfirst(SLr, current_vertex[1])            
            # insert!(SL, idx, current_vertex[4])
            # insert!(SLr, idx, current_vertex[1])
            push!(SL, current_vertex[4])
            push!(SLr, current_vertex[1])
        else
            if !isnan(intx_c[1]) && !isapprox(intx_p[1], 0)
                pushfirst!(PV, [collect(intx_p); collect(intx_c[1:2])])
            end
            pushfirst!(PV, current_vertex[[1, 2, 5, 6]])
            active_edge[1] = current_vertex[4]
            push!(active_edge_un, active_edge[1])
            #else

        end
    end
end

function intersection_pt(ed::edge, pt::Vector{T}, perp::Vector{A}) where {T,A}
    x1, y1 = pt[1], pt[2]
    pt2 = pt + perp
    x2, y2 = pt2[1], pt2[2]
    return intersection_from_pts(ed.s.x, ed.s.y, ed.e.x, ed.e.y, x1, y1, x2, y2)
end

function intersection_pt(ed::edge, pt::Vector{T}) where {T}
    intx_pt = intersection_from_pts(ed.s.x, ed.s.y, ed.e.x, ed.e.y, 0, 0, pt[1], pt[2])
    return intx_pt, (norm(intx_pt[1:2]), mod(atan(intx_pt[2], intx_pt[1]), 2π))
end

# Calculate interseciton given two points on each line segment
function intersection_from_pts(x1::Real, y1::Real, x2::Real, y2::Real,
    x3::Real, y3::Real, x4::Real, y4::Real)
    t_num = det([x1-x3 x3-x4; y1-y3 y3-y4])
    t_den = det([x1-x2 x3-x4; y1-y2 y3-y4])
    t = t_num / t_den
    px = x1 + t * (x2 - x1)
    py = y1 + t * (y2 - y1)
    return (px, py, t)
end

function check_colinear(pv_pts::Vector{Vector{T}}) where {T}
    x1, y1 = pv_pts[1][1], pv_pts[1][2]
    x2, y2 = pv_pts[2][1], pv_pts[2][2]
    x3, y3 = pv_pts[3][1], pv_pts[3][2]

    a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)
    return abs(a) < 10e-10
end

"""
Helper function to convert cartesion to polar
"""
function cart2pol(x::Real, y::Real)
    r = round(norm([x, y]), digits=5)
    θ = round(mod(atan(y, x), 2π), digits=5)
    return (r, θ)
end

"""
Helper function to conver polar to cartesion
"""
function pol2cart(r::Real, θ::Real)
    x = round(r * cos(θ), digits=5)
    y = round(r * sin(θ), digits=5)
    return (x, y)
end

"""
Sort a matrix based on a column index
"""
function sort_by_column(M::Matrix, col::Int; kwargs...)
    return M[sortperm(M[:, col]; kwargs...), :]
end