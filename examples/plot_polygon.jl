using Plots
using PolygonAlgorithms: Polygon, x_coords, y_coords, is_clockwise
using PolygonAlgorithms: Point2D, Segment2D

function plot_polygons!(canvas, polygons::AbstractVector{<:Polygon}; options...)
    for polygon in polygons
        plot_polygon!(canvas, polygon; options...)
    end
    canvas
end

function plot_polygon!(canvas, polygon::Polygon; options...)
    region = polygon.exterior
    if length(region) == 1
        scatter!(canvas, x_coords(region), y_coords(region), marker=:xcross)
        return
    end
    idxs = vcat(1:length(region), 1)
    plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]); options...)
    # Set hole colors
    options = merge((;options...), (;fill=(0, 0.7, :grey70), color=:black))
    for hole in polygon.holes
        idxs = vcat(1:length(hole), 1)
        plot!(canvas, x_coords(hole[idxs]), y_coords(hole[idxs]); options...)
    end
    canvas
end

function plot_polygons!(canvas, polygons::Vector{<:Vector{<:Point2D}}; options...)
    for polygon in polygons
        isempty(polygon) && continue
        plot_polygon!(canvas, polygon; options...)
    end
    canvas
end

function plot_polygon!(canvas, polygon::Vector{<:Point2D}; options...)
    if length(polygon) == 1
        scatter!(canvas, x_coords(polygon), y_coords(polygon), color=:black, marker=:xcross)
        return
    end
    idxs = vcat(1:length(polygon), 1)
    if is_clockwise(polygon) # assume this is a hole
        options = merge((;options...), (;fill=(0, 0.7, :grey70), color=:black))
    end
    plot!(canvas, x_coords(polygon[idxs]), y_coords(polygon[idxs]); options...)
end

function plot_segment!(canvas, segment::Segment2D; options...)
    plot!(canvas, [segment[1][1], segment[2][1]], [segment[1][2], segment[2][2]]; options...)
end