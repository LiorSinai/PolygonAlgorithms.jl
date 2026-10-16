using Plots
using PolygonAlgorithms: Polygon, x_coords, y_coords
using PolygonAlgorithms: Segment2D

function plot_polygons!(canvas, polygons::AbstractVector{<:Polygon}; options...)
    for polygon in polygons
        if length(polygon.exterior) == 1
            scatter!(canvas, x_coords(polygon.exterior), y_coords(polygon.exterior), marker=:xcross, label="")
        else 
            plot_polygon!(canvas, polygon; options...)
        end
    end
    canvas
end

function plot_polygon!(canvas, polygon::Polygon; options...)
    region = polygon.exterior
    idxs = vcat(1:length(region), 1)
    plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]); options...)
    for hole in polygon.holes
        idxs = vcat(1:length(hole), 1)
        plot!(canvas, x_coords(hole[idxs]), y_coords(hole[idxs]); fill=(0, 0.7, :grey70), color=:black, label="")
    end
    canvas
end

function plot_segment!(canvas, segment::Segment2D; options...)
    plot!(canvas, [segment[1][1], segment[2][1]], [segment[1][2], segment[2][2]]; options...)
end