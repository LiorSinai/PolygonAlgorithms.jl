using Plots
using PolygonAlgorithms
using PolygonAlgorithms: Polygon, x_coords, y_coords

function plot_polygons!(canvas, polygons::AbstractVector{<:Polygon}; options...)
    for polygon in polygons
        if length(polygon.exterior) == 1
            scatter!(canvas, x_coords(polygon.exterior), y_coords(polygon.exterior), color=:black, marker=:xcross, label="")
        else 
            plot_polygon!(canvas, polygon; label="", color=:black, options...)
        end
    end
    canvas
end

function plot_polygon!(canvas, polygon::Polygon; options...)
    region = polygon.exterior
    idxs = vcat(1:length(region), 1)
    plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]); fill=(0, 0.4, :green), options...)
    for hole in polygon.holes
        idxs = vcat(1:length(hole), 1)
        plot!(canvas, x_coords(hole[idxs]), y_coords(hole[idxs]); fill=(0, 1.0, :grey70), options...)
    end
    canvas
end

θs = 0.0:0.01:6π
rs = θs
spiral = [(r * cos(θ), r * sin(θ)) for (r, θ) in zip(rs, θs)] 
spiral = vcat(spiral, reverse([0.8 .* p for p in spiral[3:end]]))
star = [
    (0.0, 18.0), (3.0, 5.0), (15.0, 5.0), (5.0, 0.0), (10.0, -12.0), (0.0, -2.0),
    (-10.0, -12.0), (-5.0, 0.0), (-15.0, 5.0), (-3.0, 5.0)
]
polygon1 = Polygon(spiral)
polygon2 = Polygon(star)

canvas_base = plot_polygon!(plot(), polygon1, aspectratio=:equal, xlabel="base", legend=:none, fill=(0, 0.5))
plot_polygon!(canvas_base, polygon2; fill=(0, 0.3))

regions_difference12 = difference_geometry(polygon1, polygon2, PolygonAlgorithms.MartinezRuedaAlg());
regions_difference21 = difference_geometry(polygon2, polygon1, PolygonAlgorithms.MartinezRuedaAlg());
regions_intersect = intersect_geometry(polygon1, polygon2, PolygonAlgorithms.MartinezRuedaAlg());
regions_union = union_geometry(polygon1, polygon2, PolygonAlgorithms.MartinezRuedaAlg());
regions_xor = xor_geometry(polygon1, polygon2, PolygonAlgorithms.MartinezRuedaAlg());

canvas_difference12 = plot_polygons!(
    plot(xlabel="difference 1-2", xlims=xlims(canvas_base), ylims=ylims(canvas_base)),
    regions_difference12, 
)
canvas_difference21 = plot_polygons!(
    plot(aspectratio=:equal, xlabel="difference 2-1", xlims=xlims(canvas_base), ylims=ylims(canvas_base)),
    regions_difference21
)
canvas_intersect = deepcopy(canvas_base)
plot!(canvas_intersect, xlabel="intersection", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_polygons!(canvas_intersect, regions_intersect)
canvas_union = plot_polygons!(
    plot(aspectratio=:equal, xlabel="union", xlims=xlims(canvas_base), ylims=ylims(canvas_base)),
    regions_union
)
canvas_xor = plot_polygons!(
    plot(aspectratio=:equal, xlabel="xor", xlims=xlims(canvas_base), ylims=ylims(canvas_base)),
    regions_xor
)

plot(canvas_base, canvas_difference12, canvas_difference21, canvas_intersect, canvas_union, canvas_xor,
    layout = (2, 3), 
    size=(900, 600),
    margin=5Plots.mm,
)
plot!()