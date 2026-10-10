using Plots
using PolygonAlgorithms
using PolygonAlgorithms: Polygon

include("plot_polygon.jl")

θs = 0.0:0.01:6π
rs = θs
spiral = [(r * cos(θ), r * sin(θ)) for (r, θ) in zip(rs, θs)] 
spiral = vcat(spiral, reverse([0.8 .* p for p in spiral]))
star = [
    (0.0, 18.0), (3.0, 5.0), (15.0, 5.0), (5.0, 0.0), (10.0, -12.0), (0.0, -2.0),
    (-10.0, -12.0), (-5.0, 0.0), (-15.0, 5.0), (-3.0, 5.0)
]
polygon1 = Polygon(spiral)
polygon2 = Polygon(star)

colors = palette(:default)
canvas_base = plot_polygon!(
    plot(), polygon1, aspectratio=:equal,
    xlabel="base", fill=(0, 0.5), color=colors[1]
)
plot_polygon!(canvas_base, polygon2; fill=(0, 0.5), color=colors[2])

regions_difference12 = difference_geometry(polygon1, polygon2);
regions_difference21 = difference_geometry(polygon2, polygon1);
regions_intersect = intersect_geometry(polygon1, polygon2);
regions_union = union_geometry(polygon1, polygon2);
regions_xor = xor_geometry(polygon1, polygon2);

canvas_difference12 = plot_polygons!(
    plot(), regions_difference12;
    xlabel="difference 1-2", label="",
    xlims=xlims(canvas_base), ylims=ylims(canvas_base),
    fill=(0, 0.5), color=:green,
)
canvas_difference21 = plot_polygons!(
    plot(),
    regions_difference21;
    xlabel="difference 2-1", label="",
    xlims=xlims(canvas_base), ylims=ylims(canvas_base),
    fill=(0, 0.5), color=:green,
)
canvas_intersect = deepcopy(canvas_base)
plot!(canvas_intersect, xlabel="intersection", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_polygons!(canvas_intersect, regions_intersect;
    fill=(0, 0.5), color=:green, label="")
canvas_union = plot_polygons!(
    plot(),
    regions_union;
    xlabel="union", label="",
    xlims=xlims(canvas_base), ylims=ylims(canvas_base),
    fill=(0, 0.5), color=:green,
)
canvas_xor = plot_polygons!(
    plot(),
    regions_xor;
    xlabel="XOR", label="",
    xlims=xlims(canvas_base), ylims=ylims(canvas_base),
    fill=(0, 0.5), color=:green,
)

plot(canvas_base, canvas_difference12, canvas_difference21, canvas_intersect, canvas_union, canvas_xor,
    layout = (2, 3), 
    size=(900, 600),
    margin=5Plots.mm,
)