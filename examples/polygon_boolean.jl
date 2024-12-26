using Plots
using PolygonAlgorithms
using PolygonAlgorithms: x_coords, y_coords

function plot_regions!(canvas, regions; options...)
    for region in regions
        if length(region) == 1
            scatter!(canvas, x_coords(region), y_coords(region), color=:black, marker=:xcross, label="")
        else
            plot_region!(canvas, region; options...)
        end
    end
    canvas
end

function plot_region!(canvas, region; fill=(0, 0.4, :green), options...)
    idxs = vcat(1:length(region), 1)
    plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]), fill=fill, label="", color=:black, options...)
end

θs = 0.0:0.01:6π
rs = θs
spiral = [(r * cos(θ), r * sin(θ)) for (r, θ) in zip(rs, θs)] 
spiral = vcat(spiral, reverse([0.8 .* p for p in spiral[3:end]]))
star = [
    (0.0, 18.0), (3.0, 5.0), (15.0, 5.0), (5.0, 0.0), (10.0, -12.0), (0.0, -2.0),
    (-10.0, -12.0), (-5.0, 0.0), (-15.0, 5.0), (-3.0, 5.0)
]
poly1 = spiral
poly2 = star

idxs1 = vcat(1:length(poly1), 1)
idxs2 = vcat(1:length(poly2), 1)
canvas_base = plot(x_coords(poly1[idxs1]), y_coords(poly1[idxs1]),
    fill=(0, 0.5), aspectratio=:equal, xlabel="base", legend=:none)
plot!(canvas_base, x_coords(poly2[idxs2]), y_coords(poly2[idxs2]),
    fill=(0, 0.3))

regions_difference12 = difference_geometry(poly1, poly2, PolygonAlgorithms.MartinezRuedaAlg());
regions_difference21 = difference_geometry(poly2, poly1, PolygonAlgorithms.MartinezRuedaAlg());
regions_intersect = intersect_geometry(poly1, poly2, PolygonAlgorithms.MartinezRuedaAlg());
regions_union = union_geometry(poly1, poly2, PolygonAlgorithms.MartinezRuedaAlg());
regions_xor = xor_geometry(poly1, poly2, PolygonAlgorithms.MartinezRuedaAlg());

canvas_difference12 = plot(aspectratio=:equal, xlabel="difference 1-2", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_regions!(canvas_difference12, regions_difference12)
canvas_difference21 = plot(aspectratio=:equal, xlabel="difference 2-1", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_regions!(canvas_difference21, regions_difference21)
#canvas_intersect = plot(aspectratio=:equal, xlabel="intersection")
canvas_intersect = deepcopy(canvas_base)
plot!(canvas_intersect, xlabel="intersection", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_regions!(canvas_intersect, regions_intersect)
canvas_union = plot(aspectratio=:equal, xlabel="union", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_regions!(canvas_union, regions_union)
canvas_xor = plot(aspectratio=:equal, xlabel="xor", xlims=xlims(canvas_base), ylims=ylims(canvas_base))
plot_regions!(canvas_xor, regions_xor)

plot(canvas_base, canvas_difference12, canvas_difference21, canvas_intersect, canvas_union, canvas_xor,
    layout = (2, 3), 
    size=(900, 600),
    margin=5Plots.mm,
)
plot!()