using Plots
using PolygonAlgorithms
using PolygonAlgorithms: Polygon, x_coords, y_coords

include("plot_polygon.jl")

# self intersecting star
base_polygon = Polygon([(-3.0, 2.0), (3.0, 2.0), (-2.0, -2.0), (0.0, 5.0), (2.0, -2.0)])
# self intersecting 9
base_polygon = Polygon([
    (10.0, 60.0), (75.0, 60.0), (80.0, 30.0), (55.0, 30.0), (50.0, 90.0),
    (40.0, 90.0), (30.0, 10.0), (100.0, 10.0), (90.0, 120.0), (20.0, 120.0),
    (10.0, 60.0)
]);

canvas_base = plot_polygon!(
    plot(), base_polygon, aspectratio=:equal, fill=(0, 0.5),
    legend=:none, arrow=true,
)
canvas_base

canvases = []
for face_selection in [PolygonAlgorithms.SPLIT_FACES, PolygonAlgorithms.MERGE_FACES]
    for fill_rule in [PolygonAlgorithms.EVEN_ODD, PolygonAlgorithms.NON_ZERO]
        out = union_geometry(base_polygon; fill_rule=fill_rule, face_selection=face_selection);
        p = plot_polygons!(plot(), out, fill=(0, 0.5), label="",
            aspectratio=:equal,
            #xlims=xlims(canvas_base),
            #ylims=ylims(canvas_base), 
        )
        push!(canvases, p)
    end
end

plot!(canvases[1], title="EVEN_ODD", ylabel="SPLIT_FACES")
plot!(canvases[2], title="NON_ZERO")
plot!(canvases[3], ylabel="MERGE_FACES")
plot!(canvases..., layout=(2, 2), link=:y)
