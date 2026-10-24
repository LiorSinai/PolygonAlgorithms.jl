using Plots
using PolygonAlgorithms
using PolygonAlgorithms: Polygon, x_coords, y_coords

include("plot_polygon.jl")

base_polygon = Polygon([(1.0, 1.0), (9.0, 1.0), (5.0, 6.0),])

canvas_base = plot_polygon!(
    plot(), base_polygon, aspectratio=:equal, fill=(0, 0.5),
    legend=:none, arrow=true, label="subject",
    xlims=(0, 9), ylims=(0, 7)
)
canvas_base

clips = [
    Polygon([(1.0, 4.0), (5.0, 4.0), (3.0, 2.0)]), # intersect
    Polygon([(3.0, 3.5), (6.0, 3.5), (5.0, 2.0)]), # touch edge
    Polygon([(3.0, 3.5), (7.0, 3.5), (5.0, 1.0)]), # touch only edges
    Polygon([(4.0, 3.5), (6.0, 3.5), (5.0, 2.0)]), # hole
]
canvases = []
for clip in clips
    p = plot_polygon!(deepcopy(canvas_base), clip; fill=(0, 0.5), label="clip")
    push!(canvases, p)
end

for face_selection in [PolygonAlgorithms.SPLIT_FACES, PolygonAlgorithms.MERGE_FACES]
    for clip in clips
        difference12 = difference_geometry(
            base_polygon,
            clip
            ; face_selection=face_selection
        );
        p = plot_polygons!(plot(), difference12, fill=(0, 0.5), label="",
            aspectratio=:equal,
            #xlims=xlims(canvas_base),
            #ylims=ylims(canvas_base), 
        )
        push!(canvases, p)
    end
end

plot!(canvases[1], ylabel="subject & clip")
plot!(canvases[5], ylabel="SPLIT_FACES")
plot!(canvases[9], ylabel="MERGE_FACES")
plot!(
    canvases..., layout=(3, 4), link=:y,
    yticks=(0:3:7), 
    size=(900, 600), margin=4Plots.mm,
    plot_title="Face selection"
)