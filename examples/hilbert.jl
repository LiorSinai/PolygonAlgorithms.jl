using Plots
using PolygonAlgorithms
using PolygonAlgorithms: x_coords, y_coords

function hilbert_curve(
    point::Tuple{T, T},
    horiz::Tuple{T, T},
    vert::Tuple{T, T},
    level::Int
    ) where T <: AbstractFloat
    # Algorithm from http://www.fundza.com/algorithmic/space_filling/hilbert/basics/
    if level <= 0
        next = ((point[1] + (horiz[1] + vert[1])/2), (point[2] + (horiz[2] + vert[2])/2))
        return [next]
    end
    points = [
        hilbert_curve(point, (horiz[2]/2, vert[2]/2), (horiz[1]/2, vert[1]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2, point[2] + horiz[2]/2), (horiz[1]/2, horiz[2]/2), (vert[1]/2, vert[2]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2+vert[1]/2, point[2] + horiz[2]/2+vert[2]/2), (horiz[1]/2, horiz[2]/2), (vert[1]/2, vert[2]/2), level-1)...,
        hilbert_curve((point[1] + horiz[1]/2+vert[1], point[2] + horiz[2]/2+vert[2]), (-vert[1]/2, -vert[2]/2), (-horiz[1]/2, -horiz[2]/2), level-1)...
    ]
    points
end

function plot_regions!(canvas, regions)
    for region in regions
        if length(region) == 1
            scatter!(canvas, x_coords(region), y_coords(region), color=:black, marker=:xcross, label="")
        else
            idxs = vcat(1:length(region), 1)
            plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]), fill=(0, 0.4, :green), label="", color=:black)
        end
    end
end

order = 4
points = hilbert_curve((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), order);
tail = points[end]
head = points[1]
## fill inner
push!(points, (-tail[1], tail[2]))
push!(points, (-head[1], head[2]))

poly1 = reverse(points);
poly2 = PolygonAlgorithms.rotate(poly1, π/2.0, (0.5, 0.5))

idxs1 = vcat(1:length(poly1), 1)
idxs2 = vcat(1:length(poly2), 1)
canvas_base = plot(x_coords(poly1[idxs1]), y_coords(poly1[idxs1]), fill=(0, 0.5), aspectratio=:equal, legend=:none)

canvas_both = deepcopy(canvas_base)
plot!(canvas_both, x_coords(poly2[idxs2]), y_coords(poly2[idxs2]), fill=(0, 0.3))

regions_intersect = intersect_geometry(poly1, poly2, PolygonAlgorithms.WeilerAthertonAlg())
canvas_intersect = deepcopy(canvas_both)
plot_regions!(canvas_intersect, regions_intersect)

plot(canvas_base, canvas_both, canvas_intersect,
    layout = (1, 3), 
    xlims=(-0.1, 1.1),
    ylims=(-0.1, 1.1), 
    size=(900, 400)
)