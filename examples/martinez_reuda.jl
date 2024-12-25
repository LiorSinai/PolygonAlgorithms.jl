using Plots
using PolygonAlgorithms
using PolygonAlgorithms: x_coords, y_coords
using PolygonAlgorithms: SegmentEvent, event_loop!, convert_to_event_queue, add_annotated_segment!
using PolygonAlgorithms: apply_selection_criteria

function plot_segment_event!(
    canvas,
    event::SegmentEvent
    ; 
    d::Float64=0.2,
    col::Union{Colorant, Symbol},
    self_col::Union{Colorant, Symbol},
    other_col::Union{Colorant, Symbol},
    options...
    )
    segment = event.segment
    plot!(canvas, [segment[1][1], segment[2][1]], [segment[1][2], segment[2][2]]; arrow=true, label="", color=col, options...)
    midpoint = ((segment[1][1]+segment[2][1])/2, (segment[1][2]+segment[2][2])/2)
    # Vectors
    Δx = (segment[2][1] - segment[1][1])
    Δy = (segment[2][2] - segment[1][2])
    Δr = sqrt(Δx * Δx + Δy *Δy )
    n = (d * Δy/Δr,  d * Δx/Δr) # normal vector to the segment
    u = (d * Δx/Δr,  d * Δy/Δr) # unit vector to the segment
    # self annotations
    annotations = event.self_annotations
    filled = isnothing(annotations.fill_above) ? :black : (annotations.fill_above ? self_col : :white)
    scatter!(canvas, [midpoint[1] - n[1] + u[1]], [midpoint[2] + n[2] + u[2]], marker=:circle, color=filled, label="")
    filled = isnothing(annotations.fill_below) ? :black : (annotations.fill_below ? self_col : :white)
    scatter!(canvas, [midpoint[1] + n[1] + u[1]], [midpoint[2] - n[2] + u[2]], marker=:circle, color=filled, label="")
    # other annotations
    annotations = event.other_annotations
    filled = isnothing(annotations.fill_above) ? :black : (annotations.fill_above ? other_col : :white)
    scatter!(canvas, [midpoint[1] - n[1] - u[1]], [midpoint[2] + n[2] - u[2]], marker=:diamond, color=filled, label="")
    filled = isnothing(annotations.fill_below) ? :black : (annotations.fill_below ? other_col : :white)
    scatter!(canvas, [midpoint[1] + n[1] - u[1]], [midpoint[2] - n[2] - u[2]], marker=:diamond, color=filled, label="")
end

function plot_segment!(
    canvas,
    segment::Segment2D
    ;
    col::Union{Colorant, Symbol}, options...
    )
    plot!(canvas, [segment[1][1], segment[2][1]], [segment[1][2], segment[2][2]]; label="", color=col, options...)
end

function calc_annotation_distance(canvas, polygon)
    (xmin, ymin, xmax, ymax) = bounds(polygon)
    width = max(xmax-xmin, ymax-ymin)
    d = width/canvas.attr[:size][2] * 4 # 4 pixels
    d
end

self_intersect = [
    (0.0, 0.0), (2.0, 2.0), (6.0, -2.0), (11.0, 2.0), (11.0, 0.0)
]
rectangle_horiz = [
    (-1.0, 0.0), (-1.0, 3.0), (12.0, 3.0), (12.0, 0.0)
]
polygon1 = self_intersect
polygon2 = rectangle_horiz

# martinez_rueda_algorithm
event_queue1 = convert_to_event_queue(polygon1; primary=true)
annotated_segments1 = event_loop!(event_queue1, self_intersection=true)
event_queue2 = convert_to_event_queue(polygon2; primary=false)
annotated_segments2 = event_loop!(event_queue2, self_intersection=true)
queue = SegmentEvent{Float64}[]
for ev in vcat(annotated_segments1, annotated_segments2)
    add_annotated_segment!(queue, ev)
end
annotated_segments3 = event_loop!(deepcopy(queue), self_intersection=false)
# for consistent reporting, swap annotations so that self annotations are always the primary
for ev in annotated_segments3
    if !ev.primary
        temp = ev.self_annotations
        ev.self_annotations = ev.other_annotations
        ev.other_annotations = temp
    end
end

# plot
colors = palette(:default)
idxs1 = vcat(1:length(polygon1), 1)
canvas_shapes = plot(x_coords(polygon1[idxs1]), y_coords(polygon1[idxs1]), aspectratio=:equal, arrow=true, fill=(0, 0.5))
idxs2 = vcat(1:length(polygon2), 1)
plot!(canvas_shapes, x_coords(polygon2[idxs2]), y_coords(polygon2[idxs2]), arrow=true, fill=(0, 0.5))

canvas_annotations = plot(aspect_ratio=:equal)
d = calc_annotation_distance(canvas_annotations, polygon2)
# for (i, event) in enumerate(annotated_segments1)
#     plot_segment_event!(canvas_annotations, event; col=colors[1], d=d, self_col=colors[1], other_col=colors[2])
# end
# for (i, event) in enumerate(annotated_segments2)
#     plot_segment_event!(canvas_annotations, event; col=colors[2], d=d, self_col=colors[2], other_col=colors[2])
# end
for (i, event) in enumerate(annotated_segments3)
    col = event.primary ? colors[1] : colors[2]
    plot_segment_event!(canvas_annotations, event; col=col, d=d, self_col=colors[1], other_col=colors[2])
end
canvas_annotations

## Selections
selected = apply_selection_criteria(annotated_segments3, PolygonAlgorithms.INTERSECTION_CRITERIA)

canvas_selected = deepcopy(canvas_shapes)
for (i, event) in enumerate(selected)
    plot_segment_event!(canvas_selected, event; col=:red, self_col=:red, other_col=colors[2], d=d)
end
canvas_selected

plot(canvas_annotations, canvas_selected)
