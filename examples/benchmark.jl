using PolygonAlgorithms
using PolygonAlgorithms: x_coords, y_coords
using PolygonAlgorithms: ChasingEdgesAlg, PointSearchAlg, MartinezRuedaAlg, WeilerAthertonAlg
using Random
using Plots
using Printf

include("polygons.jl")

function plot_regions!(canvas, regions; options...)
    for idx in 1:length(regions)
        region = regions[idx]
        isempty(region) && continue
        if length(region) == 1
            scatter!(canvas, x_coords(region), y_coords(region), color=:black, marker=:xcross, label="")
        else 
            plot_region!(canvas, region; options...)
        end
    end
    canvas
end

function plot_region!(canvas, region; options...)
    idxs = vcat(1:length(region), 1)
    plot!(canvas, x_coords(region[idxs]), y_coords(region[idxs]); options...)
end

function intersection_area_circles(r1::AbstractFloat, origin1::Tuple, r2::AbstractFloat, origin2::Tuple)
    # rotate and translate x-axis to x' so that origin is on origin1 and passes through origin2
    d = sqrt((origin2[1] - origin1[1])^2 + (origin2[2] - origin1[2])^2)
    #     x'^2 + y'^2 = r1^2
    # (x'-d)^2 + y'^2 = r2^2
    x = (d^2 - r2^2 + r1^2) / 2d
    y = sqrt(r1^2 - x^2)
    #area = sector - 2 * triangles
    #     = 2 * 0.5 * r^2 * 2θ - 2 * 0.5 * h * y 
    a1 = r1^2 * acos(x/r1) - x * y 
    a2 = r2^2 * acos((d-x)/r2) - (d-x) * y
    a1 + a2
end

function calculate_metrics(alg, fn, polygon1, polygon2, expected_area
        ; ntrials::Int=5)
    stats = @timed fn(alg, polygon1, polygon2; atol=1e-12, rtol=1e-12)
    t = stats.time
    bytes = stats.bytes
    for i in 2:ntrials
        si = @timed fn(alg, polygon1, polygon2; atol=1e-12, rtol=1e-12)
        t += si.time
        bytes += si.bytes
    end
    value = fn == intersect_convex ? [stats.value] : stats.value
    area = sum(area_polygon, value; init=0.0)
    Dict(
        "time" => t/ntrials,
        "bytes" => bytes/ntrials,
        "value" => value,
        "n1" => length(polygon1),
        "n2" => length(polygon2),
        "area" => area,
        "ratio" => area / expected_area,
        "count" => length(value),
    )
end

rng = MersenneTwister(42)
fork = make_forks(50; step_size=5.0, height=500.0, base=5.0)
cases = [
    (   
        # worst case convex algorithm
        "circles",
        make_circle(1000), make_circle(1000; origin=(0.5, 0.4)),
        1.883178, # intersection area
    ),
    (
        # worst case convex algorithm
        "circles no overlap",
        make_circle(1000), make_circle(1000; origin=(2.0, 0.4)),
        0.0,
    ),
    (
        # worst case concave polygon
        "forks",
        fork, [(y, x) for (x, y) in fork],
        64950.0, # = (n-1)^2*step_size^2 + (2(height-step_size)*step_size - step_size^2)
    ),
    (
        # fractal intersection with many regions of intersection
        "hilbert 30°",
        make_hilbert_polygon(6), PolygonAlgorithms.rotate(make_hilbert_polygon(6), π/3.0, (0.5, 0.5)),
        0.216820,
    ),
    (
        # fractal intersection with many overlapping segments
        "hilbert 90°",
        make_hilbert_polygon(6), PolygonAlgorithms.rotate(make_hilbert_polygon(6), π/2.0, (0.5, 0.5)),
        0.276855,
    ),
    (
        # random shape with spiked points
        "random", random_spiked_balls(rng, 300), random_spiked_balls(rng, 300, (0.3, 0.0)),
        0.175942,
    ),
    (
        # random shape with self-intersections
        "random self-intersecting",
        random_self_intersecting(rng, 100), random_self_intersecting(rng, 100, (0.3, 0.0)),
        0.136359
    ),
]

alg = MartinezRuedaAlg() # ChasingEdgesAlg, PointSearchAlg, MartinezRuedaAlg, WeilerAthertonAlg
alg_name = typeof(alg).name.name
fn = intersect_geometry # intersect_convex, intersect_geometry, difference_geometry, union_geometry, xor_geometry

# warm up
region = fn(alg, [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)], [(0.0, 0.0), (4.0, 1.0), (2.0, 0.0)])

println("Start benchmarking for $alg_name")
metrics = Dict{String, Dict}()
for (case, polygon1, polygon2, area) in cases
    metrics[case] = calculate_metrics(alg, fn, polygon1, polygon2, area)
end
println("Done benchmarking")
println("$alg_name | $fn")
@printf("%-25s | %4s | %4s | %8s | %11s | %6s | %6s \n" ,
    "case", "n1", "n2", "time (s)", "memory (MB)", "area", "ratio") 
for (case, _, _) in cases
    case_metrics = metrics[case]
    @printf("%-25s | %4s | %4s | %1.6f | %11d | %1.4f | %1.4f\n" , 
        case, case_metrics["n1"], case_metrics["n2"],
        case_metrics["time"],
        case_metrics["bytes"]/1024/1024,
        case_metrics["area"],
        case_metrics["ratio"],
    )
end

## Plotting
println("\nPlotting ...")
canvases = Dict{String, Any}()
for (case, polygon1, polygon2) in cases
    println(" $case")
    idxs1 = vcat(1:length(polygon1), 1)
    canvas_base = plot(
        x_coords(polygon1[idxs1]), y_coords(polygon1[idxs1]), 
        ; aspectratio=:equal, fill=true, alpha=0.5, legend=:none,
        )
    idxs2 = vcat(1:length(polygon2), 1)
    plot!(canvas_base, x_coords(polygon2[idxs2]), y_coords(polygon2[idxs2]), fill=true, alpha=0.5)
    regions = metrics[case]["value"]
    canvas_algo = plot_regions!(
        deepcopy(canvas_base),
        #plot(),
        regions,
        ; fill=(0, 0.4, :green), color=:black, legend=:none,
        aspectratio=:equal, xlims=xlims(canvas_base), ylims=ylims(canvas_base),
    )
    canvases[case] = [canvas_base, canvas_algo]
end

println(first.(cases))
case, _ = cases[4]
plot(canvases[case]..., layout=(1, 2), plot_title="$(alg_name)")

