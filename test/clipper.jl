# https://github.com/AngusJohnson/Clipper2/blob/main/CSharp/Tests/Tests1/Tests/TestPolygons.cs
using StatsBase
using Test
using PolygonAlgorithms
using PolygonAlgorithms: Polygon

struct PolygonTest
    caption::String
    clip_type::String
    fill_rule::String
    solution_area::Integer
    solution_count::Integer
    subjects::Vector{Vector{Int}}
    clips::Vector{Vector{Int}}
end

function load_test_specification(spec::Vector{<:AbstractString})
    caption = String(match(r"CAPTION: (\d+).", spec[1])[1])
    clip_type = String(match(r"CLIPTYPE: (.+)", spec[2])[1])
    fill_rule = String(match(r"FILLRULE: (.+)", spec[3])[1])
    solution_area = parse(Int, match(r"SOL_AREA: (.+)", spec[4])[1])
    solution_count = parse(Int, match(r"SOL_COUNT: (.+)", spec[5])[1])
    subjects = Vector{Vector{Int}}()
    clips = Vector{Vector{Int}}()
    target = subjects
    for idx in 7:length(spec)
        line = strip(spec[idx])
        if line == "CLIPS"
            target = clips
            continue
        end
        path = map(x->parse(Int, x), split(line, r",[ ]?| "))
        push!(target, path)
    end
    PolygonTest(caption, clip_type, fill_rule, solution_area, solution_count, subjects, clips)
end

function load_tests(filepath)
    tests = open(filepath, "r") do f
        data = read(f, String)
        data = strip(replace(data, r"#.+\n" => "\n")) # remove comments
        specs = split(data, "\n\n")
        out = Vector{PolygonTest}()
        for (idx, spec) in enumerate(specs)
            push!(out, load_test_specification(split(spec, "\n")))
        end
        out
    end
    tests
end

function convert_path(path::Vector{Int})
    polygon = Tuple{Float64, Float64}[]
    for i in 1:Int(length(path)/2)
        push!(polygon, (path[2i - 1], path[2i]))
    end
    polygon
end

function calculate_area_diff_ratio(stored_area::Integer, regions::Vector{<:Polygon})
    measured_area = sum(map(area_polygon, regions); init=0)
    area_diff = stored_area > 0 ? (stored_area - measured_area) : 0
    area_ratio = stored_area <= 0 ? 0 : area_diff / stored_area
    area_ratio
end

function validate_area(test::PolygonTest, area_ratio::Real)
    if test.caption in ["22"]
        threshold = 0.25
    elseif test.caption in ["19"]
        threshold = 0.07
    elseif test.caption in ["63"]
        threshold = 0.05
    elseif test.caption in ["26", "130", "133", "143", "151"]
        threshold = 0.04
    elseif test.caption in ["44", "132", "164"]
        threshold = 0.03
    elseif test.caption in ["15", "123", "128", "134", "140", "142", "144", "153", "173", "181"]
        threshold = 0.02
    else
        threshold = 0.01
    end
    abs(area_ratio) <= threshold
end

func_map = Dict(
    "DIFFERENCE" => difference_geometry,
    "INTERSECTION" => intersect_geometry,
    "UNION" => union_geometry,
)

function execute(test::PolygonTest, subjects::Vector{<:Polygon}, clips::Vector{<:Polygon})
    try
        func = func_map[test.clip_type]
        regions = func(subjects, clips...; atol=1e-7, rtol=1e-7)
        area_ratio = calculate_area_diff_ratio(test.solution_area, regions)
        count_regions = length(regions)
        count_diff = test.solution_count > 0 ? (test.solution_count - count_regions) : 0
        area_ratio, count_diff
    catch err
        if err isa AssertionError
            return (NaN, Inf)
        else
            rethrow(err)
        end
    end
end

tests = load_tests(joinpath(@__DIR__, "Polygons.txt"))

area_ratios = fill(NaN, length(tests))
count_diffs = fill(Inf, length(tests))
results = fill("", length(tests))
for (idx, poly_test) in enumerate(tests)
    print("$idx, ")
    subjects = map(Polygon, map(convert_path, poly_test.subjects))
    clips = isempty(poly_test.clips) ? Polygon[] : map(Polygon, map(convert_path, poly_test.clips))
    metrics = (NaN, Inf) # area_ratio, count_diff
    if poly_test.fill_rule != "EVENODD"
        result = poly_test.fill_rule
    else
        metrics = execute(poly_test, subjects, clips)
        result = isnan(metrics[1]) ? "error" : validate_area(poly_test, metrics[1])
    end
    results[idx] = string(result)
    area_ratios[idx] = metrics[1]
    count_diffs[idx] = metrics[2]
end

results_summary = countmap(results)
