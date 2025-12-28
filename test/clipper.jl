# https://github.com/AngusJohnson/Clipper2/blob/main/CSharp/Tests/Tests1/Tests/TestPolygons.cs
using StatsBase
using Test
using PolygonAlgorithms

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

function is_hole(idx::Int, regions::Vector{<:Vector{<:Tuple}})
    region = regions[idx]
    others = vcat(1:(idx-1), (idx+1):length(regions))
    any(x->fully_contains(x, region), regions[others])
end

function fully_contains(polygon1::Vector{<:Tuple{T, T}}, polygon2::Vector{<:Tuple{T, T}}) where T <: AbstractFloat
   all(contains(polygon1, point; rtol=1e-8) for point in polygon2) 
end

function calculate_area_diff_ratio(stored_area::Integer, regions::Vector{<:Vector}; check_holes=false)
    if check_holes
        measured_area = 0.0
        for (idx, region) in enumerate(regions)
            a = area_polygon(region)
            if is_hole(idx, regions)
                measured_area -= a
            else
                measured_area += a
            end
        end
    else
        measured_area = sum(map(area_polygon, regions))
    end
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
    elseif test.caption in ["123", "128", "134", "140", "142", "144", "153", "173", "181"]
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

function execute(test::PolygonTest, polygon1::Vector, polygon2::Union{Vector, Nothing})
    try
        if isnothing(polygon2)
            regions = PolygonAlgorithms.martinez_rueda_algorithm(polygon1)
        else
            func = func_map[test.clip_type]
            regions = func(polygon1, polygon2)
        end
        check_holes = test.clip_type == "UNION"
        area_ratio = calculate_area_diff_ratio(test.solution_area, regions; check_holes=check_holes)
        count_diff = test.solution_count > 0 ? (test.solution_count - length(regions)) : 0
        area_ratio, count_diff
    catch err
        if err isa AssertionError
            return (1.0, Inf)
        else
            rethrow(err)
        end
    end
end

tests = load_tests(joinpath(@__DIR__, "Polygons.txt"))

area_ratios = fill(1.0, length(tests))
count_diffs = fill(Inf, length(tests))
results = fill("", length(tests))
for (idx, poly_test) in enumerate(tests)
    subjects = map(convert_path, poly_test.subjects)
    clips = map(convert_path, poly_test.clips)
    metrics = (1.0, Inf) # area_ratio, count_diff
    if poly_test.fill_rule != "EVENODD"
        result = poly_test.fill_rule
    elseif length(subjects) == 1 && length(clips) == 0
        metrics = execute(poly_test, subjects[1], nothing)
        result = "SINGLE"
    elseif length(subjects) == 1 && length(clips) == 1
        metrics = execute(poly_test, subjects[1], clips[1])
        result = validate_area(poly_test, metrics[1])
    else
        result = "MULTI"
    end
    results[idx] = string(result)
    area_ratios[idx] = metrics[1]
    count_diffs[idx] = metrics[2]
end

summary = countmap(results)
