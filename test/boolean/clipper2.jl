# https://github.com/AngusJohnson/Clipper2/blob/main/CSharp/Tests/Tests1/Tests/TestPolygons.cs
using Test
using PolygonAlgorithms
using PolygonAlgorithms: Polygon

struct PolygonTest
    caption::Int
    clip_type::String
    fill_rule::String
    solution_area::Integer
    solution_count::Integer
    subjects::Vector{Vector{Int}}
    clips::Vector{Vector{Int}}
end

function load_test_specification(spec::Vector{<:AbstractString})
    caption = parse(Int, match(r"CAPTION: (\d+).", spec[1])[1])
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
        data = strip(replace(data, r"#.+\r?\n" => "\n")) # remove comments
        specs = split(data, r"\r?\n\r?\n")
        out = Vector{PolygonTest}()
        for (idx, spec) in enumerate(specs)
            push!(out, load_test_specification(split(spec, r"\r?\n")))
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

function calculate_area_diff_ratio(stored_area::Integer, regions::Vector)
    measured_area = sum(map(first_moment, regions); init=0.0) # first moment is signed
    area_diff = stored_area > 0 ? (stored_area - measured_area) : 0
    area_ratio = stored_area <= 0 ? 0 : area_diff / stored_area
    area_ratio
end

function area_threshold(test::PolygonTest)
    if test.caption in [22]
        threshold = 0.25
    elseif test.caption in [193]
        threshold = 0.17
    elseif test.caption in [19]
        threshold = 0.07
    elseif test.caption in [63]
        threshold = 0.07
    elseif test.caption in [26]
        threshold = 0.04
    elseif test.caption in [44]
        threshold = 0.03
    elseif test.caption in [15, 184]
        threshold = 0.02
    elseif test.caption in [25, 27, 33, 52, 58, 64, 66, 78, 103, 117, 185]
        threshold = 0.01
    else
        threshold = 0.005
    end
    threshold
end

func_map = Dict(
    "DIFFERENCE" => difference_geometry,
    "INTERSECTION" => intersect_geometry,
    "UNION" => union_geometry,
)

fill_map = Dict(
    "EVENODD" => PolygonAlgorithms.EVEN_ODD,
    "NONZERO" => PolygonAlgorithms.NON_ZERO,
    "POSITIVE" => PolygonAlgorithms.POSITIVE,
    "NEGATIVE" => PolygonAlgorithms.NEGATIVE,
)

Base.@kwdef mutable struct ClipperResult
    caption::Int
    fill_rule::String
    result::String = "error"
    area_ratio::Float64 = NaN
    count_difference::Float64 = Inf
end

function execute_clipper_test(test::PolygonTest, subjects::Vector, clips::Vector)
    func = func_map[test.clip_type]
    fill_rule = fill_map[test.fill_rule]
    regions = func(subjects, clips...; fill_rule=fill_rule)
    area_ratio = calculate_area_diff_ratio(test.solution_area, regions)
    count_regions = length(regions)
    count_diff = test.solution_count > 0 ? (test.solution_count - count_regions) : 0
    area_ratio, count_diff
end

@testset "Clipper2" begin
    tests = load_tests(joinpath(@__DIR__, "clipper.txt"))
    results = ClipperResult[]
    for (idx, poly_test) in enumerate(tests)
        @testset "caption $(poly_test.caption)" begin
            record = ClipperResult(caption=poly_test.caption, fill_rule=poly_test.fill_rule)
            push!(results, record)
            subjects = map(Polygon, map(convert_path, poly_test.subjects))
            clips = isempty(poly_test.clips) ? empty(subjects) : map(Polygon, map(convert_path, poly_test.clips))
            area_ratio, count_diff = execute_clipper_test(poly_test, subjects, clips)
            threshold = area_threshold(poly_test)
            @test abs(area_ratio) <= threshold
            record.result = string(abs(area_ratio) <= threshold)
            record.area_ratio = area_ratio
            record.count_difference = count_diff
        end
    end
end
