module PolygonAlgorithms

import Base: contains, getindex, insert!, iterate, merge!, push!, reverse, ==
import Base: length, show

include("definitions.jl")
include("utils.jl")

# Data structures
include("data_structures/linked_list.jl")
include("data_structures/point_set.jl")
include("data_structures/segment_event.jl")

# Foundation algorithms
include("bounds.jl")
include("moments.jl")
include("orientation.jl")
include("intersect.jl")
include("point_in_polygon.jl")

# Complex algorithms
include("convex_hull.jl")
include("line_sweep.jl")
include("data_structures/polygon.jl")

include("graphs/face_computation.jl")
include("mappings.jl")

include("boolean.jl")

include("deprecations.jl")

export get_orientation, Orientation, on_segment
export bounds
export area_polygon, first_moment, centroid_polygon, is_counter_clockwise, is_clockwise
export do_intersect, intersect_geometry, intersect_edges
export convex_hull
export any_intersect
export difference_geometry, union_geometry, xor_geometry, intersect_convex

end