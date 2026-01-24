module PolygonAlgorithms

import Base: contains, insert!, iterate, merge!, push!, ==
import Base: length, show

include("definitions.jl")
include("utils.jl")
include("data_structures/linked_list.jl")
include("data_structures/point_set.jl")
include("data_structures/segment_event.jl")

include("orientation.jl")
include("bounds.jl")
include("convex_hull.jl")
include("intersect.jl")
include("moments.jl")

include("point_in_polygon.jl")
include("line_sweep.jl")

include("data_structures/polygon.jl")

include("polygon_boolean.jl")
include("deprecations.jl")

export get_orientation, Orientation, on_segment
export bounds, convex_hull
export area_polygon, first_moment, centroid_polygon, is_counter_clockwise, is_clockwise
export do_intersect, intersect_geometry, intersect_edges
export difference_geometry, union_geometry, xor_geometry 
export intersect_convex

end