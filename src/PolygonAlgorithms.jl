module PolygonAlgorithms

import Base: push!, insert!, iterate, merge!, contains
import Base: length, show

include("definitions.jl")
include("orientation.jl")
include("linked_list.jl")
include("point_set.jl")

include("moments.jl")
include("point_in_polygon.jl")
include("intersect.jl")

include("intersect_poly.jl")
include("intersect_convex.jl")
include("convex_hull.jl")

include("deprecations.jl")

export Point2D, Polygon2D, Segment2D, Line2D
export get_orientation, Orientation, on_segment
export area_polygon, first_moment, centroid_polygon, is_counter_clockwise, is_clockwise
export do_intersect, intersect_geometry, intersect_edges
export intersect_convex
export convex_hull

end