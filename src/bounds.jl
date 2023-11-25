"""
    bounds(polygon)

A rectangle which bounds the polygon, given as `(xmin, ymin, xmax, ymax)`.
"""
function bounds(polygon::Polygon2D)
    start = polygon[1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for (x, y) in polygon
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    (xmin, ymin, xmax, ymax)
end

"""
    bounds(polygons)

A rectangle which bounds all the polygons, given as `(xmin, ymin, xmax, ymax)`.
"""
function bounds(polygons::Vector{<:Polygon2D})
    start = polygons[1][1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for polygon in polygons
        x0, y0, x1, y1 = bounds(polygon)
        xmin = min(xmin, x0, x1)
        xmax = max(xmax, x0, x1)
        ymin = min(ymin, y0, y1)
        ymax = max(ymax, y0, y1)
    end
   (xmin, ymin, xmax, ymax)
end

function box(xmin::T, ymin::T, xmax::T, ymax::T) where T
    [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)]
end