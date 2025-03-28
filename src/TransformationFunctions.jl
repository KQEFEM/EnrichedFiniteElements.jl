module TransformationFunctions
using PolygonOps
export all
"""
Corrects the orientation of triangle nodes to ensure that triangle normals point in a consistent direction.

This function calculates a 2D cross product based on the coordinates of the triangle nodes.
If the cross product is negative, it indicates an incorrect orientation of the triangle, and the function performs the following adjustment:
1. Permutes the order of the second and third nodes in `triangle_nodes` to correct the normal direction.
2. Applies the same permutation to `connect_triangle`.

The function is designed to work with triangle nodes represented as a matrix where each row corresponds to a node and columns represent coordinates. The first two columns of `triangle_nodes` are used for the cross product calculation, assumed to be x and y coordinates, respectively.  `connect_triangle` is assumed to be a related indexing structure that needs to be updated in sync with the `triangle_nodes` permutation.

# Arguments
- `triangle_nodes`: A matrix representing the nodes of a triangle. It is expected to have at least 3 rows and 2 columns. The first column is used as the x-coordinate and the second column as the y-coordinate for the normal direction check.  This function modifies the order of rows within this matrix if the orientation is corrected.
- `connect_triangle`: An array (vector or matrix) related to triangle connectivity. Its elements are reordered based on the same permutation applied to `triangle_nodes` when the orientation is corrected. The exact structure and meaning depend on the broader application, but it should have a compatible row structure to be indexed by `[1, 3, 2]`.

# Returns
- `triangle_nodes`: The (potentially) reoriented `triangle_nodes` matrix. If the cross product was non-negative, the returned matrix will be the same as the input. If the orientation was corrected, the rows corresponding to the triangle nodes will be permuted within the returned matrix.
- `connect_triangle`: The (potentially) reordered `connect_triangle` array, adjusted in sync with `triangle_nodes` if reorientation occurred.

# Notes
- The function uses a 2D cross product calculation based on the first two columns of `triangle_nodes` to determine triangle orientation in 2D.
- The `!` at the end of the function name is a Julia convention indicating that the function may modify its arguments, although technically in this implementation, it returns new arrays with reordered data rather than modifying the original arrays in-place in the caller's scope.
- It is crucial to understand the role of `connect_triangle` within your specific application to ensure correct usage of this function.
"""
function correct_triangle_orientation!(triangle_nodes, connect_triangle)
    cross_product =
        (triangle_nodes[2, 1] - triangle_nodes[1, 1]) *
        (triangle_nodes[3, 2] - triangle_nodes[1, 2]) -
        (triangle_nodes[2, 2] - triangle_nodes[1, 2]) *
        (triangle_nodes[3, 1] - triangle_nodes[1, 1])

    if cross_product < 0  # Ensures all of the normals point in the correct direction
        ABC = triangle_nodes[[1, 3, 2], :]
        connect_temp = connect_triangle[[1, 3, 2]] # Assuming connect_triangle is a matrix or vector that can be indexed like this
        connect_triangle = connect_temp
        triangle_nodes = ABC # Update triangle_nodes to the permuted ABC (which was based on original triangle_nodes)
    end
    return triangle_nodes, connect_triangle
end


function Gradients_Larson(x, y)
    # Ensure x and y are vectors of length 3 (triangle vertices)
    if length(x) != 3 || length(y) != 3
        error("x and y must have exactly 3 elements (triangle vertices).")
    end

    # Compute the area of the triangle using the shoelace formula
    tri_area = 0.5 * abs(x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2]))

    # Check if the area is zero (degenerate triangle)
    if tri_area == 0
        error("The triangle area is zero (degenerate triangle).")
    end

    # Compute gradients b and c
    b = reshape([y[2] - y[3]; y[3] - y[1]; y[1] - y[2]] / (2 * tri_area), 3, 1)
    c = reshape([x[3] - x[2]; x[1] - x[3]; x[2] - x[1]] / (2 * tri_area), 3, 1)
    #! The reshape is to ensure it is a matrix so it dits with the numerical integration
    gradients = (b, c)
    return tri_area, gradients
end

function arb_triangle_to_ref(nodes)
    x1, y1 = nodes[1, 1], nodes[1, 2]
    x2, y2 = nodes[2, 1], nodes[2, 2]
    x3, y3 = nodes[3, 1], nodes[3, 2]

    X21 = x2 - x1
    X31 = x3 - x1 
    Y21 = y2 - y1
    Y31 = y3 - y1

    return X21,X31,Y21,Y31,x1,y1
end 

end
