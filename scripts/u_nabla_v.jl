using LinearAlgebra
using SparseArrays
using Revise
using EnrichedFiniteElements

const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator

#! Set up for the matrix creation
domain = ((0, 1), (0, 1))
num_nodes = 10

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges
function compute_u_nabla_v(nodes, elements, u, v)
    # nodes: N×d matrix (N nodes, d dimensions)
    # elements: connectivity matrix (triangles in 2D, tetrahedra in 3D)
    # u, v: vectors of function values at nodes

    N, d = size(nodes)
    ∇v = zeros(N, d)  # Gradient of v at each node

    # Compute gradients using finite element shape functions
    for elem in eachrow(elements)
        coords = nodes[elem, :]
        v_local = v[elem]

        # Compute gradient of shape functions
        A = hcat(coords[2:end, :] .- coords[1, :])  # Transformation matrix
        B = ones(d + 1)
        B[2:end] = 0
        grad_N = A \ B  # Solve for gradients

        # Compute gradient of v
        ∇v_local = grad_N' * v_local
        ∇v[elem, :] .+= ∇v_local / length(elem)  # Distribute contribution
    end

    # Compute u * ∇v
    u_nabla_v = u .* ∇v  # Element-wise multiplication

    return u_nabla_v
end

compute_u_nabla_v(nodes, connectivity)
