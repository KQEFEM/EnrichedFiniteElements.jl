module MeshCreation

using Gmsh, FileIO
using Makie, CairoMakie
using GeometryBasics

export MeshData
# Define a structure to hold the mesh data
struct MeshData
    nodes::Matrix{Float64}        # Node coordinates (nx2 matrix)
    elem_nodes::Vector{Vector{Int}} # Element nodes (list of node indices for each element)
    elem_types::Vector{Int}       # Element types (Gmsh element type codes)
    connectivity::Matrix{Int}     # Triangle connectivity (mx3 matrix)
    boundary_idx::Vector{Int}   # Unique boundary node indices
    boundary_edges::Matrix{Int}   # Edge connectivity (bx2 matrix)
end

function rectangle_domain(
    domain;
    mesh_size::Real = 0.5,
    filename::String = "./mesh/rectangle_mesh.msh",
)
    # Create the directory if it doesn't exist
    dir_path = dirname(filename)  # Get the directory part of the filename
    if !isdir(dir_path)
        mkpath(dir_path)  # Create the directory and any necessary parent directories
    end

    (x_min, x_max), (y_min, y_max) = domain

    # Initialize Gmsh
    Gmsh.initialize()
    gmsh.model.add(filename)

    # Define the rectangle's corners
    p1 = gmsh.model.geo.addPoint(x_min, y_min, 0.0)
    p2 = gmsh.model.geo.addPoint(x_max, y_min, 0.0)
    p3 = gmsh.model.geo.addPoint(x_max, y_max, 0.0)
    p4 = gmsh.model.geo.addPoint(x_min, y_max, 0.0)

    # Add lines to form the rectangle
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Create a curve loop and a surface
    curve_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surface = gmsh.model.geo.addPlaneSurface([curve_loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Set mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    # Generate the mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh to a file
    gmsh.write(filename)

    # Get node information
    _, node_coords, node_tags = gmsh.model.mesh.getNodes()
    nodes = reshape(node_coords, 3, :)'[:, 1:2]  # Extract (x, y) coordinates

    # Get element connectivity
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(2)  # Get 2D elements

    # Extract triangle connectivity
    triangle_list = []
    for (etype, conn) in zip(elem_types, elem_nodes)
        if etype == 2  # Triangle element type
            push!(triangle_list, reshape(conn, 3, :)')
        end
    end
    connectivity = vcat(triangle_list...)  # Convert list of arrays to a matrix

    # Extract boundary edges
    boundary_edges = []
    edge_entities = gmsh.model.getEntities(1)  # Get 1D elements (edges)

    for (dim, tag) in edge_entities
        elem_type, _, elem_nodes = gmsh.model.mesh.getElements(dim, tag)
        for (etype, conn) in zip(elem_type, elem_nodes)
            if etype == 1  # Edge element type
                push!(boundary_edges, reshape(conn, 2, :)')
            end
        end
    end

    boundary_edges = vcat(boundary_edges...)  # Convert list to a matrix
    boundary_nodes = unique(vcat(boundary_edges...))  # Get unique boundary node indices

    # Finalize Gmsh
    Gmsh.finalize()

    return MeshData(
        nodes,
        elem_nodes,
        elem_types,
        connectivity,
        boundary_nodes,
        boundary_edges,
    )
end



function extract_mesh_info(filename::String = "rectangle_mesh.msh")
    Gmsh.initialize()

    # Read the generated mesh file (assumed to be "rectangle_mesh.msh")
    gmsh.open(filename)

    # Get node information
    _, node_coords, _ = gmsh.model.mesh.getNodes()
    nodes = reshape(node_coords, 3, :)'[:, 1:2]  # Reshape and extract (x, y)

    # Get element connectivity
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(2)  # 2D elements

    # Extract triangle elements if they exist (type 2 is typically for triangles)
    triangles = []
    for (etype, connectivity) in zip(elem_types, elem_nodes)
        if etype == 2  # Triangle element type
            push!(triangles, reshape(connectivity, 3, :)')
        end
    end
    # Get boundary edges (Corrected)
    boundary_entities = gmsh.model.mesh.getBoundary([2], true, false, true) # Get boundary edges of the surface
    boundary_edges = Matrix{Int}(undef, length(boundary_entities), 2)
    boundary_nodes_indices = unique(Int[])
    for (i, entity) in enumerate(boundary_entities)
        _, edge_nodes, _ = gmsh.model.mesh.getNodes(1, entity[2]) #Get nodes of the edge
        boundary_edges[i, :] = edge_nodes
        push!(boundary_nodes_indices, edge_nodes...)
    end
    boundary_nodes = nodes[unique(boundary_nodes_indices), :] #Get coordinates of the boundary nodes


    Gmsh.finalize()

    return nodes, triangles
end

function plot_mesh(mesh)

    nodes = mesh.nodes
    connectivity = mesh.connectivity
    boundary_index = mesh.boundary_idx
    boundary_edges = mesh.boundary_edges

    # Create a figure and an axis *together*
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect()) #Aspect ratio is important

    # Plot the elements (triangles)
    n_triangles = length(connectivity) ÷ 3 #Number of triangles. Each triangle has 3 nodes
    for i = 1:n_triangles
        triangle_nodes = nodes[connectivity[i, :], :]
        # triangle_nodes = @view connectivity[(i-1)*3+1:i*3] #Get the 3 nodes for the triangle
        x_vals = triangle_nodes[:, 1] #Correct indexing
        y_vals = triangle_nodes[:, 2] #Correct indexing
        lines!(ax, [x_vals..., x_vals[1]], [y_vals..., y_vals[1]], color = :black) # Plot into the axis
    end

    # Plot the boundary edges (highlighted) - Corrected
    for edge in eachrow(boundary_edges)
        x_vals = nodes[edge, 1]
        y_vals = nodes[edge, 2]
        lines!(ax, x_vals, y_vals, color = :red, linewidth = 2, label = "Boundary Edges") # Plot the boundary edges
    end


    # Legend(fig[1,2], ax) #Add legend

    return fig #Display the figure
end

end  # module MeshCreation
