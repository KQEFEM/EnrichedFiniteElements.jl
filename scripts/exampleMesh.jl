using Revise
using EnrichedFiniteElements
using Makie, CairoMakie
using GeometryBasics

domain = ((0, 1), (0, 1))
num_nodes = 10
const mesh_create = EnrichedFiniteElements.MeshCreation

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges

f = mesh_create.plot_mesh(mesh)