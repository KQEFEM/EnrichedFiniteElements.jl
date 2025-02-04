function combine_wavenumber_with_all_nodes(wavenumbers, nodes)
    n_wavenumbers = size(wavenumbers, 1)
    n_nodes = size(nodes, 1)

    result = Array{Tuple{Vector{eltype(wavenumbers)}, Vector{eltype(nodes)}}}(undef, n_wavenumbers * n_nodes)

    for k_iter in 1:n_wavenumbers
        wavenumber = wavenumbers[k_iter, :]  # Current wavenumber row
        for node_iter in 1:n_nodes
          result[(k_iter-1)*n_nodes + node_iter] = (wavenumber, nodes[node_iter, :])
        end
    end

    return result
end


# Example usage:
n_wavenumbers = 3
n_nodes = 5
wavenumbers = rand(n_wavenumbers, 3)
nodes = rand(n_nodes, 2)
wavenumbers = [-1,0,1]; nodes =[-1,0,1];
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)
c = combine_wavenumber_with_all_nodes(combined, combined)

# Print or process the combinations as needed
for (wavenumber, node) in combined
    println("Wavenumber: $wavenumber, Node: $node")
end

println("Number of combinations: ", length(combined)) #Should be n_wavenumbers * n_nodes

# Accessing a specific combination:
wavenumber_2_node_3 = combined[(2-1)*n_nodes + 3] #Wavenumber 2 with node 3
println("Wavenumber 2, Node 3: Wavenumber $(wavenumber_2_node_3[1]), Node $(wavenumber_2_node_3[2])")

wavenumber_idx = 2
node_idx = 4
wavenumber_node_combination = combined[(wavenumber_idx-1)*n_nodes + node_idx]
println("Wavenumber $wavenumber_idx, Node $node_idx: Wavenumber $(wavenumber_node_combination[1]), Node $(wavenumber_node_combination[2])")


wavenumbers = collect(Iterators.product(-2:1, -1:1))
wavenumbers = reshape(wavenumbers,:,1)
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)


