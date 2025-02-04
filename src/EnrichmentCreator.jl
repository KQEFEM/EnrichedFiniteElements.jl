module EnrichmentCreator


function create_wavenumbers(x_enrichments::Real, y_enrichments::Real)
    """ Creates the wavenumber matrix to use """
    wavenumbers = collect(Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments))
    return reshape(wavenumbers,:,1)
    
end 

function wavenumber_creation(wavenumbers_ansatz::Matrix{Real}, y_enrichments::Matrix{Real})
    


end 

function combine_wavenumber_with_all_nodes(matrix_1, matrix_2)
    n_wavenumbers = size(matrix_1, 1)
    n_nodes = size(matrix_2, 1)

    result = Array{Tuple{Vector{eltype(wavenumbers)}, Vector{eltype(nodes)}}}(undef, n_wavenumbers * n_nodes)

    for k_iter in 1:n_wavenumbers
        wavenumber = wavenumbers[k_iter, :]  # Current wavenumber row
        for node_iter in 1:n_nodes
          result[(k_iter-1)*n_nodes + node_iter] = (wavenumber, nodes[node_iter, :])
        end
    end

    return result
end

end 