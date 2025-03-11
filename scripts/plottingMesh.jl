using Plots

function plot_mesh(nodes, connectivity)
    plot(legend = false, aspect_ratio = :equal, grid = false, framestyle = :none)

    # Plot each triangle
    for element in eachrow(connectivity)
        x = [
            nodes[element[1], 1],
            nodes[element[2], 1],
            nodes[element[3], 1],
            nodes[element[1], 1],
        ]
        y = [
            nodes[element[1], 2],
            nodes[element[2], 2],
            nodes[element[3], 2],
            nodes[element[1], 2],
        ]
        plot!(x, y, lw = 1, linecolor = :black)
    end

    # Plot nodes and indices
    scatter!(nodes[:, 1], nodes[:, 2], markersize = 4, color = :red)
    for (i, (x, y)) in enumerate(eachrow(nodes))
        annotate!(x, y, text(string(i), 8, :blue))
    end

    display(plot!)
end
