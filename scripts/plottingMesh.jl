using Plots
function plot_mesh(nodes, connectivity)
    plt = plot(legend = false, aspect_ratio = :equal, grid = false, framestyle = :none)

    # Plot each triangle
    for element in eachrow(connectivity)
        x = [
            nodes[element[1], 1],
            nodes[element[2], 1],
            nodes[element[3], 1],
            nodes[element[1], 1],  # Closing the triangle
        ]
        y = [
            nodes[element[1], 2],
            nodes[element[2], 2],
            nodes[element[3], 2],
            nodes[element[1], 2],
        ]
        plot!(plt, x, y, lw = 1, linecolor = :black)
    end

    # Plot nodes
    scatter!(plt, nodes[:, 1], nodes[:, 2], markersize = 4, color = :red)

    # Annotate node indices with an offset for readability
    dx, dy = 0.02, 0.02  # Adjust the offsets based on your mesh scale
    for (i, (x, y)) in enumerate(eachrow(nodes))
        annotate!(plt, x + dx, y + dy, text(string(i), 8, :blue))
    end

    display(plt)
end

plot_mesh(nodes, connectivity)
