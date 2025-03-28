
using DelimitedFiles
using SparseArrays

# Load the data from the text file
data = readdlm("test/testdata/pdtq_space_time.txt")

# Extract row indices, column indices, and values
rows = Int.(data[:, 1])
cols = Int.(data[:, 2])
vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts

# Reconstruct the sparse matrix
exact = sparse(rows, cols, vals);

println(norm((final_matrix - conj(Mass_matlab))))

