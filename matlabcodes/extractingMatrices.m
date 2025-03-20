Mass_DIAG_test = sparse([1, 2, 3], [1, 2, 3], [1 + 2i, 3 + 4i, 5 + 6i], 3, 3);

% Extract row indices, column indices, and values
[row, col, val] = find(cell2mat(Mass_DIAG));
% Separate real and imaginary parts
val_real = real(val);
val_imag = imag(val);
% Combine into a matrix for export
export_data = [row, col, val_real, val_imag];
% Export to a text file
dlmwrite('E:/Git/EnrichedFiniteElements.jl/test/testdata/MassMatrixEnriched_noFrequencies.txt', export_data, 'delimiter', ' ', 'precision', '%.15f');


% dlmwrite('E:/Git/EnrichedFiniteElements.jl/test/testdata/Mass_DIAG.txt', full(cell2mat(Mass_DIAG(end,end-1))), 'delimiter', ' ', 'precision', '%.15f');