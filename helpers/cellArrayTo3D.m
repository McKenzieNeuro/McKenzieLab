function B = cellArrayTo3D(A)

A = A(:)';
 B = cell2mat(permute(A,[3 1 2]));
 
end