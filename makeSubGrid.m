function sg = makeSubGrid(grid, sel)
% A 1 x K "grid" holding only pixels sel of grid, for fast backprojection
% onto a subset of pixels with bpNav (keeps the stripmap beam prefilter).
[~, col] = ind2sub(size(grid.X), sel);
sg = grid;
sg.X = grid.X(sel)'; sg.Y = grid.Y(sel)'; sg.Z = grid.Z(sel)';
sg.along_vec = reshape(grid.along_vec(col), 1, []);
sg.dpx = abs(mean(diff(grid.along_vec)));
sg.corners = [grid.X([1 end 1 end]); grid.Y([1 end 1 end]); grid.Z([1 end 1 end])];
end
