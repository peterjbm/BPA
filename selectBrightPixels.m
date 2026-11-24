function sel = selectBrightPixels(pw, hits, K, ntiles)
% Linear indices of the K brightest pixels of power image pw (pixels with
% hits == 0 excluded), the same number from each of ntiles along-track
% (column) tiles so bright features along a whole strip are represented.
[~, Nx] = size(pw);
pw(hits == 0) = 0;
edges = round(linspace(0, Nx, ntiles+1));
per = ceil(K/ntiles);
sel = [];
for t = 1:ntiles
    cols = edges(t)+1:edges(t+1);
    tile = pw(:, cols);
    [~, ord] = sort(tile(:), 'descend');
    ord = ord(1:min(per, nnz(tile)));
    [r, c] = ind2sub(size(tile), ord);
    sel = [sel; sub2ind(size(pw), r(:), reshape(cols(c),[],1))]; %#ok<AGROW>
end
end
