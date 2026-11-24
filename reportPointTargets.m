function rep = reportPointTargets(img, grid, truth)
% Prints peak position error, -3 dB widths and peak level for each
% simulated point target (truth from simulateNavSARData).

A = abs(img.im);
Amax = max(A(:));
Tenu = nav_geodetic2enu(truth.targets_lat, truth.targets_lon, truth.targets_h, grid.ref);
Tp = Tenu - repmat(grid.center_enu, 1, size(Tenu,2));
t_al = grid.u_along' * Tp;  t_cr = grid.u_cross' * Tp;
da = mean(diff(grid.along_vec)); dc = mean(diff(grid.cross_vec));
fprintf('\n tgt   true along/cross (m)   error along/cross (m)   -3dB width along/cross (m)   peak (dB)\n');
rep = zeros(numel(t_al), 5);
for t = 1:numel(t_al)
    [~, ix] = min(abs(grid.along_vec - t_al(t)));
    [~, iy] = min(abs(grid.cross_vec - t_cr(t)));
    wx = max(1,ix-1):min(numel(grid.along_vec),ix+1);
    wy = max(1,iy-1):min(numel(grid.cross_vec),iy+1);
    sub = A(wy,wx); [m, k] = max(sub(:)); [r, q] = ind2sub(size(sub), k);
    py = wy(r); px = wx(q);
    % sub-pixel peak by parabolic fit
    ox = 0; oy = 0;
    if px > 1 && px < size(A,2)
        y3 = A(py, px-1:px+1); ox = 0.5*(y3(1)-y3(3))/(y3(1)-2*y3(2)+y3(3));
    end
    if py > 1 && py < size(A,1)
        y3 = A(py-1:py+1, px); oy = 0.5*(y3(1)-y3(3))/(y3(1)-2*y3(2)+y3(3));
    end
    ea = grid.along_vec(px) + ox*da - t_al(t);
    ec = grid.cross_vec(py) + oy*dc - t_cr(t);
    row = A(py, max(1,px-20):min(end,px+20)); col = A(max(1,py-20):min(end,py+20), px);
    wa = sum(row >= m/sqrt(2))*da; wc = sum(col >= m/sqrt(2))*dc;
    pk = 20*log10(m/Amax);
    fprintf(' %2d   %8.2f / %7.2f      %7.3f / %7.3f          %6.2f / %6.2f          %6.1f\n', ...
        t, t_al(t), t_cr(t), ea, ec, wa, wc, pk);
    rep(t,:) = [ea ec wa wc pk];
end
end
