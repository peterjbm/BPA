function y = nav_sinc(x)
% sin(pi x)/(pi x), toolbox-free
y = ones(size(x)); nz = x ~= 0;
y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end
