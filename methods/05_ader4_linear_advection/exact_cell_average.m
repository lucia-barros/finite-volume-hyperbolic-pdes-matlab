function q_exact = exact_cell_average(x, dx, t, lambda)
% EXACT_CELL_AVERAGE
% Media celular exacta de:
%
% q(x,t) = sin(2*pi*(x-lambda*t))
%
% sobre celdas de tamaño dx centradas en x.

    factor = sin(pi * dx) / (pi * dx);

    q_exact = factor .* ...
        sin(2 * pi * (x - lambda * t));
end