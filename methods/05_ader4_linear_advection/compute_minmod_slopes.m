function slope = compute_minmod_slopes(q, dx)
% COMPUTE_MINMOD_SLOPES
% Calcula pendientes limitadas con minmod.
%
% q     : medias celulares
% dx    : tamaño de celda
% slope : pendiente reconstruida en cada celda

    q_left  = circshift(q, 1);
    q_right = circshift(q, -1);

    delta_left  = (q - q_left) / dx;
    delta_right = (q_right - q) / dx;

    slope = zeros(size(q));

    same_sign = delta_left .* delta_right > 0;

    slope(same_sign) = sign(delta_left(same_sign)) .* ...
        min(abs(delta_left(same_sign)), ...
            abs(delta_right(same_sign)));
end