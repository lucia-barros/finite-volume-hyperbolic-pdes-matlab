function slope = compute_slopes(q, dx, method)
% COMPUTE_SLOPES Calcula las pendientes de reconstrucción periódicas.
%
% q      : vector fila con las medias celulares
% dx     : tamaño de celda
% method : 'eno', 'minmod' o 'centered'
%
% slope  : pendiente reconstruida en cada celda

    arguments
        q (1,:) double
        dx (1,1) double {mustBePositive}
        method char = 'eno'
    end

    % Valores vecinos usando condiciones periódicas
    q_left  = circshift(q, 1);
    q_right = circshift(q, -1);

    % Gradientes a izquierda y derecha
    delta_left  = (q - q_left) / dx;
    delta_right = (q_right - q) / dx;

    switch lower(method)

        case 'eno'
            % Reconstrucción ENO lineal de Toro, ecuación (20.19):
            % se elige el gradiente de menor valor absoluto.
            slope = delta_left;

            use_right = abs(delta_right) < abs(delta_left);
            slope(use_right) = delta_right(use_right);

        case 'minmod'
            % Limitador minmod:
            % si las pendientes cambian de signo, se toma pendiente cero.
            slope = zeros(size(q));

            same_sign = delta_left .* delta_right > 0;

            slope(same_sign) = sign(delta_left(same_sign)) .* ...
                min(abs(delta_left(same_sign)), ...
                    abs(delta_right(same_sign)));

        case 'centered'
            % Pendiente centrada, adecuada para soluciones suaves.
            % No es no-oscilatoria cerca de discontinuidades.
            slope = 0.5 * (delta_left + delta_right);

        otherwise
            error(['Método de reconstrucción desconocido. ', ...
                   'Usa ''eno'', ''minmod'' o ''centered''.']);
    end
end