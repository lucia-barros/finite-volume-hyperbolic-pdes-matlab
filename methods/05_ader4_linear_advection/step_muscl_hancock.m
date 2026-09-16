function q_new = step_muscl_hancock(q, dx, dt, lambda)
% STEP_MUSCL_HANCOCK
% Un paso temporal MUSCL-Hancock para:
%
%       q_t + lambda*q_x = 0
%
% con condiciones periódicas.

    % 1. Pendientes limitadas
    slope = compute_minmod_slopes(q, dx);

    % 2. Reconstrucción lineal dentro de cada celda
    q_left  = q - 0.5 * dx * slope;
    q_right = q + 0.5 * dx * slope;

    % 3. Predictor de Hancock hasta t^{n+1/2}
    %
    % q_t = -lambda*q_x
    %
    q_left_half = q_left - 0.5 * dt * lambda * slope;
    q_right_half = q_right - 0.5 * dt * lambda * slope;

    % 4. Problema de Riemann upwind en cada interfaz
    if lambda >= 0
        q_interface = q_right_half;
    else
        q_interface = circshift(q_left_half, -1);
    end

    flux_right = lambda * q_interface;
    flux_left = circshift(flux_right, 1);

    % 5. Actualización de volúmenes finitos
    q_new = q ...
        - (dt / dx) * (flux_right - flux_left);
end