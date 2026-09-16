function q_new = step_ader4(q, dx, dt, lambda)
% STEP_ADER4
% ADER de cuarto orden para la ecuación:
%
%       q_t + lambda*q_x = 0
%
% con lambda > 0 y condiciones periódicas.
%
% El método emplea:
%   - reconstrucción cúbica a partir de cuatro medias celulares;
%   - derivadas espaciales en la interfaz;
%   - Cauchy-Kowalewski;
%   - integración temporal del flujo hasta cuarto orden.
%
% IMPORTANTE:
% Esta implementación es para soluciones suaves. La reconstrucción
% cúbica no está limitada y no debe aplicarse directamente a shocks.

    if lambda <= 0
        error(['Esta versión didáctica de ADER4 está implementada ', ...
               'para lambda > 0. Usa lambda = 1.']);
    end

    % ===============================================================
    % 1. Valores de las cuatro celdas del stencil
    %
    % Para la interfaz x_{i+1/2}, la reconstrucción se realiza
    % utilizando las medias:
    %
    % q_{i-1}, q_i, q_{i+1}, q_{i+2}
    % ===============================================================

    q_im1 = circshift(q, 1);
    q_i   = q;
    q_ip1 = circshift(q, -1);
    q_ip2 = circshift(q, -2);

    % ===============================================================
    % 2. Reconstrucción cúbica en la interfaz x_{i+1/2}
    %
    % Se obtiene:
    %   q*
    %   q_x*
    %   q_xx*
    %   q_xxx*
    %
    % Las expresiones se derivan imponiendo que el polinomio cúbico
    % reproduzca las medias celulares de las cuatro celdas del stencil.
    % ===============================================================

    q_star = ...
        (-q_im1 + 7*q_i + 7*q_ip1 - q_ip2) / 12;

    qx_star = ...
        (q_im1 - 15*q_i + 15*q_ip1 - q_ip2) ...
        / (12 * dx);

    qxx_star = ...
        (q_im1 - q_i - q_ip1 + q_ip2) ...
        / (2 * dx^2);

    qxxx_star = ...
        (-q_im1 + 3*q_i - 3*q_ip1 + q_ip2) ...
        / dx^3;

    % ===============================================================
    % 3. Procedimiento de Cauchy-Kowalewski
    %
    % Para:
    %
    % q_t + lambda*q_x = 0
    %
    % se tiene:
    %
    % q_t   = -lambda*q_x
    % q_tt  =  lambda^2*q_xx
    % q_ttt = -lambda^3*q_xxx
    % ===============================================================

    qt_star = -lambda * qx_star;

    qtt_star = lambda^2 * qxx_star;

    qttt_star = -lambda^3 * qxxx_star;

    % ===============================================================
    % 4. Flujo medio temporal
    %
    % q(tau) ≈ q* + tau*q_t*
    %             + tau^2/2*q_tt*
    %             + tau^3/6*q_ttt*
    %
    % Promediando entre 0 y dt:
    %
    % q_bar = q*
    %       + dt/2*q_t*
    %       + dt^2/6*q_tt*
    %       + dt^3/24*q_ttt*
    % ===============================================================

    q_time_average = ...
        q_star ...
        + 0.5 * dt * qt_star ...
        + (dt^2 / 6) * qtt_star ...
        + (dt^3 / 24) * qttt_star;

    % Para advección lineal:
    %
    % f(q) = lambda*q
    %
    flux_right = lambda * q_time_average;

    % F_{i-1/2}
    flux_left = circshift(flux_right, 1);

    % ===============================================================
    % 5. Actualización conservativa
    % ===============================================================

    q_new = q ...
        - (dt / dx) * (flux_right - flux_left);
end