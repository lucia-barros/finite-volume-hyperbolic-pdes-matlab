function F = godunov_flux_burgers(uL, uR)
%GODUNOV_FLUX_BURGERS Flujo exacto de Godunov para Burgers.
%
%   La ecuación es:
%
%       u_t + (u^2/2)_x = 0
%
%   uL y uR pueden ser escalares o vectores del mismo tamaño.

    if ~isequal(size(uL), size(uR))
        error('uL y uR deben tener el mismo tamaño.');
    end

    flux = @(u) 0.5*u.^2;

    F = zeros(size(uL));

    %% Caso 1: rarefacción, uL <= uR

    rarefaction = (uL <= uR);

    % Rarefacción situada completamente a la derecha de cero
    index = rarefaction & (uL >= 0);
    F(index) = flux(uL(index));

    % Rarefacción situada completamente a la izquierda de cero
    index = rarefaction & (uR <= 0);
    F(index) = flux(uR(index));

    % Cuando uL < 0 < uR, el flujo es f(0) = 0.
    % F ya está inicializado a cero.

    %% Caso 2: shock, uL > uR

    shock = ~rarefaction;

    % Velocidad de Rankine-Hugoniot:
    % s = (f(uL)-f(uR))/(uL-uR) = (uL+uR)/2

    shockSpeed = 0.5*(uL + uR);

    % Shock que se desplaza hacia la derecha
    index = shock & (shockSpeed >= 0);
    F(index) = flux(uL(index));

    % Shock que se desplaza hacia la izquierda
    index = shock & (shockSpeed < 0);
    F(index) = flux(uR(index));

end