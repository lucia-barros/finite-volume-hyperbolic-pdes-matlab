function [x, u, t, diagnostics] = solve_burgers_godunov( ...
    N, xmin, xmax, tFinal, CFL, caseName, boundaryName)
%SOLVE_BURGERS_GODUNOV Resuelve Burgers con Godunov de primer orden.
%
%   Se utilizan:
%       - volúmenes finitos;
%       - flujo exacto de Godunov;
%       - paso temporal adaptativo mediante CFL;
%       - condiciones de contorno transmisivas.

    arguments
        N (1,1) double {mustBeInteger, mustBePositive}
        xmin (1,1) double
        xmax (1,1) double
        tFinal (1,1) double {mustBeNonnegative}
        CFL (1,1) double {mustBePositive}
        caseName (1,1) string
        boundaryName (1,1) string = "transmissive"
    end

    if xmax <= xmin
        error('Debe cumplirse xmax > xmin.');
    end

    if CFL > 1
        warning('Para Godunov se recomienda CFL <= 1.');
    end

    %% Malla de volúmenes finitos

    dx = (xmax - xmin)/N;

    % Centros de las celdas
    x = linspace( ...
        xmin + dx/2, ...
        xmax - dx/2, ...
        N).';

    %% Condición inicial

    u = initial_condition(x, caseName);

    t = 0;

    %% Inicialización de diagnósticos

    diagnostics.time = 0;
    diagnostics.mass = dx*sum(u);
    diagnostics.expectedMass = diagnostics.mass;
    diagnostics.TV = compute_total_variation(u, boundaryName);

    expectedMass = diagnostics.mass;

    %% Bucle temporal

    while t < tFinal

        % Máxima velocidad característica para Burgers:
        % f'(u) = u

        maximumSpeed = max(abs(u));

        if maximumSpeed < 1e-14
            dt = tFinal - t;
        else
            dt = CFL*dx/maximumSpeed;
            dt = min(dt, tFinal - t);
        end

        %% Condiciones de contorno transmisivas

        % Añadimos una celda fantasma a cada lado.
        % En contorno transmisivo se copia el valor más próximo.

        switch lower(boundaryName)
        
            case "transmissive"
        
                uGhost = [
                    u(1);
                    u;
                    u(end)
                ];
        
            case "periodic"
        
                uGhost = [
                    u(end);
                    u;
                    u(1)
                ];
        
            otherwise
        
                error('Condición de contorno desconocida: %s', ...
                    boundaryName);
        
        end

        %% Flujos en todas las interfaces

        uLeft = uGhost(1:end-1);
        uRight = uGhost(2:end);

        F = godunov_flux_burgers(uLeft, uRight);

        %% Actualización conservativa

        u = u ...
            - (dt/dx)*(F(2:end) - F(1:end-1));

        %% Avance temporal

        t = t + dt;

        %% Balance conservativo esperado

        % Para un dominio abierto, la masa no tiene por qué ser constante.
        % Su variación debe estar explicada por los flujos de los extremos:
        %
        % M^{n+1} = M^n - dt*(F_derecha - F_izquierda)

        expectedMass = expectedMass ...
            - dt*(F(end) - F(1));

        %% Guardamos diagnósticos

        diagnostics.time(end+1,1) = t;
        diagnostics.mass(end+1,1) = dx*sum(u);
        diagnostics.expectedMass(end+1,1) = expectedMass;
        diagnostics.TV(end+1,1) = ...
            compute_total_variation(u, boundaryName);

    end
end