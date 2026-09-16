function [x, u, t, diagnostics] = solve_burgers_muscl( ...
    N, xmin, xmax, tFinal, CFL, ...
    caseName, limiterName, boundaryName)
%SOLVE_BURGERS_MUSCL Resuelve Burgers mediante MUSCL-Hancock.
%
%   El método utiliza:
%       - reconstrucción lineal por celda;
%       - limitadores TVD;
%       - predictor de Hancock;
%       - flujo exacto de Godunov;
%       - actualización conservativa.

    arguments
        N (1,1) double {mustBeInteger, mustBePositive}
        xmin (1,1) double
        xmax (1,1) double
        tFinal (1,1) double {mustBeNonnegative}
        CFL (1,1) double {mustBePositive}
        caseName (1,1) string
        limiterName (1,1) string
        boundaryName (1,1) string = "transmissive"
    end

    if xmax <= xmin
        error('Debe cumplirse xmax > xmin.');
    end

    if CFL > 1
        warning('Para MUSCL-Hancock se recomienda CFL <= 1.');
    end

    %% Malla

    dx = (xmax - xmin)/N;

    x = linspace( ...
        xmin + dx/2, ...
        xmax - dx/2, ...
        N).';

    %% Condición inicial

    u = initial_condition(x, caseName);

    t = 0;

    %% Diagnósticos

    initialMass = dx*sum(u);

    diagnostics.time = 0;
    diagnostics.mass = initialMass;
    diagnostics.expectedMass = initialMass;
    diagnostics.TV = compute_total_variation(u, boundaryName);

    expectedMass = initialMass;

    %% Bucle temporal

    while t < tFinal

        % Para Burgers, la velocidad característica es f'(u)=u.

        maximumSpeed = max(abs(u));

        if maximumSpeed < 1e-14
            dt = tFinal - t;
        else
            dt = CFL*dx/maximumSpeed;
            dt = min(dt, tFinal - t);
        end

        %% 1. Celdas fantasma transmisivas

        switch lower(boundaryName)

            case "transmissive"
        
                uExtended = [
                    u(1);
                    u;
                    u(end)
                ];
        
            case "periodic"
        
                uExtended = [
                    u(end);
                    u;
                    u(1)
                ];
        
            otherwise
        
                error('Condición de contorno desconocida: %s', ...
                    boundaryName);
        
        end

        %% 2. Diferencias izquierda y derecha

        deltaMinus = ...
            uExtended(2:end-1) - uExtended(1:end-2);

        deltaPlus = ...
            uExtended(3:end) - uExtended(2:end-1);

        %% 3. Pendiente limitada

        sigma = slope_limiter( ...
            deltaMinus, deltaPlus, limiterName);

        %% 4. Reconstrucción MUSCL

        uLeft = u - 0.5*sigma;
        uRight = u + 0.5*sigma;

        %% 5. Predictor de Hancock

        fluxLeft = 0.5*uLeft.^2;
        fluxRight = 0.5*uRight.^2;

        temporalCorrection = ...
            0.5*(dt/dx)*(fluxRight - fluxLeft);

        uLeftHalf = uLeft - temporalCorrection;
        uRightHalf = uRight - temporalCorrection;

        %% 6. Estados en las interfaces

        % Hay N+1 interfaces.
        %
        % En cada interfaz:
        %   estado izquierdo = borde derecho de la celda izquierda
        %   estado derecho   = borde izquierdo de la celda derecha

        switch lower(boundaryName)
        
            case "transmissive"
        
                interfaceLeftStates = [
                    u(1);
                    uRightHalf
                ];
        
                interfaceRightStates = [
                    uLeftHalf;
                    u(end)
                ];
        
            case "periodic"
        
                % La interfaz izquierda conecta la última celda
                % con la primera.
        
                interfaceLeftStates = [
                    uRightHalf(end);
                    uRightHalf
                ];
        
                interfaceRightStates = [
                    uLeftHalf;
                    uLeftHalf(1)
                ];
        
        end
        %% 7. Flujo de Godunov

        F = godunov_flux_burgers( ...
            interfaceLeftStates, interfaceRightStates);

        %% 8. Actualización conservativa

        u = u ...
            - (dt/dx)*(F(2:end) - F(1:end-1));

        %% 9. Avance temporal

        t = t + dt;

        %% 10. Balance de masa esperado

        expectedMass = expectedMass ...
            - dt*(F(end) - F(1));

        %% 11. Guardamos diagnósticos

        diagnostics.time(end+1,1) = t;
        diagnostics.mass(end+1,1) = dx*sum(u);
        diagnostics.expectedMass(end+1,1) = expectedMass;
        diagnostics.TV(end+1,1) = ...
            compute_total_variation(u, boundaryName);

    end
end