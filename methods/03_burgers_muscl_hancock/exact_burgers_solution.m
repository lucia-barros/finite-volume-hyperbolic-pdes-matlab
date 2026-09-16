function uExact = exact_burgers_solution(x, t, caseName)
%EXACT_BURGERS_SOLUTION Soluciones exactas de referencia.
%
%   Incluye shocks y rarefacciones de Riemann, y una solución suave
%   obtenida mediante características resolviendo la relación implícita
%   con el método de Newton antes de la formación del shock.

    arguments
        x (:,1) double
        t (1,1) double {mustBeNonnegative}
        caseName (1,1) string
    end

    switch lower(caseName)

        case "shock"

            % uL = 1, uR = 0
            %
            % Velocidad del shock:
            % s = (uL + uR)/2 = 0.5

            shockSpeed = 0.5;
            shockPosition = shockSpeed*t;

            uExact = double(x < shockPosition);

        case "rarefaction"

            if t == 0
                uExact = double(x >= 0);
                return;
            end

            % Para uL = 0 y uR = 1:
            %
            % u = 0,       x/t <= 0
            % u = x/t,     0 < x/t < 1
            % u = 1,       x/t >= 1

            similarityVariable = x/t;

            uExact = zeros(size(x));

            fan = similarityVariable > 0 ...
                & similarityVariable < 1;

            uExact(fan) = similarityVariable(fan);
            uExact(similarityVariable >= 1) = 1;

        case "smooth"

            % Condición inicial:
            %
            % u0(x) = 0.5 + 0.25 sin(pi x)
            %
            % Antes de la formación del shock, las características
            % satisfacen:
            %
            % x = xi + t*u0(xi)
            %
            % y la solución es:
            %
            % u(x,t) = u0(xi)
            %
            % La posición inicial xi se calcula mediante Newton.

            if t == 0
                uExact = 0.5 + 0.25*sin(pi*x);
                return;
            end

            % Aproximación inicial para Newton
            xi = x - 0.5*t;

            maximumIterations = 50;
            tolerance = 1e-13;
            converged = false;

            for iteration = 1:maximumIterations

                initialValue = ...
                    0.5 + 0.25*sin(pi*xi);

                functionValue = ...
                    xi + t*initialValue - x;

                functionDerivative = ...
                    1 + 0.25*pi*t*cos(pi*xi);

                correction = ...
                    functionValue./functionDerivative;

                xi = xi - correction;

                if max(abs(correction)) < tolerance
                    converged = true;
                    break;
                end
            end

            if ~converged
                warning(['Newton no alcanzó la tolerancia ', ...
                    'en la solución exacta suave.']);
            end

            uExact = 0.5 + 0.25*sin(pi*xi);

    end
end