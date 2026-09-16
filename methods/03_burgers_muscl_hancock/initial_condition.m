function u0 = initial_condition(x, caseName)
%INITIAL_CONDITION Define las condiciones iniciales para Burgers.
%
%   u0 = INITIAL_CONDITION(x, caseName)
%
%   caseName puede ser:
%       "shock"
%       "rarefaction"
%       "smooth"

    arguments
        x (:,1) double
        caseName (1,1) string
    end

    switch lower(caseName)

        case "shock"
            % Problema de Riemann:
            % uL = 1, uR = 0
            u0 = double(x < 0);

        case "rarefaction"
            % Problema de Riemann:
            % uL = 0, uR = 1
            u0 = double(x >= 0);

        case "smooth"
            % Condición inicial suave.
            % El término constante mantiene u positiva.
            u0 = 0.5 + 0.25*sin(pi*x);

        otherwise
            error("Caso desconocido: %s", caseName);

    end
end