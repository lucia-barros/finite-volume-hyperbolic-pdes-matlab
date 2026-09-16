function sigma = slope_limiter(deltaMinus, deltaPlus, limiterName)
%SLOPE_LIMITER Calcula pendientes limitadas para MUSCL.
%
%   deltaMinus = u_i - u_{i-1}
%   deltaPlus  = u_{i+1} - u_i
%
%   Limitadores disponibles:
%       "minmod"
%       "superbee"
%       "vanleer"
%       "vanalbada"

    if ~isequal(size(deltaMinus), size(deltaPlus))
        error('Las dos diferencias deben tener el mismo tamaño.');
    end

    limiterName = lower(string(limiterName));

    switch limiterName

        case "minmod"

            sigma = minmod_function(deltaMinus, deltaPlus);

        case "superbee"

            candidate1 = minmod_function( ...
                2*deltaMinus, deltaPlus);

            candidate2 = minmod_function( ...
                deltaMinus, 2*deltaPlus);

            sigma = maxmod_function(candidate1, candidate2);

        case "vanleer"

            sigma = zeros(size(deltaMinus));

            sameSign = deltaMinus.*deltaPlus > 0;

            sigma(sameSign) = ...
                2*deltaMinus(sameSign).*deltaPlus(sameSign) ...
                ./ ...
                (deltaMinus(sameSign) + deltaPlus(sameSign));

        case "vanalbada"

            sigma = zeros(size(deltaMinus));

            sameSign = deltaMinus.*deltaPlus > 0;

            numerator = ...
                deltaMinus(sameSign).*deltaPlus(sameSign) ...
                .* ...
                (deltaMinus(sameSign) + deltaPlus(sameSign));

            denominator = ...
                deltaMinus(sameSign).^2 ...
                + deltaPlus(sameSign).^2;

            sigma(sameSign) = numerator./denominator;

        otherwise

            error('Limitador desconocido: %s', limiterName);

    end
end


function value = minmod_function(a, b)
% Devuelve el argumento de menor módulo si ambos tienen el mismo signo.

    value = 0.5*(sign(a) + sign(b)).*min(abs(a), abs(b));
end


function value = maxmod_function(a, b)
% Devuelve el argumento de mayor módulo si ambos tienen el mismo signo.

    value = 0.5*(sign(a) + sign(b)).*max(abs(a), abs(b));
end