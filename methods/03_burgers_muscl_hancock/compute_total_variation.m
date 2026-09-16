function TV = compute_total_variation(u, boundaryName)
%COMPUTE_TOTAL_VARIATION Calcula la variación total discreta.
%
% Para un dominio periódico se incluye también la diferencia entre
% la última y la primera celda.

    arguments
        u (:,1) double
        boundaryName (1,1) string = "transmissive"
    end

    TV = sum(abs(diff(u)));

    if lower(boundaryName) == "periodic"
        TV = TV + abs(u(1) - u(end));
    end
end