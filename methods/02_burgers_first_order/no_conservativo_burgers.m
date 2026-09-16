function wn = no_conservativo_burgers(wa,lambda)
% NO_CONSERVATIVO_BURGERS
% Esquema upwind aplicado a la forma no conservativa:
%
%       w_t + w*w_x = 0
%
% Este método no conserva correctamente la velocidad
% de propagación de los choques.

N = length(wa);

wn = wa;

for j = 2:N-1

    if wa(j) >= 0
        % La información se desplaza hacia la derecha
        wn(j) = wa(j) ...
              - lambda*wa(j)*(wa(j)-wa(j-1));

    else
        % La información se desplaza hacia la izquierda
        wn(j) = wa(j) ...
              - lambda*wa(j)*(wa(j+1)-wa(j));
    end
end

% Condiciones de contorno transmisivas
wn(1) = wn(2);
wn(N) = wn(N-1);

end