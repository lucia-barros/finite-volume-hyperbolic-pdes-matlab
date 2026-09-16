function wn = godunov_burgers(wa,lambda)
% GODUNOV_BURGERS
% Un paso temporal del método de Godunov
% para la ecuación de Burgers.
%
% wa     : solución en el instante actual
% lambda : Delta t / Delta x
% wn     : solución en el siguiente instante

N = length(wa);

% Flujo numérico en las interfaces
F = zeros(1,N-1);

for j = 1:N-1

    wl = wa(j);
    wr = wa(j+1);

    if wl <= wr
        % Rarefacción

        if wl >= 0
            % Toda la rarefacción se mueve a la derecha
            F(j) = 0.5*wl^2;

        elseif wr <= 0
            % Toda la rarefacción se mueve a la izquierda
            F(j) = 0.5*wr^2;

        else
            % La rarefacción cruza la velocidad cero
            F(j) = 0;
        end

    else
        % Choque

        s = 0.5*(wl+wr);

        if s >= 0
            F(j) = 0.5*wl^2;
        else
            F(j) = 0.5*wr^2;
        end
    end
end

% Actualización conservativa
wn = wa;

for j = 2:N-1
    wn(j) = wa(j)-lambda*(F(j)-F(j-1));
end

% Condiciones de contorno transmisivas
wn(1) = wn(2);
wn(N) = wn(N-1);

end