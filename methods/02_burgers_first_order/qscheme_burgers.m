function wn = qscheme_burgers(wa,lambda)
% QSCHEME_BURGERS
% Q-scheme de Roe/van Leer para Burgers.

N = length(wa);

% Flujo numérico en las interfaces
F = zeros(1,N-1);

for j = 1:N-1

    wl = wa(j);
    wr = wa(j+1);

    % Velocidad de Roe para Burgers
    a_roe = 0.5*(wl+wr);

    % Flujo numérico de Roe
    F(j) = 0.5*(0.5*wl^2 + 0.5*wr^2) ...
         - 0.5*abs(a_roe)*(wr-wl);
end

wn = wa;

for j = 2:N-1
    wn(j) = wa(j)-lambda*(F(j)-F(j-1));
end

% Condiciones de contorno transmisivas
wn(1) = wn(2);
wn(N) = wn(N-1);

end