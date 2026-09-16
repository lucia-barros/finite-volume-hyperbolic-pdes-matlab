function wn = lax_friedrichs_burgers(wa,lambda)
% LAX_FRIEDRICHS_BURGERS
% Un paso temporal del esquema de Lax-Friedrichs
% para la ecuación de Burgers.

N = length(wa);

wn = zeros(size(wa));

% Flujo físico de Burgers
f = @(w) 0.5*w.^2;

for j = 2:N-1
    wn(j) = 0.5*(wa(j+1)+wa(j-1)) ...
          - 0.5*lambda*(f(wa(j+1))-f(wa(j-1)));
end

% Condiciones de contorno transmisivas
wn(1) = wn(2);
wn(N) = wn(N-1);

end