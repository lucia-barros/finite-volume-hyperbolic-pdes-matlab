function w = exact_burgers(x,t,xb,wl,wr)
% EXACT_BURGERS
% Solución entrópica exacta del problema de Riemann
% para la ecuación de Burgers:
%
%       w_t + (w^2/2)_x = 0
%
% Entradas:
%   x  : vector de posiciones
%   t  : tiempo
%   xb : posición inicial de la discontinuidad
%   wl : estado inicial izquierdo
%   wr : estado inicial derecho
%
% Salida:
%   w  : solución exacta evaluada en x y t

w = zeros(size(x));

% Evitamos dividir entre cero en el instante inicial
if t == 0
    w(x <= xb) = wl;
    w(x > xb)  = wr;
    return
end

if wl > wr
    % Caso 1: choque entrópico
    %
    % Velocidad de Rankine-Hugoniot:
    % s = (f(wr)-f(wl))/(wr-wl) = (wl+wr)/2

    s = 0.5*(wl + wr);

    w((x-xb)/t < s)  = wl;
    w((x-xb)/t >= s) = wr;

elseif wl < wr
    % Caso 2: onda de rarefacción

    xi = (x-xb)/t;

    w(xi < wl)  = wl;
    w(xi > wr)  = wr;

    zona_rarefaccion = (xi >= wl) & (xi <= wr);
    w(zona_rarefaccion) = xi(zona_rarefaccion);

else
    % Caso 3: solución constante
    w(:) = wl;
end

end