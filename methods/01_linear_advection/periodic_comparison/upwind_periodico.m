function w_nueva = upwind_periodico(w_actual, lambda, dt, dx)

    m = length(w_actual);

    w_derecha = zeros(size(w_actual));
    w_izquierda = zeros(size(w_actual));

    % Vecino derecho
    w_derecha(1:m-1) = w_actual(2:m);
    w_derecha(m) = w_actual(1);

    % Vecino izquierdo
    w_izquierda(2:m) = w_actual(1:m-1);
    w_izquierda(1) = w_actual(m);

    CFLpositivo = 0.5*(lambda + abs(lambda))*dt/dx;
    CFLnegativo = 0.5*(lambda - abs(lambda))*dt/dx;

    % Esquema Upwind
    w_nueva = w_actual ...
        - CFLpositivo*(w_actual-w_izquierda) ...
        - CFLnegativo*(w_derecha-w_actual);

end