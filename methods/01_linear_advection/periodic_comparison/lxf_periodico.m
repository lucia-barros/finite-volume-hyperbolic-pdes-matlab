function w_nueva = lxf_periodico(w_actual, CFL)

    m = length(w_actual);

    w_derecha = zeros(size(w_actual));
    w_izquierda = zeros(size(w_actual));

    % Vecino derecho
    w_derecha(1:m-1) = w_actual(2:m);

    % Periodicidad: el vecino derecho del ultimo nodo
    % es el primer nodo
    w_derecha(m) = w_actual(1);

    % Vecino izquierdo
    w_izquierda(2:m) = w_actual(1:m-1);

    % Periodicidad: el vecino izquierdo del primer nodo
    % es el ultimo nodo
    w_izquierda(1) = w_actual(m);

    % Esquema de Lax-Friedrichs
    w_nueva = 0.5*(w_derecha + w_izquierda) ...
            - 0.5*CFL*(w_derecha - w_izquierda);

end