%% Upwind expresado como metodo de volumenes finitos
clear;
clc;
close all;

%% 1. Parametros del problema

a = -6;
b = 6;

lambda = 1;

m = 12;                 % Numero de celdas
dx = (b-a)/m;

% Bordes de las celdas
x = a:dx:b;

% Centros de las celdas
xc = a + dx/2 : dx : b - dx/2;

CFL = 0.9;

% Paso temporal calculado a partir del CFL
dt = CFL*dx/abs(lambda);

tmax = 4;
nt = round(tmax/dt);

fprintf('dx = %.4f\n', dx);
fprintf('dt = %.4f\n', dt);
fprintf('CFL = %.4f\n', CFL);
fprintf('Numero de pasos temporales = %d\n', nt);

if CFL > 1
    warning('El esquema puede ser inestable porque CFL > 1');
end

%% 2. Condicion inicial

% Cada valor representa el promedio dentro de una celda
w0 = zeros(1,m);

% Pulso situado en la celda central
[~, indice_centro] = min(abs(xc));
w0(indice_centro) = 1;

%% 3. Inicializacion

w_actual = w0;
w_nueva = zeros(size(w_actual));

%% 4. Representacion de la condicion inicial

figure;

stairs(x, [w0 w0(end)], 'LineWidth', 1.5);

xlabel('x');
ylabel('w(x,0)');
title('Condicion inicial por volumenes finitos');
grid on;

ylim([-0.1 1.1]);

pause(1);

%% 5. Coeficientes Upwind

CFLpositivo = 0.5*(lambda + abs(lambda))*dt/dx;
CFLnegativo = 0.5*(lambda - abs(lambda))*dt/dx;

%% 6. Evolucion temporal

for n = 1:nt

    %% Actualizacion de las celdas interiores

    for i = 2:m-1

        w_nueva(i) = w_actual(i) ...
            - CFLpositivo*(w_actual(i)-w_actual(i-1)) ...
            - CFLnegativo*(w_actual(i+1)-w_actual(i));

    end

    %% Condiciones de contorno transmisivas

    if lambda > 0

        % La informacion entra por la izquierda
        w_nueva(1) = w_actual(1);

        % La informacion sale por la derecha
        w_nueva(m) = w_actual(m) ...
            - CFLpositivo*(w_actual(m)-w_actual(m-1));

    else

        % La informacion sale por la izquierda
        w_nueva(1) = w_actual(1) ...
            - CFLnegativo*(w_actual(2)-w_actual(1));

        % La informacion entra por la derecha
        w_nueva(m) = w_actual(m);

    end

    %% Tiempo actual

    t = n*dt;

    %% Solucion exacta aproximada en los centros de celda

    w_exacta = zeros(1,m);

    posicion_pulso = lambda*t;

    [~, indice_pulso] = min(abs(xc-posicion_pulso));

    if indice_pulso >= 1 && indice_pulso <= m
        w_exacta(indice_pulso) = 1;
    end

    %% Representacion

    stairs(x, [w_nueva w_nueva(end)], 'b', 'LineWidth', 1.5);
    hold on;

    stairs(x, [w_exacta w_exacta(end)], 'r--', 'LineWidth', 1.5);

    plot(xc, w_nueva, 'bo');
    plot(xc, w_exacta, 'ro');

    hold off;

    xlabel('x');
    ylabel('w(x,t)');
    title(['Upwind FV, t = ', num2str(t, '%.2f')]);
    legend('Solucion numerica', 'Solucion exacta');

    grid on;
    ylim([-0.1 1.1]);

    drawnow;
    pause(0.2);

    %% Actualizacion

    w_actual = w_nueva;

end

%% 7. Comprobacion de conservacion

masa_inicial = dx*sum(w0);
masa_final = dx*sum(w_actual);

error_conservacion = abs(masa_inicial-masa_final);

fprintf('\n');
fprintf('Comprobacion de conservacion\n');
fprintf('----------------------------\n');
fprintf('Masa inicial = %.10f\n', masa_inicial);
fprintf('Masa final   = %.10f\n', masa_final);
fprintf('Error de conservacion = %.10e\n', error_conservacion);