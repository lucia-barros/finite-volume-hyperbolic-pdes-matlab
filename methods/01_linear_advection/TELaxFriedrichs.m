%% Lax-Friedrichs para la ecuacion de transporte
clear;
clc;
close all;

%% 1. Parametros del problema

a = -3;           % Extremo izquierdo
b = 3;            % Extremo derecho
lambda = 1;       % Velocidad de propagacion

m = 60;           % Numero de intervalos espaciales
dx = (b-a)/m;     % Paso espacial
x = a:dx:b;       % Nodos de la malla

dt = 0.05;         % Paso temporal
tmax = 1;         % Tiempo final
nt = round(tmax/dt);

CFL = lambda*dt/dx;

fprintf('dx = %.4f\n', dx);
fprintf('dt = %.4f\n', dt);
fprintf('CFL = %.4f\n', CFL);

if abs(CFL) > 1
    warning('El esquema puede ser inestable porque |CFL| > 1');
end

%% 2. Condicion inicial: pulso

w0 = zeros(size(x));
w0(abs(x) <= 0.5) = 1;

%% 3. Inicializacion

w_actual = w0;
w_nueva = zeros(size(w_actual));

%% 4. Representacion de la condicion inicial

figure;

plot(x, w0, 'o-');
xlabel('x');
ylabel('w(x,0)');
title('Condicion inicial');
grid on;

%% 5. Evolucion temporal

for n = 1:nt

    % Actualizacion de los nodos interiores
    for i = 2:length(x)-1

        w_nueva(i) = 0.5*(w_actual(i+1) + w_actual(i-1)) ...
                   - 0.5*CFL*(w_actual(i+1) - w_actual(i-1));

    end

    % Condiciones de contorno transmisivas
    w_nueva(1) = w_nueva(2);
    w_nueva(end) = w_nueva(end-1);

    % Tiempo actual
    t = n*dt;

    % Solucion exacta
    w_exacta = zeros(size(x));
    w_exacta(abs(x-lambda*t) <= 0.5) = 1;
    % w_exacta = zeros(size(x));
    %posicion_pulso = lambda*t;

    %[~, indice_pulso] = min(abs(x-posicion_pulso));
    %w_exacta(indice_pulso) = 1;

    % Representacion
    plot(x, w_nueva, 'b-o', ...
         x, w_exacta, 'r--o');

    xlabel('x');
    ylabel('w(x,t)');
    title(['Lax-Friedrichs, t = ', num2str(t)]);
    legend('Solucion numerica', 'Solucion exacta');
    grid on;

    drawnow;
    pause(0.2);

    % Preparar el siguiente paso temporal
    w_actual = w_nueva;

end