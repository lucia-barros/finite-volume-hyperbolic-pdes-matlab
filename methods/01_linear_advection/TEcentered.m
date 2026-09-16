%% Esquema centrado para la ecuacion de transporte
clear;
clc;
close all;

%% 1. Parametros del problema

a = -3;
b = 3;

lambda = 1;

m = 60;
dx = (b-a)/m;
x = a:dx:b;

dt = 0.025;
tmax = 3;

nt = round(tmax/dt);

CFL = lambda*dt/dx;

fprintf('dx = %.4f\n', dx);
fprintf('dt = %.4f\n', dt);
fprintf('CFL = %.4f\n', CFL);

%% 2. Condicion inicial: pulso rectangular

w0 = sin(2*pi*x);

%% 3. Inicializacion

w_actual = w0;
w_nueva = zeros(size(w_actual));

%% 4. Representacion inicial

figure;

plot(x, w0, 'b-o');

xlabel('x');
ylabel('w(x,0)');
title('Condicion inicial');
grid on;

pause(1);

%% 5. Evolucion temporal

for n = 1:nt

    %% Esquema centrado en los nodos interiores

    for i = 2:length(x)-1

        w_nueva(i) = w_actual(i) ...
            - 0.5*CFL*(w_actual(i+1)-w_actual(i-1));

    end

    %% Condiciones de contorno transmisivas

    w_nueva(1) = w_nueva(2);
    w_nueva(end) = w_nueva(end-1);

    %% Tiempo actual

    t = n*dt;

    %% Solucion exacta

    w_exacta = sin(2*pi*(x-lambda*t));

    %% Representacion

    plot(x, w_nueva, 'b-o', ...
         x, w_exacta, 'r--o');

    xlabel('x');
    ylabel('w(x,t)');
    title(['Esquema centrado, t = ', num2str(t, '%.2f')]);

    legend('Solucion numerica', 'Solucion exacta');

    grid on;

    drawnow;
    pause(0.1);

    %% Actualizacion

    w_actual = w_nueva;

end