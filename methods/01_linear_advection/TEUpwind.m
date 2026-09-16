%% Esquema Upwind para la ecuacion de transporte
clear;
clc;
close all;

%% 1. Parametros del problema

a = -3;
b = 3;

lambda = -1;       % Velocidad de propagacion

m = 60;           % Numero de intervalos espaciales
dx = (b-a)/m;
x = a:dx:b;

dt = 0.05;
tmax = 1;
nt = round(tmax/dt);

CFL = lambda*dt/dx;

fprintf('dx = %.4f\n', dx);
fprintf('dt = %.4f\n', dt);
fprintf('CFL = %.4f\n', CFL);

if abs(CFL) > 1
    warning('El esquema puede ser inestable porque |CFL| > 1');
end

%% 2. Condicion inicial: pulso rectangular

w0 = zeros(size(x));
w0(abs(x) <= 0.5) = 1;

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

%% 5. Coeficientes Upwind

CFLpositivo = 0.5*(lambda + abs(lambda))*dt/dx;
CFLnegativo = 0.5*(lambda - abs(lambda))*dt/dx;

%% 6. Evolucion temporal

for n = 1:nt

    for i = 2:length(x)-1

        w_nueva(i) = w_actual(i) ...
            - CFLpositivo*(w_actual(i)-w_actual(i-1)) ...
            - CFLnegativo*(w_actual(i+1)-w_actual(i));

    end

    %% Condiciones de contorno transmisivas

    w_nueva(1) = w_nueva(2);
    w_nueva(end) = w_nueva(end-1);

    %% Tiempo actual

    t = n*dt;

    %% Solucion exacta

    w_exacta = zeros(size(x));
    w_exacta(abs(x-lambda*t) <= 0.5) = 1;

    %% Representacion

    plot(x, w_nueva, 'b-o', ...
         x, w_exacta, 'r--o');

    xlabel('x');
    ylabel('w(x,t)');
    title(['Upwind, t = ', num2str(t)]);
    legend('Solucion numerica', 'Solucion exacta');
    grid on;

    drawnow;
    pause(0.1);

    %% Actualizacion

    w_actual = w_nueva;

end