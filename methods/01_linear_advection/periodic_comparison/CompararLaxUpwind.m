%% Comparacion Lax-Friedrichs y Upwind
clear;
clc;
close all;

%% 1. Parametros

a = 0;
b = 1;

lambda = 1;

m = 240;
dx = (b-a)/m;

% Usamos m puntos sin repetir el extremo b
x = a:dx:b-dx;

CFL = 0.9;

dt = CFL*dx/abs(lambda);

tmax = 1;

% Ajustamos el numero de pasos para alcanzar tmax
nt = round(tmax/dt);

% Reajustamos dt para terminar exactamente en tmax
dt = tmax/nt;

% CFL real despues del reajuste
CFLreal = lambda*dt/dx;

fprintf('dx = %.6f\n', dx);
fprintf('dt = %.6f\n', dt);
fprintf('Numero de pasos = %d\n', nt);
fprintf('CFL real = %.6f\n', CFLreal);

if abs(CFLreal) > 1
    warning('El esquema puede ser inestable porque |CFL| > 1');
end

%% 2. Condicion inicial sinusoidal

w0 = sin(2*pi*x);

%% 3. Inicializacion

w_lax = w0;
w_upwind = w0;

%% 4. Grafica inicial

figure;

plot(x, w0, 'k', 'LineWidth', 1.5);

xlabel('x');
ylabel('w(x,0)');
title('Condicion inicial: seno');
grid on;

pause(1);

%% 5. Evolucion temporal

for n = 1:nt

    % Lax-Friedrichs
    w_lax_nueva = lxf_periodico(w_lax, CFLreal);

    % Upwind
    w_upwind_nueva = upwind_periodico( ...
        w_upwind, lambda, dt, dx);

    % Tiempo actual
    t = n*dt;

    % Solucion exacta
    w_exacta = sin(2*pi*(x-lambda*t));

    % Representacion
    plot(x, w_exacta, 'k', 'LineWidth', 2);
    hold on;

    plot(x, w_lax_nueva, 'b--', 'LineWidth', 1.5);
    plot(x, w_upwind_nueva, 'r-.', 'LineWidth', 1.5);

    hold off;

    xlabel('x');
    ylabel('w(x,t)');
    title(['Comparacion en t = ', num2str(t, '%.3f')]);

    legend('Solucion exacta', ...
           'Lax-Friedrichs', ...
           'Upwind');

    grid on;
    ylim([-1.2 1.2]);

    drawnow;
    pause(0.05);

    % Actualizacion
    w_lax = w_lax_nueva;
    w_upwind = w_upwind_nueva;

end

%% 6. Errores finales

error_lax_L1 = dx*sum(abs(w_lax-w_exacta));
error_upwind_L1 = dx*sum(abs(w_upwind-w_exacta));

error_lax_Linf = max(abs(w_lax-w_exacta));
error_upwind_Linf = max(abs(w_upwind-w_exacta));

fprintf('\nErrores en t = %.2f\n', tmax);
fprintf('--------------------------\n');

fprintf('Lax-Friedrichs, error L1   = %.8e\n', ...
    error_lax_L1);

fprintf('Upwind, error L1           = %.8e\n', ...
    error_upwind_L1);

fprintf('Lax-Friedrichs, error Linf = %.8e\n', ...
    error_lax_Linf);

fprintf('Upwind, error Linf         = %.8e\n', ...
    error_upwind_Linf);