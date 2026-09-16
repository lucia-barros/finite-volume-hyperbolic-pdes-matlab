%% PROYECTO MUSCL PARA LA ECUACIÓN DE BURGERS
% Comparación de métodos de volúmenes finitos
% Primera parte: Godunov de primer orden

clear;
clc;
close all;

%% Parámetros del problema

N = 200;            % Número de celdas
xmin = -1;          % Extremo izquierdo del dominio
xmax = 1;           % Extremo derecho del dominio
tFinal = 0.6;       % Tiempo final
CFL = 0.8;          % Número CFL

% Opciones disponibles:
% "shock"
% "rarefaction"
% "smooth"

caseName = "shock";

%% Resolución numérica con Godunov

[x, uNumerical, t, diagnostics] = solve_burgers_godunov( ...
    N, xmin, xmax, tFinal, CFL, caseName);

dx = (xmax - xmin)/N;

%% Condición inicial y solución exacta

uInitial = initial_condition(x, caseName);
uExact = exact_burgers_solution(x, t, caseName);

%% Cálculo de errores

if all(isfinite(uExact))
    errorL1 = dx*sum(abs(uNumerical - uExact));
    errorLinf = max(abs(uNumerical - uExact));

    fprintf('\nERROR NUMÉRICO\n');
    fprintf('Error L1      = %.6e\n', errorL1);
    fprintf('Error infinito = %.6e\n', errorLinf);
end

%% Comprobación de conservación

conservationResidual = max(abs( ...
    diagnostics.mass - diagnostics.expectedMass));

fprintf('\nDIAGNÓSTICOS\n');
fprintf('Número de celdas              = %d\n', N);
fprintf('Tiempo final alcanzado        = %.6f\n', t);
fprintf('Variación total inicial       = %.6f\n', diagnostics.TV(1));
fprintf('Variación total final         = %.6f\n', diagnostics.TV(end));
fprintf('Residuo máximo de conservación = %.6e\n', ...
    conservationResidual);

%% Figura 1: solución numérica

figure;

plot(x, uInitial, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Condición inicial');

hold on;

plot(x, uNumerical, 'o-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 3, ...
    'DisplayName', 'Godunov 1.er orden');

if all(isfinite(uExact))
    plot(x, uExact, '-', ...
        'LineWidth', 2, ...
        'DisplayName', 'Solución exacta');
end

grid on;
xlabel('x');
ylabel('u(x,t)');
title(sprintf( ...
    'Ecuación de Burgers: %s, t = %.2f, N = %d', ...
    caseName, t, N));

legend('Location', 'best');
xlim([xmin xmax]);

%% Figura 2: variación total

figure;

plot(diagnostics.time, diagnostics.TV, ...
    'LineWidth', 1.5);

grid on;
xlabel('Tiempo');
ylabel('TV(u)');
title('Evolución de la variación total');

%% Figura 3: conservación

figure;

plot(diagnostics.time, diagnostics.mass, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Masa numérica');

hold on;

plot(diagnostics.time, diagnostics.expectedMass, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Balance conservativo esperado');

grid on;
xlabel('Tiempo');
ylabel('Masa total');
title('Comprobación de conservación');
legend('Location', 'best');