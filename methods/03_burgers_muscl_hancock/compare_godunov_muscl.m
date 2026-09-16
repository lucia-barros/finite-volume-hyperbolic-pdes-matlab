%% COMPARACIÓN GODUNOV VS MUSCL-HANCOCK

clear;
clc;
close all;

%% Parámetros

N = 200;
xmin = -1;
xmax = 1;
tFinal = 0.6;
CFL = 0.8;

caseName = "shock";
limiterName = "vanleer";

dx = (xmax - xmin)/N;

%% Godunov de primer orden

[x, uGodunov, tGodunov, diagnosticGodunov] = ...
    solve_burgers_godunov( ...
    N, xmin, xmax, tFinal, CFL, caseName);

%% MUSCL-Hancock

[~, uMUSCL, tMUSCL, diagnosticMUSCL] = ...
    solve_burgers_muscl( ...
    N, xmin, xmax, tFinal, CFL, ...
    caseName, limiterName);

%% Soluciones inicial y exacta

uInitial = initial_condition(x, caseName);
uExact = exact_burgers_solution(x, tFinal, caseName);

%% Errores

errorGodunovL1 = dx*sum(abs(uGodunov - uExact));
errorMUSCLL1 = dx*sum(abs(uMUSCL - uExact));

errorGodunovLinf = max(abs(uGodunov - uExact));
errorMUSCLLinf = max(abs(uMUSCL - uExact));

%% Conservación

residualGodunov = max(abs( ...
    diagnosticGodunov.mass ...
    - diagnosticGodunov.expectedMass));

residualMUSCL = max(abs( ...
    diagnosticMUSCL.mass ...
    - diagnosticMUSCL.expectedMass));

%% Mostrar resultados

fprintf('\n============================================\n');
fprintf('COMPARACIÓN GODUNOV VS MUSCL-HANCOCK\n');
fprintf('============================================\n');

fprintf('Caso: %s\n', caseName);
fprintf('Limitador MUSCL: %s\n', limiterName);
fprintf('Número de celdas: %d\n\n', N);

fprintf('Error L1 Godunov = %.8e\n', errorGodunovL1);
fprintf('Error L1 MUSCL   = %.8e\n\n', errorMUSCLL1);

fprintf('Error Linf Godunov = %.8e\n', errorGodunovLinf);
fprintf('Error Linf MUSCL   = %.8e\n\n', errorMUSCLLinf);

fprintf('Residuo conservación Godunov = %.8e\n', ...
    residualGodunov);

fprintf('Residuo conservación MUSCL   = %.8e\n', ...
    residualMUSCL);

fprintf('TV inicial       = %.12f\n', ...
    diagnosticMUSCL.TV(1));

fprintf('TV final Godunov = %.12f\n', ...
    diagnosticGodunov.TV(end));

fprintf('TV final MUSCL   = %.12f\n', ...
    diagnosticMUSCL.TV(end));

%% Figura 1: comparación global

figure;

plot(x, uInitial, '--', ...
    'LineWidth', 1.3, ...
    'DisplayName', 'Condición inicial');

hold on;

plot(x, uExact, '-', ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Solución exacta');

plot(x, uGodunov, 'o-', ...
    'LineWidth', 1.1, ...
    'MarkerSize', 3, ...
    'DisplayName', 'Godunov 1.er orden');

plot(x, uMUSCL, 's-', ...
    'LineWidth', 1.1, ...
    'MarkerSize', 3, ...
    'DisplayName', ...
    sprintf('MUSCL-%s', limiterName));

grid on;
xlabel('x');
ylabel('u(x,t)');

title(sprintf( ...
    'Burgers: Godunov vs MUSCL, %s, t = %.2f', ...
    caseName, tFinal));

legend('Location', 'best');
xlim([xmin xmax]);

%% Figura 2: ampliación alrededor del shock

figure;

plot(x, uExact, '-', ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Solución exacta');

hold on;

plot(x, uGodunov, 'o-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5, ...
    'DisplayName', 'Godunov');

plot(x, uMUSCL, 's-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5, ...
    'DisplayName', ...
    sprintf('MUSCL-%s', limiterName));

grid on;
xlabel('x');
ylabel('u(x,t)');
title('Ampliación de la discontinuidad');
legend('Location', 'best');

xlim([0.22 0.38]);
ylim([-0.05 1.05]);

%% Figura 3: variación total

figure;

plot( ...
    diagnosticGodunov.time, ...
    diagnosticGodunov.TV, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov');

hold on;

plot( ...
    diagnosticMUSCL.time, ...
    diagnosticMUSCL.TV, ...
    'LineWidth', 1.5, ...
    'DisplayName', sprintf('MUSCL-%s', limiterName));

grid on;
xlabel('Tiempo');
ylabel('TV(u)');
title('Comparación de la variación total');
legend('Location', 'best');

ylim([0.99 1.01]);

%% Figura 4: residuo de conservación

figure;

semilogy( ...
    diagnosticGodunov.time, ...
    abs(diagnosticGodunov.mass ...
        - diagnosticGodunov.expectedMass) + eps, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov');

hold on;

semilogy( ...
    diagnosticMUSCL.time, ...
    abs(diagnosticMUSCL.mass ...
        - diagnosticMUSCL.expectedMass) + eps, ...
    'LineWidth', 1.5, ...
    'DisplayName', sprintf('MUSCL-%s', limiterName));

grid on;
xlabel('Tiempo');
ylabel('Residuo de conservación');
title('Error del balance conservativo');
legend('Location', 'best');