%% EXPERIMENTO 2: RAREFACCIÓN PARA LA ECUACIÓN DE BURGERS
%
% Problema de Riemann:
%
%       uL = 0
%       uR = 1
%
% La solución exacta contiene un abanico de rarefacción:
%
%       u(x,t) = 0,       x/t <= 0
%       u(x,t) = x/t,     0 < x/t < 1
%       u(x,t) = 1,       x/t >= 1

clear;
clc;
close all;

%% Parámetros

N = 200;
xmin = -1;
xmax = 1;
tFinal = 0.6;
CFL = 0.8;

caseName = "rarefaction";

limiters = [
    "minmod"
    "vanleer"
    "vanalbada"
    "superbee"
];

numberOfLimiters = numel(limiters);

dx = (xmax - xmin)/N;

%% Solución de Godunov de primer orden

[x, uGodunov, ~, diagnosticGodunov] = ...
    solve_burgers_godunov( ...
    N, xmin, xmax, tFinal, CFL, caseName);

%% Soluciones inicial y exacta

uInitial = initial_condition(x, caseName);
uExact = exact_burgers_solution(x, tFinal, caseName);

%% Reservamos memoria

solutions = zeros(N, numberOfLimiters);

errorL1 = zeros(numberOfLimiters, 1);
errorLinf = zeros(numberOfLimiters, 1);
finalTV = zeros(numberOfLimiters, 1);
minimumValue = zeros(numberOfLimiters, 1);
maximumValue = zeros(numberOfLimiters, 1);
conservationResidual = zeros(numberOfLimiters, 1);

diagnostics = cell(numberOfLimiters, 1);

%% Ejecutamos MUSCL para cada limitador

for k = 1:numberOfLimiters

    limiterName = limiters(k);

    [~, uMUSCL, ~, diagnostic] = ...
        solve_burgers_muscl( ...
        N, xmin, xmax, tFinal, CFL, ...
        caseName, limiterName);

    solutions(:,k) = uMUSCL;
    diagnostics{k} = diagnostic;

    errorL1(k) = dx*sum(abs(uMUSCL - uExact));
    errorLinf(k) = max(abs(uMUSCL - uExact));

    finalTV(k) = diagnostic.TV(end);

    minimumValue(k) = min(uMUSCL);
    maximumValue(k) = max(uMUSCL);

    conservationResidual(k) = max(abs( ...
        diagnostic.mass ...
        - diagnostic.expectedMass));
end

%% Resultados de Godunov

errorGodunovL1 = dx*sum(abs(uGodunov - uExact));
errorGodunovLinf = max(abs(uGodunov - uExact));

residualGodunov = max(abs( ...
    diagnosticGodunov.mass ...
    - diagnosticGodunov.expectedMass));

%% Tabla de resultados

resultsTable = table( ...
    limiters, ...
    errorL1, ...
    errorLinf, ...
    finalTV, ...
    minimumValue, ...
    maximumValue, ...
    conservationResidual, ...
    'VariableNames', { ...
        'Limitador', ...
        'Error_L1', ...
        'Error_Linf', ...
        'TV_final', ...
        'Minimo', ...
        'Maximo', ...
        'Residuo_conservacion'});

fprintf('\n========================================================\n');
fprintf('EXPERIMENTO DE RAREFACCIÓN\n');
fprintf('========================================================\n');
fprintf('N = %d, CFL = %.2f, tFinal = %.2f\n\n', ...
    N, CFL, tFinal);

fprintf('Godunov de primer orden:\n');
fprintf('  Error L1              = %.8e\n', errorGodunovL1);
fprintf('  Error Linf            = %.8e\n', errorGodunovLinf);
fprintf('  TV final              = %.12f\n', ...
    diagnosticGodunov.TV(end));
fprintf('  Residuo conservación  = %.8e\n\n', ...
    residualGodunov);

disp(resultsTable);

%% Figura 1: comparación completa

figure;

plot(x, uInitial, '--', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Condición inicial');

hold on;

plot(x, uExact, ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Solución exacta');

plot(x, uGodunov, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov 1.er orden');

for k = 1:numberOfLimiters

    plot(x, solutions(:,k), ...
        'LineWidth', 1.3, ...
        'DisplayName', ...
        sprintf('MUSCL-%s', limiters(k)));

end

grid on;
xlabel('x');
ylabel('u(x,t)');

title(sprintf( ...
    'Rarefacción de Burgers, t = %.2f, N = %d', ...
    tFinal, N));

legend('Location', 'eastoutside');
xlim([xmin xmax]);
ylim([-0.05 1.05]);

%% Figura 2: ampliación del abanico de rarefacción

figure;

plot(x, uExact, ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Solución exacta');

hold on;

plot(x, uGodunov, 'o-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 3, ...
    'DisplayName', 'Godunov');

for k = 1:numberOfLimiters

    plot(x, solutions(:,k), ...
        'LineWidth', 1.3, ...
        'DisplayName', ...
        sprintf('MUSCL-%s', limiters(k)));

end

grid on;
xlabel('x');
ylabel('u(x,t)');
title('Ampliación del abanico de rarefacción');

legend('Location', 'eastoutside');

xlim([-0.05 0.65]);
ylim([-0.05 1.05]);

%% Figura 3: error absoluto de cada método

figure;

plot(x, abs(uGodunov - uExact), ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov');

hold on;

for k = 1:numberOfLimiters

    plot(x, abs(solutions(:,k) - uExact), ...
        'LineWidth', 1.3, ...
        'DisplayName', ...
        sprintf('MUSCL-%s', limiters(k)));

end

grid on;
xlabel('x');
ylabel('|u_{num}-u_{exacta}|');
title('Distribución espacial del error absoluto');
legend('Location', 'eastoutside');

xlim([-0.1 0.7]);

%% Figura 4: comparación del error L1

figure;

allMethods = [
    "Godunov"
    "minmod"
    "vanleer"
    "vanalbada"
    "superbee"
];

allErrorsL1 = [
    errorGodunovL1
    errorL1
];

bar(categorical(allMethods), allErrorsL1);

grid on;
xlabel('Método');
ylabel('Error L^1');
title('Error L^1 en la rarefacción');

%% Figura 5: variación total

figure;

plot( ...
    diagnosticGodunov.time, ...
    diagnosticGodunov.TV, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov');

hold on;

for k = 1:numberOfLimiters

    plot( ...
        diagnostics{k}.time, ...
        diagnostics{k}.TV, ...
        'LineWidth', 1.3, ...
        'DisplayName', ...
        sprintf('MUSCL-%s', limiters(k)));

end

grid on;
xlabel('Tiempo');
ylabel('TV(u)');
title('Evolución de la variación total');
legend('Location', 'best');

ylim([0.98 1.02]);

%% Figura 6: comprobación de máximos y mínimos

figure;

barData = [minimumValue, maximumValue];

bar(categorical(limiters), barData);

hold on;

yline(0, '--', 'Mínimo exacto');
yline(1, '--', 'Máximo exacto');

grid on;
xlabel('Limitador');
ylabel('Valor numérico');
title('Comprobación de overshoots y undershoots');

legend( ...
    'Mínimo numérico', ...
    'Máximo numérico', ...
    'Location', 'best');