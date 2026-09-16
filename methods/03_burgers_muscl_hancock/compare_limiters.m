%% COMPARACIÓN DE LIMITADORES MUSCL-HANCOCK
%
% Se comparan:
%   - minmod
%   - van Leer
%   - van Albada
%   - Superbee
%
% Se calculan:
%   - error L1
%   - error infinito
%   - variación total final
%   - máximo y mínimo numéricos
%   - residuo de conservación

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

limiters = [
    "minmod"
    "vanleer"
    "vanalbada"
    "superbee"
];

numberOfLimiters = numel(limiters);

dx = (xmax - xmin)/N;

%% Solución de referencia: Godunov

[x, uGodunov, ~, diagnosticGodunov] = ...
    solve_burgers_godunov( ...
    N, xmin, xmax, tFinal, CFL, caseName);

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

%% Ejecutamos todos los limitadores

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
fprintf('COMPARACIÓN DE LIMITADORES MUSCL-HANCOCK\n');
fprintf('========================================================\n');
fprintf('Caso: %s\n', caseName);
fprintf('N = %d, CFL = %.2f, tFinal = %.2f\n\n', ...
    N, CFL, tFinal);

fprintf('Godunov:\n');
fprintf('  Error L1   = %.8e\n', errorGodunovL1);
fprintf('  Error Linf = %.8e\n\n', errorGodunovLinf);

disp(resultsTable);

%% Figura 1: solución completa

figure;

plot(x, uExact, ...
    'LineWidth', 2.4, ...
    'DisplayName', 'Solución exacta');

hold on;

plot(x, uGodunov, '--', ...
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
    'Comparación de limitadores: %s, t = %.2f', ...
    caseName, tFinal));

legend('Location', 'best');
xlim([xmin xmax]);

%% Figura 2: zoom alrededor del shock

figure;

plot(x, uExact, ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Solución exacta');

hold on;

plot(x, uGodunov, 'o--', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 4, ...
    'DisplayName', 'Godunov');

for k = 1:numberOfLimiters

    plot(x, solutions(:,k), 'o-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4, ...
        'DisplayName', ...
        sprintf('MUSCL-%s', limiters(k)));

end

grid on;
xlabel('x');
ylabel('u(x,t)');
title('Ampliación de la discontinuidad');

legend('Location', 'eastoutside');

xlim([0.27 0.33]);
ylim([-0.05 1.05]);

%% Figura 3: error L1

figure;

bar(categorical(limiters), errorL1);

hold on;

yline(errorGodunovL1, '--', ...
    'Godunov', ...
    'LineWidth', 1.5);

grid on;
xlabel('Limitador');
ylabel('Error L^1');
title('Error L^1 para cada limitador');

%% Figura 4: evolución de la variación total

figure;
hold on;

plot( ...
    diagnosticGodunov.time, ...
    diagnosticGodunov.TV, ...
    '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Godunov');

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

ylim([0.99 1.01]);

%% Figura 5: comprobación de máximos y mínimos

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
legend('Mínimo numérico', 'Máximo numérico', ...
    'Location', 'best');