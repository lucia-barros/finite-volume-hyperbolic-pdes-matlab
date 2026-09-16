clear;
clc;
close all;

%% ================================================================
% EVOLUCIÓN DEL ERROR CON EL NÚMERO DE PERIODOS
%
% Comparación:
%   - MUSCL-Hancock de segundo orden
%   - ADER4 de cuarto orden
%
% Ecuación:
%       q_t + lambda*q_x = 0
%
% Dominio periódico:
%       x in [0,1]
%
% Condición inicial:
%       q(x,0) = sin(2*pi*x)
%% ================================================================

%% 1. Parámetros del problema

xmin = 0;
xmax = 1;

N = 80;                  % Número de celdas
lambda = 1.0;            % Velocidad de advección
CFL = 0.5;               % Número CFL

% Número de vueltas que se quieren analizar
period_values = [1, 10, 25, 50];

if lambda == 0
    error('lambda debe ser distinto de cero.');
end

%% 2. Construcción de la malla

domain_length = xmax - xmin;

dx = domain_length / N;

x = xmin + ((1:N) - 0.5) * dx;

% Tiempo necesario para recorrer una vez todo el dominio
period_time = domain_length / abs(lambda);

%% 3. Condición inicial exacta como media celular

q_initial = exact_cell_average(x, dx, 0, lambda);

%% 4. Reservar memoria para los resultados

number_of_tests = length(period_values);

L1_muscl = zeros(number_of_tests, 1);
L1_ader4 = zeros(number_of_tests, 1);

L2_muscl = zeros(number_of_tests, 1);
L2_ader4 = zeros(number_of_tests, 1);

Linf_muscl = zeros(number_of_tests, 1);
Linf_ader4 = zeros(number_of_tests, 1);

relative_amplitude_muscl = zeros(number_of_tests, 1);
relative_amplitude_ader4 = zeros(number_of_tests, 1);

mass_error_muscl = zeros(number_of_tests, 1);
mass_error_ader4 = zeros(number_of_tests, 1);

%% 5. Masa inicial

initial_mass = dx * sum(q_initial);

%% 6. Bucle sobre los distintos números de periodos

for k = 1:number_of_tests

    number_of_periods = period_values(k);

    t_final = number_of_periods * period_time;

    % Reiniciar siempre desde la condición inicial
    q_muscl = q_initial;
    q_ader4 = q_initial;

    t = 0;
    iteration = 0;

    while t < t_final

        dt = CFL * dx / abs(lambda);

        % Ajustar el último paso
        if t + dt > t_final
            dt = t_final - t;
        end

        q_muscl = step_muscl_hancock( ...
            q_muscl, dx, dt, lambda);

        q_ader4 = step_ader4( ...
            q_ader4, dx, dt, lambda);

        t = t + dt;
        iteration = iteration + 1;
    end

    %% Solución exacta

    q_exact = exact_cell_average( ...
        x, dx, t_final, lambda);

    %% Errores puntuales

    error_muscl = q_muscl - q_exact;
    error_ader4 = q_ader4 - q_exact;

    %% Normas del error

    L1_muscl(k) = dx * sum(abs(error_muscl));
    L1_ader4(k) = dx * sum(abs(error_ader4));

    L2_muscl(k) = sqrt(dx * sum(error_muscl.^2));
    L2_ader4(k) = sqrt(dx * sum(error_ader4.^2));

    Linf_muscl(k) = max(abs(error_muscl));
    Linf_ader4(k) = max(abs(error_ader4));

    %% Conservación de la amplitud

    amplitude_exact = max(q_exact) - min(q_exact);

    amplitude_muscl = max(q_muscl) - min(q_muscl);
    amplitude_ader4 = max(q_ader4) - min(q_ader4);

    relative_amplitude_muscl(k) = ...
        amplitude_muscl / amplitude_exact;

    relative_amplitude_ader4(k) = ...
        amplitude_ader4 / amplitude_exact;

    %% Conservación de la masa

    mass_muscl = dx * sum(q_muscl);
    mass_ader4 = dx * sum(q_ader4);

    mass_error_muscl(k) = ...
        abs(mass_muscl - initial_mass);

    mass_error_ader4(k) = ...
        abs(mass_ader4 - initial_mass);

    %% Mostrar avance

    fprintf(['Periodos: %3d | Iteraciones: %5d | ', ...
             'L1 MUSCL: %.4e | L1 ADER4: %.4e\n'], ...
             number_of_periods, iteration, ...
             L1_muscl(k), L1_ader4(k));
end

%% 7. Cociente de errores

error_ratio = L1_muscl ./ L1_ader4;

%% 8. Tabla de resultados

results = table( ...
    period_values.', ...
    L1_muscl, ...
    L1_ader4, ...
    error_ratio, ...
    relative_amplitude_muscl, ...
    relative_amplitude_ader4, ...
    mass_error_muscl, ...
    mass_error_ader4, ...
    'VariableNames', { ...
        'Periods', ...
        'L1_MUSCL', ...
        'L1_ADER4', ...
        'ErrorRatio_MUSCL_ADER4', ...
        'RelativeAmplitude_MUSCL', ...
        'RelativeAmplitude_ADER4', ...
        'MassError_MUSCL', ...
        'MassError_ADER4' ...
    });

fprintf('\n');
fprintf('============================================================\n');
fprintf(' EVOLUCIÓN DEL ERROR CON EL NÚMERO DE PERIODOS\n');
fprintf('============================================================\n');

disp(results);

%% 9. Figura 1: error L1 frente al número de periodos

figure;

semilogy( ...
    period_values, ...
    L1_muscl, ...
    'o-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'MUSCL-Hancock');

hold on;

semilogy( ...
    period_values, ...
    L1_ader4, ...
    's-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'ADER4');

hold off;

xlabel('Número de periodos');
ylabel('Error L_1');
title('Acumulación del error durante la propagación');
legend('Location', 'best');
grid on;

%% 10. Figura 2: conservación de la amplitud

figure;

plot( ...
    period_values, ...
    relative_amplitude_muscl, ...
    'o-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'MUSCL-Hancock');

hold on;

plot( ...
    period_values, ...
    relative_amplitude_ader4, ...
    's-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'ADER4');

yline( ...
    1, ...
    '--', ...
    'Amplitud exacta', ...
    'LineWidth', 1.2);

hold off;

xlabel('Número de periodos');
ylabel('Amplitud relativa');
title('Conservación de la amplitud');
legend('Location', 'best');
grid on;

%% 11. Figura 3: ventaja relativa de ADER4

figure;

semilogy( ...
    period_values, ...
    error_ratio, ...
    'o-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7);

xlabel('Número de periodos');
ylabel('Error MUSCL / Error ADER4');
title('Diferencia de precisión entre los métodos');
grid on;

%% 12. Figura 4: conservación de la masa

figure;

semilogy( ...
    period_values, ...
    mass_error_muscl, ...
    'o-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'MUSCL-Hancock');

hold on;

semilogy( ...
    period_values, ...
    mass_error_ader4, ...
    's-', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'DisplayName', 'ADER4');

hold off;

xlabel('Número de periodos');
ylabel('Error absoluto de masa');
title('Conservación de la masa discreta');
legend('Location', 'best');
grid on;