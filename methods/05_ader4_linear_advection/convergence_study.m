clear;
clc;
close all;

%% ================================================================
% ESTUDIO DE CONVERGENCIA
%
% MUSCL-Hancock vs ADER4
%% ================================================================

xmin = 0;
xmax = 1;

lambda = 1.0;
CFL = 0.5;

% Una vuelta completa
t_final = 1.0;

N_values = [20, 40, 80, 160, 320, 640];

number_of_meshes = length(N_values);

L1_muscl = zeros(number_of_meshes, 1);
L1_ader4 = zeros(number_of_meshes, 1);

L2_muscl = zeros(number_of_meshes, 1);
L2_ader4 = zeros(number_of_meshes, 1);

Linf_muscl = zeros(number_of_meshes, 1);
Linf_ader4 = zeros(number_of_meshes, 1);

cpu_muscl = zeros(number_of_meshes, 1);
cpu_ader4 = zeros(number_of_meshes, 1);

%% Bucle sobre las mallas

for k = 1:number_of_meshes

    N = N_values(k);

    dx = (xmax - xmin) / N;

    x = xmin + ((1:N) - 0.5) * dx;

    q_initial = exact_cell_average(x, dx, 0, lambda);

    q_muscl = q_initial;
    q_ader4 = q_initial;

    %% MUSCL-Hancock

    t = 0;

    tic;

    while t < t_final

        dt = CFL * dx / abs(lambda);

        if t + dt > t_final
            dt = t_final - t;
        end

        q_muscl = step_muscl_hancock( ...
            q_muscl, dx, dt, lambda);

        t = t + dt;
    end

    cpu_muscl(k) = toc;

    %% ADER4

    t = 0;

    tic;

    while t < t_final

        dt = CFL * dx / abs(lambda);

        if t + dt > t_final
            dt = t_final - t;
        end

        q_ader4 = step_ader4( ...
            q_ader4, dx, dt, lambda);

        t = t + dt;
    end

    cpu_ader4(k) = toc;

    %% Solución exacta

    q_exact = exact_cell_average( ...
        x, dx, t_final, lambda);

    %% Errores

    error_muscl = q_muscl - q_exact;
    error_ader4 = q_ader4 - q_exact;

    L1_muscl(k) = dx * sum(abs(error_muscl));
    L1_ader4(k) = dx * sum(abs(error_ader4));

    L2_muscl(k) = sqrt(dx * sum(error_muscl.^2));
    L2_ader4(k) = sqrt(dx * sum(error_ader4.^2));

    Linf_muscl(k) = max(abs(error_muscl));
    Linf_ader4(k) = max(abs(error_ader4));
end

%% Órdenes observados

order_L1_muscl = NaN(number_of_meshes, 1);
order_L1_ader4 = NaN(number_of_meshes, 1);

order_L2_muscl = NaN(number_of_meshes, 1);
order_L2_ader4 = NaN(number_of_meshes, 1);

order_Linf_muscl = NaN(number_of_meshes, 1);
order_Linf_ader4 = NaN(number_of_meshes, 1);

for k = 2:number_of_meshes

    refinement_ratio = N_values(k) / N_values(k-1);

    order_L1_muscl(k) = ...
        log(L1_muscl(k-1) / L1_muscl(k)) ...
        / log(refinement_ratio);

    order_L1_ader4(k) = ...
        log(L1_ader4(k-1) / L1_ader4(k)) ...
        / log(refinement_ratio);

    order_L2_muscl(k) = ...
        log(L2_muscl(k-1) / L2_muscl(k)) ...
        / log(refinement_ratio);

    order_L2_ader4(k) = ...
        log(L2_ader4(k-1) / L2_ader4(k)) ...
        / log(refinement_ratio);

    order_Linf_muscl(k) = ...
        log(Linf_muscl(k-1) / Linf_muscl(k)) ...
        / log(refinement_ratio);

    order_Linf_ader4(k) = ...
        log(Linf_ader4(k-1) / Linf_ader4(k)) ...
        / log(refinement_ratio);
end

%% Tabla principal

results = table( ...
    N_values.', ...
    L1_muscl, order_L1_muscl, ...
    L1_ader4, order_L1_ader4, ...
    cpu_muscl, cpu_ader4, ...
    'VariableNames', { ...
    'N', ...
    'L1_MUSCL', ...
    'Order_MUSCL', ...
    'L1_ADER4', ...
    'Order_ADER4', ...
    'CPU_MUSCL', ...
    'CPU_ADER4'});

disp(results);

%% Tabla completa de normas

complete_results = table( ...
    N_values.', ...
    L1_muscl, order_L1_muscl, ...
    L2_muscl, order_L2_muscl, ...
    Linf_muscl, order_Linf_muscl, ...
    L1_ader4, order_L1_ader4, ...
    L2_ader4, order_L2_ader4, ...
    Linf_ader4, order_Linf_ader4);

disp(complete_results);

%% Gráfica: error frente a resolución

figure;

loglog(N_values, L1_muscl, 'o-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock');

hold on;

loglog(N_values, L1_ader4, 's-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'ADER4');

hold off;

xlabel('Número de celdas N');
ylabel('Error L_1');
title('Convergencia espacial-temporal');
legend('Location', 'best');
grid on;

%% Gráfica: error frente a tiempo de CPU

figure;

loglog(cpu_muscl, L1_muscl, 'o-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock');

hold on;

loglog(cpu_ader4, L1_ader4, 's-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'ADER4');

hold off;

xlabel('Tiempo de CPU [s]');
ylabel('Error L_1');
title('Precisión frente a coste computacional');
legend('Location', 'best');
grid on;