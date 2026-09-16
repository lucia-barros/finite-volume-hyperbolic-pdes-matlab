clear;
clc;
close all;

%% ================================================================
% MUSCL-HANCOCK VS ADER4
%
% PDE:
%
%       q_t + lambda*q_x = 0
%
% Dominio periódico:
%
%       x in [0,1]
%
% Condición inicial:
%
%       q(x,0) = sin(2*pi*x)
%% ================================================================

%% 1. Parámetros

xmin = 0;
xmax = 1;

N = 80;

lambda = 1.0;

CFL = 0.5;

% Número de vueltas completas de la onda
number_of_periods = 50;

% Como lambda = 1 y el dominio mide 1,
% una vuelta completa necesita t = 1.
t_final = number_of_periods / lambda;

%% 2. Malla

dx = (xmax - xmin) / N;

x = xmin + ((1:N) - 0.5) * dx;

%% 3. Condición inicial como media celular exacta

q_initial = exact_cell_average(x, dx, 0, lambda);

q_muscl = q_initial;
q_ader4 = q_initial;

%% 4. Integración temporal

t = 0;
iteration = 0;

tic;

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

total_cpu_time = toc;

%% 5. Solución exacta

q_exact = exact_cell_average( ...
    x, dx, t_final, lambda);

%% 6. Errores

error_muscl = q_muscl - q_exact;
error_ader4 = q_ader4 - q_exact;

L1_muscl = dx * sum(abs(error_muscl));
L1_ader4 = dx * sum(abs(error_ader4));

L2_muscl = sqrt(dx * sum(error_muscl.^2));
L2_ader4 = sqrt(dx * sum(error_ader4.^2));

Linf_muscl = max(abs(error_muscl));
Linf_ader4 = max(abs(error_ader4));

%% 7. Amplitudes

amplitude_exact = max(q_exact) - min(q_exact);

amplitude_muscl = max(q_muscl) - min(q_muscl);

amplitude_ader4 = max(q_ader4) - min(q_ader4);

relative_amplitude_muscl = ...
    amplitude_muscl / amplitude_exact;

relative_amplitude_ader4 = ...
    amplitude_ader4 / amplitude_exact;

%% 8. Mostrar resultados

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' MUSCL-HANCOCK VS ADER4\n');
fprintf('=========================================================\n');
fprintf('Número de celdas             : %d\n', N);
fprintf('dx                            : %.6e\n', dx);
fprintf('CFL                           : %.3f\n', CFL);
fprintf('Tiempo final                  : %.3f\n', t_final);
fprintf('Vueltas completas             : %d\n', number_of_periods);
fprintf('Iteraciones                   : %d\n', iteration);
fprintf('Tiempo total de CPU           : %.4f s\n', total_cpu_time);
fprintf('---------------------------------------------------------\n');
fprintf('%-18s %-13s %-13s %-13s\n', ...
    'Método', 'L1', 'L2', 'Linf');
fprintf('---------------------------------------------------------\n');
fprintf('%-18s %.5e   %.5e   %.5e\n', ...
    'MUSCL-Hancock', ...
    L1_muscl, L2_muscl, Linf_muscl);
fprintf('%-18s %.5e   %.5e   %.5e\n', ...
    'ADER4', ...
    L1_ader4, L2_ader4, Linf_ader4);
fprintf('---------------------------------------------------------\n');
fprintf('Amplitud relativa MUSCL       : %.6f\n', ...
    relative_amplitude_muscl);
fprintf('Amplitud relativa ADER4       : %.6f\n', ...
    relative_amplitude_ader4);
fprintf('=========================================================\n');

%% 9. Comparación de soluciones

figure;

plot(x, q_exact, ...
    'LineWidth', 2, ...
    'DisplayName', 'Solución exacta');

hold on;

plot(x, q_muscl, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock');

plot(x, q_ader4, ':', ...
    'LineWidth', 2, ...
    'DisplayName', 'ADER4');

hold off;

xlabel('x');
ylabel('Media celular de q');
title(sprintf( ...
    'Propagación después de %d periodos, N = %d', ...
    number_of_periods, N));

legend('Location', 'best');
grid on;
xlim([xmin, xmax]);

%% 10. Error puntual

figure;

semilogy(x, abs(error_muscl), ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Error MUSCL-Hancock');

hold on;

semilogy(x, abs(error_ader4), ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Error ADER4');

hold off;

xlabel('x');
ylabel('Error absoluto');
title('Error puntual después de la propagación');
legend('Location', 'best');
grid on;

%% 11. Diferencia directa entre métodos

figure;

plot(x, q_muscl - q_ader4, ...
    'LineWidth', 1.5);

xlabel('x');
ylabel('q_{MUSCL}-q_{ADER4}');
title('Diferencia entre MUSCL-Hancock y ADER4');
grid on;
xlim([xmin, xmax]);