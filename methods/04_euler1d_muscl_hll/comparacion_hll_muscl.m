%% Comparación entre HLL de primer orden y MUSCL-Hancock-HLL
% Problema del tubo de choque de Sod

clear;
clc;
close all;

%% 1. Comprobación de los archivos

archivoHLL   = 'resultado_hll_primer_orden.mat';
archivoMUSCL = 'resultado_muscl_hancock.mat';

if ~isfile(archivoHLL)
    error(['No se encuentra el archivo ', archivoHLL, ...
           '. Ejecuta primero hll_primer_orden_sod.m']);
end

if ~isfile(archivoMUSCL)
    error(['No se encuentra el archivo ', archivoMUSCL, ...
           '. Ejecuta primero muscl_hancock_sod.m']);
end

%% 2. Carga de los resultados

HLL   = load(archivoHLL);
MUSCL = load(archivoMUSCL);

%% 3. Comprobación de compatibilidad

if HLL.N ~= MUSCL.N
    error('Las dos simulaciones deben utilizar el mismo número de celdas.');
end

if abs(HLL.tFinal - MUSCL.tFinal) > 1e-12
    error('Las dos simulaciones deben terminar en el mismo instante.');
end

if length(HLL.x) ~= length(MUSCL.x)
    error('Las mallas de las dos simulaciones no coinciden.');
end

x = HLL.x;

%% 4. Representación comparativa

figure;

% ------------------------------------------------------------
% Densidad
% ------------------------------------------------------------
subplot(2,2,1);

plot(x, HLL.rhoFinal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'HLL de primer orden');

hold on;

plot(x, MUSCL.rhoFinal, '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock-HLL');

xlabel('x');
ylabel('\rho');
title('Densidad');
legend('Location','best');
grid on;
hold off;

% ------------------------------------------------------------
% Velocidad
% ------------------------------------------------------------
subplot(2,2,2);

plot(x, HLL.uFinal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'HLL de primer orden');

hold on;

plot(x, MUSCL.uFinal, '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock-HLL');

xlabel('x');
ylabel('u');
title('Velocidad');
legend('Location','best');
grid on;
hold off;

% ------------------------------------------------------------
% Presión
% ------------------------------------------------------------
subplot(2,2,3);

plot(x, HLL.pFinal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'HLL de primer orden');

hold on;

plot(x, MUSCL.pFinal, '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock-HLL');

xlabel('x');
ylabel('p');
title('Presión');
legend('Location','best');
grid on;
hold off;

% ------------------------------------------------------------
% Energía total
% ------------------------------------------------------------
subplot(2,2,4);

plot(x, HLL.EFinal, '--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'HLL de primer orden');

hold on;

plot(x, MUSCL.EFinal, '-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'MUSCL-Hancock-HLL');

xlabel('x');
ylabel('E');
title('Energía total');
legend('Location','best');
grid on;
hold off;

sgtitle(sprintf(...
    ['Comparación para el tubo de choque de Sod: ', ...
     'N = %d, t = %.2f'], ...
    HLL.N, HLL.tFinal));

%% 5. Guardado de la figura

exportgraphics(gcf, ...
    'comparacion_hll_vs_muscl_N200.png', ...
    'Resolution', 300);

%% 6. Diferencias máximas entre los métodos

difRho = max(abs(MUSCL.rhoFinal - HLL.rhoFinal));
difU   = max(abs(MUSCL.uFinal   - HLL.uFinal));
difP   = max(abs(MUSCL.pFinal   - HLL.pFinal));
difE   = max(abs(MUSCL.EFinal   - HLL.EFinal));

fprintf('\nComparación HLL frente a MUSCL-Hancock\n');
fprintf('--------------------------------------\n');
fprintf('Número de celdas: %d\n', HLL.N);
fprintf('Tiempo final: %.4f\n', HLL.tFinal);
fprintf('Máxima diferencia en densidad: %.6e\n', difRho);
fprintf('Máxima diferencia en velocidad: %.6e\n', difU);
fprintf('Máxima diferencia en presión: %.6e\n', difP);
fprintf('Máxima diferencia en energía: %.6e\n', difE);