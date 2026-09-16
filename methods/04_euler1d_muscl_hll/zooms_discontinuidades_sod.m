%% Ampliación del contacto y del shock en el problema de Sod

clear;
clc;
close all;

%% 1. Cargar resultados

HLL = load('resultado_hll_primer_orden.mat');
MUSCL = load('resultado_muscl_hancock.mat');
L = load('resultados_limitadores_sod.mat');

x = HLL.x;

%% 2. Figura comparativa

figure;

% ================================================================
% Contacto: primer orden frente a MUSCL-Hancock
% ================================================================

subplot(2,2,1);

plot(x,HLL.rhoFinal,'--', ...
    'LineWidth',1.6, ...
    'DisplayName','HLL de primer orden');

hold on;

plot(x,MUSCL.rhoFinal,'-', ...
    'LineWidth',1.6, ...
    'DisplayName','MUSCL-Hancock-HLL');

xlim([0.62 0.74]);
ylim([0.22 0.46]);

xlabel('x');
ylabel('\rho');
title('Discontinuidad de contacto');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Shock: primer orden frente a MUSCL-Hancock
% ================================================================

subplot(2,2,2);

plot(x,HLL.rhoFinal,'--', ...
    'LineWidth',1.6, ...
    'DisplayName','HLL de primer orden');

hold on;

plot(x,MUSCL.rhoFinal,'-', ...
    'LineWidth',1.6, ...
    'DisplayName','MUSCL-Hancock-HLL');

xlim([0.80 0.90]);
ylim([0.10 0.30]);

xlabel('x');
ylabel('\rho');
title('Shock');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Contacto: comparación de limitadores
% ================================================================

subplot(2,2,3);

plot(L.resultados.minmod.x, ...
     L.resultados.minmod.rho, ...
     '-', ...
     'LineWidth',1.5, ...
     'DisplayName','Minmod');

hold on;

plot(L.resultados.vanleer.x, ...
     L.resultados.vanleer.rho, ...
     '--', ...
     'LineWidth',1.5, ...
     'DisplayName','Van Leer');

plot(L.resultados.superbee.x, ...
     L.resultados.superbee.rho, ...
     ':', ...
     'LineWidth',2.0, ...
     'DisplayName','Superbee');

xlim([0.62 0.74]);
ylim([0.22 0.46]);

xlabel('x');
ylabel('\rho');
title('Contacto: efecto del limitador');
legend('Location','best');
grid on;
hold off;

% ================================================================
% Shock: comparación de limitadores
% ================================================================

subplot(2,2,4);

plot(L.resultados.minmod.x, ...
     L.resultados.minmod.rho, ...
     '-', ...
     'LineWidth',1.5, ...
     'DisplayName','Minmod');

hold on;

plot(L.resultados.vanleer.x, ...
     L.resultados.vanleer.rho, ...
     '--', ...
     'LineWidth',1.5, ...
     'DisplayName','Van Leer');

plot(L.resultados.superbee.x, ...
     L.resultados.superbee.rho, ...
     ':', ...
     'LineWidth',2.0, ...
     'DisplayName','Superbee');

xlim([0.80 0.90]);
ylim([0.10 0.30]);

xlabel('x');
ylabel('\rho');
title('Shock: efecto del limitador');
legend('Location','best');
grid on;
hold off;

sgtitle('Ampliación de las discontinuidades del tubo de choque de Sod');

%% 3. Guardar figura

exportgraphics(gcf, ...
    'zooms_contacto_shock_sod.png', ...
    'Resolution',300);