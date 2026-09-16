%% Método HLL de primer orden para las ecuaciones de Euler 1D
% Problema del tubo de choque de Sod
%
% Variables conservadas:
% U = [rho; rho*u; E]
%
% Variables primitivas:
% W = [rho; u; p]

clear;
clc;
close all;

%% 1. Parámetros físicos y numéricos

gamma = 1.4;

xMin = 0.0;
xMax = 1.0;

N = 200;           % Número de celdas físicas
ng = 2;            % Celdas fantasma a cada lado

CFL = 0.45;
tFinal = 0.20;

dx = (xMax - xMin)/N;

% Centros de las celdas físicas
x = linspace(xMin + dx/2, xMax - dx/2, N);

% Número total de celdas, incluyendo las fantasma
Ntot = N + 2*ng;

% Índices correspondientes a las celdas físicas
iFis = ng + 1 : ng + N;

%% 2. Inicialización del problema de Sod

W = zeros(3, Ntot);

for j = 1:N

    i = j + ng;

    if x(j) < 0.5

        % Estado inicial izquierdo
        rho = 1.0;
        u   = 0.0;
        p   = 1.0;

    else

        % Estado inicial derecho
        rho = 0.125;
        u   = 0.0;
        p   = 0.1;

    end

    W(:,i) = [rho; u; p];

end

%% 3. Conversión inicial a variables conservadas

U = primToCons(W, gamma);

%% 4. Bucle temporal

t = 0;
iter = 0;

while t < tFinal

    iter = iter + 1;

    %% 4.1 Condiciones de contorno transmisivas

    U = applyTransmissiveBC(U, ng, N);

    %% 4.2 Conversión a variables primitivas

    W = consToPrim(U, gamma);

    %% 4.3 Cálculo del paso temporal mediante la condición CFL

    rho = W(1,iFis);
    vel = W(2,iFis);
    p   = W(3,iFis);

    % Velocidad local del sonido
    a = sqrt(gamma*p./rho);

    % Mayor velocidad de propagación
    sMax = max(abs(vel) + a);

    dt = CFL*dx/sMax;

    % Ajuste del último paso temporal
    if t + dt > tFinal
        dt = tFinal - t;
    end

    lambda = dt/dx;

    %% 4.4 Flujos HLL de primer orden

    Fnum = zeros(3, Ntot-1);

    for i = 1:Ntot-1

        % En primer orden, cada celda se considera constante.
        % Por tanto, se utilizan directamente los valores medios
        % de las dos celdas adyacentes a la interfaz.

        ULeftInterface  = U(:,i);
        URightInterface = U(:,i+1);

        Fnum(:,i) = hllFlux(ULeftInterface, ...
                            URightInterface, gamma);

    end

    %% 4.5 Actualización conservativa

    Unew = U;

    for i = ng+1:ng+N

        Unew(:,i) = U(:,i) ...
                  - lambda*(Fnum(:,i) - Fnum(:,i-1));

    end

    U = Unew;

    %% 4.6 Avance temporal

    t = t + dt;

    if mod(iter,20) == 0 || abs(t-tFinal) < 1e-12

        fprintf(['Iteración: %5d | Tiempo: %.6f | ', ...
                 'dt: %.3e | Smax: %.4f\n'], ...
                 iter, t, dt, sMax);

    end

end

%% 5. Recuperación de la solución final

U = applyTransmissiveBC(U, ng, N);
Wfinal = consToPrim(U, gamma);

rhoFinal = Wfinal(1,iFis);
uFinal   = Wfinal(2,iFis);
pFinal   = Wfinal(3,iFis);

EFinal = U(3,iFis);

%% 6. Representación gráfica

figure;

subplot(2,2,1);
plot(x, rhoFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('\rho');
title('Densidad');
grid on;

subplot(2,2,2);
plot(x, uFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('Velocidad');
grid on;

subplot(2,2,3);
plot(x, pFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('p');
title('Presión');
grid on;

subplot(2,2,4);
plot(x, EFinal, 'LineWidth', 1.5);
xlabel('x');
ylabel('E');
title('Energía total');
grid on;

sgtitle(sprintf(...
    'Tubo de choque de Sod: HLL de primer orden, N = %d, t = %.2f', ...
    N, tFinal));

%% 7. Guardado de la figura

exportgraphics(gcf, ...
    'sod_hll_primer_orden_N200.png', ...
    'Resolution', 300);

%% 8. Comprobaciones básicas

fprintf('\nSimulación terminada correctamente.\n');
fprintf('Número total de iteraciones: %d\n', iter);
fprintf('Densidad mínima: %.6e\n', min(rhoFinal));
fprintf('Presión mínima: %.6e\n', min(pFinal));

%% 9. Guardado de los resultados

save('resultado_hll_primer_orden.mat', ...
    'x', 'rhoFinal', 'uFinal', 'pFinal', ...
    'EFinal', 'N', 'tFinal');

%% ================================================================
% FUNCIONES LOCALES
% ================================================================

function U = primToCons(W, gamma)
% Convierte variables primitivas:
% W = [rho; u; p]
%
% en variables conservadas:
% U = [rho; rho*u; E]

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = p/(gamma-1) + 0.5*rho.*u.^2;

    U = [rho;
         rho.*u;
         E];

end

function W = consToPrim(U, gamma)
% Convierte variables conservadas:
% U = [rho; rho*u; E]
%
% en variables primitivas:
% W = [rho; u; p]

    rho = U(1,:);
    mom = U(2,:);
    E   = U(3,:);

    % Protección frente a divisiones por cero
    rho = max(rho, 1e-12);

    u = mom./rho;

    p = (gamma-1)*(E - 0.5*rho.*u.^2);

    % Protección frente a presiones negativas por redondeo
    p = max(p, 1e-12);

    W = [rho;
         u;
         p];

end

function F = physicalFlux(U, gamma)
% Calcula el flujo físico de las ecuaciones de Euler 1D:
%
% F(U) = [rho*u;
%         rho*u^2 + p;
%         u*(E+p)]

    W = consToPrim(U, gamma);

    rho = W(1,:);
    u   = W(2,:);
    p   = W(3,:);

    E = U(3,:);

    F = [rho.*u;
         rho.*u.^2 + p;
         u.*(E + p)];

end

function Fhll = hllFlux(UL, UR, gamma)
% Calcula el flujo numérico HLL entre dos estados conservados.

    WL = consToPrim(UL, gamma);
    WR = consToPrim(UR, gamma);

    rhoL = WL(1);
    uL   = WL(2);
    pL   = WL(3);

    rhoR = WR(1);
    uR   = WR(2);
    pR   = WR(3);

    % Velocidades del sonido
    aL = sqrt(gamma*pL/rhoL);
    aR = sqrt(gamma*pR/rhoR);

    % Estimación de las velocidades extremas
    SL = min(uL-aL, uR-aR);
    SR = max(uL+aL, uR+aR);

    % Flujos físicos izquierdo y derecho
    FL = physicalFlux(UL, gamma);
    FR = physicalFlux(UR, gamma);

    if SL >= 0

        % Todas las ondas se desplazan hacia la derecha
        Fhll = FL;

    elseif SR <= 0

        % Todas las ondas se desplazan hacia la izquierda
        Fhll = FR;

    else

        % La interfaz se encuentra entre ambas ondas
        Fhll = ...
            (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL);

    end

end

function U = applyTransmissiveBC(U, ng, N)
% Aplica condiciones de contorno transmisivas mediante
% la copia del estado de la celda física más próxima.

    firstPhysical = ng + 1;
    lastPhysical  = ng + N;

    % Celdas fantasma del borde izquierdo
    U(:,1:ng) = repmat(U(:,firstPhysical), 1, ng);

    % Celdas fantasma del borde derecho
    U(:,lastPhysical+1:lastPhysical+ng) = ...
        repmat(U(:,lastPhysical), 1, ng);

end